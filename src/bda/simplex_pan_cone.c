/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_pan_cone.c
 * @brief Tangent-cone direction oracle for a deficient feasible point.
 *
 * With P={j:x_j>0}, solve mixed NNLS for C=[A;c'/scale], t=[0;-1],
 * allowing d_P to be free and requiring d_N>=0. An exact fit supplies a
 * feasible descent direction. A nonzero stationary residual supplies a dual
 * candidate. Neither outcome is trusted without checks on the original LP.
 */

#include "policy/pan.h"
#include "policy/controls.h"
#include "simplex_pan_cone.h"
#include "simplex_pan_config.h"
#include "simplex_pan_service.h"
#include "simplex_pan_basis.h"
#include "simplex_pan_phase1.h"
#include "simplex_pan_recovery.h"
#include "utils.h"

#include <lp_simplex/status.h>

#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static int pan_basis_create(struct simplex_PanBasis *basis,
	const struct simplex_PanStandard *standard,
	const struct lp_simplex_Options *options,
	const struct simplex_PanConfig *cfg)
{
	return simplex_pan_basis_create_configured(
		basis, standard, options->pivot_tolerance, cfg);
}

#define PAN_SUPPORT_FAIL \
	do { \
		failure_line = __LINE__; \
		goto restore; \
	} while (0)

int simplex_pan_cone_create(struct simplex_PanStandard *cone,
	const struct simplex_PanStandard *standard, const double scale,
	double *row_scale, int implicit)
{
	int j, k, position = 0, extra = 0;

	lp_simplex_memset(cone, 0, sizeof(*cone));
	cone->report = standard->report;
	for (j = 0; j < standard->columns; j++)
		extra += !implicit && standard->objective[j] != 0.;
	if (standard->rows == INT_MAX || standard->nonzeros > INT_MAX - extra)
		return lp_simplex_EXIT_FAILURE;
	cone->rows = standard->rows + !implicit;
	cone->columns = standard->columns;
	cone->nonzeros = standard->nonzeros + extra;
	cone->column_start = (int *)lp_simplex_malloc(
		((size_t)cone->columns + 1) * sizeof(int));
	cone->row_index =
		(int *)lp_simplex_malloc((size_t)cone->nonzeros * sizeof(int));
	cone->value = (double *)lp_simplex_malloc(
		(size_t)cone->nonzeros * sizeof(double));
	cone->normalized_value = (double *)lp_simplex_malloc(
		(size_t)cone->nonzeros * sizeof(double));
	cone->column_norm = (double *)lp_simplex_malloc(
		(size_t)cone->columns * sizeof(double));
	cone->rhs = (double *)lp_simplex_malloc(
		(size_t)cone->rows * sizeof(double));
	cone->objective = (double *)lp_simplex_malloc(
		(size_t)cone->columns * sizeof(double));
	if (!cone->column_start || !cone->row_index || !cone->value ||
		!cone->normalized_value || !cone->column_norm || !cone->rhs ||
		!cone->objective)
		return lp_simplex_EXIT_FAILURE;
	lp_simplex_memset(cone->rhs, 0, (size_t)cone->rows * sizeof(double));
	if (!implicit)
		cone->rhs[standard->rows] = -1.;
	for (j = 0; j < standard->rows; j++)
		row_scale[j] = 0.;
	for (k = 0; k < standard->nonzeros; k++) {
		int row = standard->row_index[k];

		row_scale[row] = fmax(row_scale[row], fabs(standard->value[k]));
	}
	for (j = 0; j < standard->rows; j++) {
		int exponent = 0;

		if (row_scale[j] > 0.)
			frexp(row_scale[j], &exponent);
		if (exponent < 1 - DBL_MAX_EXP)
			exponent = 1 - DBL_MAX_EXP;
		row_scale[j] = row_scale[j] > 0. ? ldexp(1., -exponent) : 1.;
	}
	for (j = 0; j < cone->columns; j++) {
		double norm = 0.;

		cone->objective[j] = standard->objective[j] / scale;
		cone->column_start[j] = position;
		for (k = standard->column_start[j];
			k < standard->column_start[j + 1]; k++) {
			cone->row_index[position] = standard->row_index[k];
			cone->value[position++] =
				row_scale[standard->row_index[k]] *
				standard->value[k];
			norm = hypot(norm, cone->value[position - 1]);
		}
		if (!implicit && standard->objective[j] != 0.) {
			double c = standard->objective[j] / scale;

			cone->row_index[position] = standard->rows;
			cone->value[position++] = c;
			norm = hypot(norm, c);
		}
		cone->column_norm[j] = norm;
		for (k = cone->column_start[j]; k < position; k++)
			cone->normalized_value[k] =
				norm > 0. ? cone->value[k] / norm : 0.;
	}
	cone->column_start[cone->columns] = position;
	simplex_pan_service_derive(cone, standard);
	return lp_simplex_EXIT_SUCCESS;
}

long double simplex_pan_cone_cost(
	const struct simplex_PanStandard *standard, const double *point)
{
	long double value = 0., correction = 0.;
	int j;

	for (j = 0; j < standard->columns; j++) {
		double product = standard->objective[j] * point[j];
		long double next = value + product;

		correction += (fabsl(value) >= fabs(product)
					      ? (value - next) + product
					      : (product - next) + value) +
			      fma(standard->objective[j], point[j], -product);
		value = next;
	}
	return value + correction;
}

static long double pan_cone_objective(
	const struct simplex_PanStandard *standard, const double *point)
{
	return standard->objective_offset +
	       simplex_pan_cone_cost(standard, point);
}

long double simplex_pan_cone_audit_cost(
	const struct simplex_PanStandard *standard, const double *point)
{
	struct simplex_PanStandard original = *standard;

	if (standard->audit_objective)
		original.objective = (double *)standard->audit_objective;
	return simplex_pan_cone_cost(&original, point);
}

static long double pan_cone_audit_objective(
	const struct simplex_PanStandard *a, const double *point)
{
	return (a->audit_objective ? a->audit_objective_offset
				   : a->objective_offset) +
	       simplex_pan_cone_audit_cost(a, point);
}

struct pan_SupportColumn {
	double activity;
	int column;
};

static int pan_support_order(const void *first, const void *second)
{
	const struct pan_SupportColumn *a =
		(const struct pan_SupportColumn *)first;
	const struct pan_SupportColumn *b =
		(const struct pan_SupportColumn *)second;

	if (a->activity > b->activity)
		return -1;
	if (a->activity < b->activity)
		return 1;
	return a->column < b->column ? -1 : (a->column > b->column);
}

int simplex_pan_support_compress_configured(const struct simplex_PanStandard *a,
	double *point, const struct lp_simplex_Options *options,
	const struct simplex_PanConfig *cfg)
{
	struct simplex_PanBasis b;
	double *storage = NULL, *saved, *zero, *coefficient, *residual, *work;
	unsigned char *remove = NULL;
	struct pan_SupportColumn *order = NULL;
	int *candidate = NULL;
	int i, j, p, initial = 0, result = -1, failure_line = 0, moves = 0;
	const double *cost =
		a->audit_objective ? a->audit_objective : a->objective;
	long double before = pan_cone_audit_objective(a, point);

	lp_simplex_memset(&b, 0, sizeof(b));
	storage = (double *)lp_simplex_malloc(
		(2 * (size_t)a->columns + 3 * (size_t)a->rows) *
		sizeof(double));
	remove = (unsigned char *)lp_simplex_malloc((size_t)a->rows);
	order = (struct pan_SupportColumn *)lp_simplex_malloc(
		(size_t)a->columns * sizeof(*order));
	candidate = (int *)lp_simplex_malloc((size_t)a->columns * sizeof(int));
	if (!storage || !remove || !order || !candidate)
		goto done;
	saved = storage;
	zero = saved + a->columns;
	coefficient = zero + a->columns;
	residual = coefficient + a->rows;
	work = residual + a->rows;
	lp_simplex_memcpy(saved, point, (size_t)a->columns * sizeof(double));
	lp_simplex_memset(zero, 0, (size_t)a->columns * sizeof(double));
	if (pan_basis_create(&b, a, options, cfg))
		PAN_SUPPORT_FAIL;
	/* Select an independent subset in bounded prefixes. Rank revelation may
	 * drop positive columns here only because the point is kept separately;
	 * no primal mass is deleted by this factor operation. */
	for (j = 0; j < a->columns; j++)
		if (point[j] > 0.) {
			order[initial].column = j;
			order[initial++].activity =
				point[j] * a->column_norm[j];
		}
	/* Select the physically strongest positive coordinates first. This
	 * changes only the independent representation, never the primal point.
	 * Every null exchange still undergoes the original-model audits. */
	qsort(order, (size_t)initial, sizeof(*order), pan_support_order);
	for (p = 0; p < initial; p++)
		candidate[p] = order[p].column;
#ifdef LP_SIMPLEX_HAVE_SPQR
	/* A mildly overcomplete positive support fits one rank-revealing QR.
	 * Avoid factoring a full prefix, then refactoring it with a few extra
	 * positive columns. Large support ratios retain bounded prefixes. */
	if ((size_t)initial <=
		(size_t)a->rows +
			(size_t)a->rows / SIMPLEX_PAN_SUPPORT_ROW_FRACTION +
			SIMPLEX_PAN_SUPPORT_ROW_SLACK) {
		if (simplex_pan_basis_select_columns(
			    &b, candidate, initial, zero))
			PAN_SUPPORT_FAIL;
	} else {
#endif
		for (p = 0; p < initial; p++) {
			j = order[p].column;
			if (b.count == a->rows &&
				simplex_pan_basis_repair_rank(&b, zero))
				PAN_SUPPORT_FAIL;
			if (b.count < a->rows && simplex_pan_basis_add(&b, j))
				PAN_SUPPORT_FAIL;
		}
		if (simplex_pan_basis_repair_rank(&b, zero))
			PAN_SUPPORT_FAIL;
#ifdef LP_SIMPLEX_HAVE_SPQR
	}
#endif
	if (initial == b.count) {
		result = 0;
		goto done;
	}
	for (p = 0; p < initial; p++) {
		j = order[p].column;
		if (point[j] > 0. && b.position[j] < 0) {
			double error, gamma = cost[j], alpha = point[j],
				      sign = -1.;
			int blocking = j;

			if (simplex_pan_basis_column_projection(
				    &b, j, coefficient, &error) ||
				!isfinite(error))
				PAN_SUPPORT_FAIL;
			if (error >
				options->pivot_tolerance * a->column_norm[j]) {
				/* Removing tied zero columns may expose an independent column
				 * previously outside the span. Admit it without changing x. */
				if (simplex_pan_basis_add(&b, j) ||
					simplex_pan_basis_factorize(&b))
					PAN_SUPPORT_FAIL;
				continue;
			}
			for (i = 0; i < b.count; i++)
				gamma -= cost[b.column[i]] * coefficient[i];
			if (!isfinite(gamma))
				PAN_SUPPORT_FAIL;
			/* d_j=sign, d_B=-sign*B^+a_j; choose nonincreasing cost. */
			if (gamma < 0.) {
				sign = 1.;
				alpha = HUGE_VAL;
				blocking = -1;
			}
			for (i = 0; i < b.count; i++)
				if (sign * coefficient[i] > 0.) {
					double ratio = point[b.column[i]] /
						       (sign * coefficient[i]);

					if (ratio < alpha) {
						alpha = ratio;
						blocking = b.column[i];
					}
				}
			if (blocking < 0 || !isfinite(alpha) || alpha <= 0. ||
				moves++ >= initial)
				PAN_SUPPORT_FAIL;
			point[j] = fma(alpha, sign, point[j]);
			for (i = 0; i < b.count; i++)
				point[b.column[i]] = fma(-alpha * sign,
					coefficient[i], point[b.column[i]]);
			point[blocking] = 0.;
			for (i = 0; i < a->columns; i++) {
				if (point[i] < 0. &&
					-point[i] * a->column_norm[i] <=
						options->primal_tolerance *
							SIMPLEX_PAN_SUPPORT_MASS_FRACTION /
							(1. + initial))
					point[i] = 0.;
				if (!isfinite(point[i]) || point[i] < 0.)
					PAN_SUPPORT_FAIL;
			}
			if (simplex_pan_residual(a, point, residual, work) >
					options->primal_tolerance ||
				!isfinite(pan_cone_audit_objective(a, point)) ||
				pan_cone_audit_objective(a, point) >
					before + 32. * DBL_EPSILON *
							 (1. + fabsl(before)))
				PAN_SUPPORT_FAIL;
			{
				for (i = 0; i < b.count; i++)
					remove[i] = point[b.column[i]] == 0.;
				simplex_pan_basis_remove_marked(&b, remove);
				if (point[j] > 0. && b.position[j] < 0 &&
					simplex_pan_basis_add(&b, j))
					PAN_SUPPORT_FAIL;
				if (simplex_pan_basis_factorize(&b) &&
					simplex_pan_basis_rebuild(&b)) {
					/* A very weak exchange may lose numerical rank. Select a new
					 * independent subset while retaining the entire primal point,
					 * then revisit excluded positive columns. Each null move has
					 * removed a positive coordinate, so this restart is finite. */
					if (simplex_pan_basis_repair_rank(
						    &b, zero))
						PAN_SUPPORT_FAIL;
					p = -1;
				}
			}
		}
	}
	result = 1;
	if (getenv(SIMPLEX_ENV_PROFILE))
		fprintf(stderr,
			"pan-bda support: positive=%d -> %d objective=%.12g\n",
			initial, b.count,
			(double)pan_cone_audit_objective(a, point));
	goto done;
restore:
	if (getenv(SIMPLEX_ENV_PROFILE))
		fprintf(stderr,
			"pan-bda support failure: line=%d positive=%d rank=%d\n",
			failure_line, initial, b.count);
	lp_simplex_memcpy(point, saved, (size_t)a->columns * sizeof(double));
done:
	simplex_pan_basis_destroy(&b);
	lp_simplex_free(storage);
	lp_simplex_free(remove);
	lp_simplex_free(order);
	lp_simplex_free(candidate);
	return result;
}

static void pan_cone_clean(const struct simplex_PanStandard *a, double *point,
	const struct lp_simplex_Options *options, long double objective,
	double *row_budget)
{
	const double *cost =
		a->audit_objective ? a->audit_objective : a->objective;
	double objective_used = 0.;
	double increase_used = 0.,
	       increase_budget = 16. * DBL_EPSILON * (1. + fabsl(objective));
	double objective_budget = options->dual_tolerance *
				  SIMPLEX_PAN_INNER_TOLERANCE_FRACTION *
				  (1. + fabsl(objective));
	double limit = options->primal_tolerance *
		       SIMPLEX_PAN_INNER_TOLERANCE_FRACTION;
	int j, k;

	lp_simplex_memset(row_budget, 0, (size_t)a->rows * sizeof(double));
	for (j = 0; j < a->columns; j++)
		if (point[j] > 0. && point[j] * a->column_norm[j] <= limit) {
			double objective_change = point[j] * fabs(cost[j]);
			int valid = 1;

			if (!isfinite(objective_change) ||
				objective_change >
					objective_budget - objective_used ||
				(cost[j] < 0. && objective_change >
							 increase_budget -
								 increase_used))
				continue;
			/* Bound the sum of absolute deleted contributions separately in each
			 * original row. Disjoint rows need not share a single global budget.
			 * Tentative accumulation also handles duplicate row entries safely. */
			for (k = a->column_start[j]; k < a->column_start[j + 1];
				k++)
				row_budget[a->row_index[k]] +=
					fabs(a->value[k] * point[j]);
			for (k = a->column_start[j]; k < a->column_start[j + 1];
				k++)
				if (row_budget[a->row_index[k]] > limit)
					valid = 0;
			if (valid) {
				objective_used += objective_change;
				if (cost[j] < 0.)
					increase_used += objective_change;
				point[j] = 0.;
			} else
				for (k = a->column_start[j];
					k < a->column_start[j + 1]; k++)
					row_budget[a->row_index[k]] -=
						fabs(a->value[k] * point[j]);
		}
}

int simplex_pan_cone_check_dual(const struct simplex_PanStandard *standard,
	const double *point, const double *multiplier, const double tolerance,
	double *violation)
{
	long double primal = pan_cone_objective(standard, point),
		    dual = standard->objective_offset;
	int i, j, valid = 1;
	*violation = 0.;
	for (i = 0; i < standard->rows; i++) {
		if (!isfinite(multiplier[i]))
			return 0;
		dual -= (long double)standard->rhs[i] * multiplier[i];
	}
	for (j = 0; j < standard->columns; j++) {
		double reduced =
			standard->objective[j] +
			simplex_pan_column_dot(standard, j, multiplier);

		if (!isfinite(reduced))
			return 0;
		*violation = fmax(*violation, -reduced);
		if (reduced < -tolerance ||
			(point[j] > 0. &&
				fabs(reduced) >
					tolerance *
						(1. + fabs(standard->objective
								      [j]))))
			valid = 0;
	}
	return valid && isfinite(primal) && isfinite(dual) &&
	       fabsl(primal - dual) <= tolerance * (1. + fabsl(primal));
}

static int pan_cone_dual(const struct simplex_PanStandard *standard,
	const double *point, const double *residual, const double denominator,
	const double scale, const double *row_scale, double *multiplier,
	const double tolerance, double *violation)
{
	int i;

	if (!isfinite(denominator) || denominator >= 0.)
		return 0;
	for (i = 0; i < standard->rows; i++)
		multiplier[i] =
			scale * (residual[i] / denominator) * row_scale[i];
	return simplex_pan_cone_check_dual(
		standard, point, multiplier, tolerance, violation);
}

/* A bounded first-order prediction seeds a coherent working set; it is never
 * accepted as a feasible direction. Mandatory free columns survive rank
 * revelation, then the ordinary mixed NNLS restores positivity and certifies
 * the residual. This avoids all-matrix completion. */
static int pan_cone_seed(struct simplex_PanBasis *basis,
	const double *prediction, const unsigned char *free_mask,
	double *required, int *iterations, int limit)
{
	int column[SIMPLEX_PAN_BULK_COLUMNS], count = 0, j, p,
					      cap = limit - *iterations;
	double score[SIMPLEX_PAN_BULK_COLUMNS];

	if (cap > SIMPLEX_PAN_BULK_COLUMNS)
		cap = SIMPLEX_PAN_BULK_COLUMNS;
	if (cap > basis->standard->rows - basis->count)
		cap = basis->standard->rows - basis->count;
	if (cap <= 0)
		return 0;
	for (j = 0; j < basis->standard->columns; j++) {
		double current =
			prediction[j] * basis->standard->column_norm[j];

		if (free_mask[j] || basis->position[j] >= 0 ||
			!isfinite(current) || current <= 0.)
			continue;
		p = count;
		if (p == cap) {
			if (current <= score[p - 1])
				continue;
			p--;
		} else
			count++;
		while (p > 0 && current > score[p - 1]) {
			score[p] = score[p - 1];
			column[p] = column[p - 1];
			p--;
		}
		score[p] = current;
		column[p] = j;
	}
	if (count == 0)
		return 0;
	for (p = 0; p < count; p++)
		if (simplex_pan_basis_add(basis, column[p]))
			return -1;
	*iterations += count;
	for (j = 0; j < basis->standard->columns; j++)
		required[j] = free_mask[j] ? 1. : 0.;
	if (simplex_pan_basis_repair_rank(basis, required)) {
		unsigned char *remove = (unsigned char *)lp_simplex_malloc(
			(size_t)basis->count);

		if (!remove)
			return -1;
		for (p = 0; p < basis->count; p++)
			remove[p] = !free_mask[basis->column[p]];
		simplex_pan_basis_remove_marked(basis, remove);
		lp_simplex_free(remove);
		if (simplex_pan_basis_rebuild(basis))
			return -1;
	}
	if (getenv(SIMPLEX_ENV_PROFILE))
		fprintf(stderr,
			"pan-bda first-order seed: selected=%d rank=%d/%d\n",
			count, basis->count, basis->standard->rows);
	return 0;
}

/* Stop the auxiliary NNLS only after certifying a useful finite step on
 * the original LP, with the first-order path's actual-step error budget. */
int simplex_pan_cone_trial_step(void *context, const double *direction)
{
	struct simplex_PanDirectionTrial *t =
		(struct simplex_PanDirectionTrial *)context;
	const struct simplex_PanStandard *a = t->a;
	double alpha = HUGE_VAL, base, error;
	long double before, after;
	int j, blocking = -1;

	if ((++t->calls % SIMPLEX_PAN_TRIAL_PERIOD) != 0 ||
		*t->iterations >= t->options->iteration_limit)
		return 0;
	(*t->iterations)++;
	t->checks++;
	base = simplex_pan_residual(a, t->point, t->residual, t->work);
	if (!isfinite(base) || base >= t->options->primal_tolerance ||
		simplex_pan_cone_audit_cost(a, direction) >= 0.)
		return 0;
	for (j = 0; j < a->columns; j++) {
		if (!isfinite(direction[j]))
			return 0;
		if (direction[j] < 0.) {
			double ratio = -t->point[j] / direction[j];

			if (ratio < alpha) {
				alpha = ratio;
				blocking = j;
			}
		}
	}
	if (blocking < 0 || !isfinite(alpha) || alpha <= 0.)
		return 0;
	for (j = 0; j < a->columns; j++) {
		t->trial[j] = j == blocking
				      ? 0.
				      : fma(alpha, direction[j], t->point[j]);
		if (!isfinite(t->trial[j]) || t->trial[j] < 0.)
			return 0;
	}
	error = simplex_pan_residual(a, t->trial, t->residual, t->work);
	if (error > base + SIMPLEX_PAN_TRIAL_MARGIN_FRACTION *
				    (t->options->primal_tolerance - base))
		return 0;
	before = pan_cone_audit_objective(a, t->point);
	after = pan_cone_audit_objective(a, t->trial);
	if (!isfinite(before) || !isfinite(after) ||
		before - after <= 64. * DBL_EPSILON * (1. + fabsl(before)))
		return 0;
	t->accepted++;
	if (getenv(SIMPLEX_ENV_PROFILE))
		fprintf(stderr,
			"pan-bda nnls trial: checks=%d alpha=%.6g gain=%.6g primal=%.6g iterations=%d\n",
			t->checks, alpha, (double)(before - after), error,
			*t->iterations);
	return 1;
}

static int pan_cone_recover_support(const struct simplex_PanStandard *a,
	double *point, const struct lp_simplex_Options *options,
	int *iterations, struct simplex_PanRecovery *recovery,
	const char *reason, const struct simplex_PanConfig *cfg)
{
	int j, count = 0, result;
	clock_t start;

	for (j = 0; j < a->columns; j++)
		count += point[j] > 0.;
	if (count > options->iteration_limit - *iterations)
		return -1;
	if (!simplex_pan_recovery_compress(recovery, point, a->columns))
		return 0;
	*iterations += count;
	start = clock();
	result =
		simplex_pan_support_compress_configured(a, point, options, cfg);
	recovery->recovery_seconds +=
		(double)(clock() - start) / CLOCKS_PER_SEC;
	if (getenv(SIMPLEX_ENV_PROFILE))
		fprintf(stderr,
			"pan-bda recovery: reason=%s action=support result=%d iterations=%d\n",
			reason, result, *iterations);
	return result;
}

static int pan_cone_phase2_backend(const struct simplex_PanStandard *standard,
	double *point, double *multiplier,
	const struct lp_simplex_Options *options, int *iterations, int *status,
	double *dual_infeasibility, int implicit, int epoch_limit,
	struct simplex_PanRecovery *recovery,
	const struct simplex_PanConfig *cfg)
{
	struct simplex_PanStandard cone;
	struct simplex_PanBasis basis;
	struct simplex_PanDirectionTrial trial;
	double *storage = NULL, *row_scale = NULL, *direction, *saved,
	       *residual, *audit, *correction;
	unsigned char *free_mask = NULL, *remove = NULL;
	double scale = 0.,
	       inner_tolerance = fmax(16. * DBL_EPSILON,
		       fmin(options->primal_tolerance,
			       options->dual_tolerance) *
			       SIMPLEX_PAN_INNER_TOLERANCE_FRACTION);
	int i, j, state, result = lp_simplex_EXIT_FAILURE, steps = 0,
			 retries = 0, keep_passive = 0, stalled = 0,
			 seed_failed = 0;
	int began = *iterations, allow_freeze = 1, frozen = 0, basis_dirty = 0,
	    try_first = 1, seed_pending = 0, warm = 0;
	double face_fraction = sqrt(DBL_EPSILON);
	int face_guard = (cfg->face_guard);
	clock_t start = clock();

	lp_simplex_memset(&cone, 0, sizeof(cone));
	lp_simplex_memset(&basis, 0, sizeof(basis));
	lp_simplex_memset(&trial, 0, sizeof(trial));
	for (j = 0; j < standard->columns; j++)
		if (!cfg->dual_center || point[j] > 0. ||
			standard->objective[j] < 0.)
			scale = fmax(scale, fabs(standard->objective[j]));
	if (scale == 0.) {
		lp_simplex_memset(
			multiplier, 0, (size_t)standard->rows * sizeof(double));
		*status = lp_simplex_Success;
		*dual_infeasibility = 0.;
		return lp_simplex_EXIT_SUCCESS;
	}
	*status = lp_simplex_MemoryAllocError;
	*dual_infeasibility = HUGE_VAL;
	if (cfg->first_order_enabled || cfg->dual_polish_enabled)
		lp_simplex_memset(
			multiplier, 0, (size_t)standard->rows * sizeof(double));
	if (simplex_pan_service_first_order(standard, point, multiplier,
		    options, iterations, dual_infeasibility, NULL, NULL, 0, cfg,
		    NULL)) {
		*status = lp_simplex_Success;
		return lp_simplex_EXIT_SUCCESS;
	}
	if (*iterations >= options->iteration_limit) {
		*status = lp_simplex_ExceedIterLimit;
		return lp_simplex_EXIT_FAILURE;
	}
	row_scale = (double *)lp_simplex_malloc(
		(size_t)(standard->rows > 0 ? standard->rows : 1) *
		sizeof(double));
	if (!row_scale || simplex_pan_cone_create(
				  &cone, standard, scale, row_scale, implicit))
		goto cleanup;
	storage = (double *)lp_simplex_malloc(
		((size_t)standard->columns * 2 +
			(size_t)(standard->rows + 1) * 3) *
		sizeof(double));
	free_mask =
		(unsigned char *)lp_simplex_malloc((size_t)standard->columns);
	remove = (unsigned char *)lp_simplex_malloc((size_t)cone.rows);
	if (!storage || !free_mask || !remove)
		goto cleanup;
	direction = storage;
	saved = direction + standard->columns;
	residual = saved + standard->columns;
	audit = residual + standard->rows + 1;
	correction = audit + standard->rows + 1;
	trial.a = standard;
	trial.options = options;
	trial.point = point;
	trial.trial = saved;
	trial.residual = audit;
	trial.work = correction;
	trial.iterations = iterations;
	lp_simplex_memcpy(
		saved, point, (size_t)standard->columns * sizeof(double));
	if (implicit)
		pan_cone_clean(standard, point, options,
			pan_cone_audit_objective(standard, point), audit);
	if (simplex_pan_residual(standard, point, audit, correction) >
		options->primal_tolerance)
		lp_simplex_memcpy(point, saved,
			(size_t)standard->columns * sizeof(double));
	if (pan_basis_create(&basis, &cone, options, cfg))
		goto cleanup;
	for (j = 0; j < standard->columns; j++)
		if (point[j] > 0. && simplex_pan_basis_add(&basis, j)) {
			seed_failed = 1;
			break;
		}
	if (seed_failed || simplex_pan_basis_rebuild(&basis)) {
		if (!face_guard ||
			pan_cone_recover_support(standard, point, options,
				iterations, recovery, "initial-rank", cfg) != 1)
			goto precision;
		simplex_pan_basis_destroy(&basis);
		if (pan_basis_create(&basis, &cone, options, cfg))
			goto cleanup;
		for (j = 0; j < standard->columns; j++)
			if (point[j] > 0. && simplex_pan_basis_add(&basis, j))
				goto precision;
		if (simplex_pan_basis_rebuild(&basis))
			goto precision;
	}
	if (cfg->cone_batch) {
		int previous = basis.count;

		if (simplex_pan_basis_complete(&basis, point)) {
			if (simplex_pan_basis_rebuild(&basis))
				goto precision;
		} else {
			int added = basis.count - previous;

			if (added > options->iteration_limit - *iterations) {
				*status = lp_simplex_ExceedIterLimit;
				goto cleanup;
			}
			if (added > 0)
				*iterations += added;
		}
		if (getenv(SIMPLEX_ENV_PROFILE))
			fprintf(stderr,
				"pan-bda cone batch: rank=%d->%d/%d seconds=%.6f\n",
				previous, basis.count, cone.rows,
				(double)(clock() - start) / CLOCKS_PER_SEC);
	}
	if (getenv(SIMPLEX_ENV_PROFILE))
		fprintf(stderr,
			"pan-bda cone start: backend=%s rank=%d/%d fill=%ld seconds=%.6f\n",
			implicit ? "implicit" : "explicit", basis.count,
			cone.rows, basis.factor_nonzeros,
			(double)(clock() - start) / CLOCKS_PER_SEC);
	for (;;) {
		double alpha = HUGE_VAL, slope = 0., error;
		int blocking = -1, from_first = 0;
		clock_t direction_start = clock();
		double direction_seconds = 0.;
		long double before;

		if (*iterations >= options->iteration_limit) {
			*status = lp_simplex_ExceedIterLimit;
			goto cleanup;
		}
		frozen = 0;
		{
			double maximum = 0.;

			if (face_guard && allow_freeze)
				for (j = 0; j < standard->columns; j++)
					maximum = fmax(maximum,
						point[j] * standard->column_norm
								   [j]);
			for (j = 0; j < standard->columns; j++) {
				free_mask[j] = point[j] > 0.;
				/* Restrict the next direction on a weak positive face instead
				 * of deleting primal mass: these coordinates may increase but
				 * may not decrease. Ax and x are entirely unchanged. A failed
				 * dual audit reopens the full tangent cone before any retry. */
				if (free_mask[j] && face_guard &&
					allow_freeze &&
					point[j] * standard->column_norm[j] <
						face_fraction * maximum) {
					free_mask[j] = 0;
					frozen++;
				}
			}
		}
		if (cfg->first_order_enabled && try_first) {
			int candidate;

			lp_simplex_memset(direction, 0,
				(size_t)standard->columns * sizeof(double));
			candidate = simplex_pan_service_first_order(standard,
				point, multiplier, options, iterations,
				dual_infeasibility, direction, free_mask, warm,
				cfg, NULL);
			direction_seconds =
				(double)(clock() - direction_start) /
				CLOCKS_PER_SEC;
			recovery->directions++;
			recovery->direction_seconds += direction_seconds;
			direction_start = clock();
			warm = 1;
			if (candidate == 1) {
				*status = lp_simplex_Success;
				result = lp_simplex_EXIT_SUCCESS;
				goto cleanup;
			}
			if (*iterations >= options->iteration_limit) {
				*status = lp_simplex_ExceedIterLimit;
				goto cleanup;
			}
			if (candidate == 2) {
				basis_dirty = 1;
				from_first = 1;
				goto have_direction;
			}
			try_first = 0;
			seed_pending = 1;
		}
		if (basis_dirty) {
			/* First-order directions may have dependent positive support.
			 * NNLS only requires the currently free face to be independent;
			 * weak positive coordinates can remain constrained outside it.
			 * Compress only if rebuilding that free face actually fails. */
			simplex_pan_basis_destroy(&basis);
			if (pan_basis_create(&basis, &cone, options, cfg))
				goto cleanup;
			basis_dirty = 0;
		}
		if (!keep_passive) {
			seed_failed = 0;
			for (i = 0; i < basis.count; i++)
				remove[i] = !free_mask[basis.column[i]];
			simplex_pan_basis_remove_marked(&basis, remove);
			for (j = 0; j < standard->columns; j++)
				if (free_mask[j] && basis.position[j] < 0. &&
					simplex_pan_basis_add(&basis, j)) {
					seed_failed = 1;
					break;
				}
			if (seed_failed ||
				simplex_pan_basis_factorize(&basis)) {
				if (!face_guard ||
					pan_cone_recover_support(standard,
						point, options, iterations,
						recovery, "face-rank",
						cfg) != 1)
					goto precision;
				simplex_pan_basis_destroy(&basis);
				if (pan_basis_create(
					    &basis, &cone, options, cfg))
					goto cleanup;
				continue;
			}
		}
		keep_passive = 0;
		if (seed_pending) {
			if (pan_cone_seed(&basis, direction, free_mask, saved,
				    iterations, options->iteration_limit))
				goto precision;
			seed_pending = 0;
		}
		state = simplex_pan_service_cone(&basis, implicit, direction,
			free_mask, residual, inner_tolerance,
			options->iteration_limit, iterations, status,
			cfg->first_order_enabled ? simplex_pan_cone_trial_step
						 : NULL,
			&trial, NULL);
		{
			double elapsed = (double)(clock() - direction_start) /
					 CLOCKS_PER_SEC;

			direction_seconds += elapsed;
			recovery->direction_seconds += elapsed;
			recovery->directions++;
			if (getenv(SIMPLEX_ENV_PROFILE) &&
				(elapsed > 1. || state < 0))
				fprintf(stderr,
					"pan-bda direction: backend=%s outcome=%d rank=%d fill=%ld seconds=%.6f iterations=%d\n",
					implicit ? "implicit" : "explicit",
					state, basis.count,
					basis.factor_nonzeros, elapsed,
					*iterations);
		}
		if (state < 0)
			goto cleanup;
		if (state == 1) {
			double norm = 0.;

			for (i = 0; i <= standard->rows; i++)
				norm = hypot(norm, residual[i]);
			if (pan_cone_dual(standard, point, residual,
				    residual[standard->rows], scale, row_scale,
				    multiplier, options->dual_tolerance,
				    dual_infeasibility) ||
				pan_cone_dual(standard, point, residual,
					-norm * norm, scale, row_scale,
					multiplier, options->dual_tolerance,
					dual_infeasibility)) {
				*status = lp_simplex_Success;
				result = lp_simplex_EXIT_SUCCESS;
				goto cleanup;
			}
			if (frozen > 0) {
				if (getenv(SIMPLEX_ENV_PROFILE))
					fprintf(stderr,
						"pan-bda face: reopen=%d steps=%d dual=%.6g\n",
						frozen, steps,
						*dual_infeasibility);
				allow_freeze = 0;
				try_first = 1;
				continue;
			}
			if (retries++ < SIMPLEX_PAN_INNER_RETRIES &&
				inner_tolerance > 16. * DBL_EPSILON) {
				inner_tolerance = fmax(16. * DBL_EPSILON,
					inner_tolerance *
						SIMPLEX_PAN_INNER_RETRY_FACTOR);
				if (simplex_pan_basis_rebuild(&basis))
					goto precision;
				keep_passive = 1;
				continue;
			}
			goto precision;
		}
have_direction:
		for (j = 0; j < standard->columns; j++) {
			if (!isfinite(direction[j]) ||
				(!free_mask[j] && direction[j] < 0.))
				goto precision;
			if (direction[j] < 0.) {
				double candidate = -point[j] / direction[j];

				if (candidate < alpha) {
					alpha = candidate;
					blocking = j;
				}
			}
		}
		slope = (double)simplex_pan_cone_audit_cost(
			standard, direction);
		if (!isfinite(slope) || slope >= 0.)
			goto precision;
		if (blocking < 0) {
			/* The mixed NNLS direction is nonnegative. Audit its unaugmented
			 * homogeneous equations before declaring a recession ray. */
			lp_simplex_memset(audit, 0,
				(size_t)standard->rows * sizeof(double));
			lp_simplex_memset(residual, 0,
				(size_t)standard->rows * sizeof(double));
			lp_simplex_memset(correction, 0,
				(size_t)standard->rows * sizeof(double));
			for (j = 0; j < standard->columns; j++) {
				int k;

				for (k = standard->column_start[j];
					k < standard->column_start[j + 1];
					k++) {
					double a = standard->value[k],
					       product = a * direction[j], next;

					i = standard->row_index[k];
					next = audit[i] + product;
					correction[i] +=
						(fabs(audit[i]) >= fabs(product)
								? (audit[i] -
									  next) +
									  product
								: (product -
									  next) +
									  audit[i]) +
						fma(a, direction[j], -product);
					audit[i] = next;
					residual[i] += fabs(product);
				}
			}
			for (i = 0; i < standard->rows; i++) {
				double e = fabs(audit[i] + correction[i]);

				if (!isfinite(e) ||
					e > options->primal_tolerance ||
					e > 32. * DBL_EPSILON *
							(standard->columns +
								2.) *
							residual[i])
					goto precision;
			}
			*status = lp_simplex_Unboundedness;
			goto cleanup;
		}
		if (!isfinite(alpha) || alpha <= 0.)
			goto precision;
		if (*iterations >= options->iteration_limit) {
			*status = lp_simplex_ExceedIterLimit;
			goto cleanup;
		}
		before = pan_cone_audit_objective(standard, point);
		lp_simplex_memcpy(saved, point,
			(size_t)standard->columns * sizeof(double));
		for (j = 0; j < standard->columns; j++) {
			point[j] = fma(alpha, direction[j], point[j]);
			if (j == blocking)
				point[j] = 0.;
			if (point[j] < 0. &&
				point[j] >= -options->primal_tolerance)
				point[j] = 0.;
			/* A QR direction can introduce positive roundoff dust on a zero
			 * face. Making those entries free on the next call creates nearly
			 * zero blocking steps. Remove only a globally bounded amount of
			 * constraint activity and objective within 0.1% of their tolerances;
			 * the original-model
			 * residual and objective audits below still decide acceptance. */
		}
		if (before - pan_cone_audit_objective(standard, point) <=
			64. * DBL_EPSILON * (1. + fabsl(before)))
			stalled++;
		else
			stalled = 0;
		/* Change the numerical face only after persistent arithmetic stalls;
		 * ordinary productive directions retain their original trajectory. */
		if (stalled >= SIMPLEX_PAN_STALLED_STEPS) {
			pan_cone_clean(standard, point, options, before, audit);
			stalled = 0;
		}
		error = simplex_pan_residual(
			standard, point, audit, correction);
		for (j = 0; j < standard->columns; j++)
			if (!isfinite(point[j]) || point[j] < 0.)
				error = HUGE_VAL;
		if (isfinite(error) &&
			error > SIMPLEX_PAN_RESTORE_TRIGGER_FRACTION *
					options->primal_tolerance) {
			clock_t repair_start = clock();
			int repaired;
			struct simplex_PanServiceTrigger trigger;

			recovery->restorations++;
			lp_simplex_memset(&trigger, 0, sizeof(trigger));
			trigger.kind = SIMPLEX_PAN_TRIGGER_PRIMAL_MARGIN;
			trigger.work = *iterations;
			trigger.primal = error;
			trigger.previous_primal = HUGE_VAL;
			repaired = simplex_pan_service_restore_standard(
				standard, point, options,
				simplex_pan_cone_audit_cost(standard, saved),
				iterations, cfg, &trigger, NULL);
			recovery->recovery_seconds +=
				(double)(clock() - repair_start) /
				CLOCKS_PER_SEC;
			if (repaired == 1) {
				recovery->restored++;
				basis_dirty = 1;
				error = simplex_pan_residual(
					standard, point, audit, correction);
			}
			if (getenv(SIMPLEX_ENV_PROFILE))
				fprintf(stderr,
					"pan-bda recovery: reason=primal-margin action=restore result=%d primal=%.6g iterations=%d\n",
					repaired, error, *iterations);
		}
		if (error > options->primal_tolerance || !isfinite(before) ||
			memcmp(point, saved,
				(size_t)standard->columns * sizeof(double)) ==
				0 ||
			!isfinite(pan_cone_audit_objective(standard, point)) ||
			pan_cone_audit_objective(standard, point) >
				before + 32. * DBL_EPSILON *
						 (1. + fabsl(before))) {
			lp_simplex_memcpy(point, saved,
				(size_t)standard->columns * sizeof(double));
			recovery->rejected++;
			if (getenv(SIMPLEX_ENV_PROFILE))
				fprintf(stderr,
					"pan-bda recovery: reason=trial-audit action=%s primal=%.6g alpha=%.6g rank=%d fill=%ld\n",
					from_first ? "nnls"
						   : "refine-or-switch",
					error, alpha, basis.count,
					basis.factor_nonzeros);
			if (from_first) {
				try_first = 0;
				seed_pending = 1;
				continue;
			}
			if (retries++ < SIMPLEX_PAN_INNER_RETRIES &&
				inner_tolerance > 16. * DBL_EPSILON) {
				inner_tolerance = fmax(16. * DBL_EPSILON,
					inner_tolerance *
						SIMPLEX_PAN_INNER_RETRY_FACTOR);
				if (simplex_pan_basis_rebuild(&basis))
					goto precision;
				keep_passive = 1;
				continue;
			}
			goto precision;
		}
		if (*iterations >= options->iteration_limit) {
			lp_simplex_memcpy(point, saved,
				(size_t)standard->columns * sizeof(double));
			*status = lp_simplex_ExceedIterLimit;
			goto cleanup;
		}
		retries = 0;
		allow_freeze = 1;
		try_first = 1;
		steps++;
		(*iterations)++;
		recovery->accepted++;
		if (simplex_pan_cone_audit_cost(standard, saved) -
				simplex_pan_cone_audit_cost(standard, point) <=
			simplex_pan_progress_threshold(standard, point,
				simplex_pan_cone_audit_cost(standard, saved),
				options->dual_tolerance,
				SIMPLEX_PAN_SMALL_GAIN_UNITS))
			recovery->low_progress++;
		else
			recovery->low_progress = 0;
		if (face_guard) {
			long double gain =
				simplex_pan_cone_audit_cost(standard, saved) -
				simplex_pan_cone_audit_cost(standard, point);

			if (gain <= simplex_pan_progress_threshold(standard,
					    point,
					    simplex_pan_cone_audit_cost(
						    standard, saved),
					    options->dual_tolerance,
					    SIMPLEX_PAN_FACE_SMALL_GAIN_UNITS))
				face_fraction =
					fmin(SIMPLEX_PAN_FACE_MAX_FRACTION,
						SIMPLEX_PAN_FACE_EXPANSION *
							face_fraction);
			else if (gain >
				 simplex_pan_progress_threshold(standard, point,
					 simplex_pan_cone_audit_cost(
						 standard, saved),
					 options->dual_tolerance,
					 SIMPLEX_PAN_FACE_LARGE_GAIN_UNITS))
				face_fraction = fmax(sqrt(DBL_EPSILON),
					face_fraction *
						SIMPLEX_PAN_FACE_CONTRACTION);
			/* Persistent small steps call for an actual feasible face change,
			 * not only a stricter restriction of the next direction. A failed
			 * compression restores x and leaves the ordinary path available. */
			if (recovery->low_progress >=
				SIMPLEX_PAN_LOW_PROGRESS_STEPS) {
				if (pan_cone_recover_support(standard, point,
					    options, iterations, recovery,
					    "low-progress", cfg) == 1)
					basis_dirty = 1;
				recovery->low_progress = 0;
			}
		}
		if (getenv(SIMPLEX_ENV_PROFILE) && (steps & (steps - 1)) == 0)
			fprintf(stderr,
				"pan-bda cone: steps=%d inner=%d rank=%d/%d objective=%.12g primal=%.3g alpha=%.3g frozen=%d source=%s seconds=%.6f\n",
				steps, *iterations - began - steps, basis.count,
				cone.rows,
				(double)pan_cone_audit_objective(
					standard, point),
				error, alpha, frozen,
				from_first ? "first-order" : "nnls",
				(double)(clock() - start) / CLOCKS_PER_SEC);
		if (getenv(SIMPLEX_ENV_PROFILE) &&
			((steps & (steps - 1)) == 0 ||
				recovery->low_progress > 0))
			fprintf(stderr,
				"pan-bda progress: blocker=%d slope=%.6g predicted=%.6g actual=%.6g low=%d rank=%d fill=%ld direction-seconds=%.6f\n",
				blocking, slope, -alpha * slope,
				(double)(before - pan_cone_audit_objective(
							  standard, point)),
				recovery->low_progress, basis.count,
				basis.factor_nonzeros, direction_seconds);
		if (epoch_limit > 0 && cfg->dual_center &&
			cfg->dual_polish_enabled && steps >= epoch_limit) {
			/* The controller survives this yield; recovery occurs after
			 * three actual low-progress steps, not every epoch boundary. */
			*status = lp_simplex_CondUnsatisfied;
			goto cleanup;
		}
	}
precision:
	*status = lp_simplex_PrecisionError;
cleanup:
	if (getenv(SIMPLEX_ENV_PROFILE))
		fprintf(stderr,
			"pan-bda cone summary: steps=%d inner=%d rank=%d/%d seconds=%.6f status=%d dual=%.6g\n",
			steps, *iterations - began - steps, basis.count,
			cone.rows, (double)(clock() - start) / CLOCKS_PER_SEC,
			*status, *dual_infeasibility);
	if (getenv(SIMPLEX_ENV_PROFILE))
		fprintf(stderr,
			"pan-bda cone qr: backend=%s factors=%ld updates=%ld fill=%ld anchors=%ld backward=%.3g orthogonal=%.3g\n",
			implicit ? "implicit" : "explicit",
			basis.factorizations, basis.updates,
			basis.factor_nonzeros, basis.anchor_rebuilds,
			basis.backward_error, basis.orthogonal_error);
	if (getenv(SIMPLEX_ENV_PROFILE))
		fprintf(stderr,
			"pan-bda cone timing: factor=%.6f least-squares=%.6f minimum-norm=%.6f projection=%.6f\n",
			basis.factor_seconds, basis.least_squares_seconds,
			basis.minimum_norm_seconds, basis.projection_seconds);
	if (getenv(SIMPLEX_ENV_PROFILE))
		fprintf(stderr,
			"pan-bda cone packing: packs=%ld packed-nodes=%ld border-solves=%ld\n",
			basis.numeric_packs, basis.numeric_packed_nodes,
			basis.border_solves);
	if (getenv(SIMPLEX_ENV_PROFILE))
		fprintf(stderr,
			"pan-bda nnls trial summary: calls=%d checks=%d accepted=%d\n",
			trial.calls, trial.checks, trial.accepted);
	simplex_pan_basis_destroy(&basis);
	simplex_pan_standard_destroy(&cone);
	lp_simplex_free(storage);
	lp_simplex_free(free_mask);
	lp_simplex_free(remove);
	lp_simplex_free(row_scale);
	return result;
}

static int pan_cone_dispatch(const struct simplex_PanStandard *standard,
	double *point, double *multiplier,
	const struct lp_simplex_Options *options, int *iterations, int *status,
	double *dual_infeasibility, int allow_epoch,
	struct simplex_PanRecovery *recovery,
	const struct simplex_PanConfig *cfg)
{
	int implicit = cfg->oracle == SIMPLEX_PAN_IMPLICIT_ORACLE;
	int epoch_limit = allow_epoch ? cfg->center_epoch_limit : 0;
	int result = pan_cone_phase2_backend(standard, point, multiplier,
		options, iterations, status, dual_infeasibility, implicit,
		epoch_limit, recovery, cfg);
	/* A rank-one solve is an optimization, not a reason to reject a valid
	 * feasible point. Fall back once to explicit augmented sparse QR, with
	 * the same point and the remaining shared iteration budget. */
	if (cfg->allow_oracle_recovery && implicit && result &&
		*status == lp_simplex_PrecisionError &&
		*iterations < options->iteration_limit) {
		if (getenv(SIMPLEX_ENV_PROFILE))
			fprintf(stderr,
				"pan-bda cone: implicit recovery -> explicit augmented QR\n");
		result = pan_cone_phase2_backend(standard, point, multiplier,
			options, iterations, status, dual_infeasibility, 0,
			epoch_limit, recovery, cfg);
	}
	return result;
}

static int pan_cone_phase2_controlled(
	const struct simplex_PanStandard *standard,
	struct simplex_PanBasis *initial_basis, double *point,
	double *multiplier, const struct lp_simplex_Options *options,
	int *iterations, int *status, double *dual_infeasibility,
	struct simplex_PanRecovery *recovery,
	const struct simplex_PanConfig *cfg)
{
	struct simplex_PanStandard centered;
	double *storage, *cost, *anchor, *right, error;
	long double shift = 0., correction = 0.;
	int i, j, result, have_anchor = 0, warm_anchor = 0, epochs = 0;

	if (!cfg->dual_center)
		return pan_cone_dispatch(standard, point, multiplier, options,
			iterations, status, dual_infeasibility, 0, recovery,
			cfg);
	storage = (double *)lp_simplex_malloc(
		((size_t)standard->columns + 2 * (size_t)standard->rows) *
		sizeof(double));
	if (!storage) {
		*status = lp_simplex_MemoryAllocError;
		return lp_simplex_EXIT_FAILURE;
	}
	cost = storage;
	anchor = cost + standard->columns;
	right = anchor + standard->rows;
	lp_simplex_memset(anchor, 0, (size_t)standard->rows * sizeof(double));
recenter:
	have_anchor = 0;
	shift = 0.;
	correction = 0.;
	if (cfg->dual_polish_enabled) {
		int began = *iterations;

		if (simplex_pan_service_first_order(standard, point, anchor,
			    options, iterations, dual_infeasibility, NULL, NULL,
			    warm_anchor, cfg, NULL) == 1) {
			lp_simplex_memcpy(multiplier, anchor,
				(size_t)standard->rows * sizeof(double));
			*status = lp_simplex_Success;
			lp_simplex_free(storage);
			return lp_simplex_EXIT_SUCCESS;
		}
		have_anchor = *iterations > began;
		for (i = 0; i < standard->rows; i++)
			if (!isfinite(anchor[i]))
				have_anchor = 0;
	}
	if (*iterations >= options->iteration_limit) {
		*status = lp_simplex_ExceedIterLimit;
		lp_simplex_free(storage);
		return lp_simplex_EXIT_FAILURE;
	}
	if (!have_anchor) {
		for (i = 0; i < initial_basis->count; i++)
			right[i] =
				-standard->objective[initial_basis->column[i]];
		if (simplex_pan_basis_minimum_norm(
			    initial_basis, right, anchor, &error))
			goto uncentered;
	}
	for (j = 0; j < standard->columns; j++) {
		cost[j] = standard->objective[j] +
			  simplex_pan_column_dot(standard, j, anchor);
		if (!isfinite(cost[j]))
			goto uncentered;
	}
	for (i = 0; i < standard->rows; i++) {
		double product = standard->rhs[i] * anchor[i];
		long double next = shift + product;

		correction += (fabsl(shift) >= fabs(product)
					      ? (shift - next) + product
					      : (product - next) + shift) +
			      fma(standard->rhs[i], anchor[i], -product);
		shift = next;
	}
	centered = *standard;
	centered.objective = cost;
	centered.audit_objective = standard->audit_objective
					   ? standard->audit_objective
					   : standard->objective;
	centered.audit_objective_offset =
		standard->audit_objective ? standard->audit_objective_offset
					  : standard->objective_offset;
	centered.objective_offset =
		(double)(standard->objective_offset - shift - correction);
	if (!isfinite(centered.objective_offset))
		goto uncentered;
	if (getenv(SIMPLEX_ENV_PROFILE))
		fprintf(stderr,
			"pan-bda cone: dual-centered objective offset=%.12g\n",
			centered.objective_offset);
	simplex_pan_service_derive(&centered, standard);
	result = pan_cone_dispatch(&centered, point, multiplier, options,
		iterations, status, dual_infeasibility, 1, recovery, cfg);
	if (result != lp_simplex_EXIT_SUCCESS &&
		*status == lp_simplex_CondUnsatisfied &&
		cfg->center_epoch_enabled && cfg->dual_polish_enabled) {
		if (*iterations >= options->iteration_limit) {
			*status = lp_simplex_ExceedIterLimit;
			lp_simplex_free(storage);
			return result;
		}
		for (i = 0; i < standard->rows; i++)
			anchor[i] += multiplier[i];
		warm_anchor = 1;
		epochs++;
		if (getenv(SIMPLEX_ENV_PROFILE))
			fprintf(stderr,
				"pan-bda center epoch: epochs=%d iterations=%d\n",
				epochs, *iterations);
		goto recenter;
	}
	if (result == lp_simplex_EXIT_SUCCESS) {
		for (i = 0; i < standard->rows; i++)
			multiplier[i] += anchor[i];
		if (!simplex_pan_cone_check_dual(standard, point, multiplier,
			    options->dual_tolerance, dual_infeasibility)) {
			*status = lp_simplex_PrecisionError;
			result = lp_simplex_EXIT_FAILURE;
		}
	}
	lp_simplex_free(storage);
	return result;
uncentered:
	lp_simplex_free(storage);
	return pan_cone_dispatch(standard, point, multiplier, options,
		iterations, status, dual_infeasibility, 0, recovery, cfg);
}

int simplex_pan_cone_phase2_configured(
	const struct simplex_PanStandard *standard,
	struct simplex_PanBasis *initial_basis, double *point,
	double *multiplier, const struct lp_simplex_Options *options,
	int *iterations, int *status, double *dual_infeasibility,
	const struct simplex_PanConfig *cfg)
{
	struct simplex_PanRecovery recovery;
	int result;

	lp_simplex_memset(&recovery, 0, sizeof(recovery));
	result = pan_cone_phase2_controlled(standard, initial_basis, point,
		multiplier, options, iterations, status, dual_infeasibility,
		&recovery, cfg);
	if (getenv(SIMPLEX_ENV_PROFILE))
		fprintf(stderr,
			"pan-bda recovery summary: directions=%ld accepted=%ld rejected=%ld support=%ld skipped=%ld restore=%ld/%ld direction-seconds=%.6f recovery-seconds=%.6f status=%d\n",
			recovery.directions, recovery.accepted,
			recovery.rejected, recovery.compressions,
			recovery.compression_skips, recovery.restored,
			recovery.restorations, recovery.direction_seconds,
			recovery.recovery_seconds, *status);
	return result;
}
