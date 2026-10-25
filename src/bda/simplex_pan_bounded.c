/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_pan_bounded.c
 * @brief Native-bound Pan/BDA solving, feasibility recovery and KKT checks.
 */

#include "policy/pan.h"
#include "policy/controls.h"
#include "simplex_pan_config.h"
#include "simplex_pan_service.h"
#include "simplex_progress.h"
#include "simplex_stage.h"
#include "simplex_pan_bounded.h"
#include "simplex_pan_storage.h"
#include "simplex_pan_basis.h"
#include "simplex_pan_phase1.h"
#include "simplex_pan_cone.h"
#include "simplex_pan_recovery.h"
#include "simplex_pan_face_recovery.h"
#include "simplex_dual.h"
#include "utils.h"

#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

int simplex_pan_bounded_create(
	struct simplex_PanBounded *b, const struct simplex_Problem *p)
{
	struct simplex_PanStandard *a = &b->a;
	int i, j, k, slacks = 0, pos = 0, column;

	lp_simplex_memset(b, 0, sizeof(*b));
	b->problem = p;
	a->report = p->report;
	for (i = 0; i < p->rows; i++)
		if (p->row_type[i] != optm_CONS_T_EQ &&
			!pan_storage_add(&slacks, 1))
			return -1;
	if (p->rows < 0 || p->columns < 0 || p->matrix.nonzeros < 0 ||
		p->columns > INT_MAX - slacks ||
		p->matrix.nonzeros > INT_MAX - slacks ||
		!pan_storage_fits(p->rows, p->columns + slacks, p->columns,
			p->matrix.nonzeros + slacks, 1, (size_t)-1))
		return -1;
	a->rows = p->rows;
	a->columns = p->columns + slacks;
	a->original_columns = p->columns;
	a->nonzeros = p->matrix.nonzeros + slacks;
	a->column_start = (int *)lp_simplex_malloc(
		((size_t)a->columns + 1) * sizeof(int));
	a->row_index = (int *)lp_simplex_malloc(
		((size_t)a->nonzeros + 1) * sizeof(int));
	a->value = (double *)lp_simplex_malloc(
		((size_t)a->nonzeros + 1) * sizeof(double));
	a->normalized_value = (double *)lp_simplex_malloc(
		((size_t)a->nonzeros + 1) * sizeof(double));
	a->column_norm = (double *)lp_simplex_malloc(
		((size_t)a->columns + 1) * sizeof(double));
	a->objective = (double *)lp_simplex_malloc(
		((size_t)a->columns + 1) * sizeof(double));
	a->rhs = (double *)lp_simplex_malloc(
		((size_t)a->rows + 1) * sizeof(double));
	b->lower = (double *)lp_simplex_malloc(
		((size_t)a->columns + 1) * sizeof(double));
	b->upper = (double *)lp_simplex_malloc(
		((size_t)a->columns + 1) * sizeof(double));
	if (!a->column_start || !a->row_index || !a->value ||
		!a->normalized_value || !a->column_norm || !a->objective ||
		!a->rhs || !b->lower || !b->upper)
		goto fail;
	lp_simplex_memcpy(a->rhs, p->rhs, (size_t)p->rows * sizeof(double));
	for (j = 0; j < p->columns; j++) {
		const struct optm_VariableBound *bound = p->bounds + j;

		b->lower[j] = (bound->b_type == optm_BOUND_T_LO ||
				      bound->b_type == optm_BOUND_T_BS)
				      ? bound->lb
				      : -HUGE_VAL;
		b->upper[j] = (bound->b_type == optm_BOUND_T_UP ||
				      bound->b_type == optm_BOUND_T_BS)
				      ? bound->ub
				      : HUGE_VAL;
		if (isnan(b->lower[j]) || isnan(b->upper[j]) ||
			b->lower[j] > b->upper[j] || b->lower[j] == HUGE_VAL ||
			b->upper[j] == -HUGE_VAL || !isfinite(p->objective[j]))
			goto fail;
		a->column_start[j] = pos;
		a->objective[j] = p->objective[j];
		for (k = p->matrix.column_start[j];
			k < p->matrix.column_start[j + 1]; k++) {
			a->row_index[pos] = p->matrix.row_index[k];
			a->value[pos++] = p->matrix.value[k];
		}
	}
	column = p->columns;
	for (i = 0; i < p->rows; i++)
		if (p->row_type[i] != optm_CONS_T_EQ) {
			a->column_start[column] = pos;
			a->row_index[pos] = i;
			a->value[pos++] =
				p->row_type[i] == optm_CONS_T_LE ? 1. : -1.;
			a->objective[column] = 0.;
			b->lower[column] = 0.;
			b->upper[column] = HUGE_VAL;
			column++;
		}
	a->column_start[a->columns] = pos;
	for (j = 0; j < a->columns; j++) {
		double norm = 0.;

		for (k = a->column_start[j]; k < a->column_start[j + 1]; k++)
			norm = hypot(norm, a->value[k]);
		a->column_norm[j] = norm;
		for (k = a->column_start[j]; k < a->column_start[j + 1]; k++)
			a->normalized_value[k] =
				norm > 0. ? a->value[k] / norm : 0.;
	}
	return 0;
fail:
	simplex_pan_bounded_destroy(b);
	return -1;
}

void simplex_pan_bounded_destroy(struct simplex_PanBounded *b)
{
	simplex_pan_standard_destroy(&b->a);
	lp_simplex_free(b->lower);
	lp_simplex_free(b->upper);
	lp_simplex_memset(b, 0, sizeof(*b));
}

double simplex_pan_bounded_primal(const struct simplex_PanBounded *b,
	const double *x, double *residual, double *work)
{
	double error = simplex_pan_residual(&b->a, x, residual, work);
	int j;

	for (j = 0; j < b->a.columns; j++) {
		if (!isfinite(x[j]))
			return HUGE_VAL;
		error = fmax(
			error, fmax(b->lower[j] - x[j], x[j] - b->upper[j]));
	}
	return error;
}

int simplex_pan_bounded_encode(const struct simplex_PanBounded *b,
	const double *x, double *point, double *residual, double *work,
	double tolerance)
{
	const struct simplex_PanStandard *a = &b->a;
	int j, k;

	for (j = 0; j < a->columns; j++)
		point[j] = 0.;
	for (j = 0; j < a->original_columns; j++) {
		if (!isfinite(x[j]) || x[j] < b->lower[j] - tolerance ||
			x[j] > b->upper[j] + tolerance)
			return -1;
		point[j] = fmin(b->upper[j], fmax(b->lower[j], x[j]));
	}
	(void)simplex_pan_residual(a, point, residual, work);
	for (j = a->original_columns; j < a->columns; j++) {
		k = a->column_start[j];
		point[j] = fmax(0., residual[a->row_index[k]] / a->value[k]);
	}
	return simplex_pan_bounded_primal(b, point, residual, work) <= tolerance
		       ? 0
		       : -1;
}

static void native_product(
	long double *sum, long double *correction, double a, double b)
{
	double product = a * b;
	long double next = *sum + product;
	*correction +=
		(fabsl(*sum) >= fabs(product) ? (*sum - next) + product
					      : (product - next) + *sum) +
		fma(a, b, -product);
	*sum = next;
}

int simplex_pan_bounded_check_dual(const struct simplex_PanBounded *b,
	const double *x, const double *y, double tolerance, double *violation)
{
	const struct simplex_PanStandard *a = &b->a;
	long double primal = simplex_pan_cone_cost(a, x), dual = 0.,
		    correction = 0.;
	int i, j, valid = 1;
	*violation = 0.;
	for (i = 0; i < a->rows; i++) {
		if (!isfinite(y[i]))
			return 0;
		native_product(&dual, &correction, -a->rhs[i], y[i]);
	}
	for (j = 0; j < a->columns; j++) {
		double g = a->objective[j] + simplex_pan_column_dot(a, j, y),
		       v = 0.;
		double budget = tolerance * (1. + fabs(a->objective[j]));

		if (!isfinite(g) || !isfinite(x[j]))
			return 0;
		if (b->lower[j] == b->upper[j])
			native_product(&dual, &correction, b->lower[j], g);
		else {
			if (x[j] == b->lower[j])
				v = fmax(0., -g);
			else if (x[j] == b->upper[j])
				v = fmax(0., g);
			else
				v = fabs(g);
			if (v > budget)
				valid = 0;
			*violation = fmax(*violation, v);
			if (g >= 0. && isfinite(b->lower[j]))
				native_product(
					&dual, &correction, b->lower[j], g);
			else if (g < 0. && isfinite(b->upper[j]))
				native_product(
					&dual, &correction, b->upper[j], g);
		}
	}
	dual += correction;
	return valid && isfinite(primal) && isfinite(dual) &&
	       fabsl(primal - dual) <= tolerance * (1. + fabsl(primal));
}

/* Search for a multiplier on the same feasible face. The gap only admits
 * this bounded repair; it never certifies optimality. Every accepted repair
 * passes the unchanged componentwise KKT and gap audit. */
int simplex_pan_bounded_repair_dual(const struct simplex_PanBounded *b,
	const double *x, double *y, double tolerance, int limit,
	int *iterations, double *violation)
{
	const struct simplex_PanStandard *a = &b->a;
	long double primal, dual = 0., correction = 0.;
	double *saved;
	int i, j, k, pass;

	if (simplex_pan_bounded_check_dual(b, x, y, tolerance, violation))
		return 1;
	if (*iterations >= limit)
		return 0;
	primal = simplex_pan_cone_cost(a, x);
	for (i = 0; i < a->rows; i++) {
		if (!isfinite(y[i]))
			return 0;
		native_product(&dual, &correction, -a->rhs[i], y[i]);
	}
	for (j = 0; j < a->columns; j++) {
		double g = a->objective[j] + simplex_pan_column_dot(a, j, y);

		if (!isfinite(g) || !isfinite(x[j]))
			return 0;
		if (b->lower[j] == b->upper[j])
			native_product(&dual, &correction, b->lower[j], g);
		else if (g >= 0. && isfinite(b->lower[j]))
			native_product(&dual, &correction, b->lower[j], g);
		else if (g < 0. && isfinite(b->upper[j]))
			native_product(&dual, &correction, b->upper[j], g);
	}
	dual += correction;
	if (!isfinite(primal) || !isfinite(dual) ||
		fabsl(primal - dual) > tolerance * (1. + fabsl(primal)))
		return 0;
	saved = (double *)lp_simplex_malloc(
		((size_t)a->rows + 1) * sizeof(double));
	if (!saved)
		return 0;
	lp_simplex_memcpy(saved, y, (size_t)a->rows * sizeof(double));
	/* One complete projection sweep is one unit of shared solver work.
	 * A failed transaction consumes work but restores the original y. */
	for (pass = 0; pass < SIMPLEX_PAN_RESTORE_PASSES && *iterations < limit;
		pass++) {
		(*iterations)++;
		for (j = 0; j < a->columns; j++) {
			double g, budget, target = 0., step;
			long double norm = 0.;

			if (b->lower[j] == b->upper[j])
				continue;
			budget = tolerance * (1. + fabs(a->objective[j]));
			g = a->objective[j] + simplex_pan_column_dot(a, j, y);
			if (x[j] == b->lower[j]) {
				if (g >= -budget)
					continue;
				target = -SIMPLEX_PAN_RESTORE_TARGET_FRACTION *
					 budget;
			} else if (x[j] == b->upper[j]) {
				if (g <= budget)
					continue;
				target = SIMPLEX_PAN_RESTORE_TARGET_FRACTION *
					 budget;
			} else if (fabs(g) <= budget)
				continue;
			for (k = a->column_start[j]; k < a->column_start[j + 1];
				k++)
				norm += (long double)a->value[k] * a->value[k];
			if (norm == 0.)
				continue;
			step = (double)(((long double)target - g) / norm);
			for (k = a->column_start[j]; k < a->column_start[j + 1];
				k++) {
				i = a->row_index[k];
				y[i] = fma(step, a->value[k], y[i]);
			}
		}
		if (simplex_pan_bounded_check_dual(
			    b, x, y, tolerance, violation)) {
			lp_simplex_free(saved);
			return 1;
		}
	}
	lp_simplex_memcpy(y, saved, (size_t)a->rows * sizeof(double));
	lp_simplex_free(saved);
	(void)simplex_pan_bounded_check_dual(b, x, y, tolerance, violation);
	return 0;
}

static int native_certificate(const struct simplex_PanBounded *b,
	const double *x, double *y, const struct lp_simplex_Options *o,
	int *iterations, double *violation, const struct simplex_PanConfig *cfg)
{
	if (cfg->certificate_dual_repair)
		return simplex_pan_bounded_repair_dual(b, x, y,
			o->dual_tolerance, o->iteration_limit, iterations,
			violation);
	return simplex_pan_bounded_check_dual(
		b, x, y, o->dual_tolerance, violation);
}

static int native_free(
	const struct simplex_PanBounded *b, const double *x, int j)
{
	return x[j] > b->lower[j] && x[j] < b->upper[j];
}

/* Recover numerically active bounds using aggregate row activity, not a
 * coordinate cutoff: a tiny coordinate with a large column can be essential.
 * The entire transaction must preserve the original feasibility and cost. */
static int native_clean(const struct simplex_PanBounded *b, double *x,
	double tolerance, double *saved, double *r, double *work)
{
	int i, j, k, changed = 0;
	long double before = simplex_pan_cone_cost(&b->a, x), cost = 0.;
	double budget = tolerance * SIMPLEX_PAN_INNER_TOLERANCE_FRACTION;

	lp_simplex_memcpy(saved, x, (size_t)b->a.columns * sizeof(double));
	for (i = 0; i < b->a.rows; i++)
		work[i] = 0.;
	for (j = 0; j < b->a.columns; j++) {
		double target = x[j], delta;

		if (isfinite(b->lower[j]) &&
			fabs(x[j] - b->lower[j]) <=
				64. * DBL_EPSILON *
					(1. + fabs(x[j]) + fabs(b->lower[j])))
			target = b->lower[j];
		if (isfinite(b->upper[j]) &&
			fabs(x[j] - b->upper[j]) < fabs(target - x[j]))
			target = b->upper[j];
		/* Upper-only variables need the same proximity test. */
		if (target == x[j] && isfinite(b->upper[j]) &&
			fabs(x[j] - b->upper[j]) <=
				64. * DBL_EPSILON *
					(1. + fabs(x[j]) + fabs(b->upper[j])))
			target = b->upper[j];
		delta = target - x[j];
		if (delta == 0.)
			continue;
		if (fabsl(cost + (long double)b->a.objective[j] * delta) >
			32. * DBL_EPSILON * (1. + fabsl(before)))
			continue;
		for (k = b->a.column_start[j]; k < b->a.column_start[j + 1];
			k++)
			if (work[b->a.row_index[k]] +
					fabs(b->a.value[k] * delta) >
				budget)
				break;
		if (k < b->a.column_start[j + 1])
			continue;
		for (k = b->a.column_start[j]; k < b->a.column_start[j + 1];
			k++)
			work[b->a.row_index[k]] += fabs(b->a.value[k] * delta);
		cost += (long double)b->a.objective[j] * delta;
		x[j] = target;
		changed++;
	}
	if (changed && (simplex_pan_bounded_primal(b, x, r, work) > tolerance ||
			       simplex_pan_cone_cost(&b->a, x) >
				       before + 32. * DBL_EPSILON *
							(1. + fabsl(before)))) {
		lp_simplex_memcpy(
			x, saved, (size_t)b->a.columns * sizeof(double));
		return 0;
	}
	return changed;
}

static double native_ratio(const struct simplex_PanBounded *b, const double *x,
	const double *d, int *blocking)
{
	double alpha = HUGE_VAL;
	int j;
	*blocking = -1;
	for (j = 0; j < b->a.columns; j++) {
		double ratio = HUGE_VAL;

		if (!isfinite(d[j]))
			return -1.;
		if (d[j] > 0. && isfinite(b->upper[j]))
			ratio = (b->upper[j] - x[j]) / d[j];
		if (d[j] < 0. && isfinite(b->lower[j]))
			ratio = (b->lower[j] - x[j]) / d[j];
		if (ratio < alpha) {
			alpha = ratio;
			*blocking = j;
		}
	}
	return alpha;
}

/* Cancellation can hide a huge movement even from compensated residuals.
 * Reject a step whose roundoff amplification exhausts the absolute LP budget. */
static int native_step_safe(const struct simplex_PanBounded *b, const double *d,
	double alpha, double tolerance, double *work)
{
	int i, j, k;

	for (i = 0; i < b->a.rows; i++)
		work[i] = 0.;
	for (j = 0; j < b->a.columns; j++)
		for (k = b->a.column_start[j]; k < b->a.column_start[j + 1];
			k++) {
			int row = b->a.row_index[k];

			work[row] += fabs(b->a.value[k]) * fabs(alpha * d[j]);
			if (!isfinite(work[row]) ||
				work[row] * 64. * DBL_EPSILON > tolerance)
				return 0;
		}
	return 1;
}

struct native_Trial {
	const struct simplex_PanBounded *b;
	const struct lp_simplex_Options *options;
	const double *point;
	const signed char *sign;
	double *trial, *r, *work;
	int *iterations;
	int calls, checks, accepted, period;
};

/* An auxiliary least-squares residual need not converge to a fixed absolute
 * tolerance when its actual box-limited LP step is already accurate enough. */
static int native_trial(void *context, const double *direction)
{
	struct native_Trial *t = (struct native_Trial *)context;
	const struct simplex_PanBounded *b = t->b;
	double alpha = HUGE_VAL, base, error;
	long double before, after;
	int j, blocking = -1;

	if ((++t->calls % t->period) != 0 ||
		*t->iterations >= t->options->iteration_limit)
		return 0;
	(*t->iterations)++;
	t->checks++;
	base = simplex_pan_bounded_primal(b, t->point, t->r, t->work);
	if (!isfinite(base) || base >= t->options->primal_tolerance)
		return 0;
	for (j = 0; j < b->a.columns; j++) {
		double d = direction[j] * t->sign[j], ratio = HUGE_VAL;

		if (!isfinite(d))
			return 0;
		if (d > 0. && isfinite(b->upper[j]))
			ratio = (b->upper[j] - t->point[j]) / d;
		if (d < 0. && isfinite(b->lower[j]))
			ratio = (b->lower[j] - t->point[j]) / d;
		if (ratio < alpha) {
			alpha = ratio;
			blocking = j;
		}
	}
	if (blocking < 0 || !isfinite(alpha) || alpha <= 0. ||
		!native_step_safe(b, direction, alpha,
			t->options->primal_tolerance, t->work))
		return 0;
	for (j = 0; j < b->a.columns; j++)
		t->trial[j] =
			fma(alpha, direction[j] * t->sign[j], t->point[j]);
	t->trial[blocking] = direction[blocking] * t->sign[blocking] > 0.
				     ? b->upper[blocking]
				     : b->lower[blocking];
	error = simplex_pan_bounded_primal(b, t->trial, t->r, t->work);
	if (!isfinite(error) ||
		error > base + SIMPLEX_PAN_TRIAL_MARGIN_FRACTION *
					(t->options->primal_tolerance - base))
		return 0;
	before = simplex_pan_cone_cost(&b->a, t->point);
	after = simplex_pan_cone_cost(&b->a, t->trial);
	if (!isfinite(before) || !isfinite(after) ||
		before - after <=
			SIMPLEX_PAN_TRIAL_GAIN_RELATIVE * (1. + fabsl(before)))
		return 0;
	t->accepted++;
	if (getenv(SIMPLEX_ENV_PROFILE))
		fprintf(stderr,
			"pan-bda native trial: gain=%.6g primal=%.6g alpha=%.6g iterations=%d\n",
			(double)(before - after), error, alpha, *t->iterations);
	return 1;
}

static int native_ray(
	const struct simplex_PanBounded *b, double *d, double *r, double *work)
{
	struct simplex_PanStandard homogeneous = b->a;
	double maximum = 0.;
	long double absolute = 0.;
	int i, j, k;

	for (j = 0; j < b->a.columns; j++) {
		if (!isfinite(d[j]) || (d[j] > 0. && isfinite(b->upper[j])) ||
			(d[j] < 0. && isfinite(b->lower[j])))
			return 0;
		maximum = fmax(maximum, fabs(d[j]));
	}
	if (maximum == 0.)
		return 0;
	for (j = 0; j < b->a.columns; j++)
		d[j] /= maximum;
	for (j = 0; j < b->a.columns; j++)
		absolute += fabsl((long double)b->a.objective[j] * d[j]);
	if (simplex_pan_cone_cost(&b->a, d) >=
		-64. * DBL_EPSILON * (1. + absolute))
		return 0;
	for (i = 0; i < b->a.rows; i++)
		r[i] = 0.;
	homogeneous.rhs = r;
	/* rhs aliases residual: use the compensated multiply's separate output. */
	(void)simplex_pan_residual(&homogeneous, d, work, r);
	for (i = 0; i < b->a.rows; i++)
		r[i] = 0.;
	for (j = 0; j < b->a.columns; j++)
		for (k = b->a.column_start[j]; k < b->a.column_start[j + 1];
			k++)
			r[b->a.row_index[k]] += fabs(b->a.value[k] * d[j]);
	for (i = 0; i < b->a.rows; i++)
		if (!isfinite(work[i]) ||
			fabs(work[i]) >
				32. * DBL_EPSILON * (b->a.columns + 2.) * r[i])
			return 0;
	return 1;
}

/* Select an independent representation of the free face. Primal x is kept
 * separately: dropping a dependent column here never deletes primal mass. */
static int native_select(struct simplex_PanBasis *basis,
	const struct simplex_PanBounded *b, const double *x, double *zero,
	int *candidate, unsigned char *remove, double *coefficient,
	const struct simplex_PanConfig *config)
{
	int j, count = 0;

	for (j = 0; j < b->a.columns; j++)
		if (native_free(b, x, j) && b->a.column_norm[j] > 0.)
			candidate[count++] = j;
	if (simplex_pan_basis_create_configured(
		    basis, &b->a, basis->rank_tolerance, config))
		return -1;
#ifdef LP_SIMPLEX_HAVE_SPQR
	if ((size_t)count <=
		(size_t)b->a.rows +
			(size_t)b->a.rows / SIMPLEX_PAN_SUPPORT_ROW_FRACTION +
			SIMPLEX_PAN_SUPPORT_ROW_SLACK)
		return simplex_pan_basis_select_columns(
			basis, candidate, count, zero);
#endif
	for (j = 0; j < count; j++) {
		int i;
		double error;

		if (basis->count) {
			if (simplex_pan_basis_column_projection(
				    basis, candidate[j], coefficient, &error))
				return -1;
			if (error <= basis->rank_tolerance *
					     b->a.column_norm[candidate[j]])
				continue;
		}
		if (basis->count == b->a.rows) {
			if (simplex_pan_basis_repair_rank(basis, zero))
				return -1;
			if (basis->count == b->a.rows)
				continue;
		}
		if (simplex_pan_basis_add(basis, candidate[j]))
			return -1;
		if (simplex_pan_basis_factorize(basis)) {
			for (i = 0; i < basis->count; i++)
				remove[i] = basis->column[i] == candidate[j];
			simplex_pan_basis_remove_marked(basis, remove);
			if (simplex_pan_basis_rebuild(basis))
				return -1;
		}
	}
	return simplex_pan_basis_repair_rank(basis, zero);
}

static int native_compress(struct simplex_PanBounded *b, double *x,
	struct simplex_PanBasis *basis, const struct lp_simplex_Options *o,
	double *d, double *saved, double *coefficient, double *r, double *work,
	unsigned char *remove, int *iterations)
{
	int j, i, moves = 0;
	long double before = simplex_pan_cone_cost(&b->a, x);

	lp_simplex_memcpy(saved, x, (size_t)b->a.columns * sizeof(double));
	for (j = 0; j < b->a.columns; j++)
		if (native_free(b, x, j) && b->a.column_norm[j] > 0. &&
			basis->position[j] < 0) {
			double error, alpha, sign = 1.;
			long double gamma, absolute = 0.;
			int blocking;

			if (*iterations >= o->iteration_limit)
				goto fail;
			(*iterations)++;
			if (simplex_pan_basis_column_projection(
				    basis, j, coefficient, &error))
				goto fail;
			if (error > o->pivot_tolerance * b->a.column_norm[j]) {
				if (simplex_pan_basis_add(basis, j) ||
					simplex_pan_basis_factorize(basis))
					goto fail;
				continue;
			}
			lp_simplex_memset(
				d, 0, (size_t)b->a.columns * sizeof(double));
			d[j] = 1.;
			for (i = 0; i < basis->count; i++)
				d[basis->column[i]] = -coefficient[i];
			gamma = simplex_pan_cone_cost(&b->a, d);
			if (!isfinite(gamma))
				goto fail;
			for (i = 0; i < b->a.columns; i++)
				absolute += fabsl(
					(long double)b->a.objective[i] * d[i]);
			if (fabsl(gamma) <= 64. * DBL_EPSILON * (1. + absolute))
				gamma = 0.;
			if (gamma > 0.)
				sign = -1.;
			for (i = 0; i < b->a.columns; i++)
				d[i] *= sign;
			alpha = native_ratio(b, x, d, &blocking);
			if (!isfinite(alpha) && gamma == 0.) {
				for (i = 0; i < b->a.columns; i++)
					d[i] = -d[i];
				alpha = native_ratio(b, x, d, &blocking);
			}
			if (!isfinite(alpha)) {
				if (gamma != 0. && native_ray(b, d, r, work))
					return 2;
				/* A cost-redundant fully free column needs no primal move. */
				if (gamma == 0.)
					continue;
				goto fail;
			}
			if (alpha <= 0. || moves++ >= b->a.columns)
				goto fail;
			if (!native_step_safe(
				    b, d, alpha, o->primal_tolerance, work))
				goto fail;
			for (i = 0; i < b->a.columns; i++)
				x[i] = fma(alpha, d[i], x[i]);
			x[blocking] = d[blocking] > 0. ? b->upper[blocking]
						       : b->lower[blocking];
			if (simplex_pan_bounded_primal(b, x, r, work) >
					o->primal_tolerance ||
				simplex_pan_cone_cost(&b->a, x) >
					before + 32. * DBL_EPSILON *
							 (1. + fabsl(before)))
				goto fail;
			for (i = 0; i < basis->count; i++)
				remove[i] =
					!native_free(b, x, basis->column[i]);
			simplex_pan_basis_remove_marked(basis, remove);
			if (native_free(b, x, j) &&
				simplex_pan_basis_add(basis, j))
				goto fail;
			if (simplex_pan_basis_factorize(basis) &&
				simplex_pan_basis_rebuild(basis))
				goto fail;
		}
	return moves > 0 ? 1 : 0;
fail:
	lp_simplex_memcpy(x, saved, (size_t)b->a.columns * sizeof(double));
	return -1;
}

/* A free-face projection can miss residual components normal to that face.
 * Allow inward bound directions in a mixed NNLS correction, then enforce the
 * opposite bounds by a finite ratio step. This is a feasibility subproblem,
 * never an infeasibility or optimality certificate for the original LP. */
static int native_restore_cone(const struct simplex_PanBounded *b,
	double *point, const struct lp_simplex_Options *o, int *iterations,
	long double ceiling, const struct simplex_PanConfig *cfg,
	int *service_stop)
{
	struct simplex_PanStandard t = b->a;
	struct simplex_PanBasis basis;
	double *storage = NULL, *saved, *zero, *d, *norm, *right, *r, *work,
	       *coefficient, *matrix = NULL;
	unsigned char *free_mask = NULL, *remove = NULL;
	int *candidate = NULL, i, j, k, count, attempt,
	    status = lp_simplex_CondUnsatisfied, state, result = -1, blocking;
	int limit = o->iteration_limit;
	double alpha, before, error;

	*service_stop = SIMPLEX_STOP_NUMERICAL;
	lp_simplex_memset(&basis, 0, sizeof(basis));
	if (limit - *iterations > SIMPLEX_PAN_RESTORE_WORK_CAP)
		limit = *iterations + SIMPLEX_PAN_RESTORE_WORK_CAP;
	storage = (double *)lp_simplex_malloc(
		(4 * ((size_t)b->a.columns + 1) + 4 * ((size_t)b->a.rows + 1)) *
		sizeof(double));
	matrix = (double *)lp_simplex_malloc(
		2 * ((size_t)b->a.nonzeros + 1) * sizeof(double));
	free_mask =
		(unsigned char *)lp_simplex_malloc((size_t)b->a.columns + 1);
	remove = (unsigned char *)lp_simplex_malloc((size_t)b->a.rows + 2);
	candidate = (int *)lp_simplex_malloc(
		((size_t)b->a.columns + 1) * sizeof(int));
	if (!storage || !matrix || !free_mask || !remove || !candidate) {
		*service_stop = SIMPLEX_STOP_RESOURCE;
		goto done;
	}
	saved = storage;
	zero = saved + b->a.columns + 1;
	d = zero + b->a.columns + 1;
	norm = d + b->a.columns + 1;
	right = norm + b->a.columns + 1;
	r = right + b->a.rows + 1;
	work = r + b->a.rows + 1;
	coefficient = work + b->a.rows + 1;
	lp_simplex_memcpy(saved, point, (size_t)b->a.columns * sizeof(double));
	lp_simplex_memset(zero, 0, ((size_t)b->a.columns + 1) * sizeof(double));
	t.rhs = right;
	t.objective = zero;
	t.column_norm = norm;
	t.value = matrix;
	t.normalized_value = matrix + b->a.nonzeros + 1;
	for (attempt = 0;
		attempt < SIMPLEX_PAN_RESTORE_ATTEMPTS && *iterations < limit;
		attempt++) {
		before = simplex_pan_bounded_primal(b, point, right, work);
		if (before <= o->primal_tolerance *
				      SIMPLEX_PAN_RESTORE_TARGET_FRACTION) {
			result = 1;
			break;
		}
		simplex_pan_basis_destroy(&basis);
		basis.rank_tolerance = o->pivot_tolerance;
		if (native_select(&basis, b, point, zero, candidate, remove,
			    coefficient, cfg))
			break;
		count = basis.count;
		for (i = 0; i < count; i++)
			candidate[i] = basis.column[i];
		simplex_pan_basis_destroy(&basis);
		for (j = 0; j < b->a.columns; j++) {
			double sign =
				b->lower[j] == b->upper[j]
					? 0.
					: (point[j] == b->upper[j] ? -1. : 1.);

			free_mask[j] = (unsigned char)native_free(b, point, j);
			norm[j] = fabs(sign) * b->a.column_norm[j];
			for (k = b->a.column_start[j];
				k < b->a.column_start[j + 1]; k++) {
				t.value[k] = sign * b->a.value[k];
				t.normalized_value[k] =
					sign * b->a.normalized_value[k];
			}
		}
		if (simplex_pan_basis_create_configured(
			    &basis, &t, o->pivot_tolerance, cfg))
			break;
		for (i = 0; i < count; i++)
			if (simplex_pan_basis_add(&basis, candidate[i]))
				break;
		if (i < count || simplex_pan_basis_rebuild(&basis))
			break;
		state = simplex_pan_phase1_mixed(&basis, d, free_mask, r,
			o->primal_tolerance *
				SIMPLEX_PAN_RESTORE_INNER_FRACTION,
			o->primal_tolerance *
				SIMPLEX_PAN_RESTORE_INNER_FRACTION,
			limit, iterations, &status);
		if (state != 0)
			break;
		for (j = 0; j < b->a.columns; j++)
			if (point[j] == b->upper[j])
				d[j] = -d[j];
		alpha = native_ratio(b, point, d, &blocking);
		if (alpha <= 0.)
			break;
		if (alpha > 1.)
			alpha = 1.;
		if (!native_step_safe(b, d, alpha, o->primal_tolerance, work))
			break;
		for (j = 0; j < b->a.columns; j++)
			point[j] = fma(alpha, d[j], point[j]);
		if (alpha < 1. && blocking >= 0)
			point[blocking] = d[blocking] > 0. ? b->upper[blocking]
							   : b->lower[blocking];
		error = simplex_pan_bounded_primal(b, point, r, work);
		if (!isfinite(error) || error >= before ||
			simplex_pan_cone_cost(&b->a, point) > ceiling)
			break;
		if (error <= o->primal_tolerance *
				     SIMPLEX_PAN_RESTORE_TARGET_FRACTION) {
			result = 1;
			break;
		}
	}
	if (result != 1 || simplex_pan_cone_cost(&b->a, point) > ceiling) {
		lp_simplex_memcpy(
			point, saved, (size_t)b->a.columns * sizeof(double));
		result = -1;
	}
	if (getenv(SIMPLEX_ENV_PROFILE))
		fprintf(stderr,
			"pan-bda native cone restore: result=%d attempts=%d rank=%d status=%d iterations=%d\n",
			result, attempt, basis.count, status, *iterations);
done:
	if (result != 1 && *iterations >= limit)
		*service_stop = SIMPLEX_STOP_BUDGET;
	simplex_pan_basis_destroy(&basis);
	lp_simplex_free(storage);
	lp_simplex_free(matrix);
	lp_simplex_free(free_mask);
	lp_simplex_free(remove);
	lp_simplex_free(candidate);
	return result;
}

/* Local projection onto the original equalities, using only independently
 * selected interior columns. Bounded steps may release further face columns;
 * failure always restores the complete input. No optimality is inferred. */
int simplex_pan_bounded_restore_configured_reported(
	const struct simplex_PanBounded *b, double *point,
	const struct lp_simplex_Options *o, int *iterations,
	long double ceiling, const struct simplex_PanConfig *cfg,
	int *service_stop)
{
	int unobserved_stop;
	struct simplex_PanBasis basis;
	double *storage = NULL, *saved, *zero, *d, *r, *work, *coefficient;
	int *candidate = NULL, i, j, blocking, attempt, result = -1;
	unsigned char *remove = NULL;
	double before, alpha, error;

	if (!service_stop)
		service_stop = &unobserved_stop;
	*service_stop = SIMPLEX_STOP_NUMERICAL;
	lp_simplex_memset(&basis, 0, sizeof(basis));
	storage = (double *)lp_simplex_malloc(
		(3 * ((size_t)b->a.columns + 1) + 3 * ((size_t)b->a.rows + 1)) *
		sizeof(double));
	candidate = (int *)lp_simplex_malloc(
		((size_t)b->a.columns + 1) * sizeof(int));
	remove = (unsigned char *)lp_simplex_malloc((size_t)b->a.rows + 2);
	if (!storage || !candidate || !remove) {
		*service_stop = SIMPLEX_STOP_RESOURCE;
		goto done;
	}
	saved = storage;
	zero = saved + b->a.columns + 1;
	d = zero + b->a.columns + 1;
	r = d + b->a.columns + 1;
	work = r + b->a.rows + 1;
	coefficient = work + b->a.rows + 1;
	lp_simplex_memcpy(saved, point, (size_t)b->a.columns * sizeof(double));
	lp_simplex_memset(zero, 0, ((size_t)b->a.columns + 1) * sizeof(double));
	for (attempt = 0; attempt < SIMPLEX_PAN_RESTORE_ATTEMPTS &&
			  *iterations < o->iteration_limit;
		attempt++) {
		before = simplex_pan_bounded_primal(b, point, r, work);
		if (before <= o->primal_tolerance *
				      SIMPLEX_PAN_RESTORE_TARGET_FRACTION) {
			result = 1;
			break;
		}
		simplex_pan_basis_destroy(&basis);
		basis.rank_tolerance = o->pivot_tolerance;
		if (native_select(&basis, b, point, zero, candidate, remove,
			    coefficient, cfg))
			break;
		(*iterations)++;
		if (simplex_pan_basis_least_squares(
			    &basis, r, coefficient, &error))
			break;
		lp_simplex_memset(d, 0, (size_t)b->a.columns * sizeof(double));
		for (i = 0; i < basis.count; i++)
			d[basis.column[i]] = coefficient[i];
		alpha = native_ratio(b, point, d, &blocking);
		if (alpha <= 0.)
			break;
		if (alpha > 1.)
			alpha = 1.;
		if (!native_step_safe(b, d, alpha, o->primal_tolerance, work))
			break;
		for (j = 0; j < b->a.columns; j++)
			point[j] = fma(alpha, d[j], point[j]);
		if (alpha < 1. && blocking >= 0)
			point[blocking] = d[blocking] > 0. ? b->upper[blocking]
							   : b->lower[blocking];
		error = simplex_pan_bounded_primal(b, point, r, work);
		if (!isfinite(error) || error >= before ||
			simplex_pan_cone_cost(&b->a, point) > ceiling)
			break;
		if (error <= o->primal_tolerance *
				     SIMPLEX_PAN_RESTORE_TARGET_FRACTION) {
			result = 1;
			break;
		}
	}
	if (result != 1 || simplex_pan_cone_cost(&b->a, point) > ceiling) {
		lp_simplex_memcpy(
			point, saved, (size_t)b->a.columns * sizeof(double));
		result = -1;
	}
done:
	simplex_pan_basis_destroy(&basis);
	lp_simplex_free(storage);
	lp_simplex_free(candidate);
	lp_simplex_free(remove);
	if (result != 1 && *iterations < o->iteration_limit)
		result = native_restore_cone(
			b, point, o, iterations, ceiling, cfg, service_stop);
	*service_stop = result == 1 ? SIMPLEX_STOP_HANDOFF
				    : (*iterations >= o->iteration_limit
						      ? SIMPLEX_STOP_BUDGET
						      : *service_stop);
	return result;
}

/* The orientation is rebuilt at every accepted bound hit. Free coordinates
 * are signed only for numerical convenience; both signs remain permitted. */
static void native_tangent(const struct simplex_PanBounded *b, const double *x,
	struct simplex_PanStandard *t, signed char *sign,
	unsigned char *free_mask)
{
	int j, k, pos = 0;

	for (j = 0; j < b->a.columns; j++) {
		double norm = 0.;

		free_mask[j] = (unsigned char)native_free(b, x, j);
		sign[j] = b->lower[j] == b->upper[j]
				  ? 0
				  : (x[j] == b->upper[j] ? -1 : 1);
		t->column_start[j] = pos;
		t->objective[j] = sign[j] * b->a.objective[j];
		if (sign[j])
			for (k = b->a.column_start[j];
				k < b->a.column_start[j + 1]; k++) {
				t->row_index[pos] = b->a.row_index[k];
				t->value[pos] = sign[j] * b->a.value[k];
				norm = hypot(norm, t->value[pos++]);
			}
		t->column_norm[j] = norm;
		for (k = t->column_start[j]; k < pos; k++)
			t->normalized_value[k] =
				norm > 0. ? t->value[k] / norm : 0.;
	}
	t->column_start[t->columns] = pos;
	t->nonzeros = pos;
	for (j = 0; j < t->rows; j++)
		t->rhs[j] = 0.;
}

int simplex_pan_bounded_solve_certified_configured(
	const struct simplex_Problem *p, const struct lp_simplex_Options *o,
	const double *start, double *out, double *row_dual,
	struct lp_simplex_Result *result, const struct simplex_PanConfig *cfg)
{
	struct simplex_PanServiceContext services;
	struct simplex_PanBounded b, tangent;
	struct simplex_PanStandard cone;
	struct simplex_PanBasis basis, free_basis;
	struct simplex_PanFaceRecovery face = {0};
	double *storage = NULL, *point, *saved, *d, *zero, *y, *r, *work,
	       *row_scale, *coefficient;
	unsigned char *free_mask = NULL, *remove = NULL;
	signed char *sign = NULL;
	int *candidate = NULL;
	int i, j, state = -1, status = lp_simplex_MemoryAllocError,
		  iterations = 0, steps = 0, have_point = 0, low_progress = 0;
	const char *feasible_source = start ? "supplied" : "nnls";
	int implicit;
	double dual = HUGE_VAL, primal = HUGE_VAL;
	clock_t begun = clock();
	struct native_Trial trial;
	int adaptive;

	if (!simplex_pan_config_valid(cfg) ||
		cfg->model_form != SIMPLEX_PAN_NATIVE_BOUNDS ||
		o->algorithm != lp_simplex_ALGORITHM_PAN_BDA) {
		result->status = lp_simplex_CondUnsatisfied;
		result->iterations = 0;
		result->objective = result->dual_infeasibility = NAN;
		result->primal_infeasibility = HUGE_VAL;
		result->dual_residual_space = lp_simplex_RESIDUAL_UNAVAILABLE;
		return lp_simplex_EXIT_FAILURE;
	}
	simplex_pan_service_context_init(&services, p);
	implicit = cfg->oracle == SIMPLEX_PAN_IMPLICIT_ORACLE;
	adaptive = cfg->adaptive;
	lp_simplex_memset(&trial, 0, sizeof(trial));
	result->objective = HUGE_VAL;
	result->dual_residual_space = lp_simplex_RESIDUAL_UNAVAILABLE;
	lp_simplex_memset(&b, 0, sizeof(b));
	lp_simplex_memset(&tangent, 0, sizeof(tangent));
	lp_simplex_memset(&cone, 0, sizeof(cone));
	lp_simplex_memset(&basis, 0, sizeof(basis));
	lp_simplex_memset(&free_basis, 0, sizeof(free_basis));
	if (simplex_pan_bounded_create(&b, p) ||
		simplex_pan_bounded_create(&tangent, p))
		goto done;
	simplex_pan_service_bind(&b.a, &services);
	simplex_pan_service_bind(&tangent.a, &services);
	storage = (double *)lp_simplex_malloc(
		(4 * ((size_t)b.a.columns + 1) + 7 * ((size_t)b.a.rows + 1)) *
		sizeof(double));
	free_mask = (unsigned char *)lp_simplex_malloc((size_t)b.a.columns + 1);
	sign = (signed char *)lp_simplex_malloc((size_t)b.a.columns + 1);
	remove = (unsigned char *)lp_simplex_malloc((size_t)b.a.rows + 2);
	candidate = (int *)lp_simplex_malloc(
		((size_t)b.a.columns + 1) * sizeof(int));
	if (!storage || !free_mask || !sign || !remove || !candidate)
		goto done;
	point = storage;
	saved = point + b.a.columns + 1;
	d = saved + b.a.columns + 1;
	zero = d + b.a.columns + 1;
	y = zero + b.a.columns + 1;
	r = y + b.a.rows + 1;
	work = r + b.a.rows + 1;
	row_scale = work + b.a.rows + 1;
	coefficient = row_scale + b.a.rows + 1;
	trial.b = &b;
	trial.options = o;
	trial.point = point;
	trial.sign = sign;
	trial.period = cfg->restore ? 1 : SIMPLEX_PAN_TRIAL_SPARSE_PERIOD;
	trial.trial = saved;
	trial.r = coefficient + b.a.rows + 1;
	trial.work = trial.r + b.a.rows + 1;
	trial.iterations = &iterations;
	lp_simplex_memset(storage, 0,
		(4 * ((size_t)b.a.columns + 1) + 7 * ((size_t)b.a.rows + 1)) *
			sizeof(double));
	if (getenv(SIMPLEX_ENV_PROFILE))
		fprintf(stderr,
			"pan-bda native model: rows=%d columns=%d structural=%d nonzeros=%d\n",
			b.a.rows, b.a.columns, p->columns, b.a.nonzeros);
	if (start) {
		if (simplex_pan_bounded_encode(
			    &b, start, point, r, work, o->primal_tolerance)) {
			status = lp_simplex_PrecisionError;
			goto done;
		}
	} else if (b.a.rows == 0) {
		feasible_source = "box";
		for (j = 0; j < b.a.columns; j++)
			point[j] = fmin(b.upper[j], fmax(b.lower[j], 0.));
	} else {
		int relaxed = 1;

		for (j = 0; j < b.a.columns; j++)
			if (b.lower[j] != 0.)
				relaxed = 0;
		if (relaxed && !cfg->dual_start) {
			int phase_limit = o->iteration_limit;
			/* A short sparse NNLS crash is effective on sparse feasible
			 * vertices. Do not spend the entire solve budget on Phase I when
			 * its support becomes dense; retain the own-dual fallback. */
			if (phase_limit >= 4)
				phase_limit /= 4;
			if (phase_limit > SIMPLEX_PAN_FEASIBILITY_PROBE_WORK)
				phase_limit =
					SIMPLEX_PAN_FEASIBILITY_PROBE_WORK;
			if (simplex_pan_basis_create_configured(
				    &basis, &b.a, o->pivot_tolerance, cfg))
				goto done;
			state = simplex_pan_service_nnls(&basis, point,
				o->primal_tolerance, o->dual_tolerance,
				phase_limit, &iterations, &status, NULL);
			simplex_pan_basis_destroy(&basis);
			if (!simplex_pan_service_context_valid(&services))
				goto precision;
			if (state && status == lp_simplex_Infeasibility)
				goto done;
			relaxed = state == 0 &&
				  simplex_pan_bounded_primal(&b, point, r,
					  work) <= o->primal_tolerance;
			if (relaxed) {
				lp_simplex_memcpy(out, point,
					(size_t)p->columns * sizeof(double));
				relaxed = simplex_pan_bounded_encode(&b, out,
						  point, r, work,
						  o->primal_tolerance) == 0;
			}
		} else
			relaxed = 0;
		if (!relaxed) {
			struct lp_simplex_Options remaining = *o;
			struct lp_simplex_Result phase;

			remaining.iteration_limit -= iterations;
			feasible_source = "dual";
			if (remaining.iteration_limit <= 0) {
				status = lp_simplex_ExceedIterLimit;
				goto done;
			}
			state = simplex_pan_service_dual_feasible(&services, p,
				o, out, &iterations, &phase, NULL);
			if (state) {
				status = phase.status;
				goto done;
			}
			if (simplex_pan_bounded_encode(&b, out, point, r, work,
				    o->primal_tolerance)) {
				status = lp_simplex_PrecisionError;
				goto done;
			}
		}
	}
	have_point = 1;
	simplex_progress_emit(
		p, point, y, 1., "native-feasible", iterations, 1);
	status = lp_simplex_CondUnsatisfied;
	primal = simplex_pan_bounded_primal(&b, point, r, work);
	if (getenv(SIMPLEX_ENV_PROFILE))
		fprintf(stderr,
			"pan-bda first feasible: source=native-%s cpu=%.6f iterations=%d primal=%.6g\n",
			feasible_source,
			(double)(clock() - begun) / CLOCKS_PER_SEC, iterations,
			primal);
	if (cfg->pdhg_iterations > 0) {
		int budget = cfg->pdhg_iterations;
		/* Leave the original cone oracle a recovery/certification budget,
		 * including on unbounded models that PDHG does not certify. */
		if (budget > (o->iteration_limit - iterations) / 2)
			budget = (o->iteration_limit - iterations) / 2;
		(void)simplex_pan_service_pdhg(
			&b, point, y, o, &iterations, budget, cfg, NULL);
	}
	if (!simplex_pan_service_context_valid(&services))
		goto precision;
	for (;;) {
		double scale = 0., alpha, error,
		       inner = fmax(16. * DBL_EPSILON,
			       fmin(o->primal_tolerance, o->dual_tolerance) *
				       SIMPLEX_PAN_INNER_TOLERANCE_FRACTION);
		int blocking = -1, compressed;
		struct simplex_PanServiceOutcome multiplier_search;
		long double before;

		simplex_progress_emit(p, point, y, 1., "native", iterations, 0);
		if (cfg->restore &&
			(steps == 0 ||
				low_progress >=
					SIMPLEX_PAN_LOW_PROGRESS_STEPS)) {
			(void)simplex_pan_service_face(&face, &b, point, y, o,
				&iterations, saved, r, work,
				steps == 0 ? SIMPLEX_PAN_FACE_INITIAL
					   : SIMPLEX_PAN_FACE_LOW_PROGRESS,
				steps, low_progress, NULL, NULL);
			low_progress = 0;
		}
		(void)native_clean(
			&b, point, o->primal_tolerance, saved, r, work);
		if (!simplex_pan_service_context_valid(&services))
			goto precision;
		if (native_certificate(
			    &b, point, y, o, &iterations, &dual, cfg)) {
			status = lp_simplex_Success;
			state = 0;
			goto done;
		}
		if (iterations >= o->iteration_limit) {
			status = lp_simplex_ExceedIterLimit;
			goto done;
		}
		/* Coordinate-only improvements do not need a factorization. */
		for (j = 0; j < b.a.columns; j++)
			if (b.a.column_norm[j] == 0. &&
				b.a.objective[j] != 0.) {
				double target = b.a.objective[j] > 0.
							? b.lower[j]
							: b.upper[j];

				if (!isfinite(target)) {
					status = lp_simplex_Unboundedness;
					goto done;
				}
				point[j] = target;
			}
		if (!simplex_pan_service_context_valid(&services))
			goto precision;
		if (native_certificate(
			    &b, point, y, o, &iterations, &dual, cfg)) {
			status = lp_simplex_Success;
			state = 0;
			goto done;
		}
		free_basis.rank_tolerance = o->pivot_tolerance;
		if (native_select(&free_basis, &b, point, zero, candidate,
			    remove, coefficient, cfg)) {
			if (getenv(SIMPLEX_ENV_PROFILE))
				fprintf(stderr,
					"pan-bda native failure: stage=select rank=%d\n",
					free_basis.count);
			if (p->report)
				p->report->failure =
					SIMPLEX_FAILURE_FACTORIZATION;
			status = lp_simplex_PrecisionError;
			goto done;
		}
		compressed = native_compress(&b, point, &free_basis, o, d,
			saved, coefficient, r, work, remove, &iterations);
		if (compressed == 2) {
			status = lp_simplex_Unboundedness;
			goto done;
		}
		if (compressed < 0) {
			if (iterations >= o->iteration_limit) {
				status = lp_simplex_ExceedIterLimit;
				goto done;
			}
			/* Compression is optional. Its transaction restored the point;
			 * rebuild the seed before trying the fully audited cone oracle. */
			if (getenv(SIMPLEX_ENV_PROFILE))
				fprintf(stderr,
					"pan-bda native recovery: compression rolled back rank=%d\n",
					free_basis.count);
			simplex_pan_basis_destroy(&free_basis);
			free_basis.rank_tolerance = o->pivot_tolerance;
			if (native_select(&free_basis, &b, point, zero,
				    candidate, remove, coefficient, cfg)) {
				status = lp_simplex_PrecisionError;
				goto done;
			}
		}
		native_tangent(&b, point, &tangent.a, sign, free_mask);
		/* A synthetic feasible free-face point lets the existing dual-only
		 * gradient enforce equality on free costs. Native KKT is rechecked. */
		for (j = 0; j < b.a.columns; j++)
			saved[j] = free_mask[j] ? 1. : 0.;
		(void)simplex_pan_residual(&tangent.a, saved, r, work);
		for (i = 0; i < b.a.rows; i++)
			tangent.a.rhs[i] = -r[i];
		simplex_pan_service_derive(&tangent.a, &b.a);
		if (adaptive && cfg->predict) {
			/* Prediction audits need distinct point storage: saved is the
			 * synthetic face point passed to the gradient routine. */
			trial.trial = zero;
			trial.period = 1;
			state = simplex_pan_service_predict(&tangent.a, saved,
				y, o, &iterations, &dual, d, free_mask,
				SIMPLEX_PAN_PREDICT_WORK_CAP, native_trial,
				&trial, cfg, &multiplier_search);
			trial.trial = saved;
			trial.period =
				cfg->restore ? 1
					     : SIMPLEX_PAN_TRIAL_SPARSE_PERIOD;
			lp_simplex_memset(zero, 0,
				((size_t)b.a.columns + 1) * sizeof(double));
		} else {
			state = 0;
			(void)simplex_pan_service_first_order(&tangent.a, saved,
				y, o, &iterations, &dual, NULL, free_mask, 1,
				cfg, &multiplier_search);
		}
		if (!simplex_pan_service_context_valid(&services))
			goto precision;
		if (native_certificate(
			    &b, point, y, o, &iterations, &dual, cfg)) {
			status = lp_simplex_Success;
			state = 0;
			goto done;
		}
		if (cfg->restore &&
			simplex_pan_service_face(&face, &b, point, y, o,
				&iterations, saved, r, work,
				SIMPLEX_PAN_FACE_AFTER_MULTIPLIER_SEARCH, steps,
				low_progress, &multiplier_search, NULL) > 0) {
			simplex_pan_basis_destroy(&free_basis);
			continue;
		}
		if (iterations >= o->iteration_limit) {
			status = lp_simplex_ExceedIterLimit;
			goto done;
		}
		for (i = 0; i < b.a.rows; i++)
			tangent.a.rhs[i] = 0.;
		if (state == 2)
			goto take_direction;
		for (j = 0; j < b.a.columns; j++) {
			tangent.a.objective[j] =
				sign[j] *
				(b.a.objective[j] +
					simplex_pan_column_dot(&b.a, j, y));
			scale = fmax(scale, fabs(tangent.a.objective[j]));
		}
		if (scale == 0.)
			scale = 1.;
		simplex_pan_service_derive(&tangent.a, &b.a);
		if (simplex_pan_cone_create(
			    &cone, &tangent.a, scale, row_scale, implicit) ||
			simplex_pan_basis_create_configured(
				&basis, &cone, o->pivot_tolerance, cfg))
			goto done;
		for (i = 0; i < free_basis.count; i++)
			if (simplex_pan_basis_add(&basis, free_basis.column[i]))
				goto precision;
		if (simplex_pan_basis_rebuild(&basis))
			goto precision;
		state = simplex_pan_service_cone(&basis, implicit, d, free_mask,
			r, inner, o->iteration_limit, &iterations, &status,
			adaptive ? native_trial : NULL, &trial, NULL);
		if (state < 0) {
			if (p->report)
				p->report->failure = basis.oracle_failure;
			goto precision;
		}
		if (state == 1) {
			double denominator = r[b.a.rows];

			if (!isfinite(denominator) || denominator >= 0.)
				goto precision;
			for (i = 0; i < b.a.rows; i++)
				y[i] += scale * (r[i] / denominator) *
					row_scale[i];
			if (native_certificate(
				    &b, point, y, o, &iterations, &dual, cfg)) {
				status = lp_simplex_Success;
				state = 0;
				goto done;
			}
			goto precision;
		}
take_direction:
		for (j = 0; j < b.a.columns; j++)
			d[j] *= sign[j];
		if (!isfinite(simplex_pan_cone_cost(&b.a, d)) ||
			simplex_pan_cone_cost(&b.a, d) >= 0.)
			goto precision;
		alpha = native_ratio(&b, point, d, &blocking);
		if (!isfinite(alpha)) {
			if (native_ray(&b, d, r, work)) {
				status = lp_simplex_Unboundedness;
				goto done;
			}
			goto precision;
		}
		if (alpha <= 0.)
			goto precision;
		if (!native_step_safe(&b, d, alpha, o->primal_tolerance, work))
			goto precision;
		if (iterations >= o->iteration_limit) {
			status = lp_simplex_ExceedIterLimit;
			goto done;
		}
		before = simplex_pan_cone_cost(&b.a, point);
		lp_simplex_memcpy(
			saved, point, (size_t)b.a.columns * sizeof(double));
		for (j = 0; j < b.a.columns; j++)
			point[j] = fma(alpha, d[j], point[j]);
		point[blocking] = d[blocking] > 0. ? b.upper[blocking]
						   : b.lower[blocking];
		error = simplex_pan_bounded_primal(&b, point, r, work);
		if (error > o->primal_tolerance ||
			!isfinite(simplex_pan_cone_cost(&b.a, point)) ||
			simplex_pan_cone_cost(&b.a, point) >
				before + 32. * DBL_EPSILON *
						 (1. + fabsl(before))) {
			lp_simplex_memcpy(point, saved,
				(size_t)b.a.columns * sizeof(double));
			goto precision;
		}
		iterations++;
		steps++;
		primal = error;
		if (before - simplex_pan_cone_cost(&b.a, point) <=
			simplex_pan_progress_threshold(&b.a, point, before,
				o->dual_tolerance,
				SIMPLEX_PAN_SMALL_GAIN_UNITS))
			low_progress++;
		else
			low_progress = 0;
		if (getenv(SIMPLEX_ENV_PROFILE))
			fprintf(stderr,
				"pan-bda native step: steps=%d objective=%.12g primal=%.6g alpha=%.6g rank=%d cpu=%.6f\n",
				steps,
				(double)simplex_pan_cone_cost(&b.a, point),
				primal, alpha, basis.count,
				(double)(clock() - begun) / CLOCKS_PER_SEC);
		simplex_pan_basis_destroy(&basis);
		simplex_pan_basis_destroy(&free_basis);
		simplex_pan_standard_destroy(&cone);
		continue;
precision:
		if (iterations >= o->iteration_limit) {
			status = lp_simplex_ExceedIterLimit;
			goto done;
		}
		/* AUTO owns cross-representation retries and their budgets. */
		if (implicit && cfg->allow_oracle_recovery) {
			implicit = 0;
			simplex_pan_basis_destroy(&basis);
			simplex_pan_basis_destroy(&free_basis);
			simplex_pan_standard_destroy(&cone);
			continue;
		}
		status = lp_simplex_PrecisionError;
		goto done;
	}
done:
	if (!simplex_pan_service_context_valid(&services))
		status = lp_simplex_PrecisionError;
	if (adaptive && getenv(SIMPLEX_ENV_PROFILE))
		fprintf(stderr,
			"pan-bda native trial summary: calls=%d checks=%d accepted=%d\n",
			trial.calls, trial.checks, trial.accepted);
	if (have_point) {
		simplex_progress_emit(
			p, point, y, 1., "native-finish", iterations, 1);
		lp_simplex_memcpy(
			out, point, (size_t)p->columns * sizeof(double));
		primal = simplex_pan_bounded_primal(&b, point, r, work);
		result->objective = (double)simplex_pan_cone_cost(&b.a, point);
		if (status == lp_simplex_Success &&
			(primal > o->primal_tolerance ||
				!simplex_pan_bounded_check_dual(&b, point, y,
					o->dual_tolerance, &dual)))
			status = lp_simplex_PrecisionError;
	}
	if (p->report) {
		struct simplex_StageReport *report = p->report;
		report->stop = simplex_stage_stop_from_status(status);
		report->point_available = have_point;
		simplex_stage_observe(p, o, have_point ? out : NULL,
			have_point ? y : NULL, 1, iterations, report);
		report->work.cpu_seconds =
			(double)(clock() - begun) / CLOCKS_PER_SEC;
	}
	result->status = status;
	result->iterations = iterations;
	result->primal_infeasibility = primal;
	result->dual_infeasibility =
		have_point ? simplex_problem_dual_infeasibility(
				     p, out, y, o->primal_tolerance, 1)
			   : NAN;
	result->dual_residual_space = have_point
					      ? lp_simplex_RESIDUAL_ORIGINAL
					      : lp_simplex_RESIDUAL_UNAVAILABLE;
	if (status == lp_simplex_Success && row_dual)
		lp_simplex_memcpy(
			row_dual, y, (size_t)p->rows * sizeof(double));
	if (getenv(SIMPLEX_ENV_PROFILE))
		fprintf(stderr,
			"pan-bda native summary: steps=%d iterations=%d status=%d primal=%.6g dual=%.6g cpu=%.6f\n",
			steps, iterations, status, primal, dual,
			(double)(clock() - begun) / CLOCKS_PER_SEC);
	simplex_pan_basis_destroy(&basis);
	simplex_pan_basis_destroy(&free_basis);
	simplex_pan_standard_destroy(&cone);
	simplex_pan_face_recovery_destroy(&face);
	simplex_pan_bounded_destroy(&b);
	simplex_pan_bounded_destroy(&tangent);
	lp_simplex_free(storage);
	lp_simplex_free(free_mask);
	lp_simplex_free(sign);
	lp_simplex_free(remove);
	lp_simplex_free(candidate);
	return status == lp_simplex_Success ? 0 : -1;
}
