/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_pan_pdhg.c
 * @brief Optional matrix-free primal-dual recovery for the native box model.
 *
 * Iterates may be infeasible. Only audited feasible improvements, or a full
 * original-model KKT certificate, may leave this routine.
 */

#include "policy/pan.h"
#include "policy/controls.h"
#include "simplex_pan_config.h"
#include "simplex_pan_service.h"
#include "simplex_pan_pdhg.h"
#include "simplex_progress.h"
#include "simplex_stage.h"
#include "simplex_pan_bounded.h"
#include "simplex_pan_cone.h"
#include "simplex_pan_recovery.h"
#include "utils.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void simplex_pan_pdhg_progress_reset(struct simplex_PdhgProgress *h)
{
	lp_simplex_memset(h, 0, sizeof(*h));
}

int simplex_pan_pdhg_repair_progress(
	double primal, double previous_primal, double primal_tolerance)
{
	return primal >
		       primal_tolerance * SIMPLEX_PAN_RESTORE_TARGET_FRACTION &&
	       primal < previous_primal * SIMPLEX_PAN_REPAIR_REDUCTION;
}

int simplex_pan_pdhg_progress_stalled(
	struct simplex_PdhgProgress *h, int work, double merit)
{
	if (!isfinite(merit) || merit < 0. || work < 0)
		return 0;
	if (!h->initialized || work < h->start) {
		h->initialized = 1;
		h->start = h->progress_at = work;
		h->best = merit;
		return 0;
	}
	/* The historical best survives rebounds; returning from a worse peak
	 * cannot postpone the next stagnation observation. */
	if (merit < h->best * (1. - SIMPLEX_PDHG_MERIT_REDUCTION)) {
		h->best = merit;
		h->progress_at = work;
	}
	return work - h->start >= SIMPLEX_PDHG_OBSERVE_WORK &&
	       work - h->progress_at >= SIMPLEX_PDHG_STAGNANT_WORK;
}

static double pdhg_merit(const struct simplex_PanBounded *b, const double *x,
	const double *y, double scale, const double *tau, const double *sigma,
	double *r, double *work)
{
	double merit = 0.;
	int i, j;

	(void)simplex_pan_residual(&b->a, x, r, work);
	for (i = 0; i < b->a.rows; i++)
		merit = hypot(merit, sqrt(sigma[i]) * r[i]);
	for (j = 0; j < b->a.columns; j++) {
		double g = b->a.objective[j] / scale +
			   simplex_pan_column_dot(&b->a, j, y);
		double projected =
			fmin(b->upper[j], fmax(b->lower[j], x[j] - tau[j] * g));

		merit = hypot(merit, (x[j] - projected) / sqrt(tau[j]));
	}
	return merit;
}

/* A handoff diagnostic only: reduced-cost signs still require the independent
 * native KKT audit. A small gap here never produces a success status. */
double simplex_pan_pdhg_handoff_gap(const struct simplex_PanBounded *b,
	const double *x, const double *y, double scale)
{
	long double primal = simplex_pan_cone_cost(&b->a, x), dual = 0.;
	int i, j;

	for (i = 0; i < b->a.rows; i++)
		dual -= (long double)b->a.rhs[i] * (scale * y[i]);
	for (j = 0; j < b->a.columns; j++) {
		double g = b->a.objective[j] +
			   scale * simplex_pan_column_dot(&b->a, j, y);

		if (!isfinite(g))
			return HUGE_VAL;
		if (g == 0.)
			continue;
		if (g >= 0. && isfinite(b->lower[j]))
			dual += (long double)b->lower[j] * g;
		else if (g < 0. && isfinite(b->upper[j]))
			dual += (long double)b->upper[j] * g;
		else
			return HUGE_VAL; /* The dual infimum is unbounded. */
	}
	if (!isfinite(primal) || !isfinite(dual))
		return HUGE_VAL;
	if (primal == dual)
		return 0.;
	return (double)(fabsl(primal - dual) /
			fmaxl(fabsl(primal), fabsl(dual)));
}

static int pdhg_audit(const struct simplex_PanBounded *b, const double *x,
	const double *y, double scale, double *point, double *multiplier,
	const struct lp_simplex_Options *o, double *physical, double *r,
	double *work, int *improvements)
{
	double primal = simplex_pan_bounded_primal(b, x, r, work), dual;
	long double cost = simplex_pan_cone_cost(&b->a, x),
		    before = simplex_pan_cone_cost(&b->a, point);
	int i, j, k, certificate = 0;

	if (!isfinite(primal) ||
		primal > o->primal_tolerance * SIMPLEX_PAN_FEASIBLE_MARGIN ||
		!isfinite(cost))
		return 0;
	/* Absolute feasibility alone can hide error in a tiny row. Include the
	 * row coefficients in its scale so a zero-RHS, vanishing-activity row
	 * does not demand an exactly zero iterate. Native projection remains
	 * available to polish points rejected by this additional audit. */
	for (i = 0; i < b->a.rows; i++)
		work[i] = fabs(b->a.rhs[i]);
	for (j = 0; j < b->a.columns; j++)
		for (k = b->a.column_start[j]; k < b->a.column_start[j + 1];
			k++)
			work[b->a.row_index[k]] +=
				fabs(b->a.value[k]) * fmax(1., fabs(x[j]));
	for (i = 0; i < b->a.rows; i++)
		if (!isfinite(work[i]) ||
			fabs(r[i]) > o->primal_tolerance * work[i])
			return 0;
	for (i = 0; i < b->a.rows; i++) {
		physical[i] = scale * y[i];
		if (!isfinite(physical[i]))
			return 0;
	}
	certificate = simplex_pan_bounded_check_dual(
		b, x, physical, o->dual_tolerance, &dual);
	if (certificate ||
		cost < before - SIMPLEX_PDHG_INCUMBENT_GAIN_RELATIVE *
					(1. + fabsl(before))) {
		lp_simplex_memcpy(
			point, x, (size_t)b->a.columns * sizeof(double));
		lp_simplex_memcpy(multiplier, physical,
			(size_t)b->a.rows * sizeof(double));
		(*improvements)++;
		if (getenv(SIMPLEX_ENV_PROFILE))
			fprintf(stderr,
				"pan-bda pdhg audit: objective=%.12g primal=%.6g dual=%.6g certificate=%d\n",
				(double)cost, primal, dual, certificate);
	}
	return certificate;
}

void simplex_pan_pdhg_columns_destroy(struct simplex_PanPdhgColumns *plan)
{
	lp_simplex_free(plan->indices);
	lp_simplex_free(plan->kind);
	lp_simplex_memset(plan, 0, sizeof(*plan));
}

void simplex_pan_pdhg_columns_create(const struct simplex_PanStandard *a,
	struct simplex_PanPdhgColumns *plan)
{
	int j, first, count;
	lp_simplex_memset(plan, 0, sizeof(*plan));
	plan->indices = (int *)lp_simplex_malloc(
		2 * ((size_t)a->columns + 1) * sizeof(int));
	plan->kind = (unsigned char *)lp_simplex_malloc((size_t)a->columns + 1);
	if (!plan->indices || !plan->kind) {
		simplex_pan_pdhg_columns_destroy(plan);
		return;
	}
	plan->first = plan->indices;
	plan->second = plan->first + (size_t)a->columns + 1;
	for (j = 0; j < a->columns; j++) {
		first = a->column_start[j];
		count = a->column_start[j + 1] - first;
		plan->kind[j] = 0;
		if ((count != 1 && count != 2) || fabs(a->value[first]) != 1. ||
			(count == 2 && fabs(a->value[first + 1]) != 1.))
			continue;
		plan->kind[j] =
			count == 1 ? PDHG_COLUMN_SINGLE : PDHG_COLUMN_PAIR;
		if (a->value[first] < 0.)
			plan->kind[j] |= PDHG_COLUMN_NEGATIVE_FIRST;
		plan->first[j] = a->row_index[first];
		if (count == 2) {
			plan->second[j] = a->row_index[first + 1];
			if (a->value[first + 1] < 0.)
				plan->kind[j] |= PDHG_COLUMN_NEGATIVE_SECOND;
		}
		plan->compiled++;
	}
}

static double pdhg_compiled_dot(const struct simplex_PanStandard *a,
	const struct simplex_PanPdhgColumns *plan, int column,
	const double *vector)
{
	unsigned char kind = plan->kind ? plan->kind[column] : 0;
	double first, second;
	if (!kind)
		return simplex_pan_column_dot(a, column, vector);
	first = vector[plan->first[column]];
	if (kind & PDHG_COLUMN_NEGATIVE_FIRST)
		first = -first;
	if (kind & PDHG_COLUMN_PAIR) {
		second = vector[plan->second[column]];
		if (kind & PDHG_COLUMN_NEGATIVE_SECOND)
			second = -second;
		first += second;
	}
	if (isfinite(first))
		return first == 0. ? 0. : first;
	return simplex_pan_column_dot(a, column, vector);
}

/* Neumaier addition is still needed for unit coefficients: exact products
 * do not prevent cancellation when contributions meet in a row. */
static void pdhg_residual_add(double *sum, double *correction,
	const double product, const double product_error)
{
	double next = *sum + product;

	*correction += (fabs(*sum) >= fabs(product) ? (*sum - next) + product
						    : (product - next) + *sum) +
		       product_error;
	*sum = next;
}

/* Form b - A*(2*next - previous) in the same column/entry order as the
 * generic residual. Compute each extrapolated value once, where it is used,
 * rather than writing and rereading a full column workspace. Audits retain
 * the independent generic residual kernel. */
void simplex_pan_pdhg_residual(const struct simplex_PanStandard *a,
	const struct simplex_PanPdhgColumns *plan, const double *next,
	const double *previous, double *residual, double *correction)
{
	int j, k, row;

	lp_simplex_memcpy(residual, a->rhs, (size_t)a->rows * sizeof(double));
	lp_simplex_memset(correction, 0, (size_t)a->rows * sizeof(double));
	for (j = 0; j < a->columns; j++) {
		unsigned char kind = plan->kind ? plan->kind[j] : 0;
		double value = 2. * next[j] - previous[j];

		if (value == 0.)
			continue;
		if (kind && isfinite(value)) {
			row = plan->first[j];
			pdhg_residual_add(residual + row, correction + row,
				kind & PDHG_COLUMN_NEGATIVE_FIRST ? value
								  : -value,
				0.);
			if (kind & PDHG_COLUMN_PAIR) {
				row = plan->second[j];
				pdhg_residual_add(residual + row,
					correction + row,
					kind & PDHG_COLUMN_NEGATIVE_SECOND
						? value
						: -value,
					0.);
			}
		} else {
			for (k = a->column_start[j]; k < a->column_start[j + 1];
				k++) {
				double coefficient = -a->value[k];
				double product = coefficient * value;

				row = a->row_index[k];
				pdhg_residual_add(residual + row,
					correction + row, product,
					fma(coefficient, value, -product));
			}
		}
	}
	for (row = 0; row < a->rows; row++) {
		residual[row] += correction[row];
		if (!isfinite(residual[row]))
			return;
	}
}

int simplex_pan_bounded_pdhg_configured_reported(
	const struct simplex_PanBounded *b, double *point, double *multiplier,
	const struct lp_simplex_Options *o, int *iterations, int budget,
	const struct simplex_PanConfig *cfg, int *service_stop,
	int *service_output)
{
	int unobserved_stop, unobserved_output;
	struct simplex_PdhgProgress progress;
	struct simplex_PanPdhgColumns column_plan;
	const struct simplex_PanStandard *a = &b->a;
	double *storage, *x, *next, *average, *tau, *anchor_x, *y, *yaverage,
		*sigma, *r, *work, *physical, *anchor_y, *scaled_cost;
	double scale = 0., checkpoint, weight_scale = 1., rhs_norm = 0.,
	       cost_norm = 0., last_restore = HUGE_VAL;
	int i, j, k, steps = 0, certificate = 0, improvements = 0, epoch = 0,
		     restarts = 0, restorations = 0, restored = 0, phase_limit,
		     adaptive_weight, calibrate_weight,
		     stop = SIMPLEX_STOP_NUMERICAL;
	size_t storage_count;
	long metric_epoch = 0;
	clock_t begun = clock();

	if (!service_stop)
		service_stop = &unobserved_stop;
	*service_stop = SIMPLEX_STOP_NUMERICAL;
	if (!service_output)
		service_output = &unobserved_output;
	*service_output = SIMPLEX_PAN_OUTPUT_NONE;
	simplex_pan_pdhg_progress_reset(&progress);
	if (budget > o->iteration_limit - *iterations)
		budget = o->iteration_limit - *iterations;
	if (budget <= 0 || a->rows == 0) {
		*service_stop = budget <= 0 ? SIMPLEX_STOP_BUDGET
					    : SIMPLEX_STOP_NOT_APPLICABLE;
		return 0;
	}
	phase_limit = *iterations + budget;
	storage_count =
		6 * ((size_t)a->columns + 1) + 7 * ((size_t)a->rows + 1);
	storage = (double *)lp_simplex_malloc(storage_count * sizeof(double));
	if (!storage) {
		*service_stop = SIMPLEX_STOP_RESOURCE;
		return 0;
	}
	simplex_pan_pdhg_columns_create(a, &column_plan);
	x = storage;
	next = x + a->columns + 1;
	average = next + a->columns + 1;
	tau = average + a->columns + 1;
	anchor_x = tau + a->columns + 1;
	scaled_cost = anchor_x + a->columns + 1;
	y = scaled_cost + a->columns + 1;
	yaverage = y + a->rows + 1;
	sigma = yaverage + a->rows + 1;
	r = sigma + a->rows + 1;
	work = r + a->rows + 1;
	physical = work + a->rows + 1;
	anchor_y = physical + a->rows + 1;
	lp_simplex_memset(storage, 0, storage_count * sizeof(double));
	lp_simplex_memcpy(x, point, (size_t)a->columns * sizeof(double));
	for (j = 0; j < a->columns; j++)
		scale = fmax(scale, fabs(a->objective[j]));
	if (scale == 0.)
		scale = 1.;
	/* This division is invariant across all PDHG epochs and restarts. Store
	 * the rounded quotient itself, preserving the original arithmetic. */
	for (j = 0; j < a->columns; j++)
		scaled_cost[j] = a->objective[j] / scale;
	for (i = 0; i < a->rows; i++)
		y[i] = multiplier[i] / scale;
	/* Schur's bound gives ||sqrt(S) A sqrt(T)|| <= .99, without
	 * assuming full row rank or forming normal equations. */
	for (j = 0; j < a->columns; j++) {
		double sum = 0.;

		for (k = a->column_start[j]; k < a->column_start[j + 1]; k++) {
			double v = fabs(a->value[k]);

			sum += v;
			sigma[a->row_index[k]] += v;
		}
		tau[j] = sum > 0. ? SIMPLEX_PDHG_STEP_CONTRACTION / sum : 1.;
	}
	for (i = 0; i < a->rows; i++)
		sigma[i] = sigma[i] > 0.
				   ? SIMPLEX_PDHG_STEP_CONTRACTION / sigma[i]
				   : 1.;
	for (i = 0; i < a->rows; i++)
		rhs_norm = hypot(rhs_norm, sqrt(sigma[i]) * a->rhs[i]);
	for (j = 0; j < a->columns; j++)
		if (b->lower[j] != b->upper[j])
			cost_norm = hypot(cost_norm,
				sqrt(tau[j]) * (a->objective[j] / scale));
	if (rhs_norm > 0. && cost_norm > 0.)
		weight_scale = fmin(SIMPLEX_PDHG_WEIGHT_MAX,
			fmax(SIMPLEX_PDHG_WEIGHT_MIN, rhs_norm / cost_norm));
	adaptive_weight = cfg->restore && !cfg->pdhg_weight_set;
	/* A neutral step already selected by diagonal norm balancing keeps its
	 * productive trajectory. If the input norms required a nonneutral
	 * preconditioner, calibrate that scaling from the observed displacement
	 * pair. This is distinct from intervening after observed stagnation. */
	calibrate_weight = weight_scale > SIMPLEX_PDHG_WEIGHT_MIN;
	if (cfg->pdhg_weight_set) {
		double requested = cfg->pdhg_weight;

		if (isfinite(requested) &&
			requested >= SIMPLEX_PDHG_MANUAL_WEIGHT_MIN &&
			requested <= SIMPLEX_PDHG_MANUAL_WEIGHT_MAX)
			weight_scale = requested;
	}
	for (j = 0; j < a->columns; j++)
		tau[j] *= weight_scale;
	for (i = 0; i < a->rows; i++)
		sigma[i] /= weight_scale;
	lp_simplex_memcpy(anchor_x, x, (size_t)a->columns * sizeof(double));
	lp_simplex_memcpy(anchor_y, y, (size_t)a->rows * sizeof(double));
	if (getenv(SIMPLEX_ENV_PROFILE))
		fprintf(stderr, "pan-bda pdhg compiled columns=%d/%d\n",
			column_plan.compiled, a->columns);
	if (getenv(SIMPLEX_ENV_PROFILE))
		fprintf(stderr,
			"pan-bda pdhg scaling: primal-step-weight=%.6g rhs=%.6g cost=%.6g\n",
			weight_scale, rhs_norm, cost_norm);
	checkpoint = pdhg_merit(b, x, y, scale, tau, sigma, r, work);
	for (steps = 0; steps < budget && *iterations < phase_limit; steps++) {
		double weight = 1. / (++epoch);

		if ((steps % SIMPLEX_PDHG_REPORT_PERIOD) == 0 && b->problem) {
			/* Retained point supplies a feasible incumbent even when x
			 * and the ergodic iterate are temporarily infeasible. */
			simplex_progress_emit(b->problem, point, y, scale,
				"pdhg", *iterations, 0);
		}

		for (j = 0; j < a->columns; j++) {
			double g = scaled_cost[j] +
				   pdhg_compiled_dot(a, &column_plan, j, y);

			if (!isfinite(g) || !isfinite(x[j] - tau[j] * g))
				goto finish;
			next[j] = fmin(b->upper[j],
				fmax(b->lower[j], x[j] - tau[j] * g));
			average[j] += weight * (next[j] - average[j]);
		}
		simplex_pan_pdhg_residual(a, &column_plan, next, x, r, work);
		for (i = 0; i < a->rows; i++) {
			y[i] -= sigma[i] * r[i];
			yaverage[i] += weight * (y[i] - yaverage[i]);
			if (!isfinite(y[i]))
				goto finish;
		}
		/* The next loop overwrites every next entry; restoration explicitly
		 * copies its input before using that scratch buffer. Swap ownership of
		 * the two local vectors instead of copying the completed iterate. */
		{
			double *previous = x;
			x = next;
			next = previous;
		}
		(*iterations)++;
		if ((steps % SIMPLEX_PDHG_AUDIT_PERIOD) == 0 ||
			steps + 1 == budget) {
			double current, averaged, best;
			int stalled;

			certificate =
				pdhg_audit(b, x, y, scale, point, multiplier, o,
					physical, r, work, &improvements);
			if (!certificate)
				certificate = pdhg_audit(b, average, yaverage,
					scale, point, multiplier, o, physical,
					r, work, &improvements);
			if (certificate) {
				steps++;
				break;
			}
			/* Observe both candidates in the same metric. Productive epochs
			 * keep their steps and remain matrix-free; the caller already
			 * bounds total work. A repair requires a full restart window
			 * without material improvement, not just elapsed time. */
			current =
				pdhg_merit(b, x, y, scale, tau, sigma, r, work);
			averaged = pdhg_merit(b, average, yaverage, scale, tau,
				sigma, r, work);
			best = fmin(current, averaged);
			stalled = simplex_pan_pdhg_progress_stalled(
				&progress, steps + 1, best);

			if (cfg->restore && stalled) {
				double primal = simplex_pan_bounded_primal(
					b, x, r, work);
				long double ceiling =
					simplex_pan_cone_cost(a, point);

				ceiling += o->dual_tolerance *
					   SIMPLEX_PAN_RESTORE_TARGET_FRACTION *
					   (1. + fabsl(ceiling));
				if (simplex_pan_pdhg_repair_progress(primal,
					    last_restore,
					    o->primal_tolerance) &&
					simplex_pan_cone_cost(a, x) < ceiling &&
					*iterations < phase_limit) {
					struct lp_simplex_Options remaining =
						*o;
					struct simplex_PanServiceTrigger
						trigger;

					remaining.iteration_limit = phase_limit;
					trigger.kind =
						SIMPLEX_PAN_TRIGGER_PDHG_STAGNATION;
					trigger.work = *iterations;
					trigger.epoch = metric_epoch;
					trigger.primal = primal;
					trigger.previous_primal = last_restore;
					last_restore = primal;
					simplex_pan_pdhg_progress_reset(
						&progress);
					restorations++;
					lp_simplex_memcpy(next, x,
						(size_t)a->columns *
							sizeof(double));
					if (simplex_pan_service_restore(b, next,
						    &remaining, iterations,
						    ceiling, cfg, &trigger,
						    NULL) == 1) {
						long double previous_cost =
							simplex_pan_cone_cost(
								a, point);
						int previous_improvements =
							improvements;

						restored++;
						certificate = pdhg_audit(b,
							next, y, scale, point,
							multiplier, o, physical,
							r, work, &improvements);
						lp_simplex_memcpy(x, next,
							(size_t)a->columns *
								sizeof(double));
						checkpoint = pdhg_merit(b, x, y,
							scale, tau, sigma, r,
							work);
						epoch = 0;
						restarts++;
						lp_simplex_memcpy(anchor_x, x,
							(size_t)a->columns *
								sizeof(double));
						lp_simplex_memcpy(anchor_y, y,
							(size_t)a->rows *
								sizeof(double));
						if (certificate) {
							steps++;
							break;
						}
						if (simplex_pan_pdhg_handoff_gap(
							    b, next, y,
							    scale) <=
								o->dual_tolerance ||
							(restored >= 2 &&
								improvements >
									previous_improvements &&
								previous_cost - simplex_pan_cone_cost(
											a,
											point) <=
									simplex_pan_progress_threshold(
										a,
										point,
										previous_cost,
										o->dual_tolerance,
										SIMPLEX_PAN_FACE_SMALL_GAIN_UNITS))) {
							if (getenv(SIMPLEX_ENV_PROFILE))
								fprintf(stderr,
									"pan-bda pdhg handoff: restored face, native KKT pending\n");
							steps++;
							stop = SIMPLEX_STOP_HANDOFF;
							goto finish;
						}
					}
					if (getenv(SIMPLEX_ENV_PROFILE))
						fprintf(stderr,
							"pan-bda pdhg restore: reason=stagnation attempts=%d restored=%d initial=%.6g iterations=%d\n",
							restorations, restored,
							primal, *iterations);
				}
			}
			/* A successful repair resets epoch, so cached merits are only
			 * reused when neither candidate has changed since observation. */
			if (epoch >= SIMPLEX_PDHG_RESTART_MIN_WORK) {
				if (isfinite(best) &&
					(best <= checkpoint * SIMPLEX_PDHG_RESTART_REDUCTION ||
						epoch >=
							SIMPLEX_PDHG_RESTART_MAX_WORK)) {
					if (averaged < current) {
						lp_simplex_memcpy(x, average,
							(size_t)a->columns *
								sizeof(double));
						lp_simplex_memcpy(y, yaverage,
							(size_t)a->rows *
								sizeof(double));
					}
					if (adaptive_weight &&
						(calibrate_weight || stalled)) {
						double dx = 0., dy = 0., target,
						       updated, factor;

						for (j = 0; j < a->columns; j++)
							dx = hypot(dx,
								(x[j] - anchor_x[j]) /
									sqrt(tau[j] /
										weight_scale));
						for (i = 0; i < a->rows; i++)
							dy = hypot(dy,
								(y[i] - anchor_y[i]) /
									sqrt(sigma[i] *
										weight_scale));
						if (dx > 0. && dy > 0. &&
							isfinite(dx) &&
							isfinite(dy)) {
							target = fmin(
								SIMPLEX_PDHG_WEIGHT_MAX,
								fmax(SIMPLEX_PDHG_WEIGHT_MIN,
									dx / dy));
							updated = sqrt(
								weight_scale *
								target);
							updated = fmin(
								weight_scale *
									SIMPLEX_PDHG_WEIGHT_CHANGE,
								fmax(weight_scale /
										SIMPLEX_PDHG_WEIGHT_CHANGE,
									updated));
							factor = updated /
								 weight_scale;
							for (j = 0;
								j < a->columns;
								j++)
								tau[j] *=
									factor;
							for (i = 0; i < a->rows;
								i++)
								sigma[i] /=
									factor;
							if (getenv(SIMPLEX_ENV_PROFILE))
								fprintf(stderr,
									"pan-bda pdhg balance: reason=%s old=%.6g new=%.6g target=%.6g\n",
									calibrate_weight
										? "calibrate-preconditioner"
										: "stagnation",
									weight_scale,
									updated,
									target);
							if (factor != 1.) {
								metric_epoch++;
								simplex_pan_pdhg_progress_reset(
									&progress);
							}
							weight_scale = updated;
						}
					}
					checkpoint = pdhg_merit(b, x, y, scale,
						tau, sigma, r, work);
					epoch = 0;
					restarts++;
					lp_simplex_memcpy(anchor_x, x,
						(size_t)a->columns *
							sizeof(double));
					lp_simplex_memcpy(anchor_y, y,
						(size_t)a->rows *
							sizeof(double));
					if (getenv(SIMPLEX_ENV_PROFILE))
						fprintf(stderr,
							"pan-bda pdhg restart: steps=%d merit=%.6g source=%s\n",
							steps + 1, best,
							averaged < current
								? "average"
								: "current");
				}
			}
			if (getenv(SIMPLEX_ENV_PROFILE) &&
				((steps % SIMPLEX_PDHG_RESTART_MAX_WORK) == 0 ||
					steps + 1 == budget)) {
				double primal = simplex_pan_bounded_primal(
					b, x, r, work);

				fprintf(stderr,
					"pan-bda pdhg progress: steps=%d objective=%.12g primal=%.6g\n",
					steps + 1,
					(double)simplex_pan_cone_cost(a, x),
					primal);
			}
		}
	}
finish:
	*service_output = certificate ? SIMPLEX_PAN_OUTPUT_NATIVE_KKT
			  : improvements > 0
				  ? SIMPLEX_PAN_OUTPUT_NATIVE_FEASIBLE_POINT
				  : SIMPLEX_PAN_OUTPUT_NONE;
	*service_stop = certificate
				? SIMPLEX_STOP_OPTIMAL
				: (steps >= budget || *iterations >= phase_limit
						  ? SIMPLEX_STOP_BUDGET
						  : stop);
	if (b->problem && b->problem->report) {
		struct simplex_StageReport *report = b->problem->report;

		report->first_order_stop =
			certificate
				? SIMPLEX_STOP_OPTIMAL
				: (steps >= budget || *iterations >= phase_limit
						  ? SIMPLEX_STOP_BUDGET
						  : stop);
		report->work.first_order_steps += steps;
		report->work.repair_attempts += restorations;
		report->work.repair_successes += restored;
		report->first_order_merit.value =
			progress.initialized ? progress.best : NAN;
		report->first_order_merit.quality =
			progress.initialized ? SIMPLEX_EVIDENCE_RECOMPUTED
					     : SIMPLEX_EVIDENCE_UNKNOWN;
		report->first_order_merit.space = SIMPLEX_SPACE_WORKING;
		report->first_order_merit.work = progress.progress_at;
		report->first_order_merit.epoch = metric_epoch;
	}
	if (b->problem)
		simplex_progress_emit(b->problem, point, y, scale,
			"pdhg-finish", *iterations, 1);
	if (getenv(SIMPLEX_ENV_PROFILE))
		fprintf(stderr,
			"pan-bda pdhg summary: steps=%d improvements=%d restarts=%d restore=%d/%d certificate=%d cpu=%.6f\n",
			steps, improvements, restarts, restored, restorations,
			certificate,
			(double)(clock() - begun) / CLOCKS_PER_SEC);
	simplex_pan_pdhg_columns_destroy(&column_plan);
	lp_simplex_free(storage);
	return certificate;
}
