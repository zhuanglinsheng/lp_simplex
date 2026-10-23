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

#include "simplex_pan_config.h"
#include "simplex_progress.h"
#include "simplex_pan_bounded.h"
#include "simplex_pan_cone.h"
#include "utils.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>


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
static double pdhg_handoff_gap(const struct simplex_PanBounded *b,
	const double *x, const double *y, double scale)
{
	long double primal = simplex_pan_cone_cost(&b->a, x), dual = 0.;
	int i, j;

	for (i = 0; i < b->a.rows; i++)
		dual -= (long double)b->a.rhs[i] * (scale * y[i]);
	for (j = 0; j < b->a.columns; j++) {
		double g = b->a.objective[j] +
			   scale * simplex_pan_column_dot(&b->a, j, y);

		if (g >= 0. && isfinite(b->lower[j]))
			dual += (long double)b->lower[j] * g;
		else if (g < 0. && isfinite(b->upper[j]))
			dual += (long double)b->upper[j] * g;
	}
	return (double)(fabsl(primal - dual) / (1. + fabsl(primal)));
}


static int pdhg_audit(const struct simplex_PanBounded *b, const double *x,
	const double *y, double scale, double *point, double *multiplier,
	const struct lp_simplex_Options *o, double *physical, double *r,
	double *work, int *improvements)
{
	double primal = simplex_pan_bounded_primal(b, x, r, work), dual;
	long double cost = simplex_pan_cone_cost(&b->a, x),
		    before = simplex_pan_cone_cost(&b->a, point);
	int i, certificate = 0;

	if (!isfinite(primal) || primal > o->primal_tolerance * .5 ||
		!isfinite(cost))
		return 0;
	for (i = 0; i < b->a.rows; i++) {
		physical[i] = scale * y[i];
		if (!isfinite(physical[i]))
			return 0;
	}
	certificate = simplex_pan_bounded_check_dual(
		b, x, physical, o->dual_tolerance, &dual);
	if (certificate || cost < before - 1e-12 * (1. + fabsl(before))) {
		lp_simplex_memcpy(
			point, x, (size_t)b->a.columns * sizeof(double));
		lp_simplex_memcpy(multiplier, physical,
			(size_t)b->a.rows * sizeof(double));
		(*improvements)++;
		if (getenv("LP_SIMPLEX_PROFILE"))
			fprintf(stderr,
				"pan-bda pdhg audit: objective=%.12g primal=%.6g dual=%.6g certificate=%d\n",
				(double)cost, primal, dual, certificate);
	}
	return certificate;
}


int simplex_pan_bounded_pdhg(const struct simplex_PanBounded *b, double *point,
	double *multiplier, const struct lp_simplex_Options *o, int *iterations,
	int budget)
{
	struct simplex_PanConfig config;
	const struct simplex_PanStandard *a = &b->a;
	double *storage, *x, *next, *average, *tau, *extrapolate, *anchor_x, *y,
		*yaverage, *sigma, *r, *work, *physical, *anchor_y;
	double scale = 0., checkpoint, weight_scale = 1., rhs_norm = 0.,
	       cost_norm = 0., last_restore = HUGE_VAL;
	int i, j, k, steps = 0, certificate = 0, improvements = 0, epoch = 0,
		     restarts = 0, restorations = 0, restored = 0, phase_limit,
		     adaptive_weight;
	clock_t begun = clock();

	simplex_pan_config_init(&config, o);
	if (budget > o->iteration_limit - *iterations)
		budget = o->iteration_limit - *iterations;
	if (budget <= 0 || a->rows == 0)
		return 0;
	phase_limit = *iterations + budget;
	storage = (double *)lp_simplex_malloc(
		(6 * ((size_t)a->columns + 1) + 7 * ((size_t)a->rows + 1)) *
		sizeof(double));
	if (!storage)
		return 0;
	x = storage;
	next = x + a->columns + 1;
	average = next + a->columns + 1;
	tau = average + a->columns + 1;
	extrapolate = tau + a->columns + 1;
	anchor_x = extrapolate + a->columns + 1;
	y = anchor_x + a->columns + 1;
	yaverage = y + a->rows + 1;
	sigma = yaverage + a->rows + 1;
	r = sigma + a->rows + 1;
	work = r + a->rows + 1;
	physical = work + a->rows + 1;
	anchor_y = physical + a->rows + 1;
	lp_simplex_memset(storage, 0,
		(6 * ((size_t)a->columns + 1) + 7 * ((size_t)a->rows + 1)) *
			sizeof(double));
	lp_simplex_memcpy(x, point, (size_t)a->columns * sizeof(double));
	for (j = 0; j < a->columns; j++)
		scale = fmax(scale, fabs(a->objective[j]));
	if (scale == 0.)
		scale = 1.;
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
		tau[j] = sum > 0. ? .99 / sum : 1.;
	}
	for (i = 0; i < a->rows; i++)
		sigma[i] = sigma[i] > 0. ? .99 / sigma[i] : 1.;
	for (i = 0; i < a->rows; i++)
		rhs_norm = hypot(rhs_norm, sqrt(sigma[i]) * a->rhs[i]);
	for (j = 0; j < a->columns; j++)
		if (b->lower[j] != b->upper[j])
			cost_norm = hypot(cost_norm,
				sqrt(tau[j]) * (a->objective[j] / scale));
	if (rhs_norm > 0. && cost_norm > 0.)
		weight_scale = fmin(1e6, fmax(1., rhs_norm / cost_norm));
	adaptive_weight =
		config.restore && !config.pdhg_weight_set && weight_scale >= 4.;
	if (config.pdhg_weight_set) {
		double requested = config.pdhg_weight;

		if (isfinite(requested) && requested >= 1e-12 &&
			requested <= 1e12)
			weight_scale = requested;
	}
	for (j = 0; j < a->columns; j++)
		tau[j] *= weight_scale;
	for (i = 0; i < a->rows; i++)
		sigma[i] /= weight_scale;
	lp_simplex_memcpy(anchor_x, x, (size_t)a->columns * sizeof(double));
	lp_simplex_memcpy(anchor_y, y, (size_t)a->rows * sizeof(double));
	if (getenv("LP_SIMPLEX_PROFILE"))
		fprintf(stderr,
			"pan-bda pdhg scaling: primal-step-weight=%.6g rhs=%.6g cost=%.6g\n",
			weight_scale, rhs_norm, cost_norm);
	checkpoint = pdhg_merit(b, x, y, scale, tau, sigma, r, work);
	for (steps = 0; steps < budget && *iterations < phase_limit; steps++) {
		double weight = 1. / (++epoch);

		if ((steps % 256) == 0 && b->problem) {
			/* Retained point supplies a feasible incumbent even when x
			 * and the ergodic iterate are temporarily infeasible. */
			simplex_progress_emit(b->problem, point, y, scale,
				"pdhg", *iterations, 0);
		}

		for (j = 0; j < a->columns; j++) {
			double g = a->objective[j] / scale +
				   simplex_pan_column_dot(a, j, y);

			if (!isfinite(g) || !isfinite(x[j] - tau[j] * g))
				goto finish;
			next[j] = fmin(b->upper[j],
				fmax(b->lower[j], x[j] - tau[j] * g));
			extrapolate[j] = 2. * next[j] - x[j];
			average[j] += weight * (next[j] - average[j]);
		}
		(void)simplex_pan_residual(a, extrapolate, r, work);
		for (i = 0; i < a->rows; i++) {
			y[i] -= sigma[i] * r[i];
			yaverage[i] += weight * (y[i] - yaverage[i]);
			if (!isfinite(y[i]))
				goto finish;
		}
		lp_simplex_memcpy(x, next, (size_t)a->columns * sizeof(double));
		(*iterations)++;
		if ((steps % 128) == 0 || steps + 1 == budget) {
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
			/* Give well-scaled matrix-free iterations time to converge before
			 * paying for a rank-revealing face repair. Poorly balanced models
			 * can request the repair earlier, after step-weight adaptation. */
			if (steps >= (weight_scale >= 4. ? 8192 : 131072) &&
				config.restore) {
				double primal = simplex_pan_bounded_primal(
					b, x, r, work);
				long double ceiling =
					simplex_pan_cone_cost(a, point);

				ceiling += o->dual_tolerance * .25 *
					   (1. + fabsl(ceiling));
				if (primal > o->primal_tolerance * .25 &&
					primal < 1e-3 &&
					primal < last_restore * .5 &&
					simplex_pan_cone_cost(a, x) < ceiling &&
					*iterations < phase_limit) {
					struct lp_simplex_Options remaining =
						*o;

					remaining.iteration_limit = phase_limit;
					last_restore = primal;
					restorations++;
					lp_simplex_memcpy(next, x,
						(size_t)a->columns *
							sizeof(double));
					if (simplex_pan_bounded_restore(b, next,
						    &remaining, iterations,
						    ceiling) == 1) {
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
						if (pdhg_handoff_gap(b, next, y,
							    scale) <=
								o->dual_tolerance ||
							(restored >= 2 &&
								improvements >
									previous_improvements &&
								previous_cost - simplex_pan_cone_cost(
											a,
											point) <=
									1e-6 * (1. + fabsl(previous_cost)))) {
							if (getenv("LP_SIMPLEX_PROFILE"))
								fprintf(stderr,
									"pan-bda pdhg handoff: restored face, native KKT pending\n");
							steps++;
							goto finish;
						}
					}
					if (getenv("LP_SIMPLEX_PROFILE"))
						fprintf(stderr,
							"pan-bda pdhg restore: attempts=%d restored=%d initial=%.6g iterations=%d\n",
							restorations, restored,
							primal, *iterations);
				}
			}
			if (epoch >= 256) {
				double current = pdhg_merit(
					b, x, y, scale, tau, sigma, r, work);
				double averaged = pdhg_merit(b, average,
					yaverage, scale, tau, sigma, r, work);
				double best = fmin(current, averaged);

				if (isfinite(best) &&
					(best <= checkpoint * .5 ||
						epoch >= 8192)) {
					if (averaged < current) {
						lp_simplex_memcpy(x, average,
							(size_t)a->columns *
								sizeof(double));
						lp_simplex_memcpy(y, yaverage,
							(size_t)a->rows *
								sizeof(double));
					}
					if (adaptive_weight) {
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
							target = fmin(1e6,
								fmax(1.,
									dx / dy));
							updated = sqrt(
								weight_scale *
								target);
							updated = fmin(
								weight_scale *
									4.,
								fmax(weight_scale *
										.25,
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
							if (getenv("LP_SIMPLEX_PROFILE"))
								fprintf(stderr,
									"pan-bda pdhg balance: old=%.6g new=%.6g target=%.6g\n",
									weight_scale,
									updated,
									target);
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
					if (getenv("LP_SIMPLEX_PROFILE"))
						fprintf(stderr,
							"pan-bda pdhg restart: steps=%d merit=%.6g source=%s\n",
							steps + 1, best,
							averaged < current
								? "average"
								: "current");
				}
			}
			if (getenv("LP_SIMPLEX_PROFILE") &&
				((steps % 8192) == 0 || steps + 1 == budget)) {
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
	if (b->problem)
		simplex_progress_emit(b->problem, point, y, scale,
			"pdhg-finish", *iterations, 1);
	if (getenv("LP_SIMPLEX_PROFILE"))
		fprintf(stderr,
			"pan-bda pdhg summary: steps=%d improvements=%d restarts=%d restore=%d/%d certificate=%d cpu=%.6f\n",
			steps, improvements, restarts, restored, restorations,
			certificate,
			(double)(clock() - begun) / CLOCKS_PER_SEC);
	lp_simplex_free(storage);
	return certificate;
}
