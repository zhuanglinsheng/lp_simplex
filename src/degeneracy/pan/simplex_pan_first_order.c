/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_pan_first_order.c
 * @brief Bounded primal-dual direction prediction and dual certification for
 * Pan/BDA.
 *
 * Sparse products generate candidates; original LP audits decide acceptance.
 */

#include "simplex_pan_config.h"
#include "simplex_pan_cone.h"
#include "simplex_pan_basis.h"
#include "utils.h"

#include <lp_simplex/status.h>

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>


/* Smooth dual feasibility objective. Positive primal coordinates impose an
 * equality, others a one-sided inequality. Column normalization and a row
 * diagonal preconditioner affect only candidate search, never acceptance.
 * No Gram matrix or QR factor is formed. */
static double pan_dual_loss(const struct simplex_PanStandard *a,
	const double *point, const unsigned char *free_mask, const double *map,
	const double *y, double scale, double *gradient)
{
	double loss = 0.;
	int i, j, k;

	if (gradient)
		lp_simplex_memset(
			gradient, 0, (size_t)a->rows * sizeof(double));
	for (j = 0; j < a->columns; j++) {
		double norm = a->column_norm[j], r;

		if (norm == 0.)
			continue;
		r = (a->objective[j] / scale) / norm;
		for (k = a->column_start[j]; k < a->column_start[j + 1]; k++) {
			i = a->row_index[k];
			r += (a->value[k] / norm) * map[i] * y[i];
		}
		if (!(free_mask ? free_mask[j] : point[j] > 0.) && r >= 0.)
			continue;
		loss += .5 * r * r;
		if (gradient)
			for (k = a->column_start[j]; k < a->column_start[j + 1];
				k++) {
				i = a->row_index[k];
				gradient[i] +=
					(a->value[k] / norm) * map[i] * r;
			}
	}
	return loss;
}


/* Matrix-free Newton correction on a fixed violated-inequality pattern.
 * The Hessian is U_active U_active^T and may be singular. Preconditioned CG
 * is only a candidate generator; a convex-loss line search and the original
 * direction/dual audits remain mandatory. */
static int pan_dual_cg(const struct simplex_PanStandard *a, const double *point,
	const unsigned char *free_mask, const double *map, double scale,
	const double *x, const double *gradient, double *delta, int budget)
{
	double *storage, *r, *p, *h, *diagonal, *z;
	unsigned char *active;
	double rz = 0., initial = 0.;
	int i, j, k, steps = 0;

	storage = (double *)lp_simplex_malloc(
		5 * (size_t)a->rows * sizeof(double));
	active = (unsigned char *)lp_simplex_malloc((size_t)a->columns);
	if (!storage || !active) {
		lp_simplex_free(storage);
		lp_simplex_free(active);
		return 0;
	}
	r = storage;
	p = r + a->rows;
	h = p + a->rows;
	diagonal = h + a->rows;
	z = diagonal + a->rows;
	lp_simplex_memset(diagonal, 0, (size_t)a->rows * sizeof(double));
	lp_simplex_memset(delta, 0, (size_t)a->rows * sizeof(double));
	for (j = 0; j < a->columns; j++) {
		double norm = a->column_norm[j], value;

		active[j] = 0;
		if (norm == 0.)
			continue;
		value = (a->objective[j] / scale) / norm;
		for (k = a->column_start[j]; k < a->column_start[j + 1]; k++)
			value += (a->value[k] / norm) * map[a->row_index[k]] *
				 x[a->row_index[k]];
		active[j] = (free_mask ? free_mask[j] : point[j] > 0.) ||
			    value < 0.;
		if (active[j])
			for (k = a->column_start[j]; k < a->column_start[j + 1];
				k++) {
				double v = (a->value[k] / norm) *
					   map[a->row_index[k]];

				diagonal[a->row_index[k]] += v * v;
			}
	}
	for (i = 0; i < a->rows; i++) {
		r[i] = -gradient[i];
		initial = hypot(initial, r[i]);
		z[i] = diagonal[i] > 0. ? r[i] / diagonal[i] : r[i];
		p[i] = z[i];
		rz += r[i] * z[i];
	}
	for (steps = 0; steps < budget && rz > 0. && isfinite(rz); steps++) {
		double denominator = 0., alpha, next_rz = 0., norm = 0., beta;

		lp_simplex_memset(h, 0, (size_t)a->rows * sizeof(double));
		for (j = 0; j < a->columns; j++)
			if (active[j]) {
				double product = 0.,
				       column_norm = a->column_norm[j];

				for (k = a->column_start[j];
					k < a->column_start[j + 1]; k++)
					product += (a->value[k] / column_norm) *
						   map[a->row_index[k]] *
						   p[a->row_index[k]];
				for (k = a->column_start[j];
					k < a->column_start[j + 1]; k++)
					h[a->row_index[k]] +=
						(a->value[k] / column_norm) *
						map[a->row_index[k]] * product;
			}
		for (i = 0; i < a->rows; i++)
			denominator += p[i] * h[i];
		if (!(denominator > 0.) || !isfinite(denominator))
			break;
		alpha = rz / denominator;
		if (!isfinite(alpha))
			break;
		for (i = 0; i < a->rows; i++) {
			delta[i] = fma(alpha, p[i], delta[i]);
			r[i] = fma(-alpha, h[i], r[i]);
			norm = hypot(norm, r[i]);
			z[i] = diagonal[i] > 0. ? r[i] / diagonal[i] : r[i];
			next_rz += r[i] * z[i];
		}
		if (norm <= fmax(1e-15, initial * 1e-10)) {
			steps++;
			break;
		}
		beta = next_rz / rz;
		for (i = 0; i < a->rows; i++)
			p[i] = z[i] + beta * p[i];
		rz = next_rz;
	}
	lp_simplex_free(storage);
	lp_simplex_free(active);
	return steps;
}


static int pan_first_order(const struct simplex_PanStandard *a,
	const double *point, double *multiplier,
	const struct lp_simplex_Options *options, int *iterations,
	double *violation, double *direction, const unsigned char *free_mask,
	int warm, int requested, int newton,
	const struct simplex_PanConfig *config,
	int (*accept)(void *, const double *), void *context)
{
	double *storage, *map, *row_sum, *x, *y, *next, *gradient;
	double scale = 0., lipschitz = 0., t = 1., loss, initial, checkpoint,
	       direction_error = HUGE_VAL, primal_error;
	int i, j, k, steps = 0, certificate = 0, limit, end_limit;
	clock_t start = clock();

	limit = requested;
	if (limit < 1)
		return 0;
	if (limit > options->iteration_limit - *iterations)
		limit = options->iteration_limit - *iterations;
	if (limit <= 0)
		return 0;
	end_limit = *iterations + limit;
	storage = (double *)lp_simplex_malloc(
		6 * (size_t)a->rows * sizeof(double));
	if (!storage)
		return 0;
	map = storage;
	row_sum = map + a->rows;
	x = row_sum + a->rows;
	y = x + a->rows;
	next = y + a->rows;
	gradient = next + a->rows;
	primal_error = simplex_pan_residual(a, point, x, gradient);
	if (primal_error > options->primal_tolerance)
		goto done;
	for (j = 0; j < a->columns; j++) {
		if (!isfinite(point[j]) || point[j] < 0.)
			goto done;
		if (!config->dual_center || point[j] > 0. ||
			a->objective[j] < 0.)
			scale = fmax(scale, fabs(a->objective[j]));
	}
	if (scale == 0.)
		scale = 1.;
	lp_simplex_memset(storage, 0, 6 * (size_t)a->rows * sizeof(double));
	for (j = 0; j < a->columns; j++)
		if (a->column_norm[j] > 0.)
			for (k = a->column_start[j]; k < a->column_start[j + 1];
				k++) {
				double v = a->value[k] / a->column_norm[j];

				map[a->row_index[k]] += v * v;
			}
	for (i = 0; i < a->rows; i++)
		map[i] = map[i] > 0. ? 1. / sqrt(map[i]) : 1.;
	/* Infinity norm of |U||U|^T bounds the spectral Lipschitz constant
	 * for every subset of currently violated columns. */
	for (j = 0; j < a->columns; j++)
		if (a->column_norm[j] > 0.) {
			double sum = 0.;

			for (k = a->column_start[j]; k < a->column_start[j + 1];
				k++)
				sum += fabs((a->value[k] / a->column_norm[j]) *
					    map[a->row_index[k]]);
			for (k = a->column_start[j]; k < a->column_start[j + 1];
				k++)
				row_sum[a->row_index[k]] +=
					sum *
					fabs((a->value[k] / a->column_norm[j]) *
						map[a->row_index[k]]);
		}
	for (i = 0; i < a->rows; i++)
		lipschitz = fmax(lipschitz, row_sum[i]);
	if (!isfinite(lipschitz) || lipschitz <= 0.)
		goto done;
	lipschitz *= 1.01;
	loss = initial = checkpoint =
		pan_dual_loss(a, point, free_mask, map, x, scale, NULL);
	if (!isfinite(loss))
		goto done;
	if (warm) {
		double candidate_loss;

		for (i = 0; i < a->rows; i++) {
			x[i] = (multiplier[i] / scale) / map[i];
			if (!isfinite(x[i]))
				break;
		}
		candidate_loss = i == a->rows
					 ? pan_dual_loss(a, point, free_mask,
						   map, x, scale, NULL)
					 : HUGE_VAL;
		if (isfinite(candidate_loss) && candidate_loss < loss) {
			loss = initial = checkpoint = candidate_loss;
			lp_simplex_memcpy(
				y, x, (size_t)a->rows * sizeof(double));
		} else
			lp_simplex_memset(
				x, 0, (size_t)a->rows * sizeof(double));
	}
	for (steps = 0; steps < limit && *iterations < end_limit; steps++) {
		double current = pan_dual_loss(
			       a, point, free_mask, map, y, scale, gradient),
		       trial, new_t, beta;

		if (!isfinite(current))
			break;
		for (i = 0; i < a->rows; i++)
			next[i] = y[i] - gradient[i] / lipschitz;
		trial = pan_dual_loss(
			a, point, free_mask, map, next, scale, NULL);
		(*iterations)++;
		if (!isfinite(trial))
			break;
		/* Monotone restart prevents extrapolation from preserving a poor
		 * active inequality pattern. It does not change the certificate. */
		if (trial > loss + 32. * DBL_EPSILON * fmax(loss, DBL_MIN)) {
			lp_simplex_memcpy(
				y, x, (size_t)a->rows * sizeof(double));
			t = 1.;
			continue;
		}
		new_t = (1. + sqrt(1. + 4. * t * t)) * .5;
		beta = (t - 1.) / new_t;
		for (i = 0; i < a->rows; i++) {
			y[i] = next[i] + beta * (next[i] - x[i]);
			x[i] = next[i];
		}
		t = new_t;
		loss = trial;
		if (direction && newton && *iterations < end_limit &&
			steps >= (accept ? 32 : 512) &&
			(steps % (accept ? 32 : 128)) == 0 &&
			limit - steps > 1) {
			double old = pan_dual_loss(a, point, free_mask, map, x,
				       scale, gradient),
			       step = 1., dot = 0.;
			int cg_budget = limit - steps - 1, used, line;

			if (cg_budget > end_limit - *iterations)
				cg_budget = end_limit - *iterations;
			if (cg_budget > 256)
				cg_budget = 256;
			if (config->newton_qr &&
				simplex_pan_dual_newton_qr(a, point, free_mask,
					map, x, scale,
					next) == lp_simplex_EXIT_SUCCESS)
				used = 1;
			else
				used = pan_dual_cg(a, point, free_mask, map,
					scale, x, gradient, next, cg_budget);
			steps += used;
			*iterations += used;
			for (i = 0; i < a->rows; i++)
				dot += gradient[i] * next[i];
			if (isfinite(dot) && dot < 0.)
				for (line = 0; line < 20; line++) {
					double candidate_loss;

					for (i = 0; i < a->rows; i++)
						y[i] = x[i] + step * next[i];
					candidate_loss = pan_dual_loss(a, point,
						free_mask, map, y, scale, NULL);
					if (isfinite(candidate_loss) &&
						candidate_loss <=
							old +
								1e-4 * step *
									dot +
								32. * DBL_EPSILON *
									fmax(old,
										DBL_MIN)) {
						lp_simplex_memcpy(x, y,
							(size_t)a->rows *
								sizeof(double));
						loss = candidate_loss;
						break;
					}
					step *= .5;
				}
			lp_simplex_memcpy(
				y, x, (size_t)a->rows * sizeof(double));
			t = 1.;
		}
		if ((steps % 16) == 0 || steps + 1 == limit) {
			for (i = 0; i < a->rows; i++)
				multiplier[i] = scale * map[i] * x[i];
			if (simplex_pan_cone_check_dual(a, point, multiplier,
				    options->dual_tolerance * .1, violation)) {
				certificate = 1;
				steps++;
				break;
			}
			if (direction && loss > DBL_EPSILON * DBL_EPSILON) {
				struct simplex_PanStandard homogeneous = *a;
				double slope, error;

				for (j = 0; j < a->columns; j++) {
					double norm = a->column_norm[j],
					       r = a->objective[j] / scale;

					for (k = a->column_start[j];
						k < a->column_start[j + 1]; k++)
						r += a->value[k] *
						     map[a->row_index[k]] *
						     x[a->row_index[k]];
					direction[j] =
						norm > 0. && (free_mask[j] ||
								     r < 0.)
							? -(r / norm) / norm
							: 0.;
				}
				slope = (double)simplex_pan_cone_audit_cost(
					a, direction);
				if (accept) {
					if (*iterations < end_limit &&
						isfinite(slope) && slope < 0. &&
						accept(context, direction)) {
						certificate = 2;
						steps++;
						break;
					}
				} else if (isfinite(slope) && slope < 0.) {
					int has_negative = 0, valid = 1;
					double alpha = HUGE_VAL,
					       inner_limit = fmax(
						       16. * DBL_EPSILON,
						       fmin(options->primal_tolerance,
							       options->dual_tolerance) *
							       1e-3);

					for (j = 0; j < a->columns; j++)
						direction[j] *= -scale / slope;
					lp_simplex_memset(row_sum, 0,
						(size_t)a->rows *
							sizeof(double));
					homogeneous.rhs = row_sum;
					error = simplex_pan_residual(
						&homogeneous, direction,
						gradient, next);
					direction_error = error;
					for (j = 0; j < a->columns; j++) {
						if (direction[j] < 0.) {
							has_negative = 1;
							alpha = fmin(alpha,
								-point[j] /
									direction
										[j]);
						}
						for (k = a->column_start[j];
							k <
							a->column_start[j + 1];
							k++)
							row_sum[a->row_index
									[k]] +=
								fabs(a->value[k] *
									direction
										[j]);
					}
					/* A candidate recession ray needs a relative cancellation
					 * check as well: a tiny nonzero column cannot prove a ray.
					 * Bounded steps additionally undergo the primal point audit. */
					if (!has_negative)
						for (i = 0; i < a->rows; i++)
							if (fabs(gradient[i]) >
								32. * DBL_EPSILON *
									(a->columns +
										2.) *
									row_sum[i])
								valid = 0;
					/* Certify error in the actual bounded step, rather than
					 * oversolving every direction to one absolute threshold.
					 * Keep a margin for rounding; the outer compensated primal
					 * audit still accepts or rolls back the complete trial point. */
					if (isfinite(error) && valid &&
						(error <= inner_limit ||
							(has_negative &&
								isfinite(
									alpha) &&
								alpha > 0. &&
								error * alpha <=
									.25 * (options->primal_tolerance -
										      primal_error)))) {
						certificate = 2;
						steps++;
						break;
					}
				}
			}
		}
		if (loss == 0.) {
			steps++;
			break;
		}
		if (!direction && steps > 0 && (steps % 128) == 0) {
			if (loss >= checkpoint * .99999) {
				steps++;
				break;
			}
			checkpoint = loss;
		}
	}
	if (getenv("LP_SIMPLEX_PROFILE")) {
		if (direction)
			fprintf(stderr,
				"pan-bda first order: steps=%d initial=%.6g final=%.6g outcome=%s direction-error=%.6g seconds=%.6f\n",
				steps, initial, loss,
				certificate == 2   ? "direction"
				: certificate == 1 ? "dual"
						   : "prediction",
				direction_error,
				(double)(clock() - start) / CLOCKS_PER_SEC);
		else
			fprintf(stderr,
				"pan-bda dual polish: steps=%d initial=%.6g final=%.6g certificate=%d seconds=%.6f\n",
				steps, initial, loss, certificate,
				(double)(clock() - start) / CLOCKS_PER_SEC);
	}
done:
	lp_simplex_free(storage);
	return certificate;
}


int simplex_pan_cone_first_order(const struct simplex_PanStandard *a,
	const double *point, double *multiplier,
	const struct lp_simplex_Options *options, int *iterations,
	double *violation, double *direction, const unsigned char *free_mask,
	int warm)
{
	struct simplex_PanConfig config;
	int budget;

	simplex_pan_config_init(&config, options);
	budget = direction ? config.first_order_iterations
			   : config.dual_polish_iterations;
	if (budget <= 0)
		return 0;
	return pan_first_order(a, point, multiplier, options, iterations,
		violation, direction, free_mask, warm, budget, config.newton_qr,
		&config, NULL, NULL);
}


int simplex_pan_cone_predict_audited(const struct simplex_PanStandard *a,
	const double *point, double *multiplier,
	const struct lp_simplex_Options *options, int *iterations,
	double *violation, double *direction, const unsigned char *free_mask,
	int budget, int (*accept)(void *, const double *), void *context)
{
	struct simplex_PanConfig config;

	simplex_pan_config_init(&config, options);
	return pan_first_order(a, point, multiplier, options, iterations,
		violation, direction, free_mask, 1, budget, 1, &config, accept,
		context);
}
