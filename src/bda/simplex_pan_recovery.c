/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_pan_recovery.c
 * @brief Standard-form BDA point repair and support recovery state.
 */

#include "policy/pan.h"
#include "simplex_pan_recovery.h"
#include "simplex_pan_config.h"
#include "simplex_pan_cone.h"
#include "simplex_pan_basis.h"
#include "simplex_stage.h"
#include "simplex_storage.h"
#include "utils.h"

#include <float.h>
#include <math.h>
#include <string.h>

long double simplex_pan_progress_threshold(const struct simplex_PanStandard *a,
	const double *point, long double reference, double tolerance,
	double units)
{
	const double *cost =
		a->audit_objective ? a->audit_objective : a->objective;
	long double scale = fabsl(reference), absolute = 0.;
	double relative = fmax(64. * DBL_EPSILON, tolerance * units);
	int j;

	for (j = 0; j < a->columns; j++)
		absolute += fabsl((long double)cost[j] * point[j]);
	if (absolute > scale)
		scale = absolute;
	return relative * scale;
}

int simplex_pan_recovery_compress(
	struct simplex_PanRecovery *s, const double *point, int columns)
{
	const unsigned char *bytes = (const unsigned char *)point;
	size_t k, n = (size_t)columns * sizeof(double);
	unsigned long first = 5381UL, second = 2166136261UL;

	for (k = 0; k < n; k++) {
		first = (first * 33UL) ^ bytes[k];
		second = (second ^ bytes[k]) * 16777619UL;
	}
	if (s->have_compressed && s->compressed_point[0] == first &&
		s->compressed_point[1] == second) {
		s->compression_skips++;
		return 0;
	}
	s->have_compressed = 1;
	s->compressed_point[0] = first;
	s->compressed_point[1] = second;
	s->compressions++;
	return 1;
}

int simplex_pan_restore_configured_reported(const struct simplex_PanStandard *a,
	double *point, const struct lp_simplex_Options *options,
	long double ceiling, int *iterations,
	const struct simplex_PanConfig *cfg, int *stop)
{
	struct simplex_PanBasis b;
	size_t bytes;
	double *storage = NULL, *saved, *trial, *zero, *rhs, *work, *delta;
	double initial, error, norm, alpha;
	int stage, j, i, cap, result = -1, blocking;
	long double cost;

	*stop = SIMPLEX_STOP_NUMERICAL;
	lp_simplex_memset(&b, 0, sizeof(b));
	if (!isfinite(ceiling) || options->primal_tolerance <= 0.)
		return -1;
	if (a->rows < 0 || a->columns < 0 ||
		!simplex_storage_bytes((size_t)a->rows + (size_t)a->columns,
			3 * sizeof(double), (size_t)-1, &bytes)) {
		*stop = SIMPLEX_STOP_RESOURCE;
		return -1;
	}
	for (j = 0; j < a->columns; j++)
		if (!isfinite(point[j]) || point[j] < 0.)
			return -1;
	storage = (double *)lp_simplex_malloc(bytes);
	if (!storage) {
		*stop = SIMPLEX_STOP_RESOURCE;
		return -1;
	}
	saved = storage;
	trial = saved + a->columns;
	zero = trial + a->columns;
	rhs = zero + a->columns;
	work = rhs + a->rows;
	delta = work + a->rows;
	lp_simplex_memcpy(saved, point, (size_t)a->columns * sizeof(double));
	lp_simplex_memset(zero, 0, (size_t)a->columns * sizeof(double));
	initial = simplex_pan_residual(a, point, rhs, work);
	if (!isfinite(initial))
		goto done;
	if (initial <= SIMPLEX_PAN_RESTORE_TARGET_FRACTION *
			       options->primal_tolerance) {
		result = 0;
		goto done;
	}
	/* Bounded local repair, not Phase I: never change zero coordinates or
	 * drop their physical mass. Grow the candidate space only on failure. */
	for (stage = 0; stage < 3; stage++) {
		cap = stage == 0 ? SIMPLEX_PAN_RESTORE_RANK_INITIAL
				 : (stage == 1 ? SIMPLEX_PAN_RESTORE_RANK_SECOND
					       : SIMPLEX_PAN_RESTORE_RANK_MAX);
		if (cap > a->rows)
			cap = a->rows;
		if (*iterations >= options->iteration_limit)
			break;
		simplex_pan_basis_destroy(&b);
		if (simplex_pan_basis_create_configured(
			    &b, a, options->pivot_tolerance, cfg)) {
			*stop = SIMPLEX_STOP_RESOURCE;
			goto done;
		}
		for (j = 0; j < a->columns && b.count < cap; j++)
			if (saved[j] > 0.) {
				if (*iterations >= options->iteration_limit)
					goto done;
				(*iterations)++;
				if (simplex_pan_basis_add(&b, j))
					goto done;
			}
		if (b.count == 0)
			break;
		if (simplex_pan_basis_repair_rank(&b, zero))
			continue;
		if (*iterations >= options->iteration_limit)
			break;
		(*iterations)++;
		(void)simplex_pan_residual(a, saved, rhs, work);
		if (simplex_pan_basis_least_squares_qr(&b, rhs, delta, &norm))
			continue;
		alpha = 1.;
		blocking = -1;
		for (i = 0; i < b.count; i++) {
			if (!isfinite(delta[i]))
				break;
			if (delta[i] < 0. &&
				-saved[b.column[i]] / delta[i] < alpha) {
				alpha = -saved[b.column[i]] / delta[i];
				blocking = b.column[i];
			}
		}
		if (i < b.count || !isfinite(alpha) || alpha <= 0.)
			continue;
		lp_simplex_memcpy(
			trial, saved, (size_t)a->columns * sizeof(double));
		for (i = 0; i < b.count; i++)
			trial[b.column[i]] =
				fma(alpha, delta[i], saved[b.column[i]]);
		if (blocking >= 0)
			trial[blocking] = 0.;
		for (j = 0; j < a->columns; j++)
			if (!isfinite(trial[j]) || trial[j] < 0.)
				break;
		if (j < a->columns)
			continue;
		error = simplex_pan_residual(a, trial, rhs, work);
		cost = simplex_pan_cone_audit_cost(a, trial);
		if (error <= SIMPLEX_PAN_FEASIBLE_MARGIN *
					options->primal_tolerance &&
			error <= SIMPLEX_PAN_REPAIR_REDUCTION * initial &&
			isfinite(cost) &&
			cost <= ceiling + 32. * DBL_EPSILON *
						  (1. + fabsl(ceiling))) {
			lp_simplex_memcpy(point, trial,
				(size_t)a->columns * sizeof(double));
			result = 1;
			break;
		}
	}
done:
	simplex_pan_basis_destroy(&b);
	lp_simplex_free(storage);
	*stop = result == 1   ? SIMPLEX_STOP_HANDOFF
		: result == 0 ? SIMPLEX_STOP_NOT_APPLICABLE
		: *iterations >= options->iteration_limit ? SIMPLEX_STOP_BUDGET
							  : *stop;
	return result;
}
