/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_pan_recovery.c
 * @brief Pan support recovery and bounded original-model corrections.
 */

#include "simplex_pan_recovery.h"
#include "simplex_pan_cone.h"
#include "simplex_pan_basis.h"
#include "utils.h"

#include <float.h>
#include <math.h>
#include <string.h>


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


int simplex_pan_restore(const struct simplex_PanStandard *a, double *point,
	const struct lp_simplex_Options *options, long double ceiling,
	int *iterations)
{
	struct simplex_PanBasis b;
	double *storage = NULL, *saved, *trial, *zero, *rhs, *work, *delta;
	double initial, error, norm, alpha;
	int stage, j, i, cap, result = -1, blocking;
	long double cost;

	lp_simplex_memset(&b, 0, sizeof(b));
	if (!isfinite(ceiling) || options->primal_tolerance <= 0.)
		return -1;
	for (j = 0; j < a->columns; j++)
		if (!isfinite(point[j]) || point[j] < 0.)
			return -1;
	storage = (double *)lp_simplex_malloc(
		(3 * (size_t)a->columns + 3 * (size_t)a->rows) *
		sizeof(double));
	if (!storage)
		return -1;
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
	if (initial <= .25 * options->primal_tolerance) {
		result = 0;
		goto done;
	}
	/* Bounded local repair, not Phase I: never change zero coordinates or
	 * drop their physical mass. Grow the candidate space only on failure. */
	for (stage = 0; stage < 3; stage++) {
		cap = stage == 0 ? 64 : (stage == 1 ? 256 : 2048);
		if (cap > a->rows)
			cap = a->rows;
		if (*iterations >= options->iteration_limit)
			break;
		simplex_pan_basis_destroy(&b);
		if (simplex_pan_basis_create(&b, a, options->pivot_tolerance))
			goto done;
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
		if (error <= .5 * options->primal_tolerance &&
			error <= .5 * initial && isfinite(cost) &&
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
	return result;
}
