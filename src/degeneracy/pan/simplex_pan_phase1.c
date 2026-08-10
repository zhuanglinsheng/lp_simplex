/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 *
 * Reliable Phase I for Pan/BDA.  This is a Lawson-Hanson active-set NNLS
 * solve of min ||Az-b|| with z >= 0.  The passive set is kept linearly
 * independent, so a zero residual is immediately a deficient feasible basis.
 */
#include "simplex_pan_phase1.h"
#include "utils.h"

#include <lp_simplex/status.h>

#include <math.h>


static double pan_phase1_norm(const double *vector, const int count)
{
	double square = 0.;
	int i;
	for (i = 0; i < count; i++)
		square += vector[i] * vector[i];
	return sqrt(square);
}


static void pan_phase1_residual(
		const struct simplex_PanBasis *basis,
		const double *value, double *residual)
{
	const struct simplex_PanStandard *standard = basis->standard;
	int i, k;
	lp_simplex_memcpy(residual, standard->rhs,
		(size_t)standard->rows * sizeof(double));
	for (i = 0; i < basis->count; i++) {
		int column = basis->column[i];
		if (value[column] == 0.)
			continue;
		for (k = standard->column_start[column];
		     k < standard->column_start[column + 1]; k++)
			residual[standard->row_index[k]] -=
				value[column] * standard->value[k];
	}
}


static double pan_phase1_column_dot(
		const struct simplex_PanStandard *standard,
		const int column, const double *vector)
{
	int k;
	double dot = 0.;
	for (k = standard->column_start[column];
	     k < standard->column_start[column + 1]; k++)
		dot += standard->value[k] * vector[standard->row_index[k]];
	return dot;
}


static int pan_phase1_entering(
		const struct simplex_PanBasis *basis,
		const double *residual, const unsigned char *blocked,
		const double tolerance)
{
	const struct simplex_PanStandard *standard = basis->standard;
	double best = tolerance;
	int entering = -1;
	int j;
	for (j = 0; j < standard->columns; j++) {
		double norm, score;
		if (basis->position[j] >= 0 || blocked[j])
			continue;
		norm = simplex_pan_column_norm(standard, j);
		if (norm == 0.)
			continue;
		score = pan_phase1_column_dot(standard, j, residual) / norm;
		if (score > best) {
			best = score;
			entering = j;
		}
	}
	return entering;
}


static int pan_phase1_positive_solution(
		struct simplex_PanBasis *basis, double *value,
		double *trial, unsigned char *remove,
		const double tolerance)
{
	double residual_norm;
	int i;
	for (;;) {
		double alpha = 1.;
		int has_nonpositive = 0;
		if (simplex_pan_basis_least_squares(basis,
			basis->standard->rhs, trial, &residual_norm) ==
		    lp_simplex_EXIT_FAILURE)
			return lp_simplex_EXIT_FAILURE;
		for (i = 0; i < basis->count; i++) {
			int column = basis->column[i];
			if (trial[i] <= tolerance) {
				double denominator = value[column] - trial[i];
				double candidate = denominator > 0.
					? value[column] / denominator : 0.;
				has_nonpositive = 1;
				if (candidate < alpha)
					alpha = candidate;
			}
		}
		if (!has_nonpositive) {
			for (i = 0; i < basis->count; i++)
				value[basis->column[i]] = trial[i];
			return lp_simplex_EXIT_SUCCESS;
		}
		lp_simplex_memset(remove, 0, (size_t)basis->count);
		for (i = 0; i < basis->count; i++) {
			int column = basis->column[i];
			value[column] += alpha * (trial[i] - value[column]);
			if (value[column] <= tolerance) {
				value[column] = 0.;
				remove[i] = 1;
			}
		}
		simplex_pan_basis_remove_marked(basis, remove);
		if (simplex_pan_basis_factorize(basis) == lp_simplex_EXIT_FAILURE)
			return lp_simplex_EXIT_FAILURE;
	}
}


int simplex_pan_phase1_run(
		struct simplex_PanBasis *basis, double *value,
		const double primal_tolerance, const double dual_tolerance,
		const int iteration_limit, int *iterations, int *status)
{
	const struct simplex_PanStandard *standard = basis->standard;
	unsigned char *blocked = NULL;
	unsigned char *remove = NULL;
	double *residual = NULL;
	double *trial = NULL;
	double rhs_norm, residual_norm;
	int entering;
	int state = lp_simplex_EXIT_FAILURE;
	blocked = (unsigned char *)lp_simplex_malloc(
		(size_t)(standard->columns > 0 ? standard->columns : 1));
	remove = (unsigned char *)lp_simplex_malloc((size_t)standard->rows);
	residual = (double *)lp_simplex_malloc(
		(size_t)standard->rows * sizeof(double));
	trial = (double *)lp_simplex_malloc(
		(size_t)standard->rows * sizeof(double));
	if (blocked == NULL || remove == NULL || residual == NULL || trial == NULL) {
		*status = lp_simplex_MemoryAllocError;
		goto cleanup;
	}
	lp_simplex_memset(blocked, 0, (size_t)standard->columns);
	lp_simplex_memset(value, 0, (size_t)standard->columns * sizeof(double));
	rhs_norm = pan_phase1_norm(standard->rhs, standard->rows);
	for (;;) {
		if (*iterations >= iteration_limit) {
			*status = lp_simplex_ExceedIterLimit;
			goto cleanup;
		}
		pan_phase1_residual(basis, value, residual);
		residual_norm = pan_phase1_norm(residual, standard->rows);
		if (residual_norm <= primal_tolerance * (1. + rhs_norm)) {
			*status = lp_simplex_CondUnsatisfied;
			state = lp_simplex_EXIT_SUCCESS;
			goto cleanup;
		}
		entering = pan_phase1_entering(basis, residual, blocked,
			dual_tolerance * (1. + residual_norm));
		if (entering < 0) {
			*status = lp_simplex_Infeasibility;
			goto cleanup;
		}
		if (simplex_pan_basis_add(basis, entering) ==
		    lp_simplex_EXIT_FAILURE) {
			*status = lp_simplex_MemoryAllocError;
			goto cleanup;
		}
		if (simplex_pan_basis_factorize(basis) == lp_simplex_EXIT_FAILURE) {
			lp_simplex_memset(remove, 0, (size_t)basis->count);
			remove[basis->count - 1] = 1;
			simplex_pan_basis_remove_marked(basis, remove);
			blocked[entering] = 1;
			continue;
		}
		(*iterations)++;
		{
			int count_before_positive_solve = basis->count;
		if (pan_phase1_positive_solution(basis, value, trial, remove,
			primal_tolerance) == lp_simplex_EXIT_FAILURE) {
			*status = lp_simplex_PrecisionError;
			goto cleanup;
		}
			/* Rank failures are relative to the current passive space.  A
			 * deletion can make an earlier column independent again, so only
			 * retain the block on a column that immediately left once more. */
			if (basis->count < count_before_positive_solve) {
				lp_simplex_memset(blocked, 0,
					(size_t)standard->columns);
				if (basis->position[entering] < 0)
					blocked[entering] = 1;
			}
		}
	}
cleanup:
	lp_simplex_free(blocked);
	lp_simplex_free(remove);
	lp_simplex_free(residual);
	lp_simplex_free(trial);
	return state;
}
