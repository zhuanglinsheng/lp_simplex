/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_pan_phase1.c
 * @brief Reliable Phase I for Pan/BDA.
 *
 * This is a Lawson-Hanson active-set NNLS solve of min ||Az-b|| with z >= 0.
 * The passive set is kept linearly independent, so a zero residual is
 * immediately a deficient feasible basis.
 */

#include "simplex_pan_phase1.h"
#include "utils.h"

#include <lp_simplex/status.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>


static double pan_phase1_norm(const double *vector, const int count)
{
	double square = 0.;
	int i;

	for (i = 0; i < count; i++)
		square = hypot(square, vector[i]);
	return square;
}


static void pan_phase1_residual(const struct simplex_PanBasis *basis,
	const double *value, double *residual, double *correction)
{
	simplex_pan_residual(basis->standard, value, residual, correction);
}


static double pan_phase1_column_dot(const struct simplex_PanStandard *standard,
	const int column, const double *vector)
{
	return simplex_pan_column_dot(standard, column, vector);
}


static int pan_phase1_entering(const struct simplex_PanBasis *basis,
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


/* A bounded shortlist from one global correlation scan. QR still checks each
 * insertion; sharing a least-squares solve never implies sharing a rank test.
 * Equal scores retain CSC order for reproducible, deterministic selection. */
int simplex_pan_nnls_shortlist(const struct simplex_PanBasis *basis,
	const double *residual, const unsigned char *blocked, double tolerance,
	int limit, int *candidate, int implicit_cost)
{
	double score[32];
	int j, count = 0;

	for (j = 0; j < basis->standard->columns; j++) {
		double norm, current;
		int p;

		if (basis->position[j] >= 0 || blocked[j])
			continue;
		norm = simplex_pan_column_norm(basis->standard, j);
		if (implicit_cost)
			norm = hypot(norm, basis->standard->objective[j]);
		if (norm == 0.)
			continue;
		current = pan_phase1_column_dot(basis->standard, j, residual) /
			  norm;
		if (implicit_cost)
			current += basis->standard->objective[j] *
				   residual[basis->standard->rows] / norm;
		if (!(current > tolerance))
			continue;
		p = count;
		if (p == limit) {
			if (current <= score[p - 1])
				continue;
			p--;
		} else
			count++;
		while (p > 0 && current > score[p - 1]) {
			score[p] = score[p - 1];
			candidate[p] = candidate[p - 1];
			p--;
		}
		score[p] = current;
		candidate[p] = j;
	}
	/* Do not fill a block with arbitrarily weak correlations. */
	while (count > 1 && score[count - 1] < .25 * score[0])
		count--;
	return count;
}


static int pan_phase1_positive_solution(struct simplex_PanBasis *basis,
	double *value, double *trial, unsigned char *remove,
	const unsigned char *free_mask)
{
	double residual_norm;
	int i;

	for (;;) {
		double alpha = 1.;
		int has_nonpositive = 0, blocking = -1;

		if (simplex_pan_basis_least_squares(basis, basis->standard->rhs,
			    trial, &residual_norm) == lp_simplex_EXIT_FAILURE)
			return lp_simplex_EXIT_FAILURE;
		for (i = 0; i < basis->count; i++) {
			int column = basis->column[i];

			if (!isfinite(trial[i]))
				return lp_simplex_EXIT_FAILURE;
			if ((free_mask == NULL || !free_mask[column]) &&
				trial[i] <= 0.) {
				double denominator = value[column] - trial[i];
				double candidate =
					denominator > 0.
						? value[column] / denominator
						: 0.;

				has_nonpositive = 1;
				if (candidate < alpha ||
					(candidate == alpha && blocking < 0)) {
					alpha = candidate;
					blocking = i;
				}
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
			if ((free_mask == NULL || !free_mask[column]) &&
				(i == blocking || value[column] <= 0.)) {
				value[column] = 0.;
				remove[i] = 1;
			}
		}
		simplex_pan_basis_remove_marked(basis, remove);
		if (simplex_pan_basis_factorize(basis) ==
			lp_simplex_EXIT_FAILURE)
			return lp_simplex_EXIT_FAILURE;
	}
}


int simplex_pan_phase1_mixed_audited(struct simplex_PanBasis *basis,
	double *value, const unsigned char *free_mask, double *residual_out,
	const double primal_tolerance, const double dual_tolerance,
	const int iteration_limit, int *iterations, int *status,
	int (*accept)(void *, const double *), void *context)
{
	const struct simplex_PanStandard *standard = basis->standard;
	unsigned char *blocked = NULL;
	unsigned char *remove = NULL;
	double *residual = NULL;
	double *trial = NULL;
	double *correction = NULL;
	double residual_norm;
	int entering;
	int block_limit = 1, block_size = 1;
	int rounds = 0, additions = 0, deletions = 0;
	int state = lp_simplex_EXIT_FAILURE;

	if (free_mask != NULL && basis->nnls_block_limit > 1) {
		block_limit = basis->nnls_block_limit;
		if (block_limit < 1)
			block_limit = 1;
		if (block_limit > 32)
			block_limit = 32;
		block_size = block_limit < 4 ? block_limit : 4;
	}
	blocked = (unsigned char *)lp_simplex_malloc(
		(size_t)(standard->columns > 0 ? standard->columns : 1));
	remove = (unsigned char *)lp_simplex_malloc((size_t)standard->rows);
	residual = (double *)lp_simplex_malloc(
		(size_t)standard->rows * sizeof(double));
	correction = (double *)lp_simplex_malloc(
		(size_t)standard->rows * sizeof(double));
	trial = (double *)lp_simplex_malloc(
		(size_t)standard->rows * sizeof(double));
	if (blocked == NULL || remove == NULL || residual == NULL ||
		trial == NULL || correction == NULL) {
		*status = lp_simplex_MemoryAllocError;
		goto cleanup;
	}
	lp_simplex_memset(blocked, 0, (size_t)standard->columns);
	lp_simplex_memset(value, 0, (size_t)standard->columns * sizeof(double));
	/* Optional batch crash: one rank-revealing QR selects an independent
	 * subset of the full matrix. Positivity restoration then removes columns
	 * along the usual NNLS segment, before ordinary correlation pricing.
	 * It is intentionally separate from mixed direction subproblems. */
	if (free_mask == NULL && basis->count == 0 && basis->block_crash &&
		simplex_pan_basis_complete(basis, value) ==
			lp_simplex_EXIT_SUCCESS) {
		if (basis->count > iteration_limit - *iterations) {
			*status = lp_simplex_ExceedIterLimit;
			goto cleanup;
		}
		*iterations += basis->count;
		if (pan_phase1_positive_solution(
			    basis, value, trial, remove, NULL)) {
			*status = lp_simplex_PrecisionError;
			goto cleanup;
		}
		if (getenv("LP_SIMPLEX_PROFILE"))
			fprintf(stderr,
				"pan-bda profile: block crash rank=%d/%d\n",
				basis->count, standard->rows);
	}
	/* Disjoint singleton columns are orthogonal and have their exact NNLS
	 * minimizer rhs[row]/a without an active-set iteration.  In particular,
	 * upper-bound slack rows need not consume thousands of Phase-I pivots.
	 * This only seeds an empty passive set; ordinary residual certification
	 * and the active-set algorithm still handle all coupled rows. */
	if (basis->count == 0 && basis->singleton_crash) {
		int j;

		lp_simplex_memset(remove, 0, (size_t)standard->rows);
		for (j = 0; j < standard->columns; j++) {
			int k = standard->column_start[j];
			int row;
			double candidate;

			if (standard->column_start[j + 1] - k != 1 ||
				standard->value[k] == 0.)
				continue;
			row = standard->row_index[k];
			if (remove[row])
				continue;
			candidate = standard->rhs[row] / standard->value[k];
			if (!isfinite(candidate) || candidate <= 0.)
				continue;
			if (simplex_pan_basis_add(basis, j) ==
				lp_simplex_EXIT_FAILURE) {
				*status = lp_simplex_MemoryAllocError;
				goto cleanup;
			}
			value[j] = candidate;
			remove[row] = 1;
		}
		if (simplex_pan_basis_rebuild(basis) ==
			lp_simplex_EXIT_FAILURE) {
			*status = lp_simplex_PrecisionError;
			goto cleanup;
		}
		if (getenv("LP_SIMPLEX_PROFILE") != NULL)
			fprintf(stderr,
				"pan-bda profile: singleton crash rank=%d/%d\n",
				basis->count, standard->rows);
	}
	if (free_mask != NULL && basis->count > 0 &&
		pan_phase1_positive_solution(basis, value, trial, remove,
			free_mask) == lp_simplex_EXIT_FAILURE) {
		*status = lp_simplex_PrecisionError;
		goto cleanup;
	}
	for (;;) {
		if (*iterations >= iteration_limit) {
			*status = lp_simplex_ExceedIterLimit;
			goto cleanup;
		}
		pan_phase1_residual(basis, value, residual, correction);
		residual_norm = pan_phase1_norm(residual, standard->rows);
		if (!isfinite(residual_norm)) {
			*status = lp_simplex_PrecisionError;
			goto cleanup;
		}
		{
			int i, feasible = 1;

			for (i = 0; i < standard->rows; i++)
				if (!isfinite(residual[i]) ||
					fabs(residual[i]) > primal_tolerance)
					feasible = 0;
			if (feasible) {
				*status = lp_simplex_CondUnsatisfied;
				state = lp_simplex_EXIT_SUCCESS;
				goto cleanup;
			}
		}
		if (free_mask != NULL && accept && accept(context, value)) {
			*status = lp_simplex_CondUnsatisfied;
			state = 2;
			goto cleanup;
		}
		if (*iterations >= iteration_limit) {
			*status = lp_simplex_ExceedIterLimit;
			goto cleanup;
		}
		entering = pan_phase1_entering(basis, residual, blocked,
			fmin(dual_tolerance * (1. + residual_norm),
				primal_tolerance /
					(1. + sqrt((double)standard->rows))));
		if (entering < 0) {
			int j, certificate = 1;

			if (free_mask != NULL) {
				*status = lp_simplex_CondUnsatisfied;
				state = 1;
				goto cleanup;
			}
			/* An exhausted/blocked numerical search is not infeasibility.
			 * Check every column, including temporarily blocked candidates. */
			for (j = 0; j < standard->columns; j++)
				if (pan_phase1_column_dot(
					    standard, j, residual) > 0.)
					certificate = 0;
			{
				long double product = 0., correction = 0.,
					    absolute = 0.;
				int i;

				for (i = 0; i < standard->rows; i++) {
					double p =
						standard->rhs[i] * residual[i];
					long double next = product + p;

					correction +=
						(fabsl(product) >= fabs(p)
								? (product -
									  next) +
									  p
								: (p - next) +
									  product) +
						fma(standard->rhs[i],
							residual[i], -p);
					absolute += fabs(p);
					product = next;
				}
				product += correction;
				if (!isfinite(product) ||
					product <= 64. * DBL_EPSILON * absolute)
					certificate = 0;
			}
			*status = certificate ? lp_simplex_Infeasibility
					      : lp_simplex_PrecisionError;
			goto cleanup;
		}
		if (block_limit > 1) {
			int candidate[32], count, p, inserted = 0,
						     previous = basis->count;
			int cap = iteration_limit - *iterations;

			if (cap > block_size)
				cap = block_size;
			count = simplex_pan_nnls_shortlist(basis, residual,
				blocked,
				fmin(dual_tolerance * (1. + residual_norm),
					primal_tolerance /
						(1. + sqrt((double)standard
								      ->rows))),
				cap, candidate, 0);
			for (p = 0; p < count && basis->count < standard->rows;
				p++) {
				int column = candidate[p];

				if (simplex_pan_basis_add(basis, column)) {
					*status = lp_simplex_MemoryAllocError;
					goto cleanup;
				}
				(*iterations)++;
				if (simplex_pan_basis_factorize(basis)) {
					int q;

					for (q = 0; q < basis->count; q++)
						remove[q] = basis->column[q] ==
							    column;
					simplex_pan_basis_remove_marked(
						basis, remove);
					if (simplex_pan_basis_rebuild(basis)) {
						*status =
							lp_simplex_PrecisionError;
						goto cleanup;
					}
					blocked[column] = 1;
				} else
					inserted++;
			}
			if (inserted == 0) {
				/* A full passive set with a nonzero residual is numerical
				 * failure, not a stationary NNLS or infeasibility witness. */
				if (basis->count == standard->rows) {
					*status = lp_simplex_PrecisionError;
					goto cleanup;
				}
				continue;
			}
			if (pan_phase1_positive_solution(
				    basis, value, trial, remove, free_mask)) {
				*status = lp_simplex_PrecisionError;
				goto cleanup;
			}
			rounds++;
			additions += inserted;
			{
				int dropped =
					previous + inserted - basis->count;

				deletions += dropped;
				if (dropped > 0) {
					lp_simplex_memset(blocked, 0,
						(size_t)standard->columns);
					/* A column rejected jointly may be improving alone. Only
					 * a single-column restoration may temporarily block it. */
					if (count == 1 &&
						basis->position[candidate[0]] <
							0)
						blocked[candidate[0]] = 1;
				}
				if (dropped * 2 >= inserted)
					block_size =
						block_size > 1
							? (block_size + 1) / 2
							: 1;
				else if (dropped == 0 &&
					 block_size < block_limit) {
					block_size *= 2;
					if (block_size > block_limit)
						block_size = block_limit;
				}
			}
			continue;
		}
		if (simplex_pan_basis_add(basis, entering) ==
			lp_simplex_EXIT_FAILURE) {
			*status = lp_simplex_MemoryAllocError;
			goto cleanup;
		}
		if (simplex_pan_basis_factorize(basis) ==
			lp_simplex_EXIT_FAILURE) {
			lp_simplex_memset(remove, 0, (size_t)basis->count);
			remove[basis->count - 1] = 1;
			simplex_pan_basis_remove_marked(basis, remove);
			blocked[entering] = 1;
			continue;
		}
		(*iterations)++;
		if (getenv("LP_SIMPLEX_PROFILE") && (*iterations % 1024) == 0)
			fprintf(stderr,
				"pan-bda nnls: iterations=%d rank=%d/%d residual=%.6g\n",
				*iterations, basis->count, standard->rows,
				residual_norm);
		{
			int count_before_positive_solve = basis->count;

			if (pan_phase1_positive_solution(basis, value, trial,
				    remove,
				    free_mask) == lp_simplex_EXIT_FAILURE) {
				*status = lp_simplex_PrecisionError;
				goto cleanup;
			}
			/* Rank failures are relative to the current passive space.  A
			 * deletion can make an earlier column independent again, so only
			 * retain the block on a column that immediately left once more. */
			if (basis->count < count_before_positive_solve) {
				lp_simplex_memset(
					blocked, 0, (size_t)standard->columns);
				if (basis->position[entering] < 0)
					blocked[entering] = 1;
			}
		}
	}
cleanup:
	if (getenv("LP_SIMPLEX_PROFILE") && additions >= 128)
		fprintf(stderr,
			"pan-bda nnls blocks: rounds=%d added=%d deleted=%d rank=%d limit=%d\n",
			rounds, additions, deletions, basis->count,
			block_limit);
	if (state >= 0 && residual_out != NULL) {
		pan_phase1_residual(basis, value, residual, correction);
		lp_simplex_memcpy(residual_out, residual,
			(size_t)standard->rows * sizeof(double));
	}
	lp_simplex_free(blocked);
	lp_simplex_free(remove);
	lp_simplex_free(residual);
	lp_simplex_free(trial);
	lp_simplex_free(correction);
	return state;
}


int simplex_pan_phase1_mixed(struct simplex_PanBasis *basis, double *value,
	const unsigned char *free_mask, double *residual_out,
	double primal_tolerance, double dual_tolerance, int iteration_limit,
	int *iterations, int *status)
{
	return simplex_pan_phase1_mixed_audited(basis, value, free_mask,
		residual_out, primal_tolerance, dual_tolerance, iteration_limit,
		iterations, status, NULL, NULL);
}


int simplex_pan_phase1_run(struct simplex_PanBasis *basis, double *value,
	const double primal_tolerance, const double dual_tolerance,
	const int iteration_limit, int *iterations, int *status)
{
	return simplex_pan_phase1_mixed(basis, value, NULL, NULL,
		primal_tolerance, dual_tolerance, iteration_limit, iterations,
		status);
}
