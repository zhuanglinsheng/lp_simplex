/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_pan_cone_implicit.c
 * @brief Mixed NNLS for [A;c^T]d=[0;-1], without factoring the objective row.
 *
 * While A_B is independent, v solves B^T v=c_B with minimum norm and
 * d_B=-B^+v/(1+v^Tv). A dependent entering column either supplies an exact
 * augmented fit or is rejected. All residuals are recomputed from CSC.
 */

#include "policy/pan.h"
#include "policy/controls.h"
#include "simplex_pan_cone.h"
#include "simplex_stage.h"
#include "simplex_pan_basis.h"
#include "utils.h"

#include <lp_simplex/status.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define CONE_PRECISION \
	do { \
		failure_line = __LINE__; \
		goto precision; \
	} while (0)

static double cone_cost(const struct simplex_PanStandard *a, const double *d)
{
	double sum = 0., correction = 0.;
	int j;

	for (j = 0; j < a->columns; j++) {
		double product = a->objective[j] * d[j], next = sum + product;

		correction +=
			(fabs(sum) >= fabs(product) ? (sum - next) + product
						    : (product - next) + sum) +
			fma(a->objective[j], d[j], -product);
		sum = next;
	}
	return sum + correction;
}

static void cone_residual(const struct simplex_PanStandard *a, const double *d,
	double *r, double *work)
{
	simplex_pan_residual(a, d, r, work);

	r[a->rows] = -1. - cone_cost(a, d);
}

static int cone_update_factor(struct simplex_PanBasis *basis)
{
	if (basis->prefer_updates &&
		simplex_pan_basis_factorize(basis) == lp_simplex_EXIT_SUCCESS)
		return 0;
	return simplex_pan_basis_rebuild(basis) ? -1 : 0;
}

/* Accept a trial or move to its first constrained zero. Free coordinates
 * are never removed. The entering column may be outside the current basis. */
static int cone_positive(struct simplex_PanBasis *b, double *d,
	const double *trial, const unsigned char *free_mask,
	unsigned char *remove, int entering, int *bulk_available)
{
	double alpha = 1.;
	int i, j, blocking = -1, positive_blockers = 0;

	for (j = 0; j < b->standard->columns; j++) {
		if (!isfinite(trial[j]))
			return -1;
		if (!free_mask[j] && trial[j] < 0.) {
			double ratio;

			if (d[j] > 0.)
				positive_blockers++;
			ratio = d[j] / (d[j] - trial[j]);
			if (ratio < alpha) {
				alpha = ratio;
				blocking = j;
			}
		}
	}
	/* One bounded working-set crash per oracle, never repeated bulk pivots:
	 * remove a large group of currently positive constrained coefficients
	 * whose unconstrained minimizer is negative. The remaining d is still
	 * NNLS-feasible; ordinary Lawson-Hanson restoration resumes afterwards.
	 * This changes no primal LP point and cannot itself return a direction. */
	if (bulk_available && *bulk_available && positive_blockers >= 8) {
		int dropped = 0;
		*bulk_available = 0;
		for (i = 0; i < b->count; i++) {
			j = b->column[i];
			remove[i] = !free_mask[j] && trial[j] <= 0.;
			if (remove[i]) {
				d[j] = 0.;
				dropped++;
			}
		}
		simplex_pan_basis_remove_marked(b, remove);
		if (getenv(SIMPLEX_ENV_PROFILE))
			fprintf(stderr,
				"pan-bda implicit bulk: removed=%d rank=%d\n",
				dropped, b->count);
		return cone_update_factor(b);
	}
	if (blocking < 0) {
		lp_simplex_memcpy(d, trial,
			(size_t)b->standard->columns * sizeof(double));
		return 1;
	}
	for (j = 0; j < b->standard->columns; j++)
		d[j] = fma(alpha, trial[j] - d[j], d[j]);
	d[blocking] = 0.;
	for (i = 0; i < b->count; i++) {
		j = b->column[i];
		remove[i] = !free_mask[j] && d[j] <= 0.;
		if (remove[i])
			d[j] = 0.;
	}
	simplex_pan_basis_remove_marked(b, remove);
	if (entering >= 0 && d[entering] > 0. && b->position[entering] < 0 &&
		simplex_pan_basis_add(b, entering))
		return -1;
	return cone_update_factor(b);
}

int simplex_pan_cone_implicit_audited(struct simplex_PanBasis *basis,
	double *direction, const unsigned char *free_mask, double *residual,
	double tolerance, int limit, int *iterations, int *status,
	int (*accept)(void *, const double *), void *context)
{
	const struct simplex_PanStandard *a = basis->standard;
	double *storage = NULL, *trial, *v, *work, *coefficient, *right;
	unsigned char *remove = NULL, *blocked = NULL;
	int state = -1, i, j, entering, failure_line = 0, rank_rejections = 0;
	int bulk_available = basis->implicit_bulk;
	unsigned long history_first[SIMPLEX_PAN_IMPLICIT_HISTORY],
		history_second[SIMPLEX_PAN_IMPLICIT_HISTORY];
	int history_rank[SIMPLEX_PAN_IMPLICIT_HISTORY], history_count = 0,
							history_next = 0;

	basis->oracle_failure = SIMPLEX_FAILURE_UNKNOWN;
	storage = (double *)lp_simplex_malloc(
		((size_t)a->columns + 4 * (size_t)(a->rows + 1)) *
		sizeof(double));
	remove = (unsigned char *)lp_simplex_malloc((size_t)(a->rows + 1));
	blocked = (unsigned char *)lp_simplex_malloc((size_t)a->columns);
	if (!storage || !remove || !blocked) {
		*status = lp_simplex_MemoryAllocError;
		goto done;
	}
	trial = storage;
	v = trial + a->columns;
	work = v + a->rows + 1;
	coefficient = work + a->rows + 1;
	right = coefficient + a->rows + 1;
	lp_simplex_memset(direction, 0, (size_t)a->columns * sizeof(double));
	lp_simplex_memset(blocked, 0, (size_t)a->columns);
	for (;;) {
		double norm = 0., denominator, error, best;
		int accepted;

		if (*iterations >= limit) {
			*status = lp_simplex_ExceedIterLimit;
			goto done;
		}
		for (i = 0; i < basis->count; i++)
			right[i] = a->objective[basis->column[i]];
		if (simplex_pan_basis_minimum_norm(basis, right, v, &error))
			CONE_PRECISION;
		for (i = 0; i < a->rows; i++)
			norm = hypot(norm, v[i]);
		denominator = 1. + norm * norm;
		if (!isfinite(denominator))
			CONE_PRECISION;
		if (simplex_pan_basis_least_squares(
			    basis, v, coefficient, &error))
			CONE_PRECISION;
		lp_simplex_memset(
			trial, 0, (size_t)a->columns * sizeof(double));
		for (i = 0; i < basis->count; i++)
			trial[basis->column[i]] = -coefficient[i] / denominator;
		{
			int before = bulk_available;

			accepted = cone_positive(basis, direction, trial,
				free_mask, remove, -1, &bulk_available);
			if (before && !bulk_available)
				(*iterations)++;
		}
		if (accepted < 0)
			CONE_PRECISION;
		if (!accepted) {
			lp_simplex_memset(blocked, 0, (size_t)a->columns);
			continue;
		}
		cone_residual(a, direction, residual, work);
		error = fabs(residual[a->rows]);
		for (i = 0; i < a->rows; i++)
			error = fmax(error, fabs(residual[i]));
		if (error <= tolerance) {
			state = 0;
			goto done;
		}
		if (accept && accept(context, direction)) {
			state = 2;
			goto done;
		}
		if (*iterations >= limit) {
			*status = lp_simplex_ExceedIterLimit;
			goto done;
		}
		{
			unsigned long first = 0UL, second = 0UL;
			int h;
			/* Restricted least squares is unique on an independent passive
			 * set. Repeating the same membership and blocked candidates after
			 * restoration is a numerical cycle, not useful NNLS progress.
			 * Ignore internal QR order; include blocked candidates so ordinary
			 * rank rejection does not look like a cycle. */
			for (j = 0; j < a->columns; j++) {
				unsigned long key = (unsigned long)j + 1UL;

				if (basis->position[j] >= 0) {
					first ^= key * 2654435761UL;
					second += key * key * 2246822519UL;
				}
				if (blocked[j]) {
					first ^= key * 3266489917UL;
					second += key * key * 668265263UL;
				}
			}
			for (h = 0; h < history_count; h++)
				if (history_rank[h] == basis->count &&
					history_first[h] == first &&
					history_second[h] == second) {
					if (getenv(SIMPLEX_ENV_PROFILE))
						fprintf(stderr,
							"pan-bda implicit cycle: rank=%d iterations=%d -> explicit recovery\n",
							basis->count,
							*iterations);
					basis->oracle_failure =
						SIMPLEX_FAILURE_IMPLICIT_ORACLE;
					CONE_PRECISION;
				}
			history_first[history_next] = first;
			history_second[history_next] = second;
			history_rank[history_next] = basis->count;
			history_next = (history_next + 1) %
				       SIMPLEX_PAN_IMPLICIT_HISTORY;
			if (history_count < SIMPLEX_PAN_IMPLICIT_HISTORY)
				history_count++;
		}
		best = tolerance / (1. + sqrt((double)a->rows + 1.));
		entering = -1;
		for (j = 0; j < a->columns; j++)
			if (basis->position[j] < 0 && !blocked[j]) {
				double length = hypot(
					a->column_norm[j], a->objective[j]);
				double score =
					length > 0.
						? (simplex_pan_column_dot(
							   a, j, residual) +
							  a->objective[j] *
								  residual[a->rows]) /
							  length
						: 0.;

				if (score > best) {
					best = score;
					entering = j;
				}
			}
		if (entering < 0) {
			state = 1;
			goto done;
		}
		if (simplex_pan_basis_column_projection(
			    basis, entering, coefficient, &error))
			CONE_PRECISION;
		(*iterations)++;
		if (getenv(SIMPLEX_ENV_PROFILE) && (*iterations % 1024) == 0)
			fprintf(stderr,
				"pan-bda implicit nnls: iterations=%d rank=%d/%d residual=%.6g\n",
				*iterations, basis->count, a->rows, error);
		if (error <= basis->rank_tolerance * a->column_norm[entering]) {
			double gamma = a->objective[entering], sum = 0.,
			       compensation = 0.;
			/* Numerically zero projection entries must not select a fictitious
			 * blocker and replace an independent column by a dependent one.
			 * Budget their total row activity, then audit the resulting ray. */
			for (i = 0; i < basis->count; i++) {
				j = basis->column[i];
				if (!free_mask[j] &&
					fabs(coefficient[i]) *
							a->column_norm[j] <=
						tolerance *
							SIMPLEX_PAN_INNER_RETRY_FACTOR /
							(1. + basis->count))
					coefficient[i] = 0.;
			}
			for (i = 0; i < basis->count; i++) {
				double cost = a->objective[basis->column[i]],
				       p = cost * coefficient[i],
				       next = sum + p;

				compensation +=
					(fabs(sum) >= fabs(p)
							? (sum - next) + p
							: (p - next) + sum) +
					fma(cost, coefficient[i], -p);
				sum = next;
			}
			gamma -= sum + compensation;
			if (!isfinite(gamma) || gamma >= 0.) {
				blocked[entering] = 1;
				continue;
			}
			lp_simplex_memset(
				trial, 0, (size_t)a->columns * sizeof(double));
			for (i = 0; i < basis->count; i++)
				trial[basis->column[i]] =
					coefficient[i] / gamma;
			trial[entering] = -1. / gamma;
			accepted = cone_positive(basis, direction, trial,
				free_mask, remove, entering, NULL);
			if (accepted < 0)
				CONE_PRECISION;
			if (accepted) {
				cone_residual(a, direction, residual, work);
				error = fabs(residual[a->rows]);
				for (i = 0; i < a->rows; i++)
					error = fmax(error, fabs(residual[i]));
				if (error > tolerance) {
					if (accept &&
						accept(context, direction)) {
						state = 2;
						goto done;
					}
					CONE_PRECISION;
				}
				state = 0;
				goto done;
			}
			lp_simplex_memset(blocked, 0, (size_t)a->columns);
		} else {
			if (simplex_pan_basis_add(basis, entering)) {
				*status = lp_simplex_MemoryAllocError;
				goto done;
			}
			if (simplex_pan_basis_factorize(basis)) {
				for (i = 0; i < basis->count; i++)
					remove[i] =
						basis->column[i] == entering;
				simplex_pan_basis_remove_marked(basis, remove);
				if (simplex_pan_basis_rebuild(basis))
					CONE_PRECISION;
				blocked[entering] = 1;
				/* A large Q-tail prediction followed by repeated QR rank
				 * rejection is a numerical disagreement, not a reason to scan
				 * thousands more columns on the same passive space. */
				if (++rank_rejections >= 4) {
					if (getenv(SIMPLEX_ENV_PROFILE))
						fprintf(stderr,
							"pan-bda implicit rank disagreement: rejects=%d -> explicit recovery\n",
							rank_rejections);
					basis->oracle_failure =
						SIMPLEX_FAILURE_IMPLICIT_ORACLE;
					CONE_PRECISION;
				}
			} else
				rank_rejections = 0;
		}
	}
precision:
	if (getenv(SIMPLEX_ENV_PROFILE))
		fprintf(stderr,
			"pan-bda implicit failure: line=%d rank=%d backward=%.3g projection=%.3g\n",
			failure_line, basis->count, basis->backward_error,
			basis->projection_solve_error);
	*status = lp_simplex_PrecisionError;
done:
	if (a->report && *status == lp_simplex_PrecisionError)
		a->report->failure = basis->oracle_failure;
	lp_simplex_free(storage);
	lp_simplex_free(remove);
	lp_simplex_free(blocked);
	return state;
}
