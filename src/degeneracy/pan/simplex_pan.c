/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 *
 * Pan's basis-deficiency-allowing generalized simplex method.  In the notation
 * of the original algorithm, the standard-form LP is the dual problem
 *
 *     max b'z  subject to Az=c, z>=0,
 *
 * with b=-objective and c=rhs.  Its complementary inequality problem supplies
 * the minimum-norm multiplier used for pricing.  The number of columns in the
 * basis changes dynamically and is never padded to the row dimension.
 */
#include "simplex_pan.h"
#include "simplex_pan_basis.h"
#include "simplex_pan_phase1.h"
#include "simplex_pan_standard.h"
#include "utils.h"

#include <lp_simplex/status.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#ifdef LP_SIMPLEX_HAVE_SPQR
#define PAN_BASIS_BACKEND_NAME "dynamic-sparse-qr+spqr"
#else
#define PAN_BASIS_BACKEND_NAME "dense-fallback"
#endif

struct pan_Workspace {
	double *value;
	double *multiplier;
	double *basis_right;
	double *basis_value;
	double *direction;
	unsigned char *remove;
};


static void pan_profile_failure(const char *stage)
{
	if (getenv("LP_SIMPLEX_PROFILE") != NULL)
		fprintf(stderr, "pan-bda precision failure: %s\n", stage);
}


static double pan_column_dot(
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


static int pan_workspace_create(
		struct pan_Workspace *work,
		const struct simplex_PanStandard *standard)
{
	lp_simplex_memset(work, 0, sizeof(*work));
	work->value = (double *)lp_simplex_malloc(
		(size_t)(standard->columns > 0 ? standard->columns : 1) *
		sizeof(double));
	work->multiplier = (double *)lp_simplex_malloc(
		(size_t)standard->rows * sizeof(double));
	work->basis_right = (double *)lp_simplex_malloc(
		(size_t)standard->rows * sizeof(double));
	work->basis_value = (double *)lp_simplex_malloc(
		(size_t)standard->rows * sizeof(double));
	work->direction = (double *)lp_simplex_malloc(
		(size_t)standard->rows * sizeof(double));
	work->remove = (unsigned char *)lp_simplex_malloc(
		(size_t)standard->rows);
	if (work->value == NULL || work->multiplier == NULL ||
	    work->basis_right == NULL || work->basis_value == NULL ||
	    work->direction == NULL || work->remove == NULL)
		return lp_simplex_EXIT_FAILURE;
	return lp_simplex_EXIT_SUCCESS;
}


static void pan_workspace_destroy(struct pan_Workspace *work)
{
	lp_simplex_free(work->value);
	lp_simplex_free(work->multiplier);
	lp_simplex_free(work->basis_right);
	lp_simplex_free(work->basis_value);
	lp_simplex_free(work->direction);
	lp_simplex_free(work->remove);
	lp_simplex_memset(work, 0, sizeof(*work));
}


static int pan_choose_violated(
		const struct simplex_PanBasis *basis,
		const double *multiplier, const double tolerance,
		const int pricing,
		double *maximum_infeasibility)
{
	const struct simplex_PanStandard *standard = basis->standard;
	double best_score = 0.;
	int entering = -1;
	int j;
	*maximum_infeasibility = 0.;
	for (j = 0; j < standard->columns; j++) {
		double residual, violation, norm, score;
		if (basis->position[j] >= 0)
			continue;
		residual = pan_column_dot(standard, j, multiplier) +
			standard->objective[j];
		violation = residual < 0. ? -residual : 0.;
		if (violation > *maximum_infeasibility)
			*maximum_infeasibility = violation;
		if (violation <= tolerance)
			continue;
		norm = simplex_pan_column_norm(standard, j);
		if (norm == 0.)
			return -2;
		score = pricing == lp_simplex_PRICING_PAN_NORMALIZED
			? violation / norm : violation;
		if (score > best_score ||
		    (score == best_score && (entering < 0 || j < entering))) {
			best_score = score;
			entering = j;
		}
	}
	return entering;
}


static double pan_standard_residual(
		const struct simplex_PanStandard *standard,
		const double *value, double *residual)
{
	double square = 0.;
	int i, j, k;
	lp_simplex_memcpy(residual, standard->rhs,
		(size_t)standard->rows * sizeof(double));
	for (j = 0; j < standard->columns; j++) {
		if (value[j] == 0.)
			continue;
		for (k = standard->column_start[j];
		     k < standard->column_start[j + 1]; k++)
			residual[standard->row_index[k]] -=
				value[j] * standard->value[k];
	}
	for (i = 0; i < standard->rows; i++)
		square += residual[i] * residual[i];
	return sqrt(square);
}


/* Certify the maintained feasible point against the original sparse columns.
 * A rank-revealing orthogonal solve is used only as a final recovery path; it
 * is deliberately not the ordinary update, because Pan's exchange itself
 * preserves feasibility and can leave a non-unique solution on a deficient
 * basis. */
static int pan_certify_feasible_values(
		struct simplex_PanBasis *basis, struct pan_Workspace *work,
		const double tolerance)
{
	double residual_norm;
	double rhs_norm = 0.;
	double feasibility_tolerance;
	int i, needs_recovery = 0;
	for (i = 0; i < basis->standard->rows; i++)
		rhs_norm += basis->standard->rhs[i] * basis->standard->rhs[i];
	rhs_norm = sqrt(rhs_norm);
	feasibility_tolerance = tolerance * (1. + rhs_norm);
	residual_norm = pan_standard_residual(basis->standard, work->value,
		work->multiplier);
	if (residual_norm > feasibility_tolerance)
		needs_recovery = 1;
	for (i = 0; i < basis->standard->columns; i++)
		if (work->value[i] < -feasibility_tolerance)
			needs_recovery = 1;
	if (!needs_recovery)
		return lp_simplex_EXIT_SUCCESS;
	if (simplex_pan_basis_least_squares_qr(basis, basis->standard->rhs,
		work->basis_value, &residual_norm) == lp_simplex_EXIT_FAILURE) {
		pan_profile_failure("feasibility QR recovery");
		return lp_simplex_EXIT_FAILURE;
	}
	if (residual_norm > feasibility_tolerance) {
		if (getenv("LP_SIMPLEX_PROFILE") != NULL)
			fprintf(stderr, "pan-bda feasibility residual %.6g > %.6g\n",
				residual_norm, feasibility_tolerance);
		return lp_simplex_EXIT_FAILURE;
	}
	lp_simplex_memset(work->value, 0,
		(size_t)basis->standard->columns * sizeof(double));
	for (i = 0; i < basis->count; i++) {
		double value = work->basis_value[i];
		if (value < -feasibility_tolerance) {
			if (getenv("LP_SIMPLEX_PROFILE") != NULL)
				fprintf(stderr,
					"pan-bda negative refreshed value %.6g < %.6g\n",
					value, -feasibility_tolerance);
			return lp_simplex_EXIT_FAILURE;
		}
		work->value[basis->column[i]] = value > 0. ? value : 0.;
	}
	residual_norm = pan_standard_residual(basis->standard, work->value,
		work->multiplier);
	if (residual_norm > feasibility_tolerance)
		return lp_simplex_EXIT_FAILURE;
	return lp_simplex_EXIT_SUCCESS;
}


static int pan_add_independent(
		struct simplex_PanBasis *basis, struct pan_Workspace *work,
		const int entering, const double tolerance)
{
	int old_count = basis->count;
	(void)tolerance;
	if (simplex_pan_basis_add(basis, entering) == lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	if (simplex_pan_basis_factorize(basis) == lp_simplex_EXIT_FAILURE) {
		lp_simplex_memset(work->remove, 0, (size_t)basis->count);
		work->remove[old_count] = 1;
		simplex_pan_basis_remove_marked(basis, work->remove);
		if (simplex_pan_basis_factorize(basis) == lp_simplex_EXIT_FAILURE)
			return lp_simplex_EXIT_FAILURE;
		/* The projection classified a numerically dependent column as
		 * independent.  Ask the caller to use the exchange path instead. */
		return 1;
	}
	work->value[entering] = 0.;
	return 0;
}


static int pan_exchange(
		struct simplex_PanBasis *basis, struct pan_Workspace *work,
		const int entering, const double primal_tolerance,
		const double pivot_tolerance, int *unbounded)
{
	double minimum = 0.;
	int bounded = 0;
	int i, leaving = -1, remove_count = 0;
	*unbounded = 0;
	for (i = 0; i < basis->count; i++) {
		if (work->direction[i] > pivot_tolerance) {
			double ratio = work->value[basis->column[i]] /
				work->direction[i];
			if (!bounded || ratio < minimum) {
				minimum = ratio;
				leaving = i;
				bounded = 1;
			}
		}
	}
	if (!bounded) {
		*unbounded = 1;
		return lp_simplex_EXIT_SUCCESS;
	}
	lp_simplex_memset(work->remove, 0, (size_t)basis->count);
	for (i = 0; i < basis->count; i++) {
		double new_value = work->value[basis->column[i]] -
			minimum * work->direction[i];
		if (work->direction[i] > pivot_tolerance) {
			double ratio = work->value[basis->column[i]] /
				work->direction[i];
			if (__lp_simplex_ABS__(ratio - minimum) <=
			    1e-12 * (1. + __lp_simplex_ABS__(minimum)))
				work->remove[i] = 1;
		}
		if (work->remove[i])
			work->value[basis->column[i]] = 0.;
		else if (new_value < -primal_tolerance) {
			pan_profile_failure("negative value after exchange");
			return lp_simplex_EXIT_FAILURE;
		} else
			work->value[basis->column[i]] = new_value > 0. ? new_value : 0.;
	}
	/* The exact minimum-ratio row must leave even if its recomputed ratio is
	 * separated by a few ulps from the comparison above. */
	if (leaving >= 0) {
		work->remove[leaving] = 1;
		work->value[basis->column[leaving]] = 0.;
	}
	work->value[entering] = minimum;
	for (i = 0; i < basis->count; i++)
		remove_count += work->remove[i] != 0;
	if (remove_count == 1) {
		if (simplex_pan_basis_replace_factorized(basis, leaving, entering) ==
		    lp_simplex_EXIT_FAILURE) {
			pan_profile_failure("exchange factor update");
			return lp_simplex_EXIT_FAILURE;
		}
		return lp_simplex_EXIT_SUCCESS;
	}
	simplex_pan_basis_remove_marked(basis, work->remove);
	if (simplex_pan_basis_add(basis, entering) == lp_simplex_EXIT_FAILURE) {
		pan_profile_failure("exchange basis insertion");
		return lp_simplex_EXIT_FAILURE;
	}
	if (simplex_pan_basis_factorize(basis) == lp_simplex_EXIT_FAILURE) {
		pan_profile_failure("exchange factorization");
		return lp_simplex_EXIT_FAILURE;
	}
	return lp_simplex_EXIT_SUCCESS;
}


static int pan_phase2(
		struct simplex_PanBasis *basis, struct pan_Workspace *work,
		const struct lp_simplex_Options *options,
		int *iterations, int *status, double *dual_infeasibility)
{
	double residual_norm;
	for (;;) {
		double column_norm, projection_residual;
		int entering, i, unbounded;
		if (*iterations >= options->iteration_limit) {
			*status = lp_simplex_ExceedIterLimit;
			return lp_simplex_EXIT_FAILURE;
		}
		for (i = 0; i < basis->count; i++)
			work->basis_right[i] =
				-basis->standard->objective[basis->column[i]];
		if (simplex_pan_basis_minimum_norm(basis, work->basis_right,
			work->multiplier, &residual_norm) == lp_simplex_EXIT_FAILURE ||
		    residual_norm > options->dual_tolerance *
			(1. + (double)basis->count)) {
			/* Dynamic factor modification is accepted only when the original
			 * basis equations certify it.  Rebuild from the sparse columns on
			 * failure; this is residual-driven, not an update-count heuristic. */
			if (simplex_pan_basis_factorize(basis) == lp_simplex_EXIT_FAILURE ||
			    simplex_pan_basis_minimum_norm(basis, work->basis_right,
				work->multiplier, &residual_norm) == lp_simplex_EXIT_FAILURE ||
			    residual_norm > options->dual_tolerance *
				(1. + (double)basis->count)) {
				pan_profile_failure("minimum-norm multiplier");
				*status = lp_simplex_PrecisionError;
				return lp_simplex_EXIT_FAILURE;
			}
		}
		entering = pan_choose_violated(basis, work->multiplier,
			options->dual_tolerance, options->pricing,
			dual_infeasibility);
		if (entering == -2) {
			*status = lp_simplex_Unboundedness;
			return lp_simplex_EXIT_FAILURE;
		}
		if (entering < 0) {
			*status = lp_simplex_Success;
			return lp_simplex_EXIT_SUCCESS;
		}
		if (simplex_pan_basis_column_projection(basis, entering,
			work->direction, &projection_residual) ==
		    lp_simplex_EXIT_FAILURE) {
			pan_profile_failure("entering-column projection");
			*status = lp_simplex_PrecisionError;
			return lp_simplex_EXIT_FAILURE;
		}
		column_norm = simplex_pan_column_norm(basis->standard, entering);
		if (projection_residual > options->pivot_tolerance *
		    (1. + column_norm)) {
			int add_state = pan_add_independent(basis, work, entering,
				options->primal_tolerance);
			if (add_state == lp_simplex_EXIT_FAILURE) {
				pan_profile_failure("independent-column expansion");
				*status = lp_simplex_PrecisionError;
				return lp_simplex_EXIT_FAILURE;
			}
			if (add_state == 1) {
				if (pan_exchange(basis, work, entering,
					options->primal_tolerance,
					options->pivot_tolerance, &unbounded) ==
				    lp_simplex_EXIT_FAILURE) {
					pan_profile_failure("rank-recovery exchange");
					*status = lp_simplex_PrecisionError;
					return lp_simplex_EXIT_FAILURE;
				}
				if (unbounded) {
					*status = lp_simplex_Unboundedness;
					return lp_simplex_EXIT_FAILURE;
				}
			}
		} else {
			if (pan_exchange(basis, work, entering,
				options->primal_tolerance, options->pivot_tolerance,
				&unbounded) == lp_simplex_EXIT_FAILURE) {
				pan_profile_failure("dependent-column exchange");
				*status = lp_simplex_PrecisionError;
				return lp_simplex_EXIT_FAILURE;
			}
			if (unbounded) {
				*status = lp_simplex_Unboundedness;
				return lp_simplex_EXIT_FAILURE;
			}
		}
		(*iterations)++;
	}
}


int simplex_pan_solve_problem(
		const struct simplex_Problem *problem,
		const struct lp_simplex_Options *options,
		double *x, struct lp_simplex_Result *result)
{
	struct simplex_PanStandard standard;
	struct simplex_PanBasis basis;
	struct pan_Workspace work;
	clock_t started = clock();
	double objective = 0.;
	int iterations = 0;
	int phase1_iterations = 0;
	long phase1_factorizations = 0;
	int status = lp_simplex_CondUnsatisfied;
	int state = lp_simplex_EXIT_FAILURE;
	int j;
	lp_simplex_memset(&standard, 0, sizeof(standard));
	lp_simplex_memset(&basis, 0, sizeof(basis));
	lp_simplex_memset(&work, 0, sizeof(work));
	if (simplex_pan_standard_create(&standard, problem) ==
	    lp_simplex_EXIT_FAILURE ||
	    simplex_pan_basis_create(&basis, &standard,
		options->pivot_tolerance) == lp_simplex_EXIT_FAILURE ||
	    pan_workspace_create(&work, &standard) == lp_simplex_EXIT_FAILURE) {
		status = lp_simplex_MemoryAllocError;
		goto cleanup;
	}
	state = simplex_pan_phase1_run(&basis, work.value,
		options->primal_tolerance, options->dual_tolerance,
		options->iteration_limit, &iterations, &status);
	if (state == lp_simplex_EXIT_FAILURE)
		goto cleanup;
	phase1_iterations = iterations;
	phase1_factorizations = basis.factorizations;
	state = pan_phase2(&basis, &work, options, &iterations, &status,
		&result->dual_infeasibility);
	if (status != lp_simplex_Success)
		goto cleanup;
	if (pan_certify_feasible_values(&basis, &work,
		options->primal_tolerance) == lp_simplex_EXIT_FAILURE) {
		status = lp_simplex_PrecisionError;
		state = lp_simplex_EXIT_FAILURE;
		goto cleanup;
	}
	simplex_pan_standard_recover(&standard, work.value, x);
	for (j = 0; j < problem->columns; j++)
		objective += problem->objective[j] * x[j];
	result->objective = objective;
cleanup:
	result->status = status;
	result->iterations = iterations;
	if (getenv("LP_SIMPLEX_PROFILE") != NULL)
		fprintf(stderr,
			"pan-bda profile: total=%.6f phase-iterations=%d "
			"rank=%d/%d backend=%s factorizations=%ld "
				"[phase1=%d/%ld] updates=%ld refinements=%ld "
				"extensions=%ld downdates=%ld rotations=%ld "
				"symbolic-updates=%ld "
				"anchors=%ld local-retriangularizations=%ld "
				"block-householders=%ld max-front=%ld fill=%ld "
				"orthogonal-error=%.3g backward-error=%.3g "
				"factor-flops=%.6g\n",
			(double)(clock() - started) / (double)CLOCKS_PER_SEC,
			iterations, basis.count, standard.rows,
			PAN_BASIS_BACKEND_NAME, basis.factorizations,
			phase1_iterations, phase1_factorizations,
				basis.updates, basis.refinements,
				basis.extensions, basis.downdates, basis.rotations,
				basis.symbolic_updates,
				basis.anchor_rebuilds,
				basis.local_retriangularizations,
				basis.block_householders,
				basis.maximum_front, basis.factor_nonzeros,
				basis.orthogonal_error, basis.backward_error,
				basis.factor_flops);
	pan_workspace_destroy(&work);
	simplex_pan_basis_destroy(&basis);
	simplex_pan_standard_destroy(&standard);
	return state;
}
