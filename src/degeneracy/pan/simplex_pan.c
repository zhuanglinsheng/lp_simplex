/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_pan.c
 * @brief Pan's basis-deficiency-allowing generalized simplex method.
 *
 * In the notation of the original algorithm, the standard-form LP is the dual
 * problem
 *
 *     max b'z  subject to Az=c, z>=0,
 *
 * with b=-objective and c=rhs.  Its complementary inequality problem supplies
 * the minimum-norm multiplier used for pricing.  The number of columns in the
 * basis changes dynamically and is never padded to the row dimension.
 */

#include "simplex_pan.h"
#include "simplex_pan_config.h"
#include "simplex_pan_basis.h"
#include "simplex_pan_phase1.h"
#include "simplex_pan_standard.h"
#include "simplex_pan_cone.h"
#include "simplex_pan_bounded.h"
#include "simplex_dual.h"
#include "simplex_degeneracy.h"
#include "utils.h"

#include <lp_simplex/status.h>

#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>


static int pan_basis_create(struct simplex_PanBasis *basis,
	const struct simplex_PanStandard *standard,
	const struct lp_simplex_Options *options)
{
	struct simplex_PanConfig config;
	int state = simplex_pan_basis_create(
		basis, standard, options->pivot_tolerance);

	if (state == lp_simplex_EXIT_SUCCESS) {
		simplex_pan_config_init(&config, options);
		simplex_pan_basis_configure(basis, &config);
	}
	return state;
}

#ifdef LP_SIMPLEX_HAVE_SPQR
#define PAN_BASIS_BACKEND_NAME "dynamic-sparse-qr+spqr"
#else
#define PAN_BASIS_BACKEND_NAME "dense-fallback"
#endif

struct pan_Workspace {
	double *value;
	double *multiplier;
	double *anchor;
	double *probe;
	long compatible_choices;
	double *basis_right;
	double *basis_value;
	double *direction;
	unsigned char *remove;
	unsigned char *state;
	struct simplex_DegeneracyControl history;
	double previous_objective;
	int has_objective;
	int cycle_stage;
	long cycles;
};


static void pan_profile_failure(const char *stage)
{
	if (getenv("LP_SIMPLEX_PROFILE") != NULL)
		fprintf(stderr, "pan-bda precision failure: %s\n", stage);
}


static double pan_column_dot(const struct simplex_PanStandard *standard,
	const int column, const double *vector)
{
	return simplex_pan_column_dot(standard, column, vector);
}


static int pan_workspace_create(struct pan_Workspace *work,
	const struct simplex_PanStandard *standard,
	const struct lp_simplex_Options *options)
{
	lp_simplex_memset(work, 0, sizeof(*work));

	work->value = (double *)lp_simplex_malloc(
		(size_t)(standard->columns > 0 ? standard->columns : 1) *
		sizeof(double));
	if (work->value != NULL)
		lp_simplex_memset(work->value, 0,
			(size_t)standard->columns * sizeof(double));
	work->multiplier = (double *)lp_simplex_malloc(
		(size_t)standard->rows * sizeof(double));
	work->anchor = (double *)lp_simplex_malloc(
		(size_t)standard->rows * sizeof(double));
	if (work->anchor != NULL)
		lp_simplex_memset(work->anchor, 0,
			(size_t)standard->rows * sizeof(double));
	if (simplex_pan_setting(options, "LP_SIMPLEX_PAN_CONE") == NULL &&
		simplex_pan_setting(
			options, "LP_SIMPLEX_PAN_COMPATIBLE_PRICING") != NULL)
		work->probe = (double *)lp_simplex_malloc(
			(size_t)standard->rows * 2 * sizeof(double));
	work->basis_right = (double *)lp_simplex_malloc(
		(size_t)standard->rows * sizeof(double));
	work->basis_value = (double *)lp_simplex_malloc(
		(size_t)standard->rows * sizeof(double));
	work->direction = (double *)lp_simplex_malloc(
		(size_t)standard->rows * sizeof(double));
	work->state = (unsigned char *)lp_simplex_malloc(
		(size_t)(standard->columns > 0 ? standard->columns : 1));
	work->remove =
		(unsigned char *)lp_simplex_malloc((size_t)standard->rows);
	if (work->value == NULL || work->multiplier == NULL ||
		work->anchor == NULL || work->basis_right == NULL ||
		work->basis_value == NULL || work->direction == NULL ||
		work->remove == NULL || work->state == NULL ||
		(simplex_pan_setting(options, "LP_SIMPLEX_PAN_CONE") == NULL &&
			simplex_pan_setting(options,
				"LP_SIMPLEX_PAN_COMPATIBLE_PRICING") != NULL &&
			work->probe == NULL))
		return lp_simplex_EXIT_FAILURE;
	return lp_simplex_EXIT_SUCCESS;
}


static void pan_workspace_destroy(struct pan_Workspace *work)
{
	lp_simplex_free(work->value);
	lp_simplex_free(work->multiplier);
	lp_simplex_free(work->anchor);
	lp_simplex_free(work->probe);
	lp_simplex_free(work->basis_right);
	lp_simplex_free(work->basis_value);
	lp_simplex_free(work->direction);
	lp_simplex_free(work->remove);
	lp_simplex_free(work->state);
	simplex_degeneracy_destroy(&work->history);
	lp_simplex_memset(work, 0, sizeof(*work));
}


static int pan_choose_violated(const struct simplex_PanBasis *basis,
	const double *multiplier, const double tolerance, const int pricing,
	const double *span_probe, const double *edge_probe,
	long *compatible_choices, double *maximum_infeasibility)
{
	const struct simplex_PanStandard *standard = basis->standard;
	double best_score = 0., preferred_score = 0.;
	int entering = -1, preferred = -1;
	int j;
	*maximum_infeasibility = 0.;
	for (j = 0; j < standard->columns; j++) {
		double residual, violation, norm, score;

		if (basis->position[j] >= 0)
			continue;
		residual = pan_column_dot(standard, j, multiplier) +
			   standard->objective[j];
		if (!isfinite(residual))
			return -3;
		violation = residual < 0. ? -residual : 0.;
		if (violation > *maximum_infeasibility)
			*maximum_infeasibility = violation;
		if (violation <= tolerance)
			continue;
		if (pricing == lp_simplex_PRICING_BLAND)
			return j;
		norm = simplex_pan_column_norm(standard, j);
		if (norm == 0.)
			return -2;
		score = pricing == lp_simplex_PRICING_PAN_NORMALIZED
				? violation / norm
				: violation;
		if (span_probe != NULL) {
			double threshold = 64. * DBL_EPSILON *
					   (1. + sqrt((double)standard->rows)) *
					   fmax(1., norm);

			if (fabs(pan_column_dot(standard, j, span_probe)) <=
					threshold &&
				pan_column_dot(standard, j, edge_probe) <=
					threshold &&
				score > preferred_score) {
				preferred_score = score;
				preferred = j;
			}
		}
		if (score > best_score ||
			(score == best_score &&
				(entering < 0 || j < entering))) {
			best_score = score;
			entering = j;
		}
	}
	if (preferred >= 0) {
		(*compatible_choices)++;
		return preferred;
	}
	return entering;
}


/* Certify the maintained feasible point against the original sparse columns.
 * A rank-revealing orthogonal solve is used only as a final recovery path; it
 * is deliberately not the ordinary update, because Pan's exchange itself
 * preserves feasibility and can leave a non-unique solution on a deficient
 * basis. */
static int pan_certify_feasible_values(struct simplex_PanBasis *basis,
	struct pan_Workspace *work, const double tolerance)
{
	double residual_norm;
	double feasibility_tolerance = tolerance;
	int i, needs_recovery = 0;

	residual_norm = simplex_pan_residual(basis->standard, work->value,
		work->multiplier, work->direction);
	if (residual_norm > feasibility_tolerance)
		needs_recovery = 1;
	for (i = 0; i < basis->standard->columns; i++)
		if (!isfinite(work->value[i]) ||
			work->value[i] < -feasibility_tolerance)
			needs_recovery = 1;
	if (!needs_recovery)
		return lp_simplex_EXIT_SUCCESS;
	if (simplex_pan_basis_least_squares_qr(basis, basis->standard->rhs,
		    work->basis_value,
		    &residual_norm) == lp_simplex_EXIT_FAILURE) {
		pan_profile_failure("feasibility QR recovery");
		return lp_simplex_EXIT_FAILURE;
	}
	lp_simplex_memset(work->value, 0,
		(size_t)basis->standard->columns * sizeof(double));
	for (i = 0; i < basis->count; i++) {
		double value = work->basis_value[i];

		if (!isfinite(value) || value < -feasibility_tolerance) {
			if (getenv("LP_SIMPLEX_PROFILE") != NULL)
				fprintf(stderr,
					"pan-bda negative refreshed value %.6g < %.6g\n",
					value, -feasibility_tolerance);
			return lp_simplex_EXIT_FAILURE;
		}
		work->value[basis->column[i]] = value > 0. ? value : 0.;
	}
	residual_norm = simplex_pan_residual(basis->standard, work->value,
		work->multiplier, work->direction);
	if (residual_norm > feasibility_tolerance)
		return lp_simplex_EXIT_FAILURE;
	return lp_simplex_EXIT_SUCCESS;
}


static int pan_add_independent(struct simplex_PanBasis *basis,
	struct pan_Workspace *work, const int entering, const double tolerance)
{
	int old_count = basis->count;

	(void)tolerance;
	if (simplex_pan_basis_add(basis, entering) == lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	if (simplex_pan_basis_factorize(basis) == lp_simplex_EXIT_FAILURE) {
		lp_simplex_memset(work->remove, 0, (size_t)basis->count);
		work->remove[old_count] = 1;
		simplex_pan_basis_remove_marked(basis, work->remove);
		if (simplex_pan_basis_rebuild(basis) ==
			lp_simplex_EXIT_FAILURE) {
			if (simplex_pan_basis_repair_rank(basis, work->value) ==
				lp_simplex_EXIT_FAILURE)
				return lp_simplex_EXIT_FAILURE;
			/* A changed active set invalidates the caller's projection and
			 * multiplier. Restart pricing from the preserved feasible point. */
			return 2;
		}
		/* The projection classified a numerically dependent column as
		 * independent.  Ask the caller to use the exchange path instead. */
		return 1;
	}
	work->value[entering] = 0.;
	return 0;
}


static int pan_exchange(struct simplex_PanBasis *basis,
	struct pan_Workspace *work, const int entering,
	const double primal_tolerance, const double pivot_tolerance,
	int *unbounded)
{
	double minimum = 0.;
	int bounded = 0;
	int i, leaving = -1, remove_count = 0;
	*unbounded = 0;
	for (i = 0; i < basis->count; i++) {
		if (!isfinite(work->direction[i]))
			return lp_simplex_EXIT_FAILURE;
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
		/* A tiny positive component or inaccurate reconstruction is not a
		 * nonnegative ray certificate, even when below the pivot threshold. */
		for (i = 0; i < basis->count; i++)
			if (work->direction[i] > 0.)
				return lp_simplex_EXIT_FAILURE;
		if (!isfinite(basis->projection_solve_error) ||
			basis->projection_solve_error >
				32. * DBL_EPSILON *
					(1. + basis->standard->rows) *
					simplex_pan_column_norm(
						basis->standard, entering))
			return lp_simplex_EXIT_FAILURE;
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
					1e-12 * (1. + __lp_simplex_ABS__(
							      minimum)) &&
				fabs(new_value) * basis->column_norm[i] <=
					primal_tolerance / (1. + basis->count))
				work->remove[i] = 1;
		}
		if (work->remove[i])
			work->value[basis->column[i]] = 0.;
		else if (new_value < -primal_tolerance) {
			pan_profile_failure("negative value after exchange");
			return lp_simplex_EXIT_FAILURE;
		} else
			work->value[basis->column[i]] =
				new_value > 0. ? new_value : 0.;
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
		if (simplex_pan_basis_replace_factorized(basis, leaving,
			    entering) == lp_simplex_EXIT_FAILURE) {
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


/* Return 1 for one fresh numeric recovery, -1 for persistent cycling.
 * Pan working-set order is irrelevant, so fingerprint membership rather
 * than QR positions; a fill-reducing permutation must not hide a cycle. */
static int pan_record_progress(
	struct simplex_PanBasis *basis, struct pan_Workspace *work)
{
	const struct simplex_PanStandard *standard = basis->standard;
	long double objective = standard->objective_offset, activity = 0.;
	long double error;
	int j, repeated;

	for (j = 0; j < standard->columns; j++) {
		long double product =
			(long double)standard->objective[j] * work->value[j];

		objective += product;
		activity += fabsl(product);
		work->state[j] = basis->position[j] >= 0 ? 1 : 0;
	}
	if (!isfinite(objective) || !isfinite(activity))
		return -1;
	error = 4. * LDBL_EPSILON * (standard->columns + 2.) * (1. + activity);
	if (work->has_objective && work->previous_objective - objective > error)
		simplex_degeneracy_reset_state(&work->history);
	work->has_objective = 1;
	work->previous_objective = (double)objective;
	work->history.fingerprint_initialized = 0;
	work->history.current_basis_hash = 0;
	work->history.current_status_hash = 0;
	repeated = simplex_degeneracy_record_state(
		&work->history, NULL, 0, work->state, standard->columns);
	if (!repeated)
		return 0;
	work->cycles++;
	if (work->cycle_stage == 0) {
		work->cycle_stage = 1;
		return 0;
	}
	if (work->cycle_stage == 1) {
		work->cycle_stage = 2;
		return 1;
	}
	return -1;
}


static int pan_phase2(struct simplex_PanBasis *basis,
	struct pan_Workspace *work, const struct lp_simplex_Options *options,
	int *iterations, int *status, double *dual_infeasibility)
{
	double residual_norm;
	int phase_steps = 0, protected_column = -1;
	int stabilize = simplex_pan_setting(
				options, "LP_SIMPLEX_PAN_STABILIZE") != NULL;
	int affine = simplex_pan_setting(options,
			     "LP_SIMPLEX_PAN_AFFINE_MULTIPLIER") != NULL;

	for (;;) {
		double column_norm, projection_residual;
		int entering, i, unbounded, probes_ready = 0;

		if (phase_steps > 0 && (phase_steps & (phase_steps - 1)) == 0 &&
			getenv("LP_SIMPLEX_PROFILE") != NULL)
			fprintf(stderr,
				"pan-bda phase2: steps=%d rank=%d/%d dual=%.6g objective=%.12g cycles=%ld guard=%d\n",
				phase_steps, basis->count,
				basis->standard->rows, *dual_infeasibility,
				work->previous_objective, work->cycles,
				work->cycle_stage);
		if (*iterations >= options->iteration_limit) {
			*status = lp_simplex_ExceedIterLimit;
			return lp_simplex_EXIT_FAILURE;
		}
		if (stabilize) {
			long drops = basis->stability_drops;

			if (simplex_pan_basis_stabilize(
				    basis, work->value, protected_column) ==
				lp_simplex_EXIT_FAILURE) {
				pan_profile_failure(
					"working-set condition estimate");
				*status = lp_simplex_PrecisionError;
				return lp_simplex_EXIT_FAILURE;
			}
			if (basis->stability_drops > drops)
				lp_simplex_memset(work->anchor, 0,
					(size_t)basis->standard->rows *
						sizeof(double));
		}
		{
			int recovery = pan_record_progress(basis, work);

			if (recovery < 0) {
				pan_profile_failure(
					"persistent working-set cycle");
				*status = lp_simplex_PrecisionError;
				return lp_simplex_EXIT_FAILURE;
			}
			if (recovery > 0) {
				if (simplex_pan_basis_rebuild(basis) ==
						lp_simplex_EXIT_FAILURE &&
					simplex_pan_basis_repair_rank(
						basis, work->value) ==
						lp_simplex_EXIT_FAILURE) {
					*status = lp_simplex_PrecisionError;
					return lp_simplex_EXIT_FAILURE;
				}
				lp_simplex_memset(work->anchor, 0,
					(size_t)basis->standard->rows *
						sizeof(double));
				simplex_degeneracy_reset_state(&work->history);
				continue;
			}
		}
		for (i = 0; i < basis->count; i++)
			work->basis_right[i] =
				-basis->standard->objective[basis->column[i]] -
				(affine ? pan_column_dot(basis->standard,
						  basis->column[i],
						  work->anchor)
					: 0.);
		if (simplex_pan_basis_minimum_norm(basis, work->basis_right,
			    work->multiplier,
			    &residual_norm) == lp_simplex_EXIT_FAILURE ||
			residual_norm > options->dual_tolerance *
						(1. + (double)basis->count)) {
			/* Dynamic factor modification is accepted only when the original
			 * basis equations certify it.  Rebuild from the sparse columns on
			 * failure; this is residual-driven, not an update-count heuristic. */
			if (simplex_pan_basis_rebuild(basis) ==
				lp_simplex_EXIT_FAILURE) {
				if (simplex_pan_basis_repair_rank(
					    basis, work->value) ==
					lp_simplex_EXIT_FAILURE) {
					*status = lp_simplex_PrecisionError;
					return lp_simplex_EXIT_FAILURE;
				}
				lp_simplex_memset(work->anchor, 0,
					(size_t)basis->standard->rows *
						sizeof(double));
			}
			for (i = 0; i < basis->count; i++)
				work->basis_right[i] =
					-basis->standard
						 ->objective[basis->column[i]] -
					(affine ? pan_column_dot(
							  basis->standard,
							  basis->column[i],
							  work->anchor)
						: 0.);
			if (simplex_pan_basis_minimum_norm(basis,
				    work->basis_right, work->multiplier,
				    &residual_norm) ==
					lp_simplex_EXIT_FAILURE ||
				residual_norm >
					options->dual_tolerance *
						(1. + (double)basis->count)) {
				pan_profile_failure("minimum-norm multiplier");
				if (getenv("LP_SIMPLEX_PROFILE") != NULL)
					fprintf(stderr,
						"pan-bda multiplier residual=%.17g backward=%.6g rank=%d\n",
						residual_norm,
						basis->backward_error,
						basis->count);
				*status = lp_simplex_PrecisionError;
				return lp_simplex_EXIT_FAILURE;
			}
		}
		if (affine)
			for (i = 0; i < basis->standard->rows; i++)
				work->multiplier[i] += work->anchor[i];
		lp_simplex_memcpy(work->anchor, work->multiplier,
			(size_t)basis->standard->rows * sizeof(double));
		if (work->probe != NULL && work->cycle_stage == 0 &&
			simplex_pan_basis_pricing_probes(basis, work->value,
				work->probe,
				work->probe + basis->standard->rows) ==
				lp_simplex_EXIT_SUCCESS)
			probes_ready = 1;
		entering = pan_choose_violated(basis, work->multiplier,
			options->dual_tolerance,
			work->cycle_stage ? lp_simplex_PRICING_BLAND
					  : options->pricing,
			probes_ready ? work->probe : NULL,
			probes_ready ? work->probe + basis->standard->rows
				     : NULL,
			&work->compatible_choices, dual_infeasibility);
		if (entering == -3) {
			*status = lp_simplex_PrecisionError;
			return lp_simplex_EXIT_FAILURE;
		}
		if (entering == -2) {
			*status = lp_simplex_Unboundedness;
			return lp_simplex_EXIT_FAILURE;
		}
		if (entering < 0) {
			*status = lp_simplex_Success;
			return lp_simplex_EXIT_SUCCESS;
		}
		if (simplex_pan_basis_column_projection(basis, entering,
			    work->direction,
			    &projection_residual) == lp_simplex_EXIT_FAILURE) {
			pan_profile_failure("entering-column projection");
			*status = lp_simplex_PrecisionError;
			return lp_simplex_EXIT_FAILURE;
		}
		column_norm =
			simplex_pan_column_norm(basis->standard, entering);
		if (projection_residual / column_norm >
			options->pivot_tolerance) {
			int add_state = pan_add_independent(basis, work,
				entering, options->primal_tolerance);

			if (add_state == lp_simplex_EXIT_FAILURE) {
				pan_profile_failure(
					"independent-column expansion");
				*status = lp_simplex_PrecisionError;
				return lp_simplex_EXIT_FAILURE;
			}
			if (add_state == 2) {
				lp_simplex_memset(work->anchor, 0,
					(size_t)basis->standard->rows *
						sizeof(double));
				(*iterations)++;
				phase_steps++;
				continue;
			}
			if (add_state == 1) {
				/* A fresh QR can reorder columns: the old direction no longer
				 * corresponds to active positions. Certify a new projection. */
				if (simplex_pan_basis_column_projection(basis,
					    entering, work->direction,
					    &projection_residual) ==
						lp_simplex_EXIT_FAILURE ||
					projection_residual >
						options->primal_tolerance) {
					*status = lp_simplex_PrecisionError;
					return lp_simplex_EXIT_FAILURE;
				}
				if (pan_exchange(basis, work, entering,
					    options->primal_tolerance,
					    options->pivot_tolerance,
					    &unbounded) ==
					lp_simplex_EXIT_FAILURE) {
					pan_profile_failure(
						"rank-recovery exchange");
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
				    options->primal_tolerance,
				    options->pivot_tolerance,
				    &unbounded) == lp_simplex_EXIT_FAILURE) {
				pan_profile_failure(
					"dependent-column exchange");
				*status = lp_simplex_PrecisionError;
				return lp_simplex_EXIT_FAILURE;
			}
			if (unbounded) {
				*status = lp_simplex_Unboundedness;
				return lp_simplex_EXIT_FAILURE;
			}
		}
		if (pan_certify_feasible_values(
			    basis, work, options->primal_tolerance) ==
			lp_simplex_EXIT_FAILURE) {
			pan_profile_failure("primal certificate after update");
			*status = lp_simplex_PrecisionError;
			return lp_simplex_EXIT_FAILURE;
		}
		protected_column = entering;
		(*iterations)++;
		phase_steps++;
	}
}


int simplex_pan_solve_from_feasible(const struct simplex_Problem *problem,
	const struct lp_simplex_Options *options, const double *start,
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

	if (simplex_pan_setting(options, "LP_SIMPLEX_PAN_NATIVE_BOUNDS"))
		return simplex_pan_bounded_solve(
			problem, options, start, x, result);
	result->dual_infeasibility = HUGE_VAL;
	result->dual_residual_space = lp_simplex_RESIDUAL_UNAVAILABLE;
	result->primal_infeasibility = HUGE_VAL;
	lp_simplex_memset(&standard, 0, sizeof(standard));
	lp_simplex_memset(&basis, 0, sizeof(basis));
	lp_simplex_memset(&work, 0, sizeof(work));
	if (simplex_pan_standard_create(&standard, problem) ==
			lp_simplex_EXIT_FAILURE ||
		pan_basis_create(&basis, &standard, options) ==
			lp_simplex_EXIT_FAILURE ||
		pan_workspace_create(&work, &standard, options) ==
			lp_simplex_EXIT_FAILURE) {
		status = lp_simplex_MemoryAllocError;
		goto cleanup;
	}
	if (options->algorithm == lp_simplex_ALGORITHM_AUTO ||
		options->algorithm == SIMPLEX_PAN_AUTO_STANDARD_IMPLICIT ||
		options->algorithm == SIMPLEX_PAN_AUTO_STANDARD_EXPLICIT) {
		struct simplex_PanConfig config;

		simplex_pan_config_init(&config, options);
		simplex_pan_basis_configure(&basis, &config);
	}
	if (start != NULL) {
		int rank_failure = 0;

		if (simplex_pan_standard_encode(&standard, start, work.value,
			    work.multiplier, options->primal_tolerance) ==
			lp_simplex_EXIT_FAILURE) {
			status = lp_simplex_PrecisionError;
			goto cleanup;
		}
		for (j = 0; j < standard.columns; j++)
			if (work.value[j] > 0. &&
				simplex_pan_basis_add(&basis, j) ==
					lp_simplex_EXIT_FAILURE) {
				rank_failure = 1;
				break;
			}
		if (rank_failure || simplex_pan_basis_rebuild(&basis) ==
					    lp_simplex_EXIT_FAILURE) {
			if (simplex_pan_setting(
				    options, "LP_SIMPLEX_PAN_CONE") == NULL ||
				(simplex_pan_setting(options,
					 "LP_SIMPLEX_PAN_FACE_GUARD") == NULL &&
					simplex_pan_setting(options,
						"LP_SIMPLEX_PAN_FIRST_ORDER") ==
						NULL) ||
				simplex_pan_support_compress(
					&standard, work.value, options) != 1) {
				status = lp_simplex_PrecisionError;
				goto cleanup;
			}
			simplex_pan_basis_destroy(&basis);
			if (pan_basis_create(&basis, &standard, options)) {
				status = lp_simplex_MemoryAllocError;
				goto cleanup;
			}
			for (j = 0; j < standard.columns; j++)
				if (work.value[j] > 0. &&
					simplex_pan_basis_add(&basis, j)) {
					status = lp_simplex_PrecisionError;
					goto cleanup;
				}
			if (simplex_pan_basis_rebuild(&basis)) {
				status = lp_simplex_PrecisionError;
				goto cleanup;
			}
		}
		if (getenv("LP_SIMPLEX_PROFILE") != NULL)
			fprintf(stderr,
				"pan-bda profile: feasible warm start rank=%d/%d\n",
				basis.count, standard.rows);
		state = lp_simplex_EXIT_SUCCESS;
	} else
		state = simplex_pan_phase1_run(&basis, work.value,
			options->primal_tolerance, options->dual_tolerance,
			options->iteration_limit, &iterations, &status);
	phase1_iterations = iterations;
	phase1_factorizations = basis.factorizations;
	if (state == lp_simplex_EXIT_FAILURE)
		goto cleanup;
	if (pan_certify_feasible_values(&basis, &work,
		    options->primal_tolerance) == lp_simplex_EXIT_FAILURE) {
		status = lp_simplex_PrecisionError;
		state = lp_simplex_EXIT_FAILURE;
		goto cleanup;
	}
	if (start == NULL && getenv("LP_SIMPLEX_PROFILE"))
		fprintf(stderr,
			"pan-bda first feasible: source=nnls cpu=%.6f iterations=%d primal=%.6g\n",
			(double)(clock() - started) / CLOCKS_PER_SEC,
			iterations,
			simplex_pan_residual(&standard, work.value,
				work.multiplier, work.direction));
	if (simplex_pan_setting(options, "LP_SIMPLEX_PAN_CONE") == NULL &&
		simplex_pan_setting(options, "LP_SIMPLEX_PAN_COMPLETE_START") !=
			NULL) {
		int old_count = basis.count;

		if (simplex_pan_basis_complete(&basis, work.value) ==
			lp_simplex_EXIT_FAILURE) {
			pan_profile_failure("initial independent completion");
			status = lp_simplex_PrecisionError;
			state = lp_simplex_EXIT_FAILURE;
			goto cleanup;
		}
		if (getenv("LP_SIMPLEX_PROFILE") != NULL)
			fprintf(stderr,
				"pan-bda profile: completed start rank=%d->%d/%d\n",
				old_count, basis.count, standard.rows);
	}
	if (simplex_pan_setting(options, "LP_SIMPLEX_PAN_CONE") != NULL)
		state = simplex_pan_cone_phase2(&standard, &basis, work.value,
			work.anchor, options, &iterations, &status,
			&result->dual_infeasibility);
	else
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
	{
		long double primal = 0., dual = 0.;

		lp_simplex_memcpy(work.multiplier, work.anchor,
			(size_t)standard.rows * sizeof(double));
		for (j = 0; j < standard.columns; j++) {
			double reduced =
				standard.objective[j] +
				pan_column_dot(&standard, j, work.multiplier);
			double budget = options->dual_tolerance *
					(1. + fabs(standard.objective[j]));

			if (!isfinite(reduced) || reduced < -budget ||
				(work.value[j] > 0. &&
					fabs(reduced) > budget)) {
				pan_profile_failure("final dual certificate");
				status = lp_simplex_PrecisionError;
				state = lp_simplex_EXIT_FAILURE;
				goto cleanup;
			}
			primal += (long double)standard.objective[j] *
				  work.value[j];
		}
		for (j = 0; j < standard.rows; j++)
			dual -= (long double)standard.rhs[j] *
				work.multiplier[j];
		if (!isfinite(primal) || !isfinite(dual) ||
			fabsl(primal - dual) >
				options->dual_tolerance *
					(1. + fabsl(primal +
						      standard.objective_offset))) {
			pan_profile_failure("final primal-dual gap");
			status = lp_simplex_PrecisionError;
			state = lp_simplex_EXIT_FAILURE;
			goto cleanup;
		}
	}
	simplex_pan_standard_recover(&standard, work.value, x);
	for (j = 0; j < problem->columns; j++)
		objective += problem->objective[j] * x[j];
	result->objective = objective;
cleanup:
	/* Preserve the actual partial point for diagnostics, without changing
	 * the failure status or claiming a certificate for it. */
	if (work.value != NULL && standard.original_shift != NULL) {
		simplex_pan_standard_recover(&standard, work.value, x);
		objective = 0.;
		for (j = 0; j < problem->columns; j++)
			objective += problem->objective[j] * x[j];
		result->objective = objective;
	}
	result->status = status;
	result->dual_infeasibility = NAN;
	result->dual_residual_space = lp_simplex_RESIDUAL_UNAVAILABLE;
	if (status == lp_simplex_Success) {
		result->dual_infeasibility = 0.;
		for (j = 0; j < standard.columns; j++) {
			double g =
				standard.objective[j] +
				pan_column_dot(&standard, j, work.multiplier);
			double v = work.value[j] > options->primal_tolerance
					   ? fabs(g)
					   : fmax(0., -g);

			if (!isfinite(g)) {
				result->dual_infeasibility = HUGE_VAL;
				break;
			}
			result->dual_infeasibility =
				fmax(result->dual_infeasibility, v);
		}
		result->dual_residual_space = lp_simplex_RESIDUAL_STANDARD;
	}
	result->iterations = iterations;
	if (getenv("LP_SIMPLEX_PROFILE") != NULL)
		fprintf(stderr,
			"pan-bda profile: total=%.6f phase-iterations=%d "
			"rank=%d/%d backend=%s factorizations=%ld "
			"[phase1=%d/%ld] updates=%ld refinements=%ld rank-repairs=%ld stability-drops=%ld rcond=%.3g rank-disagreements=%ld cycles=%ld guard=%d compatible=%ld "
			"extensions=%ld downdates=%ld rotations=%ld "
			"symbolic-updates=%ld "
			"anchors=%ld local-retriangularizations=%ld "
			"block-householders=%ld max-front=%ld fill=%ld "
			"orthogonal-error=%.3g backward-error=%.3g "
			"factor-flops=%.6g\n",
			(double)(clock() - started) / (double)CLOCKS_PER_SEC,
			iterations, basis.count, standard.rows,
			PAN_BASIS_BACKEND_NAME, basis.factorizations,
			phase1_iterations, phase1_factorizations, basis.updates,
			basis.refinements, basis.rank_repairs,
			basis.stability_drops, basis.reciprocal_condition,
			basis.rank_disagreements, work.cycles, work.cycle_stage,
			work.compatible_choices, basis.extensions,
			basis.downdates, basis.rotations,
			basis.symbolic_updates, basis.anchor_rebuilds,
			basis.local_retriangularizations,
			basis.block_householders, basis.maximum_front,
			basis.factor_nonzeros, basis.orthogonal_error,
			basis.backward_error, basis.factor_flops);
	pan_workspace_destroy(&work);
	simplex_pan_basis_destroy(&basis);
	simplex_pan_standard_destroy(&standard);
	return state;
}


int simplex_pan_solve_problem(const struct simplex_Problem *problem,
	const struct lp_simplex_Options *options, double *x,
	struct lp_simplex_Result *result)
{
	if (simplex_pan_setting(options, "LP_SIMPLEX_PAN_NATIVE_BOUNDS"))
		return simplex_pan_bounded_solve(
			problem, options, NULL, x, result);

	if (simplex_pan_setting(options, "LP_SIMPLEX_PAN_DUAL_START") != NULL) {
		struct lp_simplex_Options remaining = *options;
		struct lp_simplex_Result phase;
		double *seed = (double *)lp_simplex_malloc(
			(size_t)(problem->columns > 0 ? problem->columns : 1) *
			sizeof(double));
		int state;
		clock_t feasibility_started = clock();

		if (seed == NULL) {
			result->status = lp_simplex_MemoryAllocError;
			return lp_simplex_EXIT_FAILURE;
		}
		state = simplex_dual_find_feasible(
			problem, options, seed, &phase);
		if (state == 0 && getenv("LP_SIMPLEX_PROFILE"))
			fprintf(stderr,
				"pan-bda first feasible: source=dual cpu=%.6f iterations=%d primal=%.6g\n",
				(double)(clock() - feasibility_started) /
					CLOCKS_PER_SEC,
				phase.iterations, phase.primal_infeasibility);
		if (getenv("LP_SIMPLEX_PROFILE") != NULL)
			fprintf(stderr,
				"pan-bda profile: dual feasibility start iterations=%d status=%d primal=%.6g\n",
				phase.iterations, phase.status,
				phase.primal_infeasibility);
		if (state == lp_simplex_EXIT_FAILURE) {
			*result = phase;
			lp_simplex_free(seed);
			return state;
		}
		remaining.iteration_limit -= phase.iterations;
		state = simplex_pan_solve_from_feasible(
			problem, &remaining, seed, x, result);
		result->iterations += phase.iterations;
		lp_simplex_free(seed);
		return state;
	}
	return simplex_pan_solve_from_feasible(
		problem, options, NULL, x, result);
}
