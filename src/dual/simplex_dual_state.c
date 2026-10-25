/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_dual_state.c
 * @brief Allocation, ownership, and structural policy for a dual solve.
 */

#include "policy/factorization.h"
#include "policy/controls.h"
#include "simplex_dual_state.h"
#include "simplex_dual_storage.h"
#include "simplex_dual_internal.h"
#include "simplex_dual_pricing.h"
#include "simplex_problem.h"
#include "utils.h"

#include <lp_simplex/status.h>

#include <stdlib.h>

static int dual_matrix_is_numerically_stressed(
	const struct simplex_CscMatrix *matrix, const double tolerance)
{
	double smallest = __lp_simplex_INF__, largest = 0.;
	int k;

	for (k = 0; k < matrix->nonzeros; k++) {
		double magnitude = __lp_simplex_ABS__(matrix->value[k]);

		if (magnitude == 0.)
			continue;
		smallest = __lp_simplex_MIN__(smallest, magnitude);
		largest = __lp_simplex_MAX__(largest, magnitude);
	}
	return smallest < __lp_simplex_INF__ && largest * tolerance >= smallest;
}

static int dual_matrix_has_regular_columns(
	const struct simplex_CscMatrix *matrix, int *modal_degree)
{
	int column, best = 0;
	size_t bytes;
	int *frequency;

	if (matrix->rows < 0 || !simplex_storage_bytes((size_t)matrix->rows + 1,
					sizeof(int), (size_t)-1, &bytes))
		return 0;
	frequency = (int *)lp_simplex_malloc(bytes);

	if (frequency == NULL)
		return 0;
	lp_simplex_memset(frequency, 0, bytes);
	for (column = 0; column < matrix->columns; column++) {
		int start = matrix->column_start[column];
		int end = matrix->column_start[column + 1];
		int degree;

		if (start < 0 || end < start || end > matrix->nonzeros) {
			lp_simplex_free(frequency);
			return 0;
		}
		degree = end - start;
		if (degree > matrix->rows) {
			lp_simplex_free(frequency);
			return 0;
		}
		frequency[degree]++;
		if (frequency[degree] > best) {
			best = frequency[degree];
			*modal_degree = degree;
		}
	}
	lp_simplex_free(frequency);
	return matrix->columns > 0 &&
	       best >= matrix->columns - matrix->columns / 10;
}

static void dual_bind_workspace(struct simplex_DualState *state)
{
	int rows = state->rows;
	int variables = state->variables;

	state->basis = state->workspace.index;
	state->position = state->workspace.index != NULL
				  ? state->workspace.index + rows
				  : NULL;
	state->status = state->workspace.status;
	state->lower = state->workspace.column;
	state->upper = state->workspace.column != NULL
			       ? state->workspace.column + variables
			       : NULL;
	state->cost = state->workspace.column != NULL
			      ? state->workspace.column + (size_t)2 * variables
			      : NULL;
	state->original_cost =
		state->workspace.column != NULL
			? state->workspace.column + (size_t)3 * variables
			: NULL;
	state->value = state->workspace.column != NULL
			       ? state->workspace.column + (size_t)4 * variables
			       : NULL;
	state->reduced =
		state->workspace.column != NULL
			? state->workspace.column + (size_t)5 * variables
			: NULL;
	/* Pricing owns the shared dense ratio row; bind it after creation. */
	state->alpha = NULL;
	state->breakpoint =
		state->workspace.column != NULL
			? state->workspace.column + (size_t)6 * variables
			: NULL;
	state->basic_value = state->workspace.row;
	state->basic_lower = state->workspace.row != NULL
				     ? state->workspace.row + rows
				     : NULL;
	state->basic_upper = state->workspace.row != NULL
				     ? state->workspace.row + (size_t)2 * rows
				     : NULL;
	state->feasibility_tolerance =
		state->workspace.row != NULL
			? state->workspace.row + (size_t)3 * rows
			: NULL;
	state->pi = state->workspace.row != NULL
			    ? state->workspace.row + (size_t)4 * rows
			    : NULL;
	state->rho = state->workspace.row != NULL
			     ? state->workspace.row + (size_t)5 * rows
			     : NULL;
	state->direction = state->workspace.row != NULL
				   ? state->workspace.row + (size_t)6 * rows
				   : NULL;
	state->work = state->workspace.row != NULL
			      ? state->workspace.row + (size_t)7 * rows
			      : NULL;
	state->edge_weight = state->workspace.row != NULL
				     ? state->workspace.row + (size_t)8 * rows
				     : NULL;
}

void simplex_dual_state_destroy(struct simplex_DualState *state)
{
	if (state == NULL)
		return;
	simplex_recovery_destroy(state);
	simplex_degeneracy_destroy(&state->degeneracy);
	simplex_degeneracy_destroy(&state->crash_history);
	simplex_dual_feasibility_destroy(&state->feasibility);
	simplex_basis_destroy(&state->factor);
	simplex_csc_destroy(&state->matrix);
	lp_simplex_free(state->workspace.index);
	lp_simplex_free(state->workspace.status);
	lp_simplex_free(state->workspace.column);
	lp_simplex_free(state->workspace.row);
	simplex_dual_pricing_destroy(state);
	simplex_sparse_vector_destroy(&state->packed_direction);
	simplex_sparse_vector_destroy(&state->flip_rhs);
	lp_simplex_memset(state, 0, sizeof(*state));
}

int simplex_dual_state_create(struct simplex_DualState *state,
	const struct simplex_Problem *problem,
	const struct lp_simplex_Options *options, const int basis_backend)
{
	int variables;
	struct simplex_DualStorage storage;
	struct simplex_RunConfig local;
	const struct simplex_RunConfig *controls;

	lp_simplex_memset(state, 0, sizeof(*state));
	if (!simplex_dual_storage(
		    problem->rows, problem->columns, (size_t)-1, &storage))
		return lp_simplex_EXIT_FAILURE;
	variables = storage.variables;
	controls = simplex_run_config_get(problem->execution.context, &local);
	state->rows = problem->rows;
	state->structural = problem->columns;
	state->variables = variables;
	state->options = options;
	state->pan_deferred_row = -1;
	simplex_recovery_initialize(&state->recovery);
	state->controls = controls->dual;
	state->recovery.force_bounds = state->controls.force_bounds;
	state->recovery.disable_bounds = state->controls.disable_bounds;
	state->profile.enabled = getenv(SIMPLEX_ENV_PROFILE) != NULL;
	state->matrix = problem->matrix;
	state->matrix.owns_storage = 0;
	state->numerically_stressed = dual_matrix_is_numerically_stressed(
		&state->matrix, __lp_simplex_MAX__(options->primal_tolerance,
					options->dual_tolerance));
	/* Degeneracy and numerical stress are independent.  Pivot safeguards
	 * remain active when the matrix is ill-scaled. */
	state->pan_enabled = state->controls.pan_enabled;
	state->perturb_multi_epoch = state->controls.multi_epoch;
	state->regular_columns = dual_matrix_has_regular_columns(
		&state->matrix, &state->modal_column_degree);
	state->stable_candidate_order =
		state->rows >= SIMPLEX_DUAL_STABLE_ORDER_ROWS &&
		!state->regular_columns;
	state->workspace.index = (int *)lp_simplex_malloc(storage.index);
	state->workspace.status =
		(unsigned char *)lp_simplex_malloc(storage.status);
	state->workspace.column = (double *)lp_simplex_malloc(storage.column);
	state->workspace.row = (double *)lp_simplex_malloc(storage.row);
	dual_bind_workspace(state);
	if (state->workspace.index == NULL || state->workspace.status == NULL ||
		state->workspace.column == NULL ||
		state->workspace.row == NULL ||
		simplex_dual_pricing_create(state) == lp_simplex_EXIT_FAILURE ||
		simplex_sparse_vector_create(&state->packed_direction,
			problem->rows) == lp_simplex_EXIT_FAILURE ||
		simplex_sparse_vector_create(&state->flip_rhs, problem->rows) ==
			lp_simplex_EXIT_FAILURE ||
		simplex_dual_feasibility_create(&state->feasibility,
			problem->rows) == lp_simplex_EXIT_FAILURE ||
		simplex_basis_create_configured(&state->factor, &state->matrix,
			problem->columns, state->basis,
			state->controls.stable_klu) ==
			lp_simplex_EXIT_FAILURE) {
		simplex_dual_state_destroy(state);
		return lp_simplex_EXIT_FAILURE;
	}
	/* Sparse pricing and dense pricing share one ratio row. Packed values
	 * remain a snapshot for reduced-cost updates after the basis exchange. */
	state->alpha = variables > 0 ? state->packed_alpha.dense
				     : state->workspace.column;
	lp_simplex_memset(state->feasibility_tolerance, 0,
		(size_t)problem->rows * sizeof(double));
	simplex_basis_set_backend(&state->factor, basis_backend);
	simplex_basis_set_sparse_eta(&state->factor,
		!state->regular_columns || state->modal_column_degree == 1);
	simplex_degeneracy_initialize(&state->degeneracy);
	simplex_degeneracy_initialize(&state->crash_history);
	return lp_simplex_EXIT_SUCCESS;
}
