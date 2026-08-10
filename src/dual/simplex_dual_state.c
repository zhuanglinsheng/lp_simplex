/* Allocation, ownership, and structural policy for a dual solve. */
#include "simplex_dual_state.h"
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
	int *frequency = (int *)lp_simplex_malloc(
		(size_t)(matrix->rows + 1) * sizeof(int));
	if (frequency == NULL)
		return 0;
	lp_simplex_memset(frequency, 0,
		(size_t)(matrix->rows + 1) * sizeof(int));
	for (column = 0; column < matrix->columns; column++) {
		int degree = matrix->column_start[column + 1] -
			matrix->column_start[column];
		frequency[degree]++;
		if (frequency[degree] > best) {
			best = frequency[degree];
			*modal_degree = degree;
		}
	}
	lp_simplex_free(frequency);
	return matrix->columns > 0 &&
		(long)best * 10L >= (long)matrix->columns * 9L;
}


static void dual_bind_workspace(struct simplex_DualState *state)
{
	int rows = state->rows;
	int variables = state->variables;
	state->basis = state->workspace.index;
	state->position = state->workspace.index != NULL
		? state->workspace.index + rows : NULL;
	state->status = state->workspace.status;
	state->lower = state->workspace.column;
	state->upper = state->workspace.column != NULL
		? state->workspace.column + variables : NULL;
	state->cost = state->workspace.column != NULL
		? state->workspace.column + 2 * variables : NULL;
	state->value = state->workspace.column != NULL
		? state->workspace.column + 3 * variables : NULL;
	state->reduced = state->workspace.column != NULL
		? state->workspace.column + 4 * variables : NULL;
	state->alpha = state->workspace.column != NULL
		? state->workspace.column + 5 * variables : NULL;
	state->breakpoint = state->workspace.column != NULL
		? state->workspace.column + 6 * variables : NULL;
	state->basic_value = state->workspace.row;
	state->basic_lower = state->workspace.row != NULL
		? state->workspace.row + rows : NULL;
	state->basic_upper = state->workspace.row != NULL
		? state->workspace.row + 2 * rows : NULL;
	state->feasibility_tolerance = state->workspace.row != NULL
		? state->workspace.row + 3 * rows : NULL;
	state->pi = state->workspace.row != NULL
		? state->workspace.row + 4 * rows : NULL;
	state->rho = state->workspace.row != NULL
		? state->workspace.row + 5 * rows : NULL;
	state->direction = state->workspace.row != NULL
		? state->workspace.row + 6 * rows : NULL;
	state->work = state->workspace.row != NULL
		? state->workspace.row + 7 * rows : NULL;
	state->edge_weight = state->workspace.row != NULL
		? state->workspace.row + 8 * rows : NULL;
}


void simplex_dual_state_destroy(struct simplex_DualState *state)
{
	if (state == NULL)
		return;
	simplex_degeneracy_destroy(&state->degeneracy);
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


int simplex_dual_state_create(
		struct simplex_DualState *state, const struct simplex_Problem *problem,
		const struct lp_simplex_Options *options)
{
	int variables = problem->columns + problem->rows;
	lp_simplex_memset(state, 0, sizeof(*state));
	state->rows = problem->rows;
	state->structural = problem->columns;
	state->variables = variables;
	state->options = options;
	state->pan_deferred_row = -1;
	state->profile.enabled = getenv("LP_SIMPLEX_PROFILE") != NULL;
	state->matrix = problem->matrix;
	state->matrix.owns_storage = 0;
	state->numerically_stressed = dual_matrix_is_numerically_stressed(
		&state->matrix, __lp_simplex_MAX__(options->primal_tolerance,
			options->dual_tolerance));
	state->pan_enabled = getenv("LP_SIMPLEX_DISABLE_PAN") == NULL &&
		!state->numerically_stressed;
	state->regular_columns = dual_matrix_has_regular_columns(
		&state->matrix, &state->modal_column_degree);
	state->stable_candidate_order = state->rows >= 2048 &&
		!state->regular_columns;
	state->workspace.index = (int *)lp_simplex_malloc(
		(size_t)(problem->rows + variables) * sizeof(int));
	state->workspace.status = (unsigned char *)lp_simplex_malloc(
		(size_t)variables * sizeof(unsigned char));
	state->workspace.column = (double *)lp_simplex_malloc(
		(size_t)7 * variables * sizeof(double));
	state->workspace.row = (double *)lp_simplex_malloc(
		(size_t)9 * problem->rows * sizeof(double));
	dual_bind_workspace(state);
	if (state->workspace.index == NULL || state->workspace.status == NULL ||
	    state->workspace.column == NULL || state->workspace.row == NULL ||
	    simplex_dual_pricing_create(state) == lp_simplex_EXIT_FAILURE ||
	    simplex_sparse_vector_create(&state->packed_direction, problem->rows) ==
	    lp_simplex_EXIT_FAILURE ||
	    simplex_sparse_vector_create(&state->flip_rhs, problem->rows) ==
	    lp_simplex_EXIT_FAILURE ||
	    simplex_dual_feasibility_create(&state->feasibility, problem->rows) ==
	    lp_simplex_EXIT_FAILURE ||
	    simplex_basis_create(&state->factor, &state->matrix, problem->columns,
		state->basis) == lp_simplex_EXIT_FAILURE) {
		simplex_dual_state_destroy(state);
		return lp_simplex_EXIT_FAILURE;
	}
	lp_simplex_memset(state->alpha, 0, (size_t)variables * sizeof(double));
	lp_simplex_memset(state->feasibility_tolerance, 0,
		(size_t)problem->rows * sizeof(double));
	simplex_basis_set_sparse_eta(&state->factor,
		!state->regular_columns || state->modal_column_degree == 1);
	simplex_degeneracy_initialize(&state->degeneracy);
	return lp_simplex_EXIT_SUCCESS;
}
