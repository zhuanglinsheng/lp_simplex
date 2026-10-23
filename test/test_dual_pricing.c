/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file test_dual_pricing.c
 * @brief Dual ratio ordering and symbolic tie regressions.
 */

#include "simplex_dual_pricing.h"
#include "simplex_dual_state.h"
#include "simplex_problem.h"
#include "utils.h"

#include <stdio.h>
#include <math.h>


/* Equal real ratios must be ordered by the transformed symbolic costs,
 * rather than by the entering variable id.  Exercise both leaving bounds. */
static int test_symbolic_ratio(const int kappa, const int backend)
{
	struct simplex_Problem *problem =
		simplex_problem_create_sparse(2, 2, 3);
	struct simplex_DualState state;
	struct lp_simplex_Options options;
	double theta = -1.;
	int j, valid;

	if (problem == NULL)
		return 0;
	problem->matrix.column_start[0] = 0;
	problem->matrix.column_start[1] = 1;
	problem->matrix.column_start[2] = 3;
	problem->matrix.row_index[0] = problem->matrix.row_index[1] = 0;
	problem->matrix.row_index[2] = 1;
	problem->matrix.value[0] = problem->matrix.value[1] =
		problem->matrix.value[2] = 1.;
	problem->matrix.row_start[0] = 0;
	problem->matrix.row_start[1] = 2;
	problem->matrix.row_start[2] = 3;
	problem->matrix.column_index[0] = 0;
	problem->matrix.column_index[1] = problem->matrix.column_index[2] = 1;
	problem->matrix.row_value[0] = problem->matrix.row_value[1] =
		problem->matrix.row_value[2] = 1.;
	lp_simplex_default_options(&options, lp_simplex_ALGORITHM_DUAL_REVISED);
	if (simplex_dual_state_create(&state, problem, &options, backend) !=
		lp_simplex_EXIT_SUCCESS) {
		simplex_problem_free(problem);
		return 0;
	}
	for (j = 0; j < 4; j++) {
		state.position[j] = j < 2 ? -1 : j - 2;
		state.status[j] = j < 2 ? (kappa > 0 ? DUAL_STATUS_LOWER
						     : DUAL_STATUS_UPPER)
					: DUAL_STATUS_BASIC;
		state.reduced[j] = state.cost[j] = 0.;
	}
	state.basis[0] = 2;
	state.basis[1] = 3;
	valid = simplex_basis_factorize(&state.factor) ==
		lp_simplex_EXIT_SUCCESS;
	state.rho[0] = -1.;
	state.rho[1] = 0.;
	simplex_dual_nonbasic_initialize(&state);
	simplex_degeneracy_observe(&state.degeneracy, 0., 1e-8, 2., 2., 1e-8);
	state.degeneracy.phase = SIMPLEX_PAN_LEXICOGRAPHIC;
	valid = valid && simplex_dual_prepare_ratio_test(&state, kappa) == 2 &&
		simplex_dual_next_ratio_candidate(
			&state, &theta, 1, 2, kappa) == 1 &&
		theta == 0.;
	/* A nearby but nonzero primary ratio is not a symbolic tie. A Harris
	 * window here could choose column 1 again and violate lexicographic order. */
	state.reduced[1] = kappa > 0 ? 1e-13 : -1e-13;
	valid = valid && simplex_dual_prepare_ratio_test(&state, kappa) == 2 &&
		simplex_dual_next_ratio_candidate(
			&state, &theta, 1, 2, kappa) == 0 &&
		theta == 0.;
	state.reduced[1] = 0.;
	/* Keeping the same symbolic objective after a basis exchange requires
	 * a BTRAN for an anchor coordinate which has become basic. */
	if (valid) {
		double direction[2] = {-1., 0.};

		state.basis[0] = 0;
		state.position[0] = 0;
		state.position[2] = -1;
		state.status[0] = DUAL_STATUS_BASIC;
		state.status[2] =
			kappa > 0 ? DUAL_STATUS_UPPER : DUAL_STATUS_LOWER;
		valid = simplex_basis_update(&state.factor, 0, direction) !=
			lp_simplex_EXIT_FAILURE;
		state.rho[0] = 1.;
		state.rho[1] = 0.;
		simplex_dual_nonbasic_initialize(&state);
		valid = valid &&
			simplex_dual_prepare_ratio_test(&state, -kappa) == 2 &&
			simplex_dual_next_ratio_candidate(
				&state, &theta, 1, 3, -kappa) == 2 &&
			state.profile.lex_solves > 0;
	}
	for (j = 0; j < 4; j++)
		valid = valid && state.cost[j] == 0.;
	simplex_dual_state_destroy(&state);
	simplex_problem_free(problem);
	return valid;
}


static void initialize_logical_basis(struct simplex_DualState *state)
{
	int j;

	for (j = 0; j < state->variables; j++) {
		state->position[j] =
			j < state->structural ? -1 : j - state->structural;
		state->status[j] = j < state->structural ? DUAL_STATUS_LOWER
							 : DUAL_STATUS_BASIC;
		state->reduced[j] = 1.;
	}
	for (j = 0; j < state->rows; j++)
		state->basis[j] = state->structural + j;
	simplex_dual_nonbasic_initialize(state);
}


static int test_ordered_ratio_sequence(void)
{
	struct simplex_Problem *problem =
		simplex_problem_create_sparse(1, 3, 3);
	struct simplex_DualState state;
	struct lp_simplex_Options options;
	double theta = -1.;
	int j, valid;

	if (problem == NULL)
		return 0;
	for (j = 0; j < 3; j++) {
		problem->matrix.column_start[j] = j;
		problem->matrix.row_index[j] = 0;
		problem->matrix.column_index[j] = j;
		problem->matrix.value[j] = problem->matrix.row_value[j] =
			j == 0 ? 1e-6 : 1.;
	}
	problem->matrix.column_start[3] = 3;
	problem->matrix.row_start[0] = 0;
	problem->matrix.row_start[1] = 3;
	lp_simplex_default_options(&options, lp_simplex_ALGORITHM_DUAL_REVISED);
	if (simplex_dual_state_create(&state, problem, &options,
		    SIMPLEX_BASIS_BACKEND_DEFAULT) != lp_simplex_EXIT_SUCCESS) {
		simplex_problem_free(problem);
		return 0;
	}
	initialize_logical_basis(&state);
	state.rho[0] = -1.;
	state.reduced[0] = 0.;
	state.reduced[1] = .05;
	state.reduced[2] = .075;
	/* Harris prefers the strong later pivot; flips must consume the earlier
	 * weak breakpoint first, then reuse the heap without rebuilding it. */
	valid = simplex_dual_prepare_ratio_test(&state, 1) == 3 &&
		simplex_dual_next_ratio_candidate(&state, &theta, 0, 3, 1) == 1;
	state.ratio_ordered = 1;
	for (j = 0; j < 3; j++) {
		valid = valid && simplex_dual_next_ratio_candidate(
					 &state, &theta, 0, 3, 1) == j;
		valid = valid && theta == state.breakpoint[j];
		simplex_dual_ratio_remove(&state, j);
	}
	valid = valid && simplex_dual_next_ratio_candidate(
				 &state, &theta, 0, 3, 1) == -1;
	state.rho[0] = HUGE_VAL;
	valid = valid &&
		simplex_dual_prepare_ratio_test(&state, 1) ==
			SIMPLEX_DUAL_RATIO_ERROR &&
		state.ratio_failure_status == lp_simplex_PrecisionError &&
		simplex_dual_next_ratio_candidate(&state, &theta, 0, 3, 1) ==
			SIMPLEX_DUAL_RATIO_ERROR;
	simplex_dual_state_destroy(&state);
	simplex_problem_free(problem);
	return valid;
}


static int test_sparse_audit_updates_pricing_row(void)
{
	struct simplex_Problem *problem =
		simplex_problem_create_sparse(20, 5, 68);
	struct simplex_DualState state;
	struct lp_simplex_Options options;
	int i, j, k = 4, valid;

	if (problem == NULL)
		return 0;
	problem->matrix.column_start[0] = 0;
	problem->matrix.row_index[0] = 0;
	problem->matrix.row_index[1] = 2;
	problem->matrix.row_index[2] = 1;
	problem->matrix.row_index[3] = 3;
	problem->matrix.value[0] = -1e16;
	problem->matrix.value[1] = problem->matrix.value[3] = -1.;
	problem->matrix.value[2] = 1e16;
	for (j = 1; j < 5; j++) {
		problem->matrix.column_start[j] = k;
		for (i = 4; i < 20; i++) {
			problem->matrix.row_index[k] = i;
			problem->matrix.value[k++] = 1.;
		}
	}
	problem->matrix.column_start[5] = k;
	k = 0;
	for (i = 0; i < 20; i++) {
		problem->matrix.row_start[i] = k;
		if (i < 4) {
			problem->matrix.column_index[k] = 0;
			problem->matrix.row_value[k++] =
				i == 0 ? -1e16 : (i == 1 ? 1e16 : -1.);
		} else
			for (j = 1; j < 5; j++) {
				problem->matrix.column_index[k] = j;
				problem->matrix.row_value[k++] = 1.;
			}
	}
	problem->matrix.row_start[20] = k;
	lp_simplex_default_options(&options, lp_simplex_ALGORITHM_DUAL_REVISED);
	if (simplex_dual_state_create(&state, problem, &options,
		    SIMPLEX_BASIS_BACKEND_DEFAULT) != lp_simplex_EXIT_SUCCESS) {
		simplex_problem_free(problem);
		return 0;
	}
	initialize_logical_basis(&state);
	for (i = 0; i < 20; i++)
		state.rho[i] = i < 4 ? 1. : 0.;
	/* The row accumulation gives -2 while CSC gives -1.  The audit repairs
	 * packed_alpha; the ratio test must copy that repaired value, too. */
	valid = simplex_dual_prepare_ratio_test(&state, 1) > 0 &&
		state.packed_alpha_valid && state.alpha[0] == -1. &&
		state.breakpoint[0] == 1.;
	simplex_dual_state_destroy(&state);
	simplex_problem_free(problem);
	return valid;
}


int main(void)
{
	if (!test_ordered_ratio_sequence() ||
		!test_sparse_audit_updates_pricing_row() ||
		!test_symbolic_ratio(1, SIMPLEX_BASIS_BACKEND_DEFAULT) ||
		!test_symbolic_ratio(-1, SIMPLEX_BASIS_BACKEND_DEFAULT) ||
		!test_symbolic_ratio(1, SIMPLEX_BASIS_BACKEND_INTERNAL) ||
		!test_symbolic_ratio(-1, SIMPLEX_BASIS_BACKEND_INTERNAL)) {
		printf("dual pricing regression failed\n");
		return 1;
	}
	printf("dual pricing regression passed\n");
	return 0;
}
