/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file test_primal_pricing.c
 * @brief Regression tests for primal Devex initialization and weight updates.
 */

#include "simplex_dual_primal_pricing.h"
#include "simplex_dual_kernel.h"
#include "simplex_dual_state.h"
#include "simplex_dual_pricing.h"
#include "simplex_problem.h"

#include <lp_simplex/status.h>

#include <math.h>
#include <stdio.h>
#include <string.h>

/* Check Devex against explicit rows of B^-1, including the -I logical
 * columns, a changed reference framework, and its inaccuracy restart. */
static int test_reference_update(const int restart)
{
	struct simplex_DualState state;
	int starts[5] = {0, 1, 2, 4, 6};
	int rows[6] = {0, 1, 0, 1, 0, 1};
	double values[6] = {1., 1., 2., 1., 4., 3.};
	int basis[2] = {0, 1};
	int position[6] = {0, 1, -1, -1, -1, -1};
	double direction[2], rho[2], weight[6], ref[6];
	int valid, update;

	memset(&state, 0, sizeof(state));
	state.rows = 2;
	state.structural = 4;
	state.variables = 6;
	state.matrix.rows = 2;
	state.matrix.columns = 4;
	state.matrix.nonzeros = 6;
	state.matrix.column_start = starts;
	state.matrix.row_index = rows;
	state.matrix.value = values;
	state.basis = basis;
	state.position = position;
	state.direction = direction;
	state.rho = rho;
	if (simplex_basis_create(&state.factor, &state.matrix, 4, basis) ==
		lp_simplex_EXIT_FAILURE)
		return 0;
	if (simplex_basis_factorize(&state.factor) == lp_simplex_EXIT_FAILURE) {
		simplex_basis_destroy(&state.factor);
		return 0;
	}
	simplex_primal_devex_initialize(&state, weight, ref);
	direction[0] = 2.;
	direction[1] = 1.;
	valid = simplex_primal_devex_update(&state, weight, ref, 2, 0) ==
			lp_simplex_EXIT_SUCCESS &&
		fabs(weight[3] - 4.) < 1e-12 && weight[0] == 1. &&
		weight[4] == 1.;
	basis[0] = 2;
	position[0] = -1;
	position[2] = 0;
	update = simplex_basis_update(&state.factor, 0, direction);
	if (update == lp_simplex_EXIT_FAILURE ||
		(update == 1 && simplex_basis_factorize(&state.factor) ==
					lp_simplex_EXIT_FAILURE))
		valid = 0;
	direction[0] = 4.;
	direction[1] = 3.;
	if (simplex_basis_ftran_refined(&state.factor, direction) ==
		lp_simplex_EXIT_FAILURE)
		valid = 0;
	if (restart)
		weight[3] = 1e10;
	valid = valid && simplex_primal_devex_update(&state, weight, ref, 3,
				 1) == lp_simplex_EXIT_SUCCESS;
	if (restart)
		valid = valid && weight[0] == 1. && weight[5] == 1. &&
			ref[0] == 1. && ref[2] == 0.;
	else
		valid = valid && fabs(weight[0] - 1.25) < 1e-12 &&
			fabs(weight[4] - 1.25) < 1e-12 &&
			fabs(weight[5] - 5.) < 1e-12 &&
			fabs(weight[1] - 5.) < 1e-12;
	simplex_basis_destroy(&state.factor);
	return valid;
}

/* Start at x0=1, x1=0 for min 2*x0+x1 with x0+x1 >= 1.
 * Cleanup either exchanges x0 for x1 or flips x1 to its finite upper bound.
 * Membership and structural counts must agree on every phase exit. */
static int test_cleanup_indexes(int backend, int bounded, int budget)
{
	struct simplex_Problem *problem =
		simplex_problem_create_sparse(1, 2, 2);
	struct simplex_DualState state;
	struct lp_simplex_Options options;
	int j, result, status = lp_simplex_Success, good;

	if (problem == NULL)
		return 0;
	for (j = 0; j < 2; j++) {
		problem->matrix.column_start[j] = j;
		problem->matrix.row_index[j] = 0;
		problem->matrix.value[j] = 1.;
		problem->matrix.column_index[j] = j;
		problem->matrix.row_value[j] = 1.;
	}
	problem->matrix.column_start[2] = 2;
	problem->matrix.row_start[0] = 0;
	problem->matrix.row_start[1] = 2;
	lp_simplex_default_options(&options, lp_simplex_ALGORITHM_DUAL_REVISED);
	options.iteration_limit = budget;
	if (simplex_dual_state_create(&state, problem, &options, backend) !=
		lp_simplex_EXIT_SUCCESS) {
		simplex_problem_free(problem);
		return 0;
	}
	for (j = 0; j < 3; j++) {
		state.position[j] = j == 0 ? 0 : -1;
		state.status[j] =
			j == 0 ? DUAL_STATUS_BASIC : DUAL_STATUS_LOWER;
		state.lower[j] = j == 2 ? 1. : 0.;
		state.upper[j] = j == 1 && bounded ? .5 : HUGE_VAL;
		state.cost[j] = state.original_cost[j] = 2. - j;
		state.value[j] = j == 1 ? 0. : 1.;
		state.reduced[j] = j == 0 ? 0. : (j == 1 ? -1. : 2.);
	}
	state.basis[0] = 0;
	state.basic_value[0] = 1.;
	state.basic_lower[0] = 0.;
	state.basic_upper[0] = HUGE_VAL;
	state.pi[0] = 2.;
	good = simplex_basis_factorize(&state.factor) ==
	       lp_simplex_EXIT_SUCCESS;
	if (good) {
		result = simplex_dual_primal_cleanup(&state, &status);
		good = simplex_dual_nonbasic_is_consistent(&state) &&
		       state.structural_basic == 1 &&
		       state.candidate_count == 0 && state.nonbasic_count == 2;
		if (budget >= 2)
			good = good && result == lp_simplex_EXIT_SUCCESS &&
			       fabs(state.value[1] - (bounded ? .5 : 1.)) <
				       1e-12 &&
			       fabs(state.value[0] + state.value[1] - 1.) <
				       1e-12;
		else
			good = good && result == lp_simplex_EXIT_FAILURE &&
			       status == lp_simplex_ExceedIterLimit;
	}
	simplex_dual_state_destroy(&state);
	simplex_problem_free(problem);
	return good;
}

int main(void)
{
	if (!test_reference_update(0) || !test_reference_update(1) ||
		!test_cleanup_indexes(SIMPLEX_BASIS_BACKEND_DEFAULT, 0, 2) ||
		!test_cleanup_indexes(SIMPLEX_BASIS_BACKEND_DEFAULT, 1, 2) ||
		!test_cleanup_indexes(SIMPLEX_BASIS_BACKEND_INTERNAL, 0, 2) ||
		!test_cleanup_indexes(SIMPLEX_BASIS_BACKEND_INTERNAL, 1, 2) ||
		!test_cleanup_indexes(SIMPLEX_BASIS_BACKEND_DEFAULT, 0, 0) ||
		!test_cleanup_indexes(SIMPLEX_BASIS_BACKEND_DEFAULT, 0, 1)) {
		printf("primal Devex reference update failed\n");
		return 1;
	}
	printf("primal Devex reference update passed\n");
	return 0;
}
