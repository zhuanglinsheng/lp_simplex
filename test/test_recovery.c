/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file test_recovery.c
 * @brief Regression tests for simplex recovery stages and anti-stalling
 * transitions.
 */

#include "simplex_recovery.h"
#include "simplex_dual_kernel.h"
#include "simplex_dual_internal.h"
#include "simplex_dual_state.h"
#include "simplex_problem.h"
#include "simplex_dual_pricing.h"

#include <lp_simplex/status.h>

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

static int test_policy_priority(void)
{
	struct simplex_Recovery recovery;
	int i, action = -1, valid;

	simplex_recovery_initialize(&recovery);
	valid = simplex_recovery_cost_changes_allowed(&recovery);
	simplex_recovery_begin_primal(&recovery);
	valid = valid && !simplex_recovery_cost_changes_allowed(&recovery);
	for (i = 0; i < SIMPLEX_PRIMAL_STALL_LIMIT; i++)
		action = simplex_recovery_observe_primal(
			&recovery, 1e-20, 1e-10, 0);
	valid = valid && action == SIMPLEX_PRIMAL_PERTURB_BOUNDS;
	recovery.bounds_attempted = 1;
	action = simplex_recovery_observe_primal(&recovery, 1e-20, 1e-10, 0);
	valid = valid && action == SIMPLEX_PRIMAL_BLAND && recovery.bland;
	for (i = 0; i < SIMPLEX_PRIMAL_STALL_LIMIT; i++)
		action = simplex_recovery_observe_primal(
			&recovery, 0., 1e-10, 0);
	valid = valid && action == SIMPLEX_PRIMAL_CONTINUE && !recovery.bland &&
		!recovery.cycle_guard;
	/* A real decrease releases a stagnation guard, but never a cycle guard. */
	action = simplex_recovery_observe_primal(&recovery, 1., 1e-10, 0);
	valid = valid && action == SIMPLEX_PRIMAL_CONTINUE && !recovery.bland;
	action = simplex_recovery_observe_primal(&recovery, 0., 1e-10, 1);
	valid = valid && action == SIMPLEX_PRIMAL_BLAND && recovery.cycle_guard;
	action = simplex_recovery_observe_primal(&recovery, 10., 1e-10, 0);
	valid = valid && action == SIMPLEX_PRIMAL_BLAND;
	action = simplex_recovery_observe_primal(&recovery, 0., 1e-10, 1);
	valid = valid && action == SIMPLEX_PRIMAL_REBUILD &&
		recovery.primal_rebuilds == 1 && recovery.cycle_guard;
	action = simplex_recovery_observe_primal(&recovery, 0., 1e-10, 1);
	valid = valid && action == SIMPLEX_PRIMAL_RECOVERY_FAILED;
	simplex_recovery_begin_primal(&recovery);
	valid = valid && recovery.bland && recovery.bounds_attempted;
	valid = valid && !simplex_recovery_progress_certified(HUGE_VAL, 1.) &&
		!simplex_recovery_progress_certified(1., -1.);
	simplex_degeneracy_destroy(&recovery.primal_history);
	return valid;
}

static int test_bound_snapshot_and_fixed_restore(void)
{
	struct simplex_DualState state;
	struct lp_simplex_Options options;
	int basis[2] = {0, 1}, position[5] = {0, 1, -1, -1, -1};
	unsigned char status[5] = {DUAL_STATUS_BASIC, DUAL_STATUS_BASIC,
		DUAL_STATUS_FIXED, DUAL_STATUS_FREE, DUAL_STATUS_LOWER};
	double lower[5] = {0., 2., 3., -HUGE_VAL, -DBL_MAX};
	double upper[5] = {HUGE_VAL, 2., 3., HUGE_VAL, DBL_MAX};
	double original_lower[5], original_upper[5];
	double value[5] = {0., 2., 3., 0., -DBL_MAX};
	double basic_value[2] = {0., 2.}, basic_lower[2], basic_upper[2];
	int valid;

	memset(&state, 0, sizeof(state));
	lp_simplex_default_options(&options, lp_simplex_ALGORITHM_DUAL_REVISED);
	state.options = &options;
	state.rows = 2;
	state.variables = 5;
	state.basis = basis;
	state.position = position;
	state.status = status;
	state.lower = lower;
	state.upper = upper;
	state.value = value;
	state.basic_value = basic_value;
	state.basic_lower = basic_lower;
	state.basic_upper = basic_upper;
	memcpy(original_lower, lower, sizeof(lower));
	memcpy(original_upper, upper, sizeof(upper));
	simplex_recovery_begin_primal(&state.recovery);
	valid = simplex_recovery_perturb_bounds(&state) == 1 && lower[0] < 0. &&
		lower[1] < 2. && upper[1] > 2. && lower[2] == 3. &&
		upper[2] == 3. && lower[3] == -HUGE_VAL &&
		upper[3] == HUGE_VAL && lower[4] == -DBL_MAX &&
		upper[4] == DBL_MAX && value[0] == 0. && value[1] == 2.;
	valid = valid && simplex_recovery_perturb_bounds(&state) == 0;
	/* A formerly fixed basic variable leaves at a perturbed upper bound. */
	basis[1] = 4;
	position[1] = -1;
	position[4] = 1;
	status[1] = DUAL_STATUS_UPPER;
	status[4] = DUAL_STATUS_BASIC;
	value[1] = upper[1];
	valid = valid && simplex_recovery_restore_bounds(&state) == 1 &&
		memcmp(lower, original_lower, sizeof(lower)) == 0 &&
		memcmp(upper, original_upper, sizeof(upper)) == 0 &&
		status[1] == DUAL_STATUS_FIXED && value[1] == 2. &&
		state.recovery.stage == SIMPLEX_RECOVERY_ORIGINAL_DUAL &&
		!simplex_recovery_cost_changes_allowed(&state.recovery) &&
		state.recovery.bound_epochs == 1 &&
		state.recovery.bound_restores == 1 &&
		state.recovery.saved_bounds == NULL &&
		!state.recovery.bounds_active;
	valid = valid && simplex_recovery_restore_bounds(&state) == 0 &&
		simplex_recovery_perturb_bounds(&state) == 0;
	simplex_recovery_destroy(&state);
	return valid;
}

/* Optimize -x1 subject to x0+x1=0, x>=0 from the primal-feasible x0 basis.
 * Perturbed cleanup can improve its working objective below zero; restoring
 * the bound must return the real solution x=0 before handing control back.
 * A zero pivot budget tests the same rollback on an early failure exit. */
static int test_cleanup_transaction(const int limit, const int backend)
{
	struct simplex_Problem *problem =
		simplex_problem_create_sparse(1, 2, 2);
	struct simplex_DualState state;
	struct lp_simplex_Options options;
	double original_lower[3], original_upper[3], original_cost[3];
	int j, valid, result, status = lp_simplex_CondUnsatisfied;

	if (problem == NULL)
		return 0;
	problem->matrix.column_start[0] = 0;
	problem->matrix.column_start[1] = 1;
	problem->matrix.column_start[2] = 2;
	problem->matrix.row_start[0] = 0;
	problem->matrix.row_start[1] = 2;
	for (j = 0; j < 2; j++) {
		problem->matrix.row_index[j] = 0;
		problem->matrix.column_index[j] = j;
		problem->matrix.value[j] = problem->matrix.row_value[j] = 1.;
	}
	lp_simplex_default_options(&options, lp_simplex_ALGORITHM_DUAL_REVISED);
	options.iteration_limit = limit;
	if (simplex_dual_state_create(&state, problem, &options, backend) ==
		lp_simplex_EXIT_FAILURE) {
		simplex_problem_free(problem);
		return 0;
	}
	for (j = 0; j < 3; j++) {
		state.position[j] = j == 0 ? 0 : -1;
		state.status[j] = j == 0 ? DUAL_STATUS_BASIC
					 : (j == 1 ? DUAL_STATUS_LOWER
						   : DUAL_STATUS_FIXED);
		state.lower[j] = 0.;
		state.upper[j] = j == 2 ? 0. : HUGE_VAL;
		state.value[j] = 0.;
		state.cost[j] = state.original_cost[j] = j == 1 ? -1. : 0.;
		state.reduced[j] = state.cost[j];
		state.feasibility_tolerance[0] = 0.;
	}
	state.basis[0] = 0;
	state.basic_value[0] = state.basic_lower[0] = 0.;
	state.basic_upper[0] = HUGE_VAL;
	state.recovery.force_bounds = 1;
	state.recovery.disable_bounds = 0;
	memcpy(original_lower, state.lower, sizeof(original_lower));
	memcpy(original_upper, state.upper, sizeof(original_upper));
	memcpy(original_cost, state.cost, sizeof(original_cost));
	valid = simplex_basis_factorize(&state.factor) ==
		lp_simplex_EXIT_SUCCESS;
	result = valid ? simplex_dual_primal_cleanup(&state, &status) : -1;
	valid = valid &&
		memcmp(state.lower, original_lower, sizeof(original_lower)) ==
			0 &&
		memcmp(state.upper, original_upper, sizeof(original_upper)) ==
			0 &&
		memcmp(state.cost, original_cost, sizeof(original_cost)) == 0 &&
		state.recovery.bound_epochs == 1 &&
		state.recovery.bound_restores == 1 &&
		state.recovery.bounds_attempted &&
		!state.recovery.bounds_active &&
		!state.recovery.primal_tracking &&
		!state.primal_cleanup_active &&
		state.recovery.original_cost_locked &&
		state.recovery.saved_bounds == NULL &&
		simplex_dual_nonbasic_is_consistent(&state) &&
		fabs(state.value[0]) < 1e-12 && fabs(state.value[1]) < 1e-12;
	if (limit == 0)
		valid = valid && result == lp_simplex_EXIT_FAILURE &&
			status == lp_simplex_ExceedIterLimit &&
			state.iterations == 0;
	else
		valid = valid && result == lp_simplex_EXIT_SUCCESS &&
			state.iterations > 0;
	simplex_dual_state_destroy(&state);
	simplex_problem_free(problem);
	return valid;
}

/* A near-degenerate original vertex needs an actual dual repair after the
 * relaxed optimum is snapped back to its original nonbasic bounds. */
static int test_primal_to_dual_handoff(const int backend, const int exhaust)
{
	const double coefficients[3][6] = {{1., 0., 0., -1., 1., -2.},
		{0., 1., 0., -1., 0., 2.}, {0., 0., 1., -1., -1., -1.}};
	const double rhs[3] = {1., 1., 1. - 3e-7};
	const double costs[6] = {0., 0., 0., -4., -4., -3.};
	struct simplex_Problem *problem =
		simplex_problem_create_sparse(3, 6, 11);
	struct simplex_DualState state;
	struct lp_simplex_Options options;
	double lower[9], upper[9], objective = 0.;
	int i, j, k = 0, valid, result, before_repair;
	int status = lp_simplex_CondUnsatisfied;

	if (problem == NULL)
		return 0;
	for (j = 0; j < 6; j++) {
		problem->matrix.column_start[j] = k;
		for (i = 0; i < 3; i++)
			if (coefficients[i][j] != 0.) {
				problem->matrix.row_index[k] = i;
				problem->matrix.value[k++] = coefficients[i][j];
			}
	}
	problem->matrix.column_start[6] = k;
	k = 0;
	for (i = 0; i < 3; i++) {
		problem->matrix.row_start[i] = k;
		for (j = 0; j < 6; j++)
			if (coefficients[i][j] != 0.) {
				problem->matrix.column_index[k] = j;
				problem->matrix.row_value[k++] =
					coefficients[i][j];
			}
	}
	problem->matrix.row_start[3] = k;
	lp_simplex_default_options(&options, lp_simplex_ALGORITHM_DUAL_REVISED);
	options.iteration_limit = 200;
	if (simplex_dual_state_create(&state, problem, &options, backend) ==
		lp_simplex_EXIT_FAILURE) {
		simplex_problem_free(problem);
		return 0;
	}
	for (i = 0; i < 3; i++) {
		state.basis[i] = i;
		state.pi[i] = 0.;
		state.basic_value[i] = rhs[i];
		state.basic_lower[i] = 0.;
		state.basic_upper[i] = 1.;
		state.feasibility_tolerance[i] = 0.;
	}
	for (j = 0; j < 9; j++) {
		state.position[j] = j < 3 ? j : -1;
		state.status[j] =
			j < 3 ? DUAL_STATUS_BASIC
			      : (j < 6 ? DUAL_STATUS_LOWER : DUAL_STATUS_FIXED);
		state.lower[j] = lower[j] = j < 6 ? 0. : rhs[j - 6];
		state.upper[j] = upper[j] = j < 6 ? 1. : rhs[j - 6];
		state.value[j] = j < 3 ? rhs[j] : (j < 6 ? 0. : rhs[j - 6]);
		state.cost[j] = state.original_cost[j] = j < 6 ? costs[j] : 0.;
		state.reduced[j] = state.cost[j];
	}
	state.recovery.force_bounds = 1;
	state.recovery.disable_bounds = 0;
	valid = simplex_basis_factorize(&state.factor) ==
		lp_simplex_EXIT_SUCCESS;
	result = valid ? simplex_dual_primal_cleanup(&state, &status) : -1;
	valid = valid && result == lp_simplex_EXIT_SUCCESS &&
		state.recovery.dual_repairs == 1 &&
		!state.recovery.bounds_active &&
		state.recovery.stage == SIMPLEX_RECOVERY_ORIGINAL_DUAL &&
		state.recovery.original_cost_locked;
	before_repair = state.iterations;
	if (exhaust)
		options.iteration_limit = before_repair;
	result = valid ? simplex_dual_continue(&state, &status) : -1;
	valid = valid && memcmp(state.lower, lower, sizeof(lower)) == 0 &&
		memcmp(state.upper, upper, sizeof(upper)) == 0 &&
		memcmp(state.cost, state.original_cost, 9 * sizeof(double)) ==
			0 &&
		state.recovery.bound_epochs == 1 &&
		state.recovery.bound_restores == 1;
	if (exhaust)
		valid = valid && result == lp_simplex_EXIT_FAILURE &&
			status == lp_simplex_ExceedIterLimit &&
			state.iterations == before_repair;
	else {
		valid = valid && result == lp_simplex_EXIT_SUCCESS &&
			status == lp_simplex_Success &&
			state.iterations > before_repair;
		for (j = 0; j < 6; j++) {
			objective += costs[j] * state.value[j];
			valid = valid && state.value[j] >= -1e-7 &&
				state.value[j] <= 1. + 1e-7;
		}
		valid = valid && fabs(objective + 1.2e-6) < 1e-10;
		for (i = 0; i < 3; i++) {
			double residual = -rhs[i];

			for (j = 0; j < 6; j++)
				residual += coefficients[i][j] * state.value[j];
			valid = valid && fabs(residual) < 1e-10;
		}
	}
	simplex_dual_state_destroy(&state);
	simplex_problem_free(problem);
	return valid;
}

/* A tolerance-sized wrong-sign cached reduced cost used to change c[q] when
 * a zero ratio was accepted. Original-domain dual repair must preserve c and
 * recompute the real reduced costs after exchanging the basis. */
static int test_sealed_objective_ratio_clip(const int backend)
{
	struct simplex_Problem *problem =
		simplex_problem_create_sparse(1, 1, 1);
	struct simplex_DualState state;
	struct lp_simplex_Options options;
	int valid, result, status = lp_simplex_CondUnsatisfied;

	if (problem == NULL)
		return 0;
	problem->matrix.column_start[0] = problem->matrix.row_start[0] = 0;
	problem->matrix.column_start[1] = problem->matrix.row_start[1] = 1;
	problem->matrix.row_index[0] = problem->matrix.column_index[0] = 0;
	problem->matrix.value[0] = problem->matrix.row_value[0] = 1.;
	lp_simplex_default_options(&options, lp_simplex_ALGORITHM_DUAL_REVISED);
	if (simplex_dual_state_create(&state, problem, &options, backend) ==
		lp_simplex_EXIT_FAILURE) {
		simplex_problem_free(problem);
		return 0;
	}
	state.basis[0] = 1;
	state.position[0] = -1;
	state.position[1] = 0;
	state.status[0] = DUAL_STATUS_LOWER;
	state.status[1] = DUAL_STATUS_BASIC;
	state.lower[0] = 0.;
	state.upper[0] = HUGE_VAL;
	state.lower[1] = state.upper[1] = 1.;
	state.value[0] = state.value[1] = state.basic_value[0] = 0.;
	state.basic_lower[0] = state.basic_upper[0] = 1.;
	state.feasibility_tolerance[0] = state.pi[0] = 0.;
	state.cost[0] = state.cost[1] = 0.;
	state.original_cost[0] = state.original_cost[1] = 0.;
	state.reduced[0] = -.5 * options.dual_tolerance;
	state.reduced[1] = 0.;
	state.edge_weight[0] = 1.;
	simplex_recovery_seal_objective(&state.recovery);
	simplex_dual_nonbasic_initialize(&state);
	simplex_dual_feasibility_rebuild(&state);
	valid = simplex_basis_factorize(&state.factor) ==
		lp_simplex_EXIT_SUCCESS;
	result = valid ? simplex_dual_continue(&state, &status) : -1;
	valid = valid && result == lp_simplex_EXIT_SUCCESS &&
		status == lp_simplex_Success &&
		state.recovery.protected_costs > 0 && state.cost[0] == 0. &&
		state.cost[1] == 0. && !state.cost_shift_active &&
		!state.original_dual_recompute &&
		fabs(state.value[0] - 1.) < 1e-12 && state.pi[0] == 0.;
	simplex_dual_state_destroy(&state);
	simplex_problem_free(problem);
	return valid;
}

int main(void)
{
	if (!test_policy_priority() ||
		!test_bound_snapshot_and_fixed_restore() ||
		!test_cleanup_transaction(0, SIMPLEX_BASIS_BACKEND_DEFAULT) ||
		!test_cleanup_transaction(100, SIMPLEX_BASIS_BACKEND_DEFAULT) ||
		!test_cleanup_transaction(0, SIMPLEX_BASIS_BACKEND_INTERNAL) ||
		!test_cleanup_transaction(
			100, SIMPLEX_BASIS_BACKEND_INTERNAL) ||
		!test_primal_to_dual_handoff(
			SIMPLEX_BASIS_BACKEND_DEFAULT, 0) ||
		!test_primal_to_dual_handoff(
			SIMPLEX_BASIS_BACKEND_DEFAULT, 1) ||
		!test_primal_to_dual_handoff(
			SIMPLEX_BASIS_BACKEND_INTERNAL, 0) ||
		!test_primal_to_dual_handoff(
			SIMPLEX_BASIS_BACKEND_INTERNAL, 1) ||
		!test_sealed_objective_ratio_clip(
			SIMPLEX_BASIS_BACKEND_DEFAULT) ||
		!test_sealed_objective_ratio_clip(
			SIMPLEX_BASIS_BACKEND_INTERNAL)) {
		printf("anti-degeneracy recovery regression failed\n");
		return 1;
	}
	printf("anti-degeneracy recovery regression passed\n");
	return 0;
}
