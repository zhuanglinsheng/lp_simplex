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
#include "simplex_dual_feasibility.h"
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

/* A tolerated negative reduced cost has already spent part of Harris'
 * allowance. Clipping its primary breakpoint to zero must not replenish that
 * allowance. Test both bounds and coefficient scales independently of an LP. */
static int test_harris_remaining_allowance(void)
{
	struct simplex_Problem *problem =
		simplex_problem_create_sparse(1, 2, 2);
	struct simplex_DualState state;
	struct lp_simplex_Options options;
	const double scales[] = {1e-3, 1., 1e3};
	const double tolerances[] = {1e-5, 1e-7, 1e-10};
	int i, j, bound, valid = 1;
	double theta;
	if (problem == NULL)
		return 0;
	problem->matrix.column_start[0] = 0;
	problem->matrix.column_start[1] = 1;
	problem->matrix.column_start[2] = 2;
	problem->matrix.row_index[0] = problem->matrix.row_index[1] = 0;
	problem->matrix.row_start[0] = 0;
	problem->matrix.row_start[1] = 2;
	problem->matrix.column_index[0] = 0;
	problem->matrix.column_index[1] = 1;
	lp_simplex_default_options(&options, lp_simplex_ALGORITHM_DUAL_REVISED);
	if (simplex_dual_state_create(&state, problem, &options,
		    SIMPLEX_BASIS_BACKEND_DEFAULT) != lp_simplex_EXIT_SUCCESS) {
		simplex_problem_free(problem);
		return 0;
	}
	initialize_logical_basis(&state);
	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++)
			for (bound = 0; bound < 2; bound++) {
				int sign = bound ? -1 : 1;
				options.dual_tolerance = tolerances[j];
				state.rho[0] = -scales[i];
				state.matrix.value[0] =
					state.matrix.row_value[0] = 1.;
				state.matrix.value[1] =
					state.matrix.row_value[1] = 2.;
				state.status[0] = state.status[1] =
					bound ? DUAL_STATUS_UPPER
					      : DUAL_STATUS_LOWER;
				state.reduced[0] = -sign * .9 * tolerances[j];
				state.reduced[1] = sign * tolerances[j];
				valid = valid &&
					simplex_dual_prepare_ratio_test(
						&state, sign) == 2 &&
					simplex_dual_next_ratio_candidate(
						&state, &theta, 0, 2, sign) ==
						0 &&
					theta == 0.;
			}
	simplex_dual_state_destroy(&state);
	simplex_problem_free(problem);
	return valid;
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

static int test_wide_dense_pricing(void)
{
	struct simplex_Problem *problem =
		simplex_problem_create_sparse(1, 40, 40);
	struct simplex_DualState state;
	struct lp_simplex_Options options;
	double theta;
	int j, valid;

	if (!problem)
		return 0;
	for (j = 0; j < 40; j++) {
		problem->matrix.column_start[j] = j;
		problem->matrix.row_index[j] = 0;
		problem->matrix.column_index[j] = j;
		problem->matrix.value[j] = problem->matrix.row_value[j] = 1.;
	}
	problem->matrix.column_start[40] = 40;
	problem->matrix.row_start[0] = 0;
	problem->matrix.row_start[1] = 40;
	lp_simplex_default_options(&options, lp_simplex_ALGORITHM_DUAL_REVISED);
	if (simplex_dual_state_create(
		    &state, problem, &options, SIMPLEX_BASIS_BACKEND_DEFAULT)) {
		simplex_problem_free(problem);
		return 0;
	}
	initialize_logical_basis(&state);
	state.rho[0] = -1.;
	/* Swap deletion scrambles the list; fixed columns must still be skipped,
	 * and both Harris and heap ties must choose the smallest admissible id. */
	state.status[0] = DUAL_STATUS_FIXED;
	simplex_dual_nonbasic_remove(&state, 0);
	valid = simplex_dual_prepare_ratio_test(&state, 1) == 39 &&
		!state.packed_alpha_valid &&
		simplex_dual_next_ratio_candidate(&state, &theta, 0, 40, 1) ==
			1 &&
		theta == 1.;
	state.ratio_ordered = 1;
	for (j = 1; j < 40; j++) {
		valid = valid && simplex_dual_next_ratio_candidate(
					 &state, &theta, 0, 40, 1) == j;
		simplex_dual_ratio_remove(&state, j);
	}
	/* Removing enough columns selects the compact-list path. */
	for (j = 2; j < 40; j++) {
		state.status[j] = DUAL_STATUS_FIXED;
		simplex_dual_nonbasic_remove(&state, j);
	}
	valid = valid && simplex_dual_prepare_ratio_test(&state, 1) == 1 &&
		simplex_dual_next_ratio_candidate(&state, &theta, 0, 40, 1) ==
			1;
	simplex_dual_state_destroy(&state);
	simplex_problem_free(problem);
	return valid;
}

/* Compare batched scatter to scalar operations across cancellation and reuse.
 * Disable the CSC audit: its different summation order is tested separately. */
static int test_sparse_scatter_generations(struct simplex_DualState *state)
{
	struct simplex_SparseVector reference;
	struct simplex_SparseVector *alpha = &state->packed_alpha;
	const struct simplex_CscMatrix *matrix = &state->matrix;
	int pass, i, j, k, valid = 1;

	if (simplex_sparse_vector_create(&reference, state->variables) !=
		lp_simplex_EXIT_SUCCESS)
		return 0;
	for (pass = 0; pass < 4; pass++) {
		/* Dense writes must participate in the shared clearing contract. */
		for (i = 0; i < state->rows; i++)
			state->rho[i] = 1.;
		simplex_dual_prepare_ratio_test(state, 1);
		valid = valid && !state->packed_alpha_valid;
		lp_simplex_memset(
			state->rho, 0, (size_t)state->rows * sizeof(double));
		simplex_sparse_vector_clear(&reference);
		for (i = 0; i < 4; i++) {
			double rho = (double)(i + 1);

			if (pass == 0)
				rho = 1.;
			else if (i < 2)
				rho = (double)(pass + 1);
			else if (pass == i)
				rho = 0.;
			state->rho[i] = rho;
			if (rho == 0.)
				continue;
			for (k = matrix->row_start[i];
				k < matrix->row_start[i + 1]; k++)
				simplex_sparse_vector_add(&reference,
					matrix->column_index[k],
					rho * matrix->row_value[k]);
			simplex_sparse_vector_set(
				&reference, state->structural + i, -rho);
		}
		state->pricing_validation_countdown = 127;
		simplex_sparse_vector_compact(&reference);
		simplex_dual_prepare_ratio_test(state, 1);
		valid = valid && state->packed_alpha_valid &&
			alpha->count == reference.count;
		/* Untouched coordinates must also be zero: dense pricing shares this
		 * storage, and stale unregistered entries evade sparse clear. */
		for (j = 0; j < state->variables; j++)
			valid = valid && alpha->dense[j] == reference.dense[j];
		for (j = 0; j < reference.count; j++) {
			int index = reference.index[j];

			valid = valid && alpha->index[j] == index &&
				alpha->slot[index] == j &&
				alpha->value[j] == reference.value[j] &&
				state->alpha[index] == reference.value[j] &&
				simplex_sparse_vector_get(alpha, index) ==
					reference.value[j];
		}
		/* Consuming the entering coefficient must not affect the next row. */
		if (reference.count)
			state->alpha[reference.index[0]] = 0.;
	}
	simplex_sparse_vector_destroy(&reference);
	return valid;
}

static int test_sparse_audit_updates_pricing_row(int cancellation)
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
	if (cancellation)
		problem->matrix.value[1] = 1.;
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
				i == 0 ? -1e16
				       : (i == 1 ? 1e16
						 : (cancellation && i == 2
								   ? 1.
								   : -1.));
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
	/* CSR gives -2, or exact zero in the cancellation case; CSC gives -1.
	 * The full audit must repair or reinsert that column after compaction. */
	valid = simplex_dual_prepare_ratio_test(&state, 1) > 0 &&
		state.alpha == state.packed_alpha.dense &&
		state.packed_alpha_valid && state.alpha[0] == -1. &&
		state.breakpoint[0] == 1.;
	if (cancellation) {
		state.pricing_validation_countdown = 127;
		valid = valid &&
			simplex_dual_prepare_ratio_test(&state, 1) == 1 &&
			!state.packed_alpha_valid && state.alpha[0] == -1.;
		state.pricing_validation_countdown = 0;
		valid = valid &&
			simplex_dual_prepare_ratio_test(&state, 1) == 1 &&
			state.packed_alpha_valid && state.alpha[0] == -1. &&
			simplex_sparse_vector_get(&state.packed_alpha, 0) ==
				-1.;
	} else
		valid = valid && test_sparse_scatter_generations(&state);
	simplex_dual_state_destroy(&state);
	simplex_problem_free(problem);
	return valid;
}

static int test_dse_reprice(const int upper)
{
	struct simplex_Problem *p = simplex_problem_create_sparse(2, 1, 2);
	struct simplex_DualState s;
	struct lp_simplex_Options o;
	double target, maximum;
	int kappa, valid, row;
	if (!p)
		return 0;
	p->matrix.column_start[0] = 0;
	p->matrix.column_start[1] = 2;
	p->matrix.row_index[0] = 0;
	p->matrix.row_index[1] = 1;
	p->matrix.value[0] = p->matrix.value[1] = 1.;
	p->matrix.row_start[0] = 0;
	p->matrix.row_start[1] = 1;
	p->matrix.row_start[2] = 2;
	p->matrix.column_index[0] = p->matrix.column_index[1] = 0;
	p->matrix.row_value[0] = p->matrix.row_value[1] = 1.;
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_DUAL_REVISED);
	if (simplex_dual_state_create(
		    &s, p, &o, SIMPLEX_BASIS_BACKEND_INTERNAL)) {
		simplex_problem_free(p);
		return 0;
	}
	for (row = 0; row < 2; row++) {
		s.basis[row] = 1 + row;
		s.basic_lower[row] = 0.;
		s.basic_upper[row] = 10.;
		s.basic_value[row] =
			upper ? 10. + (row ? 4. : 10.) : (row ? -4. : -10.);
		s.edge_weight[row] = 1.;
		s.feasibility_tolerance[row] = 0.;
	}
	simplex_dual_feasibility_rebuild(&s);
	valid = simplex_dual_feasibility_choose(
			&s, &target, &kappa, &maximum, 0, 0, -1) == 0;
	valid = valid && simplex_dual_feasibility_reprice(&s, 0, 100.) == 1;
	valid = valid && simplex_dual_feasibility_choose(
				 &s, &target, &kappa, &maximum, 0, 0, -1) == 1;
	/* Corrected selection follows true scores 10^2/100 < 4^2/1,
	 * without changing the primal point or feasibility tolerances. */
	valid = valid && s.basic_value[0] == (upper ? 20. : -10.) &&
		s.basic_value[1] == (upper ? 14. : -4.) &&
		simplex_dual_feasibility_merit(&s) == 14.;
	s.edge_weight[0] = 25.;
	valid = valid && !simplex_dual_feasibility_reprice(&s, 0, 100.);
	valid = valid &&
		simplex_dual_feasibility_reprice(&s, 0, HUGE_VAL) == -1 &&
		s.edge_weight[0] == 100.;
	valid = valid && simplex_dual_feasibility_reprice(&s, 0, 0.) == -1 &&
		s.edge_weight[0] == 100.;
	simplex_dual_state_destroy(&s);
	simplex_problem_free(p);
	return valid;
}

int main(void)
{
	if (!test_dse_reprice(0) || !test_dse_reprice(1) ||
		!test_harris_remaining_allowance() ||
		!test_ordered_ratio_sequence() || !test_wide_dense_pricing() ||
		!test_sparse_audit_updates_pricing_row(0) ||
		!test_sparse_audit_updates_pricing_row(1) ||
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
