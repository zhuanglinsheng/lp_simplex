/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file test_presolve.c
 * @brief Regression tests for presolve reductions and original-variable
 * reconstruction.
 */

#include <lp_simplex/lp_simplex.h>

#include "simplex_presolve_queue.h"
#include "simplex_presolve.h"
#include "simplex_dual.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

static double absolute_value(const double value)
{
	return value < 0. ? -value : value;
}

static int exercise_presolve_work_queue(void)
{
	struct simplex_PresolveQueue queue;
	int item;

	if (simplex_presolve_queue_init(&queue, 3) == lp_simplex_EXIT_FAILURE)
		return 0;
	if (simplex_presolve_queue_push(&queue, 0) == lp_simplex_EXIT_FAILURE ||
		simplex_presolve_queue_push(&queue, 1) ==
			lp_simplex_EXIT_FAILURE ||
		simplex_presolve_queue_push(&queue, 2) ==
			lp_simplex_EXIT_FAILURE ||
		simplex_presolve_queue_push(&queue, 1) ==
			lp_simplex_EXIT_FAILURE ||
		!simplex_presolve_queue_pop(&queue, &item) || item != 0 ||
		simplex_presolve_queue_push(&queue, 0) ==
			lp_simplex_EXIT_FAILURE ||
		!simplex_presolve_queue_pop(&queue, &item) || item != 1 ||
		!simplex_presolve_queue_pop(&queue, &item) || item != 2 ||
		!simplex_presolve_queue_pop(&queue, &item) || item != 0 ||
		simplex_presolve_queue_pop(&queue, &item)) {
		simplex_presolve_queue_destroy(&queue);
		return 0;
	}
	simplex_presolve_queue_destroy(&queue);
	return 1;
}

static int solve_fixed_and_empty_column(void)
{
	struct lp_Model *model = lp_model_create(1, 3);
	struct lp_simplex_Options options;
	struct lp_simplex_Result result;
	double x[3];
	int state;

	if (model == NULL)
		return 0;
	/* x0 is fixed, x2 is an empty bounded column selected by its cost. */
	model->objective[1] = 1.;
	model->objective[2] = -1.;
	model->bounds[0].b_type = optm_BOUND_T_BS;
	model->bounds[0].lb = model->bounds[0].ub = 2.;
	model->bounds[2].b_type = optm_BOUND_T_BS;
	model->bounds[2].lb = 0.;
	model->bounds[2].ub = 4.;
	model->constraints[0].type = optm_CONS_T_GE;
	model->constraints[0].rhs = 5.;
	model->constraints[0].coef[0] = 1.;
	model->constraints[0].coef[1] = 1.;
	lp_model_build_sparse(model);
	lp_simplex_default_options(&options, lp_simplex_ALGORITHM_DUAL_REVISED);
	state = lp_simplex_solve(model, &options, x, &result);
	lp_model_free(model);
	return state == lp_simplex_EXIT_SUCCESS &&
	       result.status == lp_simplex_Success &&
	       absolute_value(x[0] - 2.) <= 1e-9 &&
	       absolute_value(x[1] - 3.) <= 1e-9 &&
	       absolute_value(x[2] - 4.) <= 1e-9 &&
	       absolute_value(result.objective + 1.) <= 1e-9;
}

static int solve_with_redundant_empty_row(void)
{
	struct lp_Model *model = lp_model_create(2, 1);
	struct lp_simplex_Options options;
	struct lp_simplex_Result result;
	double x[1];
	int state;

	if (model == NULL)
		return 0;
	model->objective[0] = 1.;
	model->constraints[0].type = optm_CONS_T_GE;
	model->constraints[0].rhs = 1.;
	model->constraints[0].coef[0] = 1.;
	model->constraints[1].type = optm_CONS_T_LE;
	model->constraints[1].rhs = 5.;
	lp_model_build_sparse(model);
	lp_simplex_default_options(&options, lp_simplex_ALGORITHM_DUAL_REVISED);
	state = lp_simplex_solve(model, &options, x, &result);
	lp_model_free(model);
	return state == lp_simplex_EXIT_SUCCESS &&
	       result.status == lp_simplex_Success &&
	       absolute_value(x[0] - 1.) <= 1e-9;
}

static int detect_infeasible_empty_row(void)
{
	struct lp_Model *model = lp_model_create(2, 1);
	struct lp_simplex_Options options;
	struct lp_simplex_Result result;
	double x[1];
	int state;

	if (model == NULL)
		return 0;
	model->constraints[0].type = optm_CONS_T_GE;
	model->constraints[0].rhs = 0.;
	model->constraints[0].coef[0] = 1.;
	model->constraints[1].type = optm_CONS_T_GE;
	model->constraints[1].rhs = 1.;
	lp_model_build_sparse(model);
	lp_simplex_default_options(&options, lp_simplex_ALGORITHM_DUAL_REVISED);
	state = lp_simplex_solve(model, &options, x, &result);
	lp_model_free(model);
	return state == lp_simplex_EXIT_FAILURE &&
	       result.status == lp_simplex_Infeasibility;
}

static int solve_singleton_inequality_chain(void)
{
	struct lp_Model *model = lp_model_create(3, 2);
	struct lp_simplex_Options options;
	struct lp_simplex_Result result;
	double x[2];
	int state;

	if (model == NULL)
		return 0;
	/* x0 >= 2 tightens its lower bound.  The second singleton tightens x1;
	 * the final row is then redundant from the two implied bounds. */
	model->objective[0] = 1.;
	model->objective[1] = 1.;
	model->constraints[0].type = optm_CONS_T_GE;
	model->constraints[0].rhs = 2.;
	model->constraints[0].coef[0] = 1.;
	model->constraints[1].type = optm_CONS_T_LE;
	model->constraints[1].rhs = -3.;
	model->constraints[1].coef[1] = -1.;
	model->constraints[2].type = optm_CONS_T_GE;
	model->constraints[2].rhs = 4.;
	model->constraints[2].coef[0] = 1.;
	model->constraints[2].coef[1] = 1.;
	lp_model_build_sparse(model);
	lp_simplex_default_options(&options, lp_simplex_ALGORITHM_DUAL_REVISED);
	state = lp_simplex_solve(model, &options, x, &result);
	lp_model_free(model);
	return state == lp_simplex_EXIT_SUCCESS &&
	       result.status == lp_simplex_Success &&
	       absolute_value(x[0] - 2.) <= 1e-9 &&
	       absolute_value(x[1] - 3.) <= 1e-9 &&
	       absolute_value(result.objective - 5.) <= 1e-9;
}

static int detect_infeasible_activity(void)
{
	struct lp_Model *model = lp_model_create(1, 2);
	struct lp_simplex_Options options;
	struct lp_simplex_Result result;
	double x[2];
	int state;

	if (model == NULL)
		return 0;
	model->bounds[0].b_type = optm_BOUND_T_BS;
	model->bounds[0].lb = 0.;
	model->bounds[0].ub = 1.;
	model->bounds[1].b_type = optm_BOUND_T_BS;
	model->bounds[1].lb = 0.;
	model->bounds[1].ub = 1.;
	model->constraints[0].type = optm_CONS_T_GE;
	model->constraints[0].rhs = 3.;
	model->constraints[0].coef[0] = 1.;
	model->constraints[0].coef[1] = 1.;
	lp_model_build_sparse(model);
	lp_simplex_default_options(&options, lp_simplex_ALGORITHM_DUAL_REVISED);
	state = lp_simplex_solve(model, &options, x, &result);
	lp_model_free(model);
	return state == lp_simplex_EXIT_FAILURE &&
	       result.status == lp_simplex_Infeasibility;
}

static int solve_implied_bound_propagation(void)
{
	struct lp_Model *model = lp_model_create(3, 2);
	struct lp_simplex_Options options;
	struct lp_simplex_Result result;
	double x[2];
	int state;

	if (model == NULL)
		return 0;
	/* y <= 3 and x + y >= 10 imply x >= 7.  Together with x <= 7,
	 * propagation fixes both variables before the simplex kernel. */
	model->objective[0] = 1.;
	model->objective[1] = 1.;
	model->constraints[0].type = optm_CONS_T_GE;
	model->constraints[0].rhs = 10.;
	model->constraints[0].coef[0] = 1.;
	model->constraints[0].coef[1] = 1.;
	model->constraints[1].type = optm_CONS_T_LE;
	model->constraints[1].rhs = 3.;
	model->constraints[1].coef[1] = 1.;
	model->constraints[2].type = optm_CONS_T_LE;
	model->constraints[2].rhs = 7.;
	model->constraints[2].coef[0] = 1.;
	lp_model_build_sparse(model);
	lp_simplex_default_options(&options, lp_simplex_ALGORITHM_DUAL_REVISED);
	state = lp_simplex_solve(model, &options, x, &result);
	lp_model_free(model);
	return state == lp_simplex_EXIT_SUCCESS &&
	       result.status == lp_simplex_Success &&
	       absolute_value(x[0] - 7.) <= 1e-6 &&
	       absolute_value(x[1] - 3.) <= 1e-6 &&
	       absolute_value(result.objective - 10.) <= 1e-6;
}

static int solve_dominated_parallel_row(void)
{
	struct lp_Model *model = lp_model_create(2, 2);
	struct lp_simplex_Options options;
	struct lp_simplex_Result result;
	double x[2];
	int state;

	if (model == NULL)
		return 0;
	model->objective[0] = 1.;
	model->objective[1] = 1.;
	model->constraints[0].type = optm_CONS_T_GE;
	model->constraints[0].rhs = 2.;
	model->constraints[0].coef[0] = 1.;
	model->constraints[0].coef[1] = 1.;
	model->constraints[1].type = optm_CONS_T_GE;
	model->constraints[1].rhs = 6.;
	model->constraints[1].coef[0] = 2.;
	model->constraints[1].coef[1] = 2.;
	lp_model_build_sparse(model);
	lp_simplex_default_options(&options, lp_simplex_ALGORITHM_DUAL_REVISED);
	state = lp_simplex_solve(model, &options, x, &result);
	lp_model_free(model);
	return state == lp_simplex_EXIT_SUCCESS &&
	       result.status == lp_simplex_Success &&
	       absolute_value(result.objective - 3.) <= 1e-7;
}

static int detect_conflicting_parallel_equalities(void)
{
	struct lp_Model *model = lp_model_create(2, 2);
	struct lp_simplex_Options options;
	struct lp_simplex_Result result;
	double x[2];
	int state;

	if (model == NULL)
		return 0;
	model->constraints[0].type = optm_CONS_T_EQ;
	model->constraints[0].rhs = 1.;
	model->constraints[0].coef[0] = 1.;
	model->constraints[0].coef[1] = 1.;
	model->constraints[1].type = optm_CONS_T_EQ;
	model->constraints[1].rhs = 3.;
	model->constraints[1].coef[0] = 2.;
	model->constraints[1].coef[1] = 2.;
	lp_model_build_sparse(model);
	lp_simplex_default_options(&options, lp_simplex_ALGORITHM_DUAL_REVISED);
	state = lp_simplex_solve(model, &options, x, &result);
	lp_model_free(model);
	return state == lp_simplex_EXIT_FAILURE &&
	       result.status == lp_simplex_Infeasibility;
}

static int solve_singleton_column_substitution(void)
{
	struct lp_Model *model = lp_model_create(2, 3);
	struct lp_simplex_Options options;
	struct lp_simplex_Result result;
	double x[3];
	int state;

	if (model == NULL)
		return 0;
	/* x0 occurs only in x0 + x1 = 5.  Since 0 <= x1 <= 5 implies x0 >= 0,
	 * eliminate x0 and scatter its objective through the equality. */
	model->objective[0] = 1.;
	model->bounds[1].b_type = optm_BOUND_T_BS;
	model->bounds[1].lb = 0.;
	model->bounds[1].ub = 5.;
	model->constraints[0].type = optm_CONS_T_EQ;
	model->constraints[0].rhs = 5.;
	model->constraints[0].coef[0] = 1.;
	model->constraints[0].coef[1] = 1.;
	model->constraints[1].type = optm_CONS_T_GE;
	model->constraints[1].rhs = 0.;
	model->constraints[1].coef[2] = 1.;
	lp_model_build_sparse(model);
	lp_simplex_default_options(&options, lp_simplex_ALGORITHM_DUAL_REVISED);
	state = lp_simplex_solve(model, &options, x, &result);
	lp_model_free(model);
	return state == lp_simplex_EXIT_SUCCESS &&
	       result.status == lp_simplex_Success &&
	       absolute_value(x[0]) <= 1e-8 &&
	       absolute_value(x[1] - 5.) <= 1e-8 &&
	       absolute_value(x[2]) <= 1e-8 &&
	       absolute_value(result.objective) <= 1e-8;
}

static int solve_forcing_row(void)
{
	struct lp_Model *model = lp_model_create(1, 2);
	struct lp_simplex_Options options;
	struct lp_simplex_Result result;
	double x[2];
	int state;

	if (model == NULL)
		return 0;
	/* The minimum activity of x - y is -2.  Requiring x - y <= -2
	 * therefore forces both variables to the bounds attaining that value. */
	model->objective[0] = 1.;
	model->objective[1] = -1.;
	model->bounds[0].b_type = optm_BOUND_T_BS;
	model->bounds[0].lb = 0.;
	model->bounds[0].ub = 5.;
	model->bounds[1].b_type = optm_BOUND_T_BS;
	model->bounds[1].lb = 0.;
	model->bounds[1].ub = 2.;
	model->constraints[0].type = optm_CONS_T_LE;
	model->constraints[0].rhs = -2.;
	model->constraints[0].coef[0] = 1.;
	model->constraints[0].coef[1] = -1.;
	lp_model_build_sparse(model);
	lp_simplex_default_options(&options, lp_simplex_ALGORITHM_DUAL_REVISED);
	state = lp_simplex_solve(model, &options, x, &result);
	lp_model_free(model);
	return state == lp_simplex_EXIT_SUCCESS &&
	       result.status == lp_simplex_Success &&
	       absolute_value(x[0]) <= 1e-9 &&
	       absolute_value(x[1] - 2.) <= 1e-9 &&
	       absolute_value(result.objective + 2.) <= 1e-9;
}

static int solve_zero_cost_singleton_inequality_column(void)
{
	struct lp_Model *model = lp_model_create(2, 2);
	struct lp_simplex_Options options;
	struct lp_simplex_Result result;
	double x[2];
	int state;

	if (model == NULL)
		return 0;
	/* x0 appears only in x0 + x1 <= 5 and has zero cost.  Its lower bound
	 * relaxes the row, so x0 can be fixed there before solving for x1. */
	model->objective[1] = -1.;
	model->bounds[1].b_type = optm_BOUND_T_BS;
	model->bounds[1].lb = 0.;
	model->bounds[1].ub = 10.;
	model->constraints[0].type = optm_CONS_T_LE;
	model->constraints[0].rhs = 5.;
	model->constraints[0].coef[0] = 1.;
	model->constraints[0].coef[1] = 1.;
	model->constraints[1].type = optm_CONS_T_GE;
	model->constraints[1].rhs = 0.;
	model->constraints[1].coef[1] = 1.;
	lp_model_build_sparse(model);
	lp_simplex_default_options(&options, lp_simplex_ALGORITHM_DUAL_REVISED);
	state = lp_simplex_solve(model, &options, x, &result);
	lp_model_free(model);
	return state == lp_simplex_EXIT_SUCCESS &&
	       result.status == lp_simplex_Success &&
	       absolute_value(x[0]) <= 1e-9 &&
	       absolute_value(x[1] - 5.) <= 1e-9 &&
	       absolute_value(result.objective + 5.) <= 1e-9;
}

static int solve_substitution_fixed_point_cascade(void)
{
	struct lp_Model *model = lp_model_create(2, 3);
	struct lp_simplex_Options options;
	struct lp_simplex_Result result;
	double x[3];
	int state, j;

	if (model == NULL)
		return 0;
	for (j = 0; j < 3; j++) {
		model->bounds[j].b_type = optm_BOUND_T_BS;
		model->bounds[j].lb = 0.;
		model->bounds[j].ub = 3.;
	}
	model->objective[1] = 1.;
	/* Eliminating x0 from the doubleton leaves x2 = 1.  The resulting fixed
	 * column and the empty x1 column must be closed before entering simplex. */
	model->constraints[0].type = optm_CONS_T_EQ;
	model->constraints[0].rhs = 2.;
	model->constraints[0].coef[0] = 1.;
	model->constraints[0].coef[1] = 1.;
	model->constraints[1].type = optm_CONS_T_EQ;
	model->constraints[1].rhs = 3.;
	model->constraints[1].coef[0] = 1.;
	model->constraints[1].coef[1] = 1.;
	model->constraints[1].coef[2] = 1.;
	lp_model_build_sparse(model);
	lp_simplex_default_options(&options, lp_simplex_ALGORITHM_DUAL_REVISED);
	state = lp_simplex_solve(model, &options, x, &result);
	lp_model_free(model);
	return state == lp_simplex_EXIT_SUCCESS &&
	       result.status == lp_simplex_Success &&
	       absolute_value(x[0] - 2.) <= 1e-9 &&
	       absolute_value(x[1]) <= 1e-9 &&
	       absolute_value(x[2] - 1.) <= 1e-9 &&
	       absolute_value(result.objective) <= 1e-9;
}

static int solve_singleton_equality_projection(void)
{
	struct lp_Model *model = lp_model_create(2, 4);
	struct lp_simplex_Options options;
	struct lp_simplex_Result result;
	double x[4];
	int state;

	if (model == NULL)
		return 0;
	/* x0 = 5 - x1 - x2 and x0 >= 0 project to x1 + x2 <= 5.  Separate
	 * implied bounds x1 <= 5 and x2 <= 5 cannot replace this aggregate row. */
	model->objective[1] = -1.;
	model->objective[2] = -1.;
	model->bounds[1].b_type = optm_BOUND_T_BS;
	model->bounds[1].lb = 0.;
	model->bounds[1].ub = 10.;
	model->bounds[2].b_type = optm_BOUND_T_BS;
	model->bounds[2].lb = 0.;
	model->bounds[2].ub = 10.;
	model->constraints[0].type = optm_CONS_T_EQ;
	model->constraints[0].rhs = 5.;
	model->constraints[0].coef[0] = 1.;
	model->constraints[0].coef[1] = 1.;
	model->constraints[0].coef[2] = 1.;
	model->constraints[1].type = optm_CONS_T_GE;
	model->constraints[1].rhs = 0.;
	model->constraints[1].coef[3] = 1.;
	lp_model_build_sparse(model);
	lp_simplex_default_options(&options, lp_simplex_ALGORITHM_DUAL_REVISED);
	state = lp_simplex_solve(model, &options, x, &result);
	lp_model_free(model);
	return state == lp_simplex_EXIT_SUCCESS &&
	       result.status == lp_simplex_Success &&
	       absolute_value(x[0]) <= 1e-9 &&
	       absolute_value(x[1] + x[2] - 5.) <= 1e-9 &&
	       absolute_value(result.objective + 5.) <= 1e-9;
}

static int solve_degree_two_implied_free_substitution(void)
{
	struct lp_Model *model = lp_model_create(2, 3);
	struct lp_simplex_Options options;
	struct lp_simplex_Result result;
	double x[3];
	int state;

	if (model == NULL)
		return 0;
	/* x0 has degree two.  The first equality and 0 <= x1 <= 5 imply the
	 * lower bound of x0, so eliminating x0 aggregates the two rows without
	 * increasing the matrix nonzero count. */
	model->objective[1] = 1.;
	model->bounds[1].b_type = optm_BOUND_T_BS;
	model->bounds[1].lb = 0.;
	model->bounds[1].ub = 5.;
	model->constraints[0].type = optm_CONS_T_EQ;
	model->constraints[0].rhs = 5.;
	model->constraints[0].coef[0] = 1.;
	model->constraints[0].coef[1] = 1.;
	model->constraints[1].type = optm_CONS_T_EQ;
	model->constraints[1].rhs = 4.;
	model->constraints[1].coef[0] = 1.;
	model->constraints[1].coef[2] = 1.;
	lp_model_build_sparse(model);
	lp_simplex_default_options(&options, lp_simplex_ALGORITHM_DUAL_REVISED);
	state = lp_simplex_solve(model, &options, x, &result);
	lp_model_free(model);
	return state == lp_simplex_EXIT_SUCCESS &&
	       result.status == lp_simplex_Success &&
	       absolute_value(x[0] - 4.) <= 1e-9 &&
	       absolute_value(x[1] - 1.) <= 1e-9 &&
	       absolute_value(x[2]) <= 1e-9 &&
	       absolute_value(result.objective - 1.) <= 1e-9;
}

static int solve_with_objective_offset(void)
{
	struct lp_Model *model = lp_model_create(1, 1);
	struct lp_simplex_Options options;
	struct lp_simplex_Result result;
	double x[1];
	int state;

	if (model == NULL)
		return 0;
	model->objective[0] = 2.;
	model->objective_offset = -7.;
	model->constraints[0].type = optm_CONS_T_GE;
	model->constraints[0].rhs = 3.;
	model->constraints[0].coef[0] = 1.;
	lp_model_build_sparse(model);
	lp_simplex_default_options(&options, lp_simplex_ALGORITHM_DUAL_REVISED);
	state = lp_simplex_solve(model, &options, x, &result);
	lp_model_free(model);
	return state == lp_simplex_EXIT_SUCCESS &&
	       result.status == lp_simplex_Success &&
	       absolute_value(x[0] - 3.) <= 1e-9 &&
	       absolute_value(result.objective + 1.) <= 1e-9;
}

static int solve_dense_presolve_equivalence(void)
{
	struct lp_Model *model = lp_model_create(4, 4);
	struct lp_simplex_Options options;
	struct lp_simplex_Result with_presolve, without_presolve;
	double reduced_x[4], raw_x[4];
	int first, second, j;

	if (model == NULL)
		return 0;
	for (j = 0; j < 4; j++) {
		model->bounds[j].b_type = optm_BOUND_T_BS;
		model->bounds[j].lb = 0.;
		model->bounds[j].ub = 10.;
	}
	model->objective[0] = 2.;
	model->objective[1] = -1.;
	model->objective[2] = .5;
	model->objective[3] = 1.;
	model->constraints[0].type = optm_CONS_T_GE;
	model->constraints[0].rhs = 3.;
	model->constraints[0].coef[0] = 1.;
	model->constraints[0].coef[1] = 1.;
	model->constraints[1].type = optm_CONS_T_LE;
	model->constraints[1].rhs = 7.;
	model->constraints[1].coef[1] = 1.;
	model->constraints[1].coef[2] = 1.;
	model->constraints[2].type = optm_CONS_T_EQ;
	model->constraints[2].rhs = 1.;
	model->constraints[2].coef[0] = 1.;
	model->constraints[2].coef[2] = -1.;
	model->constraints[3].type = optm_CONS_T_GE;
	model->constraints[3].rhs = 2.;
	model->constraints[3].coef[3] = 1.;
	lp_simplex_default_options(&options, lp_simplex_ALGORITHM_DUAL_REVISED);
	first = lp_simplex_solve(model, &options, reduced_x, &with_presolve);
	options.presolve = 0;
	second = lp_simplex_solve(model, &options, raw_x, &without_presolve);
	lp_model_free(model);
	return first == lp_simplex_EXIT_SUCCESS &&
	       second == lp_simplex_EXIT_SUCCESS &&
	       with_presolve.status == lp_simplex_Success &&
	       without_presolve.status == lp_simplex_Success &&
	       absolute_value(with_presolve.objective -
			      without_presolve.objective) <= 1e-7;
}

/* Outward-rounded working bounds must not justify dropping x0 >= 0.
 * Check reduced solve + postsolve directly, so raw-model recovery cannot
 * conceal a bad reduction.  The former relative implication test admits
 * x1 = 1000.00005 and reconstructs x0 = -0.00005. */
static int preserve_nearly_implied_singleton_bound(const int upper_bound)
{
	struct lp_Model *model = lp_model_create(2, 4);
	struct simplex_Presolve presolve;
	struct lp_simplex_Options options;
	struct lp_simplex_Result result = {0};
	double x[4], *reduced_x = NULL;
	int state, passed = 0;

	if (model == NULL)
		return 0;
	model->objective[1] = -1.;
	if (upper_bound) {
		model->bounds[0].b_type = optm_BOUND_T_UP;
		model->bounds[0].ub = 0.;
	}
	model->bounds[1].b_type = optm_BOUND_T_BS;
	model->bounds[1].ub = 1000.00005;
	model->bounds[2].b_type = optm_BOUND_T_BS;
	model->bounds[2].ub = 0.00001;
	model->constraints[0].rhs = 1000.;
	model->constraints[0].coef[0] = upper_bound ? -1. : 1.;
	model->constraints[0].coef[1] = 1.;
	model->constraints[0].coef[2] = 1.;
	model->constraints[1].type = optm_CONS_T_GE;
	model->constraints[1].rhs = 500.;
	model->constraints[1].coef[1] = 1.;
	model->constraints[1].coef[2] = 1.;
	model->constraints[1].coef[3] = -1.;
	lp_simplex_default_options(&options, lp_simplex_ALGORITHM_DUAL_REVISED);
	state = simplex_presolve_run(
		&presolve, model, options.primal_tolerance);
	if (state == lp_simplex_EXIT_SUCCESS &&
		presolve.status == lp_simplex_CondUnsatisfied &&
		presolve.reduced != NULL) {
		reduced_x = (double *)calloc(
			(size_t)presolve.reduced->columns, sizeof(double));
		if (reduced_x != NULL) {
			result.status = lp_simplex_CondUnsatisfied;
			result.iterations = 0;
			state = simplex_dual_solve_problem(presolve.reduced,
				&options, reduced_x, NULL, &result, 1);
			simplex_presolve_postsolve(&presolve, reduced_x, x);
			passed = state == lp_simplex_EXIT_SUCCESS &&
				 result.status == lp_simplex_Success &&
				 (upper_bound ? x[0] <= 1e-7 : x[0] >= -1e-7) &&
				 x[2] >= -1e-7 && x[2] <= 0.0000101 &&
				 absolute_value(
					 model->constraints[0].coef[0] * x[0] +
					 x[1] + x[2] - 1000.) <= 1e-7 &&
				 absolute_value(x[1] - 1000.) <= 1e-7;
		}
	}
	free(reduced_x);
	simplex_presolve_destroy(&presolve);
	lp_model_free(model);
	return passed;
}

/* These rows are close, not parallel. The exact binary witness (-2^40,
 * 2^40) satisfies both, so approximate support matching cannot prove I. */
static int preserve_nearly_parallel_rows(void)
{
	struct lp_Model *m;
	struct simplex_Presolve pre;
	struct lp_simplex_Options o;
	struct lp_simplex_Result r;
	double x[2], witness = ldexp(1., 40);
	int cached, algorithm, valid;
	for (cached = 0; cached < 2; cached++) {
		m = lp_model_create(2, 2);
		if (!m)
			return 0;
		m->constraints[0].type = m->constraints[1].type =
			optm_CONS_T_EQ;
		m->constraints[0].coef[0] = m->constraints[0].coef[1] = 1.;
		m->constraints[1].coef[0] = 1.;
		m->constraints[1].coef[1] = 1. + ldexp(1., -40);
		m->constraints[1].rhs = 1.;
		m->bounds[0].b_type = m->bounds[1].b_type = optm_BOUND_T_FR;
		valid = -witness + m->constraints[1].coef[1] * witness == 1.;
		if (cached && lp_model_build_sparse(m)) {
			lp_model_free(m);
			return 0;
		}
		if (simplex_presolve_run(&pre, m, 1e-7)) {
			lp_model_free(m);
			return 0;
		}
		if (pre.status == lp_simplex_Infeasibility)
			fprintf(stderr,
				"feasible nearly parallel rows classified infeasible: cached=%d\n",
				cached);
		valid = valid && pre.status != lp_simplex_Infeasibility &&
			pre.status != lp_simplex_Unboundedness;
		simplex_presolve_destroy(&pre);
		for (algorithm = lp_simplex_ALGORITHM_DUAL_REVISED;
			valid && algorithm <= lp_simplex_ALGORITHM_AUTO;
			algorithm++) {
			lp_simplex_default_options(&o, algorithm);
			(void)lp_simplex_solve(m, &o, x, &r);
			valid = r.status != lp_simplex_Infeasibility &&
				r.status != lp_simplex_Unboundedness;
		}
		lp_model_free(m);
		if (!valid)
			return 0;
	}
	return 1;
}

int main(void)
{
	int passed = preserve_nearly_parallel_rows() +
		     solve_fixed_and_empty_column() +
		     solve_with_redundant_empty_row() +
		     detect_infeasible_empty_row() +
		     solve_singleton_inequality_chain() +
		     detect_infeasible_activity() +
		     solve_implied_bound_propagation() +
		     solve_dominated_parallel_row() +
		     detect_conflicting_parallel_equalities() +
		     solve_singleton_column_substitution() +
		     solve_forcing_row() +
		     solve_zero_cost_singleton_inequality_column() +
		     solve_substitution_fixed_point_cascade() +
		     solve_singleton_equality_projection() +
		     solve_degree_two_implied_free_substitution() +
		     solve_with_objective_offset() +
		     solve_dense_presolve_equivalence() +
		     exercise_presolve_work_queue() +
		     preserve_nearly_implied_singleton_bound(0) +
		     preserve_nearly_implied_singleton_bound(1);

	printf("%d/20 presolve reductions passed.\n", passed);
	return passed == 20 ? 0 : 1;
}
