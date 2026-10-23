/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file test_scaling.c
 * @brief Regression tests for row and column scaling and original-variable
 * recovery.
 */

#include "simplex_problem.h"
#include "simplex_scaling.h"
#include "simplex_dual.h"
#include "utils.h"

#include <stdio.h>


static int nearly_equal(const double left, const double right)
{
	return __lp_simplex_ABS__(left - right) <=
	       1e-12 * (1. + __lp_simplex_ABS__(right));
}


static struct simplex_Problem *create_problem(void)
{
	struct simplex_Problem *problem =
		simplex_problem_create_sparse(2, 2, 3);

	if (problem == NULL)
		return NULL;
	problem->matrix.column_start[0] = 0;
	problem->matrix.column_start[1] = 1;
	problem->matrix.column_start[2] = 3;
	problem->matrix.row_index[0] = 0;
	problem->matrix.row_index[1] = 0;
	problem->matrix.row_index[2] = 1;
	problem->matrix.value[0] = 1e-8;
	problem->matrix.value[1] = 2.;
	problem->matrix.value[2] = 1e8;
	problem->matrix.row_start[0] = 0;
	problem->matrix.row_start[1] = 2;
	problem->matrix.row_start[2] = 3;
	problem->matrix.column_index[0] = 0;
	problem->matrix.column_index[1] = 1;
	problem->matrix.column_index[2] = 1;
	problem->matrix.row_value[0] = 1e-8;
	problem->matrix.row_value[1] = 2.;
	problem->matrix.row_value[2] = 1e8;
	problem->objective[0] = 3.;
	problem->objective[1] = -5.;
	problem->bounds[0].b_type = optm_BOUND_T_BS;
	problem->bounds[0].v_type = optm_VAR_T_REAL;
	problem->bounds[0].lb = -2.;
	problem->bounds[0].ub = 8.;
	problem->bounds[1].b_type = optm_BOUND_T_LO;
	problem->bounds[1].v_type = optm_VAR_T_REAL;
	problem->bounds[1].lb = 0.;
	problem->bounds[1].ub = __lp_simplex_INF__;
	problem->rhs[0] = 4.;
	problem->rhs[1] = 1e8;
	problem->row_type[0] = optm_CONS_T_EQ;
	problem->row_type[1] = optm_CONS_T_LE;
	return problem;
}


static int test_reversible_scaling(void)
{
	struct simplex_Problem *problem = create_problem();
	struct simplex_Problem *scaled;
	struct simplex_Scaling scaling;
	double original_x[2] = {2., 1.};
	double x[2], dual[2] = {3., -4.};
	int i, j, k, valid = problem != NULL;

	if (!valid)
		return 0;
	scaled = simplex_scaling_create_problem(problem, &scaling);
	if (scaled == NULL || !scaling.active) {
		simplex_problem_free(scaled);
		simplex_scaling_destroy(&scaling);
		simplex_problem_free(problem);
		return 0;
	}
	valid = valid &&
		scaled->matrix.column_start == problem->matrix.column_start &&
		scaled->matrix.row_start == problem->matrix.row_start &&
		scaled->matrix.row_index == problem->matrix.row_index &&
		scaled->matrix.column_index == problem->matrix.column_index &&
		scaled->matrix.value != problem->matrix.value &&
		scaled->matrix.row_value != problem->matrix.row_value &&
		scaled->matrix.row_value ==
			scaled->matrix.value + problem->matrix.nonzeros;
	for (j = 0; j < problem->columns; j++) {
		x[j] = original_x[j] / scaling.column[j];
		valid = valid &&
			nearly_equal(scaled->objective[j],
				problem->objective[j] * scaling.column[j]);
	}
	valid = valid &&
		nearly_equal(scaled->bounds[0].lb,
			problem->bounds[0].lb / scaling.column[0]) &&
		nearly_equal(scaled->bounds[0].ub,
			problem->bounds[0].ub / scaling.column[0]) &&
		nearly_equal(scaled->bounds[1].lb,
			problem->bounds[1].lb / scaling.column[1]);
	for (i = 0; i < problem->rows; i++)
		valid = valid && nearly_equal(scaled->rhs[i],
					 problem->rhs[i] * scaling.row[i]);
	for (j = 0; j < problem->columns; j++)
		for (k = problem->matrix.column_start[j];
			k < problem->matrix.column_start[j + 1]; k++) {
			i = problem->matrix.row_index[k];
			valid = valid && nearly_equal(scaled->matrix.value[k],
						 problem->matrix.value[k] *
							 scaling.row[i] *
							 scaling.column[j]);
		}
	simplex_scaling_recover_primal(&scaling, x);
	valid = valid && nearly_equal(x[0], original_x[0]) &&
		nearly_equal(x[1], original_x[1]);
	simplex_scaling_recover_dual(&scaling, dual);
	valid = valid && nearly_equal(dual[0], 3. * scaling.row[0]) &&
		nearly_equal(dual[1], -4. * scaling.row[1]);
	simplex_problem_free(scaled);
	simplex_scaling_destroy(&scaling);
	/* Destroying the view must leave both source values and indices usable. */
	valid = valid && problem->matrix.column_start[2] == 3 &&
		problem->matrix.row_index[2] == 1 &&
		problem->matrix.value[0] == 1e-8;
	simplex_problem_free(problem);
	return valid;
}


static int test_identity_scaling_borrows_storage(void)
{
	struct simplex_Problem *problem = create_problem();
	struct simplex_Problem *scaled;
	struct simplex_Scaling scaling;
	int i, valid = problem != NULL;

	if (!valid)
		return 0;
	for (i = 0; i < 3; i++)
		problem->matrix.value[i] = problem->matrix.row_value[i] = 1.;
	scaled = simplex_scaling_create_problem(problem, &scaling);
	valid = scaled != NULL && !scaling.active;
	if (valid)
		valid = scaled->matrix.value == problem->matrix.value &&
			scaled->matrix.row_value == problem->matrix.row_value &&
			scaled->objective == problem->objective &&
			scaled->bounds == problem->bounds &&
			scaled->rhs == problem->rhs &&
			scaled->row_type == problem->row_type;
	simplex_problem_free(scaled);
	simplex_scaling_destroy(&scaling);
	valid = valid && problem->matrix.value[0] == 1. &&
		problem->objective[0] == 3. && problem->bounds[0].lb == -2.;
	simplex_problem_free(problem);
	return valid;
}


static int test_optional_dual_output(void)
{
	struct lp_Model *model = lp_model_create(2, 2);
	struct simplex_Problem problem;
	struct lp_simplex_Options options;
	struct lp_simplex_Result without_dual = {0}, with_dual = {0};
	double x[2], reference_x[2], dual[2] = {12345., 12345.};
	int state, reference_state, valid;

	if (model == NULL)
		return 0;
	model->objective[0] = 3.;
	model->objective[1] = 5.;
	model->constraints[0].type = model->constraints[1].type =
		optm_CONS_T_GE;
	model->constraints[0].coef[0] = 1e-8;
	model->constraints[0].rhs = 2e-8;
	model->constraints[1].coef[1] = 1e8;
	model->constraints[1].rhs = 3e8;
	if (simplex_problem_from_model(&problem, model) !=
		lp_simplex_EXIT_SUCCESS) {
		lp_model_free(model);
		return 0;
	}
	lp_simplex_default_options(&options, lp_simplex_ALGORITHM_DUAL_REVISED);
	state = simplex_dual_solve_problem(
		&problem, &options, x, NULL, &without_dual, 0);
	reference_state = simplex_dual_solve_problem(
		&problem, &options, reference_x, dual, &with_dual, 0);
	valid = state == lp_simplex_EXIT_SUCCESS && reference_state == state &&
		without_dual.status == lp_simplex_Success &&
		with_dual.status == without_dual.status &&
		without_dual.iterations == with_dual.iterations &&
		nearly_equal(without_dual.objective, 21.) &&
		nearly_equal(with_dual.objective, without_dual.objective) &&
		nearly_equal(x[0], 2.) && nearly_equal(x[1], 3.) &&
		nearly_equal(x[0], reference_x[0]) &&
		nearly_equal(x[1], reference_x[1]) &&
		nearly_equal(1e-8 * dual[0], 3.) &&
		nearly_equal(1e8 * dual[1], 5.);
	simplex_problem_destroy(&problem);
	lp_model_free(model);
	return valid;
}


int main(void)
{
	if (!test_reversible_scaling() ||
		!test_identity_scaling_borrows_storage() ||
		!test_optional_dual_output()) {
		printf("scaling regression failed\n");
		return 1;
	}
	printf("scaling regression passed\n");
	return 0;
}
