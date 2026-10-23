/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file test_pan_warmstart.c
 * @brief Regression tests for Pan starts from supplied feasible points.
 */

#include "simplex_pan.h"
#include "simplex_pan_standard.h"
#include "simplex_dual.h"
#include "simplex_problem.h"

#include <lp_simplex/lp_simplex.h>

#include <math.h>
#include <stdio.h>


static int check_warm_point(int kind)
{
	struct lp_Model *model = lp_model_create(2, 3);
	struct simplex_Problem problem;
	struct lp_simplex_Options options;
	struct lp_simplex_Result result;
	double start[3] = {1., 0., 1.}, x[3] = {0., 0., 0.};
	int j, state, valid;

	if (model == NULL)
		return 0;
	model->constraints[0].coef[0] = 1.;
	model->constraints[0].coef[1] = 1.;
	model->constraints[1].coef[1] = 1.;
	model->constraints[1].coef[2] = 1.;
	for (j = 0; j < 2; j++) {
		model->constraints[j].rhs = 1.;
		model->constraints[j].type = optm_CONS_T_EQ;
	}
	for (j = 0; j < 3; j++) {
		model->bounds[j].b_type = optm_BOUND_T_BS;
		model->bounds[j].lb = 0.;
		model->bounds[j].ub = 1.;
	}
	model->objective[1] = -1.;
	if (simplex_problem_from_model(&problem, model))
		return 0;
	lp_simplex_default_options(&options, lp_simplex_ALGORITHM_PAN_BDA);
	if (kind == 1)
		start[0] = 2.;
	if (kind == 2) {
		struct lp_simplex_Result phase;

		state = simplex_dual_find_feasible(
			&problem, &options, start, &phase);
		if (state != lp_simplex_EXIT_SUCCESS ||
			phase.status != lp_simplex_CondUnsatisfied)
			return 0;
	}
	state = simplex_pan_solve_from_feasible(
		&problem, &options, start, x, &result);
	valid = kind == 1 ? state == lp_simplex_EXIT_FAILURE &&
				    result.status == lp_simplex_PrecisionError
			  : state == lp_simplex_EXIT_SUCCESS &&
				    result.status == lp_simplex_Success &&
				    fabs(result.objective + 1.) < 1e-8 &&
				    fabs(x[1] - 1.) < 1e-8 &&
				    fabs(x[0]) < 1e-8 && fabs(x[2]) < 1e-8;
	simplex_problem_destroy(&problem);
	lp_model_free(model);
	return valid;
}


static int check_compensated_dot(void)
{
	struct simplex_PanStandard standard;
	int start[2] = {0, 3}, row[3] = {0, 1, 2};
	double entry[3] = {1., 1., 1.}, vector[3] = {1e16, 1., -1e16};

	standard.column_start = start;
	standard.row_index = row;
	standard.value = entry;
	return simplex_pan_column_dot(&standard, 0, vector) == 1.;
}


static int check_tiny_nonzero_column(void)
{
	struct lp_Model *model = lp_model_create(1, 1);
	struct lp_simplex_Options options;
	struct lp_simplex_Result result;
	double x[1] = {0.};
	int state, valid;

	if (!model)
		return 0;
	model->constraints[0].coef[0] = 1e-12;
	model->constraints[0].type = optm_CONS_T_EQ;
	model->objective[0] = -1.;
	lp_simplex_default_options(&options, lp_simplex_ALGORITHM_PAN_BDA);
	options.presolve = 0;
	state = lp_simplex_solve(model, &options, x, &result);
	valid = state == lp_simplex_EXIT_SUCCESS &&
		result.status == lp_simplex_Success && fabs(x[0]) < 1e-12;
	if (!valid)
		printf("tiny nonzero constraint status=%d value=%.6g\n",
			result.status, x[0]);
	lp_model_free(model);
	return valid;
}


static int check_required_tiny_value(void)
{
	struct lp_Model *model = lp_model_create(1, 1);
	struct lp_simplex_Options options;
	struct lp_simplex_Result result;
	double x[1] = {0.};
	int state, valid;

	if (!model)
		return 0;
	model->constraints[0].coef[0] = 1e12;
	model->constraints[0].rhs = 1.;
	model->constraints[0].type = optm_CONS_T_EQ;
	model->objective[0] = 1.;
	lp_simplex_default_options(&options, lp_simplex_ALGORITHM_PAN_BDA);
	options.presolve = 0;
	state = lp_simplex_solve(model, &options, x, &result);
	valid = state == lp_simplex_EXIT_SUCCESS &&
		result.status == lp_simplex_Success &&
		fabs(1e12 * x[0] - 1.) <= options.primal_tolerance;
	lp_model_free(model);
	return valid;
}


int main(void)
{
	int valid = check_required_tiny_value() &&
		    check_tiny_nonzero_column() && check_compensated_dot() &&
		    check_warm_point(0) && check_warm_point(1) &&
		    check_warm_point(2);

	printf("Pan feasible-start regression %s\n",
		valid ? "passed" : "failed");
	return valid ? 0 : 1;
}
