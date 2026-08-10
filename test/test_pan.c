/* Regression tests for the complete Pan basis-deficiency-allowing solver. */
#include <lp_simplex/model.h>
#include <lp_simplex/solve.h>
#include <lp_simplex/status.h>

#include <math.h>
#include <stdio.h>


static int nearly_equal(const double left, const double right)
{
	return fabs(left - right) <= 1e-7 * (1. + fabs(right));
}


static int solve_pan(
		struct lp_Model *model, double *x, struct lp_simplex_Result *result)
{
	struct lp_simplex_Options options;
	lp_simplex_default_options(&options, lp_simplex_ALGORITHM_PAN_BDA);
	options.presolve = 0;
	options.iteration_limit = 10000;
	return lp_simplex_solve(model, &options, x, result);
}


static int solve_pan_normalized(
		struct lp_Model *model, double *x, struct lp_simplex_Result *result)
{
	struct lp_simplex_Options options;
	lp_simplex_default_options(&options, lp_simplex_ALGORITHM_PAN_BDA);
	options.pricing = lp_simplex_PRICING_PAN_NORMALIZED;
	options.presolve = 0;
	options.iteration_limit = 10000;
	return lp_simplex_solve(model, &options, x, result);
}


static int test_degenerate_optimum(void)
{
	struct lp_simplex_Result result;
	struct lp_Model *model = lp_model_create(1, 2);
	double x[2];
	int state;
	if (model == NULL)
		return 0;
	model->objective[0] = -1.;
	model->objective[1] = -1.;
	model->constraints[0].coef[0] = 1.;
	model->constraints[0].coef[1] = 1.;
	model->constraints[0].rhs = 1.;
	model->constraints[0].type = optm_CONS_T_LE;
	state = solve_pan(model, x, &result);
	lp_model_free(model);
	return state == lp_simplex_EXIT_SUCCESS &&
		result.status == lp_simplex_Success &&
		nearly_equal(result.objective, -1.) && nearly_equal(x[0] + x[1], 1.);
}


static int test_equality_minimum(void)
{
	struct lp_simplex_Result result;
	struct lp_Model *model = lp_model_create(1, 2);
	double x[2];
	int state;
	if (model == NULL)
		return 0;
	model->objective[0] = 1.;
	model->objective[1] = 2.;
	model->constraints[0].coef[0] = 1.;
	model->constraints[0].coef[1] = 1.;
	model->constraints[0].rhs = 1.;
	model->constraints[0].type = optm_CONS_T_EQ;
	state = solve_pan(model, x, &result);
	lp_model_free(model);
	return state == lp_simplex_EXIT_SUCCESS &&
		result.status == lp_simplex_Success &&
		nearly_equal(result.objective, 1.) &&
		nearly_equal(x[0], 1.) && nearly_equal(x[1], 0.);
}


static int test_normalized_pricing(void)
{
	struct lp_simplex_Result result;
	struct lp_Model *model = lp_model_create(1, 2);
	double x[2];
	int state;
	if (model == NULL)
		return 0;
	model->objective[0] = 1.;
	model->objective[1] = 2.;
	model->constraints[0].coef[0] = 1.;
	model->constraints[0].coef[1] = 1.;
	model->constraints[0].rhs = 1.;
	model->constraints[0].type = optm_CONS_T_EQ;
	state = solve_pan_normalized(model, x, &result);
	lp_model_free(model);
	return state == lp_simplex_EXIT_SUCCESS &&
		result.status == lp_simplex_Success &&
		nearly_equal(result.objective, 1.);
}


static int test_infeasible(void)
{
	struct lp_simplex_Result result;
	struct lp_Model *model = lp_model_create(2, 1);
	double x[1];
	int state;
	if (model == NULL)
		return 0;
	model->constraints[0].coef[0] = 1.;
	model->constraints[0].rhs = 2.;
	model->constraints[0].type = optm_CONS_T_GE;
	model->constraints[1].coef[0] = 1.;
	model->constraints[1].rhs = 1.;
	model->constraints[1].type = optm_CONS_T_LE;
	state = solve_pan(model, x, &result);
	lp_model_free(model);
	return state == lp_simplex_EXIT_FAILURE &&
		result.status == lp_simplex_Infeasibility;
}


static int test_unbounded_zero_column(void)
{
	struct lp_simplex_Result result;
	struct lp_Model *model = lp_model_create(1, 1);
	double x[1];
	int state;
	if (model == NULL)
		return 0;
	model->objective[0] = -1.;
	model->constraints[0].rhs = 0.;
	model->constraints[0].type = optm_CONS_T_EQ;
	state = solve_pan(model, x, &result);
	lp_model_free(model);
	return state == lp_simplex_EXIT_FAILURE &&
		result.status == lp_simplex_Unboundedness;
}


static int test_unbounded_dependent_direction(void)
{
	struct lp_simplex_Result result;
	struct lp_Model *model = lp_model_create(1, 2);
	double x[2];
	int state;
	if (model == NULL)
		return 0;
	model->objective[0] = -1.;
	model->objective[1] = -1.;
	model->constraints[0].coef[0] = 1.;
	model->constraints[0].coef[1] = -1.;
	model->constraints[0].rhs = 0.;
	model->constraints[0].type = optm_CONS_T_EQ;
	state = solve_pan(model, x, &result);
	lp_model_free(model);
	return state == lp_simplex_EXIT_FAILURE &&
		result.status == lp_simplex_Unboundedness;
}


static int test_bound_transformations(void)
{
	struct lp_simplex_Result result;
	struct lp_Model *model = lp_model_create(1, 3);
	double x[3];
	int state;
	if (model == NULL)
		return 0;
	model->objective[0] = 1.;
	model->objective[1] = -1.;
	model->objective[2] = -1.;
	model->bounds[0].b_type = optm_BOUND_T_FR;
	model->bounds[1].b_type = optm_BOUND_T_BS;
	model->bounds[1].lb = 1.;
	model->bounds[1].ub = 3.;
	model->bounds[2].b_type = optm_BOUND_T_UP;
	model->bounds[2].ub = 4.;
	model->constraints[0].coef[0] = 1.;
	model->constraints[0].rhs = -2.;
	model->constraints[0].type = optm_CONS_T_EQ;
	state = solve_pan(model, x, &result);
	lp_model_free(model);
	return state == lp_simplex_EXIT_SUCCESS &&
		result.status == lp_simplex_Success &&
		nearly_equal(x[0], -2.) && nearly_equal(x[1], 3.) &&
		nearly_equal(x[2], 4.) && nearly_equal(result.objective, -9.);
}


static int test_all_fixed_feasible(void)
{
	struct lp_simplex_Result result;
	struct lp_Model *model = lp_model_create(1, 1);
	double x[1];
	int state;
	if (model == NULL)
		return 0;
	model->objective[0] = 3.;
	model->bounds[0].b_type = optm_BOUND_T_BS;
	model->bounds[0].lb = 2.;
	model->bounds[0].ub = 2.;
	model->constraints[0].coef[0] = 1.;
	model->constraints[0].rhs = 2.;
	model->constraints[0].type = optm_CONS_T_EQ;
	state = solve_pan(model, x, &result);
	lp_model_free(model);
	return state == lp_simplex_EXIT_SUCCESS &&
		result.status == lp_simplex_Success &&
		nearly_equal(x[0], 2.) && nearly_equal(result.objective, 6.);
}


static int test_all_fixed_infeasible(void)
{
	struct lp_simplex_Result result;
	struct lp_Model *model = lp_model_create(1, 1);
	double x[1];
	int state;
	if (model == NULL)
		return 0;
	model->bounds[0].b_type = optm_BOUND_T_BS;
	model->bounds[0].lb = 2.;
	model->bounds[0].ub = 2.;
	model->constraints[0].coef[0] = 1.;
	model->constraints[0].rhs = 3.;
	model->constraints[0].type = optm_CONS_T_EQ;
	state = solve_pan(model, x, &result);
	lp_model_free(model);
	return state == lp_simplex_EXIT_FAILURE &&
		result.status == lp_simplex_Infeasibility;
}


int main(void)
{
	int failures = 0;
	if (!test_degenerate_optimum()) {
		printf("Pan degenerate optimum regression failed\n");
		failures++;
	}
	if (!test_equality_minimum()) {
		printf("Pan equality regression failed\n");
		failures++;
	}
	if (!test_normalized_pricing()) {
		printf("Pan normalized-pricing regression failed\n");
		failures++;
	}
	if (!test_infeasible()) {
		printf("Pan infeasibility regression failed\n");
		failures++;
	}
	if (!test_unbounded_zero_column()) {
		printf("Pan unboundedness regression failed\n");
		failures++;
	}
	if (!test_unbounded_dependent_direction()) {
		printf("Pan dependent-direction unboundedness regression failed\n");
		failures++;
	}
	if (!test_bound_transformations()) {
		printf("Pan bound transformation regression failed\n");
		failures++;
	}
	if (!test_all_fixed_feasible()) {
		printf("Pan all-fixed feasible regression failed\n");
		failures++;
	}
	if (!test_all_fixed_infeasible()) {
		printf("Pan all-fixed infeasible regression failed\n");
		failures++;
	}
	if (failures == 0)
		printf("Pan BDA regressions passed\n");
	return failures == 0 ? 0 : 1;
}
