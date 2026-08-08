#include <lp_simplex/lp_simplex.h>
#include <stdio.h>


static double absolute_value(const double value)
{
	return value < 0. ? -value : value;
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


int main(void)
{
	int passed = solve_fixed_and_empty_column() +
		solve_with_redundant_empty_row() + detect_infeasible_empty_row();
	printf("%d/3 presolve reductions passed.\n", passed);
	return passed == 3 ? 0 : 1;
}
