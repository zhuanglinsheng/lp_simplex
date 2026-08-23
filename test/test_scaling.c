#include "simplex_problem.h"
#include "simplex_scaling.h"
#include "utils.h"

#include <stdio.h>


static int nearly_equal(const double left, const double right)
{
	return __lp_simplex_ABS__(left - right) <=
		1e-12 * (1. + __lp_simplex_ABS__(right));
}


static struct simplex_Problem *create_problem(void)
{
	struct simplex_Problem *problem = simplex_problem_create_sparse(2, 2, 3);
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
		simplex_problem_free(problem);
		return 0;
	}
	for (j = 0; j < problem->columns; j++) {
		x[j] = original_x[j] / scaling.column[j];
		valid = valid && nearly_equal(scaled->objective[j],
			problem->objective[j] * scaling.column[j]);
	}
	valid = valid && nearly_equal(scaled->bounds[0].lb,
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
				problem->matrix.value[k] * scaling.row[i] *
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
	simplex_problem_free(problem);
	return valid;
}


int main(void)
{
	if (!test_reversible_scaling()) {
		printf("scaling regression failed\n");
		return 1;
	}
	printf("scaling regression passed\n");
	return 0;
}
