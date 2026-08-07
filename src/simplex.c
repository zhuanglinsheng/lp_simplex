/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#include <lp_simplex/solve.h>
#include "simplex_phase.h"
#include "simplex_transform.h"
#include "utils.h"


static int simplex_solve_arrays(
		const double *objective,
		const struct optm_LinearConstraint *constraints,
		const struct optm_VariableBound *bounds,
		int m, int n, const char *criteria, int iteration_limit,
		double *x, double *value, int *status)
{
	int transformed_m, transformed_n;
	double *transformed_objective;
	double *transformed_x;
	double *transformed_coefficients;
	double transformed_value = 0.;
	double objective_offset = 0.;
	struct optm_LinearConstraint *transformed_constraints;
	int i, j;

	if (status == NULL)
		return lp_simplex_EXIT_FAILURE;
	if (objective == NULL || constraints == NULL || x == NULL || value == NULL ||
	    m <= 0 || n <= 0 || iteration_limit <= 0) {
		*status = lp_simplex_CondUnsatisfied;
		return lp_simplex_EXIT_FAILURE;
	}
	if (criteria == NULL)
		criteria = "";
	for (i = 0; i < m; i++) {
		if (constraints[i].coef == NULL || constraints[i].type < optm_CONS_T_EQ ||
		    constraints[i].type > optm_CONS_T_LE) {
			*status = lp_simplex_CondUnsatisfied;
			return lp_simplex_EXIT_FAILURE;
		}
	}
	if (bounds == NULL)
		return simplex_solve_standard(objective, constraints, m, n, criteria,
					      iteration_limit, x, value, status);
	for (j = 0; j < n; j++) {
		if (bounds[j].b_type < optm_BOUND_T_FR ||
		    bounds[j].b_type > optm_BOUND_T_BS ||
		    bounds[j].v_type != optm_VAR_T_REAL ||
		    (bounds[j].b_type == optm_BOUND_T_BS &&
		     bounds[j].lb > bounds[j].ub)) {
			*status = lp_simplex_CondUnsatisfied;
			return lp_simplex_EXIT_FAILURE;
		}
	}

	simplex_transform_size(bounds, m, n, &transformed_m, &transformed_n);
	if (simplex_transform_alloc(transformed_m, transformed_n,
				    &transformed_objective, &transformed_x,
				    &transformed_coefficients,
				    &transformed_constraints) == lp_simplex_EXIT_FAILURE) {
		*status = lp_simplex_MemoryAllocError;
		return lp_simplex_EXIT_FAILURE;
	}
	lp_simplex_memset(transformed_coefficients, 0,
			  (size_t)transformed_m * transformed_n * sizeof(double));
	simplex_transform_problem(objective, constraints, bounds, m, n,
				  transformed_n, transformed_objective,
				  &objective_offset, transformed_coefficients,
				  transformed_constraints);
	if (simplex_solve_standard(transformed_objective, transformed_constraints,
				   transformed_m, transformed_n, criteria,
				   iteration_limit, transformed_x,
				   &transformed_value, status) == lp_simplex_EXIT_SUCCESS) {
		simplex_transform_recover(bounds, n, transformed_x, transformed_value,
					  objective_offset, x, value);
		simplex_transform_free(transformed_objective, transformed_x,
				       transformed_coefficients,
				       transformed_constraints);
		*status = lp_simplex_Success;
		return lp_simplex_EXIT_SUCCESS;
	}
	simplex_transform_free(transformed_objective, transformed_x,
			       transformed_coefficients, transformed_constraints);
	return lp_simplex_EXIT_FAILURE;
}


int lp_simplex_solve(const struct lp_Model *model, const char *criteria,
		int iteration_limit, double *x, double *objective,
		int *status)
{
	if (model == NULL) {
		if (status != NULL)
			*status = lp_simplex_CondUnsatisfied;
		return lp_simplex_EXIT_FAILURE;
	}
	return simplex_solve_arrays(model->objective, model->constraints, model->bounds,
				    model->m, model->n, criteria, iteration_limit,
				    x, objective, status);
}
