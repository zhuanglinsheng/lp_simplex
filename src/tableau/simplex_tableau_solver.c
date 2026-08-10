/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#include "simplex_tableau_solver.h"
#include "simplex_phase.h"
#include "simplex_transform.h"
#include "utils.h"

#include <lp_simplex/status.h>


int simplex_tableau_solve_model(
		const struct lp_Model *model,
		const struct lp_simplex_Options *options,
		double *x, double *value, int *status, int *iterations)
{
	int transformed_m, transformed_n;
	double *transformed_objective;
	double *transformed_x;
	double *transformed_coefficients;
	double transformed_value = 0.;
	double objective_offset = 0.;
	struct optm_LinearConstraint *transformed_constraints;

	if (model->bounds == NULL)
		return simplex_solve_standard(model->objective, model->constraints,
				model->m, model->n, options,
				x, value, status, iterations);

	simplex_transform_size(model->bounds, model->m, model->n,
			       &transformed_m, &transformed_n);
	if (simplex_transform_alloc(transformed_m, transformed_n,
				    &transformed_objective, &transformed_x,
				    &transformed_coefficients,
				    &transformed_constraints) == lp_simplex_EXIT_FAILURE) {
		*status = lp_simplex_MemoryAllocError;
		return lp_simplex_EXIT_FAILURE;
	}
	lp_simplex_memset(transformed_coefficients, 0,
			  (size_t)transformed_m * transformed_n * sizeof(double));
	simplex_transform_problem(model->objective, model->constraints,
				  model->bounds, model->m, model->n,
				  transformed_n, transformed_objective,
				  &objective_offset, transformed_coefficients,
				  transformed_constraints);
	if (simplex_solve_standard(transformed_objective, transformed_constraints,
				   transformed_m, transformed_n, options,
				   transformed_x, &transformed_value,
				   status, iterations) == lp_simplex_EXIT_SUCCESS) {
		simplex_transform_recover(model->bounds, model->n, transformed_x,
					  transformed_value, objective_offset,
					  x, value);
		simplex_transform_free(transformed_objective, transformed_x,
				       transformed_coefficients,
				       transformed_constraints);
		return lp_simplex_EXIT_SUCCESS;
	}
	simplex_transform_free(transformed_objective, transformed_x,
			       transformed_coefficients, transformed_constraints);
	return lp_simplex_EXIT_FAILURE;
}
