/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#include "simplex_dual.h"
#include "simplex_presolve.h"
#include "simplex_singleton_dual.h"
#include "simplex_tableau_solver.h"
#include "utils.h"

#include <lp_simplex/solve.h>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>


#define SIMPLEX_DEFAULT_ITERATION_LIMIT 100000
#define SIMPLEX_DEFAULT_PRIMAL_TOLERANCE 1e-7
#define SIMPLEX_DEFAULT_DUAL_TOLERANCE 1e-7
#define SIMPLEX_DEFAULT_PIVOT_TOLERANCE 1e-9


void lp_simplex_default_options(
		struct lp_simplex_Options *options, const int algorithm)
{
	if (options == NULL)
		return;
	options->algorithm = algorithm;
	options->pricing = algorithm == lp_simplex_ALGORITHM_DUAL_REVISED
		? lp_simplex_PRICING_DUAL_STEEPEST_EDGE
		: lp_simplex_PRICING_BLAND;
	options->iteration_limit = SIMPLEX_DEFAULT_ITERATION_LIMIT;
	options->presolve = 1;
	options->primal_tolerance = SIMPLEX_DEFAULT_PRIMAL_TOLERANCE;
	options->dual_tolerance = SIMPLEX_DEFAULT_DUAL_TOLERANCE;
	options->pivot_tolerance = SIMPLEX_DEFAULT_PIVOT_TOLERANCE;
}


static int simplex_validate_sparse_storage(const struct lp_Model *model)
{
	int i, k;
	if (model->nnz < 0 || model->column_start == NULL ||
	    model->row_start == NULL ||
	    (model->nnz > 0 && (model->row_index == NULL || model->value == NULL ||
	     model->column_index == NULL || model->row_value == NULL)) ||
	    model->column_start[0] != 0 ||
	    model->column_start[model->n] != model->nnz ||
	    model->row_start[0] != 0 || model->row_start[model->m] != model->nnz)
		return 0;
	for (i = 0; i < model->n; i++)
		if (model->column_start[i] > model->column_start[i + 1])
			return 0;
	for (i = 0; i < model->m; i++)
		if (model->row_start[i] > model->row_start[i + 1])
			return 0;
	for (k = 0; k < model->nnz; k++)
		if (model->row_index[k] < 0 || model->row_index[k] >= model->m ||
		    model->column_index[k] < 0 ||
		    model->column_index[k] >= model->n)
			return 0;
	return 1;
}


static int simplex_validate_model(const struct lp_Model *model)
{
	int i, j;
	int has_sparse, sparse = 0;
	if (model == NULL || model->objective == NULL || model->constraints == NULL ||
	    model->m <= 0 || model->n <= 0)
		return 0;
	has_sparse = model->column_start != NULL || model->row_start != NULL ||
		model->row_index != NULL || model->value != NULL ||
		model->column_index != NULL || model->row_value != NULL ||
		model->nnz != 0;
	if (has_sparse) {
		if (!simplex_validate_sparse_storage(model))
			return 0;
		sparse = 1;
	}
	for (i = 0; i < model->m; i++) {
		if ((model->constraints[i].coef == NULL && !sparse) ||
		    model->constraints[i].type < optm_CONS_T_EQ ||
		    model->constraints[i].type > optm_CONS_T_LE)
			return 0;
	}
	if (model->bounds == NULL)
		return 1;
	for (j = 0; j < model->n; j++) {
		if (model->bounds[j].b_type < optm_BOUND_T_FR ||
		    model->bounds[j].b_type > optm_BOUND_T_BS ||
		    model->bounds[j].v_type != optm_VAR_T_REAL ||
		    (model->bounds[j].b_type == optm_BOUND_T_BS &&
		     model->bounds[j].lb > model->bounds[j].ub))
			return 0;
	}
	return 1;
}


static struct lp_Model *simplex_materialize_dense(
		const struct lp_Model *model)
{
	struct lp_Model *dense = lp_model_create(model->m, model->n);
	int i, j, k;
	if (dense == NULL)
		return NULL;
	lp_simplex_memcpy(dense->objective, model->objective,
		(size_t)model->n * sizeof(double));
	dense->objective_offset = model->objective_offset;
	if (model->bounds != NULL)
		lp_simplex_memcpy(dense->bounds, model->bounds,
			(size_t)model->n * sizeof(*dense->bounds));
	for (i = 0; i < model->m; i++) {
		lp_simplex_memcpy(dense->constraints[i].name,
			model->constraints[i].name,
			sizeof(dense->constraints[i].name));
		dense->constraints[i].rhs = model->constraints[i].rhs;
		dense->constraints[i].type = model->constraints[i].type;
	}
	if (model->column_start != NULL)
		for (j = 0; j < model->n; j++)
			for (k = model->column_start[j];
			     k < model->column_start[j + 1]; k++)
				dense->constraints[model->row_index[k]].coef[j] =
					model->value[k];
	else
		for (i = 0; i < model->m; i++)
			lp_simplex_memcpy(dense->constraints[i].coef,
				model->constraints[i].coef,
				(size_t)model->n * sizeof(double));
	return dense;
}


static int simplex_validate_options(const struct lp_simplex_Options *options)
{
	if (options == NULL || options->iteration_limit <= 0 ||
	    (options->presolve != 0 && options->presolve != 1) ||
	    options->primal_tolerance <= 0. || options->dual_tolerance <= 0. ||
	    options->pivot_tolerance <= 0.)
		return 0;
	if (options->algorithm == lp_simplex_ALGORITHM_TABLEAU)
		return options->pricing == lp_simplex_PRICING_BLAND ||
			options->pricing == lp_simplex_PRICING_DANTZIG;
	if (options->algorithm == lp_simplex_ALGORITHM_DUAL_REVISED)
		return options->pricing == lp_simplex_PRICING_DUAL_STEEPEST_EDGE;
	return 0;
}


static int simplex_solve_dual_raw(
		const struct lp_Model *model,
		const struct lp_simplex_Options *options,
		double *x, struct lp_simplex_Result *result)
{
	int applicable = 0, state;
	if (model->coefficients != NULL) {
		state = simplex_singleton_dual_solve(
			model, options, x, result, &applicable);
		if (applicable)
			return state;
	}
	return simplex_dual_solve(model, options, x, NULL, result, 1);
}


static int simplex_solve_dual_presolved(
		const struct lp_Model *model,
		const struct lp_simplex_Options *options,
		double *x, struct lp_simplex_Result *result)
{
	struct simplex_Presolve presolve;
	double *reduced_x = NULL;
	double objective = model->objective_offset;
	clock_t presolve_started;
	int j, state;
	if (!options->presolve || getenv("LP_SIMPLEX_DISABLE_PRESOLVE") != NULL ||
	    model->bounds == NULL) {
		state = simplex_solve_dual_raw(model, options, x, result);
		if (result->status == lp_simplex_Success)
			result->objective += model->objective_offset;
		return state;
	}
	presolve_started = clock();
	if (simplex_presolve_run(&presolve, model,
			options->primal_tolerance) == lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	if (getenv("LP_SIMPLEX_PROFILE") != NULL)
		simplex_presolve_print_profile(&presolve,
			(double)(clock() - presolve_started) /
			(double)CLOCKS_PER_SEC);
	if (presolve.status != lp_simplex_CondUnsatisfied) {
		simplex_presolve_postsolve(&presolve, NULL, x);
		result->status = presolve.status;
		state = presolve.status == lp_simplex_Success
			? lp_simplex_EXIT_SUCCESS : lp_simplex_EXIT_FAILURE;
	} else {
		reduced_x = (double *)lp_simplex_malloc(
			(size_t)presolve.reduced->columns * sizeof(double));
		if (reduced_x == NULL) {
			simplex_presolve_destroy(&presolve);
			return lp_simplex_EXIT_FAILURE;
		}
		/* Presolve uses outward-rounded bounds for safe elimination.  The crash
		 * phase then derives exact working bounds on the reduced matrix; these
		 * bounds guide initialization but do not authorize more elimination. */
		state = simplex_dual_solve_problem(
			presolve.reduced, options, reduced_x, NULL, result, 1);
		simplex_presolve_postsolve(&presolve, reduced_x, x);
	}
	if (result->status == lp_simplex_Success) {
		for (j = 0; j < model->n; j++)
			objective += model->objective[j] * x[j];
		result->objective = objective;
	}
	lp_simplex_free(reduced_x);
	simplex_presolve_destroy(&presolve);
	return state;
}


int lp_simplex_solve(
		const struct lp_Model *model,
		const struct lp_simplex_Options *options,
		double *x, struct lp_simplex_Result *result)
{
	const char *criteria;
	int state;
	int status = lp_simplex_CondUnsatisfied;
	double objective = 0.;

	if (result == NULL)
		return lp_simplex_EXIT_FAILURE;
	result->status = lp_simplex_CondUnsatisfied;
	result->iterations = 0;
	result->objective = 0.;
	result->primal_infeasibility = 0.;
	result->dual_infeasibility = 0.;
	if (x == NULL || !simplex_validate_model(model) ||
	    !simplex_validate_options(options))
		return lp_simplex_EXIT_FAILURE;

	if (options->algorithm == lp_simplex_ALGORITHM_DUAL_REVISED) {
		return simplex_solve_dual_presolved(model, options, x, result);
	}
	if (model->coefficients == NULL) {
		struct lp_Model *dense = simplex_materialize_dense(model);
		if (dense == NULL) {
			result->status = lp_simplex_MemoryAllocError;
			return lp_simplex_EXIT_FAILURE;
		}
		state = lp_simplex_solve(dense, options, x, result);
		lp_model_free(dense);
		return state;
	}

	criteria = options->pricing == lp_simplex_PRICING_DANTZIG
		? "dantzig" : "bland";
	state = simplex_tableau_solve_model(model, criteria,
					    options->iteration_limit,
					    x, &objective, &status);
	result->status = status;
	result->objective = objective + model->objective_offset;
	return state;
}
