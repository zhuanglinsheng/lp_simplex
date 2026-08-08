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
	options->primal_tolerance = SIMPLEX_DEFAULT_PRIMAL_TOLERANCE;
	options->dual_tolerance = SIMPLEX_DEFAULT_DUAL_TOLERANCE;
	options->pivot_tolerance = SIMPLEX_DEFAULT_PIVOT_TOLERANCE;
}


static int simplex_validate_model(const struct lp_Model *model)
{
	int i, j;
	int sparse;
	if (model == NULL || model->objective == NULL || model->constraints == NULL ||
	    model->m <= 0 || model->n <= 0)
		return 0;
	sparse = model->column_start != NULL && model->row_start != NULL &&
		(model->nnz == 0 || (model->row_index != NULL &&
		 model->value != NULL && model->column_index != NULL &&
		 model->row_value != NULL));
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
	double objective = 0.;
	clock_t presolve_started;
	int j, state;
	if (model->bounds == NULL)
		return simplex_solve_dual_raw(model, options, x, result);
	presolve_started = clock();
	if (simplex_presolve_run(&presolve, model,
			options->primal_tolerance) == lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	if (getenv("LP_SIMPLEX_PROFILE") != NULL)
		fprintf(stderr, "presolve: time=%.6f seconds\n",
			(double)(clock() - presolve_started) /
			(double)CLOCKS_PER_SEC);
	if (getenv("LP_SIMPLEX_PROFILE") != NULL &&
	    (presolve.removed_rows != 0 || presolve.removed_columns != 0 ||
	     presolve.tightened_bounds != 0))
		fprintf(stderr, "presolve: removed rows=%d columns=%d "
			"[fixed=%d empty=%d singleton-col=%d singleton-ineq-col=%d "
			"singleton-row=%d "
			"redundant-row=%d duplicate-row=%d forcing-row=%d forced-column=%d "
			"doubleton=%d singleton-column=%d implied-free-column=%d "
			"singleton-projection=%d "
			"tightened=%d passes=%d singleton-column-candidates=%d "
			"free-singleton=%d equality-singleton=%d exact-propagation=%d] "
			"remaining=%d/%d\n",
			presolve.removed_rows, presolve.removed_columns,
			presolve.fixed_columns, presolve.empty_columns,
			presolve.singleton_columns,
			presolve.singleton_inequality_columns,
			presolve.singleton_rows, presolve.redundant_rows,
			presolve.duplicate_rows,
			presolve.forcing_rows, presolve.forced_columns,
			presolve.doubleton_rows,
			presolve.singleton_column_rows,
			presolve.implied_free_columns,
			presolve.singleton_projection_columns,
			presolve.tightened_bounds, presolve.passes,
			presolve.singleton_column_candidates,
			presolve.free_singleton_columns,
			presolve.equality_singleton_columns,
			presolve.exact_bound_propagation,
			presolve.reduced != NULL ? presolve.reduced->rows :
				(presolve.terminal ? 0 : model->m),
			presolve.reduced != NULL ? presolve.reduced->columns :
				(presolve.terminal ? 0 : model->n));
	if (presolve.terminal) {
		simplex_presolve_postsolve(&presolve, NULL, x);
		result->status = presolve.terminal_status;
		state = presolve.terminal_status == lp_simplex_Success
			? lp_simplex_EXIT_SUCCESS : lp_simplex_EXIT_FAILURE;
	} else if (presolve.reduced == NULL) {
		state = simplex_solve_dual_raw(model, options, x, result);
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
	result->objective = objective;
	return state;
}
