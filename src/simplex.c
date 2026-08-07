/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#include <lp_simplex/solve.h>
#include "simplex_dual.h"
#include "simplex_singleton_dual.h"
#include "simplex_tableau_solver.h"
#include "utils.h"

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
	if (model == NULL || model->objective == NULL || model->constraints == NULL ||
	    model->m <= 0 || model->n <= 0)
		return 0;
	for (i = 0; i < model->m; i++) {
		if (model->constraints[i].coef == NULL ||
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
		int applicable = 0;
		state = simplex_singleton_dual_solve(
			model, options, x, result, &applicable);
		if (applicable)
			return state;
		return simplex_dual_solve(model, options, x, NULL, result);
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
