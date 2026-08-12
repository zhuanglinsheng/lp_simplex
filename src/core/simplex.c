/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#include "simplex_dual.h"
#include "simplex_pan.h"
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
		: algorithm == lp_simplex_ALGORITHM_PAN_BDA
		? lp_simplex_PRICING_DANTZIG : lp_simplex_PRICING_BLAND;
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
	if (options->algorithm == lp_simplex_ALGORITHM_PAN_BDA)
		return options->pricing == lp_simplex_PRICING_DANTZIG ||
			options->pricing == lp_simplex_PRICING_PAN_NORMALIZED;
	return 0;
}


static double simplex_primal_infeasibility(
		const struct lp_Model *model, const double *x)
{
	double maximum = 0.;
	int i, j, k;
	for (j = 0; j < model->n; j++) {
		double violation = 0.;
		if (model->bounds == NULL) {
			if (x[j] < 0.)
				violation = -x[j];
		} else {
			int type = model->bounds[j].b_type;
			if ((type == optm_BOUND_T_LO || type == optm_BOUND_T_BS) &&
			    x[j] < model->bounds[j].lb)
				violation = model->bounds[j].lb - x[j];
			if ((type == optm_BOUND_T_UP || type == optm_BOUND_T_BS) &&
			    x[j] > model->bounds[j].ub &&
			    x[j] - model->bounds[j].ub > violation)
				violation = x[j] - model->bounds[j].ub;
		}
		if (violation > maximum)
			maximum = violation;
	}
	for (i = 0; i < model->m; i++) {
		double activity = 0.;
		double violation;
		if (model->row_start != NULL) {
			for (k = model->row_start[i]; k < model->row_start[i + 1]; k++)
				activity += model->row_value[k] * x[model->column_index[k]];
		} else {
			for (j = 0; j < model->n; j++)
				activity += model->constraints[i].coef[j] * x[j];
		}
		violation = activity - model->constraints[i].rhs;
		if (model->constraints[i].type == optm_CONS_T_EQ)
			violation = __lp_simplex_ABS__(violation);
		else if (model->constraints[i].type == optm_CONS_T_GE)
			violation = violation < 0. ? -violation : 0.;
		else
			violation = violation > 0. ? violation : 0.;
		if (violation > maximum)
			maximum = violation;
	}
	return maximum;
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
	/* Raw solves are also the numerical fallback for a rejected presolve
	 * result.  Do not derive additional working bounds here: on unscaled models
	 * with extreme coefficient ranges (ETAMACRO/PILOT) repeated activity
	 * cancellation can make those optional crash bounds less reliable than the
	 * original model.  Presolved problems retain propagation below. */
	return simplex_dual_solve(model, options, x, NULL, result, 0);
}


static int simplex_solve_pan_raw(
		const struct lp_Model *model,
		const struct lp_simplex_Options *options,
		double *x, struct lp_simplex_Result *result)
{
	struct simplex_Problem problem;
	int state;
	if (simplex_problem_from_model(&problem, model) == lp_simplex_EXIT_FAILURE) {
		result->status = lp_simplex_MemoryAllocError;
		return lp_simplex_EXIT_FAILURE;
	}
	state = simplex_pan_solve_problem(&problem, options, x, result);
	simplex_problem_destroy(&problem);
	return state;
}


static int simplex_solve_sparse_raw(
		const struct lp_Model *model,
		const struct lp_simplex_Options *options,
		double *x, struct lp_simplex_Result *result)
{
	return options->algorithm == lp_simplex_ALGORITHM_PAN_BDA
		? simplex_solve_pan_raw(model, options, x, result)
		: simplex_solve_dual_raw(model, options, x, result);
}


static int simplex_solve_sparse_presolved(
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
		state = simplex_solve_sparse_raw(model, options, x, result);
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
		state = options->algorithm == lp_simplex_ALGORITHM_PAN_BDA
			? simplex_pan_solve_problem(
				presolve.reduced, options, reduced_x, result)
			: simplex_dual_solve_problem(
				presolve.reduced, options, reduced_x, NULL, result, 1);
		simplex_presolve_postsolve(&presolve, reduced_x, x);
	}
	/* A reduced model is only an optimization, never a weaker correctness
	 * contract.  Outward-rounded substitutions can still magnify cancellation
	 * during postsolve on badly scaled models.  Retry the original sparse model
	 * when postsolve violates the requested primal tolerance, or when a crash on
	 * the reduced model failed before making a simplex iteration.  The latter
	 * is a cheap and useful discriminator: continuing failures are not hidden,
	 * while models such as SHARE1B recover without disabling presolve globally. */
	{
		double postsolve_infeasibility = result->status == lp_simplex_Success
			? simplex_primal_infeasibility(model, x) : 0.;
		int invalid_postsolve = result->status == lp_simplex_Success &&
			postsolve_infeasibility > options->primal_tolerance;
		int repairable_postsolve = invalid_postsolve &&
			(postsolve_infeasibility > 10. * options->primal_tolerance ||
			 (model->m <= 1024 && postsolve_infeasibility <=
			  10. * options->primal_tolerance));
		if (invalid_postsolve && getenv("LP_SIMPLEX_PROFILE") != NULL)
			fprintf(stderr, "presolve: rejected postsolve primal=%.6g\n",
				postsolve_infeasibility);
		if (repairable_postsolve ||
		    (result->status == lp_simplex_PrecisionError &&
		     (result->iterations == 0 ||
		      result->primal_infeasibility <= options->primal_tolerance))) {
			double *saved_x = (double *)lp_simplex_malloc(
				(size_t)model->n * sizeof(double));
			struct lp_simplex_Result saved_result = *result;
			struct lp_simplex_Result retry_result;
			int retry_state;
			if (saved_x != NULL) {
				lp_simplex_memcpy(saved_x, x,
					(size_t)model->n * sizeof(double));
				retry_state = simplex_solve_sparse_raw(
					model, options, x, &retry_result);
				if (retry_result.status == lp_simplex_Success) {
					retry_result.objective += model->objective_offset;
					*result = retry_result;
					state = retry_state;
				} else {
					lp_simplex_memcpy(x, saved_x,
						(size_t)model->n * sizeof(double));
					*result = saved_result;
				}
				lp_simplex_free(saved_x);
			}
		}
		/* Report feasibility for the original model, not only the reduced
		 * problem.  A materially infeasible postsolve vector must never retain a
		 * Success status merely because the objective happens to match. */
		if (result->status == lp_simplex_Success) {
			double original_infeasibility =
				simplex_primal_infeasibility(model, x);
			result->primal_infeasibility = original_infeasibility;
			if (original_infeasibility >
			    10. * options->primal_tolerance) {
				result->status = lp_simplex_PrecisionError;
				state = lp_simplex_EXIT_FAILURE;
			}
		}
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

	if (options->algorithm == lp_simplex_ALGORITHM_DUAL_REVISED ||
	    options->algorithm == lp_simplex_ALGORITHM_PAN_BDA) {
		state = simplex_solve_sparse_presolved(model, options, x, result);
		if (result->status == lp_simplex_Success)
			result->primal_infeasibility =
				simplex_primal_infeasibility(model, x);
		return state;
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

	state = simplex_tableau_solve_model(
		model, options, x, &objective, &status, &result->iterations);
	result->status = status;
	result->objective = objective + model->objective_offset;
	if (status == lp_simplex_Success)
		result->primal_infeasibility =
			simplex_primal_infeasibility(model, x);
	return state;
}
