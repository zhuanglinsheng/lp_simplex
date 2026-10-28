/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex.c
 * @brief Public solve dispatch, input validation and original-model result
 * checks.
 */

#include "simplex_solve_context.h"
#include "simplex_execution.h"
#include "policy/defaults.h"
#include "policy/controls.h"
#include "simplex_dual.h"
#include "simplex_progress.h"
#include "simplex_pan.h"
#include "simplex_auto.h"
#include "simplex_model_audit.h"
#include "simplex_original_recovery.h"
#include "simplex_auto_routes.h"
#include "simplex_pan_standard.h"
#include "simplex_presolve.h"
#include "simplex_singleton_dual.h"
#include "simplex_tableau_solver.h"
#include "utils.h"

#include <lp_simplex/solve.h>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

void lp_simplex_default_options(
	struct lp_simplex_Options *options, const int algorithm)
{
	if (options == NULL)
		return;

	options->algorithm = algorithm;
	options->pricing = (algorithm == lp_simplex_ALGORITHM_DUAL_REVISED ||
				   algorithm == lp_simplex_ALGORITHM_AUTO)
				   ? lp_simplex_PRICING_DUAL_STEEPEST_EDGE
			   : algorithm == lp_simplex_ALGORITHM_PAN_BDA
				   ? lp_simplex_PRICING_DANTZIG
				   : lp_simplex_PRICING_BLAND;
	options->iteration_limit =
		algorithm == lp_simplex_ALGORITHM_AUTO
			? LP_SIMPLEX_DEFAULT_AUTO_ITERATION_LIMIT
			: SIMPLEX_DEFAULT_ITERATION_LIMIT;
	options->presolve = 1;
	options->progress_callback = NULL;
	options->progress_context = NULL;
	options->progress_interval = SIMPLEX_DEFAULT_PROGRESS_INTERVAL;
	options->primal_tolerance = SIMPLEX_DEFAULT_PRIMAL_TOLERANCE;
	options->dual_tolerance = SIMPLEX_DEFAULT_DUAL_TOLERANCE;
	options->pivot_tolerance = SIMPLEX_DEFAULT_PIVOT_TOLERANCE;
}

static int simplex_validate_sparse_storage(const struct lp_Model *model)
{
	int i, k;

	if (model->nnz < 0 || model->column_start == NULL ||
		model->row_start == NULL ||
		(model->nnz > 0 &&
			(model->row_index == NULL || model->value == NULL ||
				model->column_index == NULL ||
				model->row_value == NULL)) ||
		model->column_start[0] != 0 ||
		model->column_start[model->n] != model->nnz ||
		model->row_start[0] != 0 ||
		model->row_start[model->m] != model->nnz)
		return 0;
	for (i = 0; i < model->n; i++)
		if (model->column_start[i] > model->column_start[i + 1])
			return 0;
	for (i = 0; i < model->m; i++)
		if (model->row_start[i] > model->row_start[i + 1])
			return 0;
	for (k = 0; k < model->nnz; k++)
		if (model->row_index[k] < 0 ||
			model->row_index[k] >= model->m ||
			model->column_index[k] < 0 ||
			model->column_index[k] >= model->n ||
			!isfinite(model->value[k]) ||
			!isfinite(model->row_value[k]))
			return 0;
	return 1;
}

static int simplex_validate_model(const struct lp_Model *model)
{
	int i, j;
	int has_sparse, sparse = 0;

	if (model == NULL || model->objective == NULL ||
		model->constraints == NULL || model->m <= 0 || model->n <= 0 ||
		!isfinite(model->objective_offset) ||
		(model->objective_sense != 1. && model->objective_sense != -1.))
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
			!isfinite(model->constraints[i].rhs) ||
			model->constraints[i].type < optm_CONS_T_EQ ||
			model->constraints[i].type > optm_CONS_T_LE)
			return 0;
		if (model->constraints[i].coef)
			for (j = 0; j < model->n; j++)
				if (!isfinite(model->constraints[i].coef[j]))
					return 0;
	}
	for (j = 0; j < model->n; j++)
		if (!isfinite(model->objective[j]))
			return 0;
	if (model->bounds == NULL)
		return 1;
	for (j = 0; j < model->n; j++) {
		if (model->bounds[j].b_type < optm_BOUND_T_FR ||
			model->bounds[j].b_type > optm_BOUND_T_BS ||
			model->bounds[j].v_type != optm_VAR_T_REAL ||
			((model->bounds[j].b_type == optm_BOUND_T_LO ||
				 model->bounds[j].b_type == optm_BOUND_T_BS) &&
				!isfinite(model->bounds[j].lb)) ||
			((model->bounds[j].b_type == optm_BOUND_T_UP ||
				 model->bounds[j].b_type == optm_BOUND_T_BS) &&
				!isfinite(model->bounds[j].ub)) ||
			(model->bounds[j].b_type == optm_BOUND_T_BS &&
				model->bounds[j].lb > model->bounds[j].ub))
			return 0;
	}
	return 1;
}

static struct lp_Model *simplex_materialize_dense(const struct lp_Model *model)
{
	struct lp_Model *dense = lp_model_create(model->m, model->n);
	int i, j, k;

	if (dense == NULL)
		return NULL;
	lp_simplex_memcpy(dense->objective, model->objective,
		(size_t)model->n * sizeof(double));
	dense->objective_offset = model->objective_offset;
	dense->objective_sense = model->objective_sense;
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
				dense->constraints[model->row_index[k]]
					.coef[j] = model->value[k];
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
		!isfinite(options->primal_tolerance) ||
		!isfinite(options->dual_tolerance) ||
		!isfinite(options->pivot_tolerance) ||
		options->primal_tolerance <= 0. ||
		options->dual_tolerance <= 0. ||
		options->pivot_tolerance <= 0. ||
		(options->progress_callback &&
			(!isfinite(options->progress_interval) ||
				options->progress_interval <= 0.)))
		return 0;
	if (options->algorithm == lp_simplex_ALGORITHM_TABLEAU)
		return options->pricing == lp_simplex_PRICING_BLAND ||
		       options->pricing == lp_simplex_PRICING_DANTZIG;
	if (options->algorithm == lp_simplex_ALGORITHM_DUAL_REVISED ||
		options->algorithm == lp_simplex_ALGORITHM_AUTO)
		return options->pricing ==
		       lp_simplex_PRICING_DUAL_STEEPEST_EDGE;
	if (options->algorithm == lp_simplex_ALGORITHM_PAN_BDA)
		return options->pricing == lp_simplex_PRICING_DANTZIG ||
		       options->pricing == lp_simplex_PRICING_PAN_NORMALIZED;
	return 0;
}

static int simplex_solve_dual_raw(const struct lp_Model *model,
	const struct lp_simplex_Options *options, double *x,
	struct lp_simplex_Result *result,
	const struct simplex_ProgressSink *progress,
	const struct simplex_ExecutionScope *execution)
{
	int applicable = 0, state;

	if (model->coefficients != NULL) {
		state = simplex_singleton_dual_solve_reported(model, options, x,
			result, &applicable, NULL, execution);
		if (applicable)
			return state;
	}
	/* Raw solves are also the numerical fallback for a rejected presolve
	 * result.  Do not derive additional working bounds here: on unscaled models
	 * with extreme coefficient ranges (ETAMACRO/PILOT) repeated activity
	 * cancellation can make those optional crash bounds less reliable than the
	 * original model.  Presolved problems retain propagation below. */
	if (progress || execution->context || execution->history) {
		struct simplex_Problem problem;
		if (simplex_problem_from_model(&problem, model)) {
			result->status = lp_simplex_MemoryAllocError;
			return lp_simplex_EXIT_FAILURE;
		}
		problem.progress = progress;
		problem.execution = *execution;
		state = simplex_dual_solve_problem(
			&problem, options, x, NULL, result, 0);
		simplex_problem_destroy(&problem);
		return state;
	}
	return simplex_dual_solve(model, options, x, NULL, result, 0);
}

static int simplex_solve_pan_raw(const struct lp_Model *model,
	const struct lp_simplex_Options *options, double *x,
	struct lp_simplex_Result *result,
	const struct simplex_ProgressSink *progress,
	const struct simplex_ExecutionScope *execution)
{
	struct simplex_Problem problem;
	int state;

	if (simplex_problem_from_model(&problem, model) ==
		lp_simplex_EXIT_FAILURE) {
		result->status = lp_simplex_MemoryAllocError;
		return lp_simplex_EXIT_FAILURE;
	}
	problem.progress = progress;
	problem.execution = *execution;
	state = options->algorithm == lp_simplex_ALGORITHM_AUTO
			? simplex_auto_solve(&problem, options, x, result, 0,
				  NULL, NULL, NULL)
			: simplex_pan_solve_problem(
				  &problem, options, x, result);
	simplex_problem_destroy(&problem);
	return state;
}

static int simplex_solve_sparse_raw(const struct lp_Model *model,
	const struct lp_simplex_Options *options, double *x,
	struct lp_simplex_Result *result,
	const struct simplex_ProgressSink *progress,
	const struct simplex_ExecutionScope *execution)
{
	return (options->algorithm == lp_simplex_ALGORITHM_PAN_BDA ||
		       options->algorithm == lp_simplex_ALGORITHM_AUTO)
		       ? simplex_solve_pan_raw(
				 model, options, x, result, progress, execution)
		       : simplex_solve_dual_raw(model, options, x, result,
				 progress, execution);
}

static int simplex_solve_sparse_presolved(const struct lp_Model *model,
	const struct lp_simplex_Options *options, double *x,
	struct lp_simplex_Result *result,
	const struct simplex_ProgressSink *progress,
	const struct simplex_ExecutionScope *execution)
{
	struct simplex_Presolve presolve;
	struct simplex_ProgressBridge bridge;
	double *reduced_x = NULL;
	clock_t presolve_started;
	int j, state, bridge_state, reduced_event = 0;
	struct simplex_ExecutionRef reduced_actor = {0};
	struct simplex_ExecutionRef parent_actor = execution->parent_actor;
	struct simplex_ExecutionScope original_scope = *execution;
	struct simplex_RunConfig local_controls;
	const struct simplex_RunConfig *controls =
		simplex_run_config_get(execution->context, &local_controls);

	if (!options->presolve || controls->disable_presolve ||
		model->bounds == NULL) {
		state = simplex_solve_sparse_raw(
			model, options, x, result, progress, execution);
		return state;
	}
	lp_simplex_memset(&bridge, 0, sizeof(bridge));
	presolve_started = clock();
	if (simplex_presolve_run_configured(&presolve, model,
		    options->primal_tolerance,
		    controls->presolve_dependencies) ==
		lp_simplex_EXIT_FAILURE) {
		result->status = presolve.status;
		return lp_simplex_EXIT_FAILURE;
	}
	if (getenv(SIMPLEX_ENV_PROFILE) != NULL)
		simplex_presolve_print_profile(
			&presolve, (double)(clock() - presolve_started) /
					   (double)CLOCKS_PER_SEC);
	if (presolve.status != lp_simplex_CondUnsatisfied) {
		simplex_presolve_postsolve(&presolve, NULL, x);
		result->status = presolve.status;
		state = presolve.status == lp_simplex_Success
				? lp_simplex_EXIT_SUCCESS
				: lp_simplex_EXIT_FAILURE;
	} else if (options->algorithm == lp_simplex_ALGORITHM_PAN_BDA &&
		   controls->pan.model_form != SIMPLEX_PAN_NATIVE_BOUNDS &&
		   ((presolve.reduced->rows == model->m &&
			    presolve.reduced->columns == model->n) ||
			   simplex_pan_original_form_smaller(
				   model, presolve.reduced))) {
		/* With no structural reductions, implied bound tightening only adds
		 * redundant width rows to Pan's equality form. Keep the original LP. */
		state = simplex_solve_pan_raw(
			model, options, x, result, progress, execution);
	} else {
		reduced_x = (double *)lp_simplex_malloc(
			(size_t)presolve.reduced->columns * sizeof(double));
		if (reduced_x == NULL) {
			result->status = lp_simplex_MemoryAllocError;
			simplex_presolve_destroy(&presolve);
			return lp_simplex_EXIT_FAILURE;
		}
		for (j = 0; j < presolve.reduced->columns; j++)
			reduced_x[j] = NAN;
		presolve.reduced->execution = simplex_solve_context_derive(
			execution, SIMPLEX_EXECUTION_PRESOLVE);
		bridge_state = simplex_progress_bridge(&bridge, progress,
			&presolve, NULL, model->m, model->n, presolve.reduced,
			execution, NULL);
		if (bridge_state) {
			lp_simplex_free(reduced_x);
			simplex_presolve_destroy(&presolve);
			result->status =
				bridge_state == SIMPLEX_PROGRESS_BRIDGE_MEMORY
					? lp_simplex_MemoryAllocError
					: lp_simplex_PrecisionError;
			return lp_simplex_EXIT_FAILURE;
		}
		presolve.reduced->progress = progress ? &bridge.sink : NULL;
		reduced_event =
			simplex_execution_begin(&presolve.reduced->execution,
				SIMPLEX_EXECUTION_REDUCED_SOLVE,
				options->algorithm, SIMPLEX_REASON_NONE, 0);
		reduced_actor = simplex_execution_ref(
			&presolve.reduced->execution, reduced_event);
		simplex_execution_set_parent(
			&presolve.reduced->execution, reduced_event);
		/* Presolve uses outward-rounded bounds for safe elimination.  The crash
		 * phase then derives exact working bounds on the reduced matrix; these
		 * bounds guide initialization but do not authorize more elimination. */
		if (options->algorithm == lp_simplex_ALGORITHM_AUTO) {
			struct simplex_AutoOriginalInitializer context;
			struct simplex_AutoInitializer initializer;

			context.model = model;
			context.binding = *model;
			context.presolve = &presolve;
			context.original = x;
			context.controls = &controls->pan;
			context.scope = *execution;
			context.target =
				simplex_problem_identity(presolve.reduced);
			initializer.context = &context;
			initializer.run =
				simplex_auto_routes_initialize_original;
			state = simplex_auto_solve(presolve.reduced, options,
				reduced_x, result, 1, NULL, &initializer, NULL);
		} else if (options->algorithm == lp_simplex_ALGORITHM_PAN_BDA)
			state = simplex_pan_solve_problem(
				presolve.reduced, options, reduced_x, result);
		else
			state = simplex_dual_solve_problem(presolve.reduced,
				options, reduced_x, NULL, result, 1);
		reduced_event =
			simplex_execution_ref_id(&reduced_actor, execution, 0);
		simplex_execution_end(
			execution->history, reduced_event, result, NULL);
		if (result->dual_residual_space == lp_simplex_RESIDUAL_ORIGINAL)
			result->dual_residual_space =
				lp_simplex_RESIDUAL_PRESOLVED;
		simplex_presolve_postsolve(&presolve, reduced_x, x);
	}
	simplex_execution_set_parent(&original_scope,
		simplex_execution_ref_id(&parent_actor, execution, 0));
	state = simplex_original_recovery_run(model, options, x, result, state,
		progress, &original_scope, &reduced_actor,
		simplex_solve_sparse_raw);

	simplex_progress_bridge_destroy(&bridge);
	lp_simplex_free(reduced_x);
	simplex_presolve_destroy(&presolve);
	return state;
}

static int simplex_solve_impl(const struct lp_Model *model,
	const struct lp_simplex_Options *options, double *x,
	struct lp_simplex_Result *result,
	const struct simplex_ProgressSink *progress,
	const struct simplex_ExecutionScope *execution)
{
	int state, j;
	int status = lp_simplex_CondUnsatisfied;
	double objective = 0.;

	if (result == NULL)
		return lp_simplex_EXIT_FAILURE;
	result->status = lp_simplex_CondUnsatisfied;
	result->iterations = 0;
	result->objective = NAN;
	result->primal_infeasibility = HUGE_VAL;
	result->dual_infeasibility = NAN;
	result->dual_residual_space = lp_simplex_RESIDUAL_UNAVAILABLE;
	if (x == NULL || !simplex_validate_model(model) ||
		!simplex_validate_options(options)) {
		/*
		 * Integer variables are parsed (MPS INTORG/MARKER, BV/LI/UI) but not
		 * supported.  Say so, instead of failing with a bare status code that
		 * looks like a malformed model.
		 */
		if (x != NULL && model != NULL && model->bounds != NULL) {
			int j;

			for (j = 0; j < model->n; j++)
				if (model->bounds[j].v_type !=
					optm_VAR_T_REAL) {
					fprintf(stderr,
						"lp_simplex: model declares "
						"integer or binary variables; "
						"only continuous LPs are supported\n");
					break;
				}
		}
		return lp_simplex_EXIT_FAILURE;
	}

	for (j = 0; j < model->n; j++)
		x[j] = NAN;
	if (options->algorithm == lp_simplex_ALGORITHM_DUAL_REVISED ||
		options->algorithm == lp_simplex_ALGORITHM_PAN_BDA ||
		options->algorithm == lp_simplex_ALGORITHM_AUTO) {
		state = simplex_solve_sparse_presolved(
			model, options, x, result, progress, execution);
		return simplex_finalize_result(
			model, options, x, result, state);
	}
	if (model->coefficients == NULL) {
		struct lp_Model *dense = simplex_materialize_dense(model);

		if (dense == NULL) {
			result->status = lp_simplex_MemoryAllocError;
			return lp_simplex_EXIT_FAILURE;
		}
		state = simplex_solve_impl(
			dense, options, x, result, progress, execution);
		lp_model_free(dense);
		return state;
	}

	state = simplex_tableau_solve_model_reported(model, options, x,
		&objective, &status, &result->iterations, NULL, execution);
	result->status = status;
	result->objective = objective + model->objective_offset;
	return simplex_finalize_result(model, options, x, result, state);
}

static int simplex_solve_run(const struct lp_Model *model,
	const struct lp_simplex_Options *options, double *x,
	struct lp_simplex_Result *result, struct simplex_Execution *history,
	int experimental)
{
	struct simplex_ProgressTracker tracker;
	struct simplex_SolveContext context;
	struct simplex_Execution local_history;
	struct simplex_ExecutionScope execution;
	struct simplex_ExecutionRef root_actor;
	int state, root, observed = 0;

	if (!history)
		history = &local_history;
	simplex_solve_context_init(&context, history, &execution);
	simplex_solve_context_configure(&context, experimental);
	root = simplex_execution_begin(&execution, SIMPLEX_EXECUTION_SOLVE,
		options ? options->algorithm : 0, SIMPLEX_REASON_NONE, 0);
	root_actor = simplex_execution_ref(&execution, root);
	simplex_execution_set_parent(&execution, root);
	history->events[root - 1].grant =
		options ? options->iteration_limit : 0;
	history->events[root - 1].experimental_environment =
		context.controls.experimental_environment;
	/* The archive is optional; the causal history does not depend on its
	 * allocation or on whether the user requested progress notifications. */
	if (options && result && x && simplex_validate_model(model) &&
		simplex_validate_options(options) &&
		(options->progress_callback ||
			options->algorithm == lp_simplex_ALGORITHM_AUTO))
		observed = !simplex_progress_init_reported(
			&tracker, model, options, &execution);
	simplex_execution_set_parent(&execution,
		simplex_execution_ref_id(&root_actor, &execution, 0));
	state = simplex_solve_impl(model, options, x, result,
		observed ? &tracker.sink : NULL, &execution);
	if (observed) {
		simplex_progress_finish(&tracker, x, result->iterations);
		/* Original-model observations survive every reduced/scaled retry. */
		if (result->status != lp_simplex_Success &&
			tracker.have_primal &&
			simplex_prefer_point(
				model, options, tracker.point, x)) {
			lp_simplex_memcpy(x, tracker.point,
				(size_t)model->n * sizeof(double));
			result->objective =
				simplex_original_objective(model, x);
			result->primal_infeasibility = tracker.primal_error;
			result->dual_infeasibility = NAN;
			result->dual_residual_space =
				lp_simplex_RESIDUAL_UNAVAILABLE;
			if (result->status == lp_simplex_Infeasibility)
				result->status = lp_simplex_PrecisionError;
			simplex_execution_set_parent(&execution,
				simplex_execution_ref_id(
					&root_actor, &execution, 0));
			(void)simplex_progress_record_restore(
				&tracker, &execution, result);
		}
		simplex_progress_destroy(&tracker);
	}
	root = simplex_execution_ref_id(&root_actor, &execution, 0);
	simplex_execution_end(history, root, result, NULL);
	if (getenv(SIMPLEX_ENV_PROFILE)) {
		simplex_execution_print(history);
		if (context.retries.unrecorded)
			fprintf(stderr, "execution: retry-unrecorded=%d\n",
				context.retries.unrecorded);
	}
	simplex_solve_context_destroy(&context);
	return state;
}

int simplex_solve_reported(const struct lp_Model *model,
	const struct lp_simplex_Options *options, double *x,
	struct lp_simplex_Result *result, struct simplex_Execution *history)
{
	return simplex_solve_run(model, options, x, result, history, 0);
}

int simplex_solve_reported_experimental(const struct lp_Model *model,
	const struct lp_simplex_Options *options, double *x,
	struct lp_simplex_Result *result, struct simplex_Execution *history)
{
	return simplex_solve_run(model, options, x, result, history, 1);
}

int lp_simplex_solve(const struct lp_Model *model,
	const struct lp_simplex_Options *options, double *x,
	struct lp_simplex_Result *result)
{
	return simplex_solve_run(model, options, x, result, NULL, 0);
}

int lp_simplex_solve_experimental(const struct lp_Model *model,
	const struct lp_simplex_Options *options, double *x,
	struct lp_simplex_Result *result)
{
	return simplex_solve_run(model, options, x, result, NULL, 1);
}
