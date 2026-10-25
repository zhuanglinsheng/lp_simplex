/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_tableau_solver.c
 * @brief Model conversion and solve orchestration for tableau simplex.
 */

#include "simplex_solve_context.h"
#include "simplex_tableau_solver.h"
#include "simplex_phase.h"
#include "simplex_transform.h"
#include "utils.h"
#include "simplex_stage.h"
#include "simplex_execution.h"
#include <math.h>

#include <lp_simplex/status.h>

/* Raw producer facts are separate from the caller's independent audit. */
static void tableau_finish(const struct simplex_ExecutionScope *execution,
	int event, int state, int status, int iterations,
	struct simplex_StageReport *report)
{
	struct lp_simplex_Result raw;
	report->stop = status == lp_simplex_CondUnsatisfied
			       ? SIMPLEX_STOP_NOT_APPLICABLE
			       : simplex_stage_stop_from_status(status);
	report->point_available = state == lp_simplex_EXIT_SUCCESS &&
				  status == lp_simplex_Success &&
				  report->tableau.point_available;
	lp_simplex_memset(&raw, 0, sizeof(raw));
	raw.status = status;
	raw.iterations = iterations;
	raw.objective = raw.dual_infeasibility = NAN;
	raw.primal_infeasibility = HUGE_VAL;
	raw.dual_residual_space = lp_simplex_RESIDUAL_UNAVAILABLE;
	simplex_execution_end(execution->history, event, &raw, report);
}

static int tableau_standard_solve(const double *objective,
	const struct optm_LinearConstraint *constraints, int m, int n,
	const struct lp_simplex_Options *options, double *x, double *value,
	int *status, int *iterations, struct simplex_StageReport *report,
	const struct simplex_ExecutionScope *execution, int input_model,
	int converted, int *producer)
{
	int state;
	struct simplex_ExecutionRef actor;
	*producer = simplex_execution_begin(execution,
		SIMPLEX_EXECUTION_TABLEAU_STANDARD, options->pricing,
		SIMPLEX_REASON_NONE, 0);
	actor = simplex_execution_ref(execution, *producer);
	if (*producer)
		execution->history->events[*producer - 1].grant =
			options->iteration_limit;
	state = simplex_solve_standard_reported(objective, constraints, m, n,
		options, x, value, status, iterations, &report->tableau);
	report->tableau.input_model = input_model;
	report->tableau.standard_model = execution->model;
	report->tableau.standard_rows = m;
	report->tableau.standard_columns = n;
	report->tableau.converted = converted;
	*producer = simplex_execution_ref_id(&actor, execution, 0);
	tableau_finish(
		execution, *producer, state, *status, *iterations, report);
	return state;
}

static int tableau_solve_model(const struct lp_Model *model,
	const struct lp_simplex_Options *options, double *x, double *value,
	int *status, int *iterations, struct simplex_StageReport *report,
	const struct simplex_ExecutionScope *execution, int *producer)
{
	int transformed_m, transformed_n, before;
	struct simplex_ExecutionScope standard;
	double *transformed_objective;
	double *transformed_x;
	double *transformed_coefficients;
	double transformed_value = 0.;
	double objective_offset = 0.;
	struct optm_LinearConstraint *transformed_constraints;

	if (model->bounds == NULL)
		return tableau_standard_solve(model->objective,
			model->constraints, model->m, model->n, options, x,
			value, status, iterations, report, execution,
			execution->model, 0, producer);

	if (simplex_transform_size(model->bounds, model->m, model->n,
		    &transformed_m,
		    &transformed_n) != lp_simplex_EXIT_SUCCESS) {
		*status = lp_simplex_MemoryAllocError;
		return lp_simplex_EXIT_FAILURE;
	}
	if (simplex_transform_alloc(transformed_m, transformed_n,
		    &transformed_objective, &transformed_x,
		    &transformed_coefficients,
		    &transformed_constraints) == lp_simplex_EXIT_FAILURE) {
		*status = lp_simplex_MemoryAllocError;
		return lp_simplex_EXIT_FAILURE;
	}
	lp_simplex_memset(transformed_coefficients, 0,
		(size_t)transformed_m * (size_t)transformed_n * sizeof(double));
	simplex_transform_problem(model->objective, model->constraints,
		model->bounds, model->m, model->n, transformed_n,
		transformed_objective, &objective_offset,
		transformed_coefficients, transformed_constraints);
	/* The converted buffers are owned by this invocation. Derivation is
	 * runtime identity even without a journal; dimensions are not identity. */
	before = execution->history ? execution->history->count : 0;
	standard = simplex_solve_context_derive(
		execution, SIMPLEX_EXECUTION_TABLEAU_MODEL);
	if (execution->history && execution->history->count > before) {
		struct simplex_ExecutionEvent *e =
			execution->history->events + before;
		e->model_rows = transformed_m;
		e->model_columns = transformed_n;
		simplex_execution_set_parent(&standard, before + 1);
	}
	if (tableau_standard_solve(transformed_objective,
		    transformed_constraints, transformed_m, transformed_n,
		    options, transformed_x, &transformed_value, status,
		    iterations, report, &standard, execution->model, 1,
		    producer) == lp_simplex_EXIT_SUCCESS) {
		simplex_transform_recover(model->bounds, model->n,
			transformed_x, transformed_value, objective_offset, x,
			value);
		report->tableau.recovered_point = 1;
		simplex_transform_free(transformed_objective, transformed_x,
			transformed_coefficients, transformed_constraints);
		return lp_simplex_EXIT_SUCCESS;
	}
	simplex_transform_free(transformed_objective, transformed_x,
		transformed_coefficients, transformed_constraints);
	return lp_simplex_EXIT_FAILURE;
}

int simplex_tableau_solve_model_reported(const struct lp_Model *model,
	const struct lp_simplex_Options *options, double *x, double *value,
	int *status, int *iterations, struct simplex_StageReport *report,
	const struct simplex_ExecutionScope *execution)
{
	struct simplex_StageReport local;
	struct simplex_ExecutionScope scope;
	int state = lp_simplex_EXIT_FAILURE, event, producer = 0;
	if (!report)
		report = &local;
	simplex_stage_reset(report);
	report->tableau.available = 1;
	if (!status || !iterations) {
		report->stop = SIMPLEX_STOP_NOT_APPLICABLE;
		return state;
	}
	*iterations = 0;
	*status = lp_simplex_CondUnsatisfied;
	lp_simplex_memset(&scope, 0, sizeof(scope));
	if (execution)
		scope = *execution;
	report->tableau.input_model = scope.model;
	event = simplex_execution_begin(&scope, SIMPLEX_EXECUTION_TABLEAU,
		options ? options->pricing : 0, SIMPLEX_REASON_NONE, 0);
	if (event) {
		scope.history->events[event - 1].grant =
			options ? options->iteration_limit : 0;
		simplex_execution_set_parent(&scope, event);
	}
	if (model && options && x && value)
		state = tableau_solve_model(model, options, x, value, status,
			iterations, report, &scope, &producer);
	if (event)
		event = simplex_execution_parent(&scope);
	tableau_finish(&scope, event, state, *status, *iterations, report);
	if (event && report->point_available) {
		scope.history->events[event - 1].point_source = producer;
		scope.history->events[event - 1].source_model =
			report->tableau.standard_model;
	}
	return state;
}

int simplex_tableau_solve_model(const struct lp_Model *model,
	const struct lp_simplex_Options *options, double *x, double *value,
	int *status, int *iterations)
{
	return simplex_tableau_solve_model_reported(
		model, options, x, value, status, iterations, NULL, NULL);
}
