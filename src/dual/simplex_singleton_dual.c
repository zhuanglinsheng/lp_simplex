/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_singleton_dual.c
 * @brief Dualize equality models whose row residuals use paired singleton
 * columns.
 */

#include "simplex_stage.h"
#include "simplex_solve_context.h"
#include "simplex_singleton_dual.h"
#include "simplex_dual.h"
#include "utils.h"
#include "simplex_execution.h"
#include <math.h>

#include <lp_simplex/status.h>

static int singleton_dual_structure(const struct lp_Model *model, int *positive,
	int *negative, int *core_index, int *core_count)
{
	int i, j;
	for (i = 0; i < model->m; i++) {
		if (model->constraints[i].type != optm_CONS_T_EQ ||
			!model->constraints[i].coef ||
			!isfinite(model->constraints[i].rhs))
			return 0;
		positive[i] = negative[i] = -1;
	}
	*core_count = 0;
	for (j = 0; j < model->n; j++) {
		int count = 0, row = -1;
		if (!isfinite(model->objective[j]) ||
			model->bounds[j].b_type != optm_BOUND_T_LO ||
			model->bounds[j].lb != 0. ||
			model->bounds[j].v_type != optm_VAR_T_REAL)
			return 0;
		for (i = 0; i < model->m; i++) {
			if (!isfinite(model->constraints[i].coef[j]))
				return 0;
			if (model->constraints[i].coef[j] != 0.) {
				count++;
				row = i;
			}
		}
		if (count == 1 && model->objective[j] > 0.) {
			if (model->constraints[row].coef[j] > 0.) {
				if (positive[row] >= 0)
					return 0;
				positive[row] = j;
			} else {
				if (negative[row] >= 0)
					return 0;
				negative[row] = j;
			}
			core_index[j] = -1;
		} else
			core_index[j] = (*core_count)++;
	}
	for (i = 0; i < model->m; i++)
		if (positive[i] < 0 || negative[i] < 0)
			return 0;
	return 1;
}

static void singleton_result_reset(struct lp_simplex_Result *r)
{
	lp_simplex_memset(r, 0, sizeof(*r));
	r->status = lp_simplex_CondUnsatisfied;
	r->objective = r->dual_infeasibility = NAN;
	r->primal_infeasibility = HUGE_VAL;
	r->dual_residual_space = lp_simplex_RESIDUAL_UNAVAILABLE;
}

/* The paired columns are exclusive to each row, so with core variables zero
 * this is an exact structural feasibility construction over the reals.
 * Finite representation and input-model residuals are checked before use. */
static void singleton_origin(const struct lp_Model *model, const int *positive,
	const int *negative, double *x)
{
	int i, column;
	lp_simplex_memset(x, 0, (size_t)model->n * sizeof(double));
	for (i = 0; i < model->m; i++) {
		column = model->constraints[i].rhs >= 0. ? positive[i]
							 : negative[i];
		x[column] = model->constraints[i].rhs /
			    model->constraints[i].coef[column];
	}
}

int simplex_singleton_dual_solve_reported(const struct lp_Model *model,
	const struct lp_simplex_Options *options, double *x,
	struct lp_simplex_Result *result, int *applicable,
	struct simplex_StageReport *report,
	const struct simplex_ExecutionScope *execution)
{
	int i, j, core_count = 0, state = lp_simplex_EXIT_FAILURE, event;
	int *positive = NULL, *negative = NULL, *core_index = NULL;
	double *dual_x = NULL, *row_dual = NULL;
	struct lp_Model *dual = NULL;
	struct simplex_Problem auxiliary, original;
	struct simplex_StageReport local, auxiliary_report;
	struct lp_simplex_Result auxiliary_result;
	struct simplex_ExecutionScope scope;
	int auxiliary_live = 0, original_live = 0, before;
	if (!report)
		report = &local;
	simplex_stage_reset(report);
	report->singleton.available = 1;
	singleton_result_reset(result);
	*applicable = 0;
	lp_simplex_memset(&scope, 0, sizeof(scope));
	if (execution)
		scope = *execution;
	report->singleton.input_model = scope.model;
	event = simplex_execution_begin(
		&scope, SIMPLEX_EXECUTION_SINGLETON, 0, SIMPLEX_REASON_NONE, 0);
	if (event) {
		scope.history->events[event - 1].grant =
			options ? options->iteration_limit : 0;
		simplex_execution_set_parent(&scope, event);
	}
	if (!model || !options || !x || model->m <= 0 || model->n <= 0 ||
		!model->bounds || !model->constraints || !model->objective ||
		model->column_start)
		goto NOT_APPLICABLE;
	positive = (int *)lp_simplex_malloc((size_t)model->m * sizeof(int));
	negative = (int *)lp_simplex_malloc((size_t)model->m * sizeof(int));
	core_index = (int *)lp_simplex_malloc((size_t)model->n * sizeof(int));
	if (!positive || !negative || !core_index)
		goto MEMORY_FAILURE;
	if (!singleton_dual_structure(
		    model, positive, negative, core_index, &core_count))
		goto NOT_APPLICABLE;
	*applicable = report->singleton.structure_verified = 1;
	dual_x = (double *)lp_simplex_malloc((size_t)model->m * sizeof(double));
	if (!dual_x)
		goto MEMORY_FAILURE;
	if (!core_count) {
		/* No coupling remains: each positive-cost pair has an analytic optimum. */
		singleton_origin(model, positive, negative, x);
		for (i = 0; i < model->m; i++) {
			j = model->constraints[i].rhs >= 0. ? positive[i]
							    : negative[i];
			dual_x[i] = model->objective[j] /
				    model->constraints[i].coef[j];
		}
		result->status = lp_simplex_Success;
		goto AUDIT_INPUT;
	}
	dual = lp_model_create(core_count, model->m);
	row_dual = (double *)lp_simplex_malloc(
		(size_t)core_count * sizeof(double));
	if (!dual || !row_dual)
		goto MEMORY_FAILURE;
	for (i = 0; i < model->m; i++) {
		int p = positive[i], n = negative[i];
		double ap = model->constraints[i].coef[p],
		       an = model->constraints[i].coef[n];
		dual->objective[i] = -model->constraints[i].rhs;
		dual->bounds[i].b_type = optm_BOUND_T_BS;
		/* Enclose the exact ratios. An inward-rounded box could create a
		 * false dual contradiction, especially at degenerate boundaries. */
		dual->bounds[i].lb =
			nextafter(model->objective[n] / an, -HUGE_VAL);
		dual->bounds[i].ub =
			nextafter(model->objective[p] / ap, HUGE_VAL);
		/* Compact finite dual box is part of the classification proof. */
		if (!isfinite(dual->bounds[i].lb) ||
			!isfinite(dual->bounds[i].ub) ||
			dual->bounds[i].lb > dual->bounds[i].ub)
			goto NUMERICAL_FAILURE;
	}
	for (j = 0; j < model->n; j++) {
		int row = core_index[j];
		if (row < 0)
			continue;
		dual->constraints[row].type = optm_CONS_T_LE;
		dual->constraints[row].rhs = model->objective[j];
		for (i = 0; i < model->m; i++)
			dual->constraints[row].coef[i] =
				model->constraints[i].coef[j];
	}
	if (simplex_problem_from_model(&auxiliary, dual) !=
		lp_simplex_EXIT_SUCCESS)
		goto MEMORY_FAILURE;
	auxiliary_live = 1;
	before = scope.history ? scope.history->count : 0;
	auxiliary.execution = simplex_solve_context_derive(
		&scope, SIMPLEX_EXECUTION_SINGLETON_MODEL);
	if (scope.history && scope.history->count > before) {
		simplex_execution_set_parent(&auxiliary.execution, before + 1);
		scope.history->events[before].model_rows = auxiliary.rows;
		scope.history->events[before].model_columns = auxiliary.columns;
	}
	report->singleton.auxiliary_model = auxiliary.execution.model;
	simplex_stage_reset(&auxiliary_report);
	auxiliary.report = &auxiliary_report;
	singleton_result_reset(&auxiliary_result);
	state = simplex_dual_solve_problem(
		&auxiliary, options, dual_x, row_dual, &auxiliary_result, 0);
	report->singleton.auxiliary_available = 1;
	report->singleton.auxiliary_status = auxiliary_result.status;
	result->iterations = auxiliary_result.iterations;
	result->status = auxiliary_result.status;
	simplex_problem_destroy(&auxiliary);
	auxiliary_live = 0;
	if (result->iterations < 0 ||
		result->iterations > options->iteration_limit ||
		result->status < lp_simplex_Success ||
		result->status > lp_simplex_PrecisionError ||
		(state != lp_simplex_EXIT_SUCCESS &&
			state != lp_simplex_EXIT_FAILURE) ||
		((result->status == lp_simplex_Success) !=
			(state == lp_simplex_EXIT_SUCCESS))) {
		report->failure = SIMPLEX_FAILURE_STAGE_CONTRACT;
		goto NUMERICAL_FAILURE;
	}
	if (auxiliary_result.status == lp_simplex_Infeasibility) {
		/* Input feasibility follows from the exclusive +/- columns. With
		 * that fact and a certified infeasible dual, LP duality establishes
		 * input unboundedness. Never export dual infeasibility as primal I. */
		singleton_origin(model, positive, negative, x);
		result->status = lp_simplex_Unboundedness;
		goto AUDIT_INPUT;
	}
	if (auxiliary_result.status == lp_simplex_Unboundedness)
		goto NUMERICAL_FAILURE;
	if (state != lp_simplex_EXIT_SUCCESS ||
		result->status != lp_simplex_Success)
		goto FINISH;
	lp_simplex_memset(x, 0, (size_t)model->n * sizeof(double));
	for (j = 0; j < model->n; j++)
		if (core_index[j] >= 0) {
			x[j] = -row_dual[core_index[j]];
			if (x[j] < 0. && x[j] >= -options->primal_tolerance)
				x[j] = 0.;
		}
	for (i = 0; i < model->m; i++) {
		double residual = model->constraints[i].rhs;

		for (j = 0; j < model->n; j++)
			if (core_index[j] >= 0)
				residual -=
					model->constraints[i].coef[j] * x[j];
		if (residual >= 0.)
			x[positive[i]] =
				residual /
				model->constraints[i].coef[positive[i]];
		else
			x[negative[i]] =
				residual /
				model->constraints[i].coef[negative[i]];
	}

AUDIT_INPUT:
	if (simplex_problem_from_model(&original, model) !=
		lp_simplex_EXIT_SUCCESS)
		goto MEMORY_FAILURE;
	original_live = 1;
	simplex_stage_observe(&original, options, x,
		result->status == lp_simplex_Success ? dual_x : NULL, -1,
		result->iterations, report);
	/* No new lower-bound/basis handoff is exposed by this adapter. */
	report->lower_bound.quality = SIMPLEX_EVIDENCE_UNKNOWN;
	if (report->feasible_objective.quality < SIMPLEX_EVIDENCE_RECOMPUTED)
		goto NUMERICAL_FAILURE;
	if (result->status == lp_simplex_Success &&
		(!isfinite(report->dual.value) ||
			report->dual.value > options->dual_tolerance))
		goto NUMERICAL_FAILURE;
	report->singleton.input_point_checked = report->point_available = 1;
	result->objective = report->feasible_objective.value;
	result->primal_infeasibility = report->primal.value;
	if (result->status == lp_simplex_Success) {
		result->dual_infeasibility = report->dual.value;
		result->dual_residual_space = lp_simplex_RESIDUAL_ORIGINAL;
		report->singleton.mapping = SIMPLEX_SINGLETON_OPTIMUM;
		state = lp_simplex_EXIT_SUCCESS;
	} else {
		report->singleton.mapping = SIMPLEX_SINGLETON_UNBOUNDED;
		state = lp_simplex_EXIT_FAILURE;
	}
	goto FINISH;
MEMORY_FAILURE:
	result->status = lp_simplex_MemoryAllocError;
	state = lp_simplex_EXIT_FAILURE;
	goto FINISH;
NOT_APPLICABLE:
	*applicable = 0;
	result->status = lp_simplex_CondUnsatisfied;
	state = lp_simplex_EXIT_FAILURE;
	goto FINISH;
NUMERICAL_FAILURE:
	result->status = lp_simplex_PrecisionError;
	report->singleton.mapping = SIMPLEX_SINGLETON_REJECTED;
	state = lp_simplex_EXIT_FAILURE;
FINISH:
	report->stop = result->status == lp_simplex_CondUnsatisfied
			       ? SIMPLEX_STOP_NOT_APPLICABLE
			       : simplex_stage_stop_from_status(result->status);
	if (auxiliary_live)
		simplex_problem_destroy(&auxiliary);
	if (original_live)
		simplex_problem_destroy(&original);
	if (event)
		event = simplex_execution_parent(&scope);
	simplex_execution_end(scope.history, event, result, report);
	if (event && report->singleton.input_point_checked &&
		report->singleton.auxiliary_model)
		scope.history->events[event - 1].source_model =
			report->singleton.auxiliary_model;
	lp_model_free(dual);
	lp_simplex_free(dual_x);
	lp_simplex_free(row_dual);
	lp_simplex_free(positive);
	lp_simplex_free(negative);
	lp_simplex_free(core_index);
	return state;
}

int simplex_singleton_dual_solve(const struct lp_Model *model,
	const struct lp_simplex_Options *options, double *x,
	struct lp_simplex_Result *result, int *applicable)
{
	return simplex_singleton_dual_solve_reported(
		model, options, x, result, applicable, NULL, NULL);
}
