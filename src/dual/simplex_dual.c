/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_dual.c
 * @brief Dual solve entries, scaling and input-space result assembly.
 */

#include "simplex_dual.h"
#include "simplex_dual_backend.h"
#include "simplex_dual_storage.h"
#include "simplex_basis.h"
#include "simplex_scaling.h"
#include "simplex_stage.h"
#include "simplex_progress.h"
#include "simplex_run_config.h"
#include "utils.h"

#include <lp_simplex/status.h>

#include <math.h>

/* A scaled reduced cost is not an original-units residual. Exported pi
 * uses c-A^T*pi; retain that convention and recompute after recovery. */
static void dual_report_input_residual(const struct simplex_Problem *p,
	const struct lp_simplex_Options *o, const double *x, const double *dual,
	struct lp_simplex_Result *r)
{
	int i;
	/* Only the feasibility-only service omits exported multipliers. */
	if (!dual) {
		r->dual_residual_space = lp_simplex_RESIDUAL_AUXILIARY;
		return;
	}
	for (i = 0; i < p->rows; i++)
		if (!isfinite(dual[i])) {
			r->dual_infeasibility = NAN;
			r->dual_residual_space =
				lp_simplex_RESIDUAL_UNAVAILABLE;
			return;
		}
	r->dual_infeasibility = simplex_problem_dual_infeasibility(
		p, x, dual, o->primal_tolerance, -1);
	r->dual_residual_space = lp_simplex_RESIDUAL_ORIGINAL;
	if (p->report) {
		simplex_stage_observe(
			p, o, x, dual, -1, r->iterations, p->report);
		r->dual_infeasibility = p->report->dual.value;
	}
}

static int dual_solve_problem_mode(const struct simplex_Problem *problem,
	const struct lp_simplex_Options *options, double *x, double *row_dual,
	struct lp_simplex_Result *result, const int propagate_bounds,
	const int feasibility_only, struct simplex_DualBasisSeed *seed,
	const struct simplex_Problem *retry_input)
{
	struct simplex_Scaling scaling;
	struct simplex_DualStorage storage;
	struct simplex_RunConfig local_controls;
	const struct simplex_RunConfig *controls = simplex_run_config_get(
		problem->execution.context, &local_controls);
	struct simplex_DualRetry retry;
	struct simplex_ProgressBridge progress_bridge;
	struct simplex_Problem *scaled;
	double *dual = row_dual;
	int own_dual = 0;
	int i, j, state, bridge_state;

	simplex_dual_seed_begin(seed);
	if (problem->report)
		simplex_stage_reset(problem->report);
	simplex_dual_retry_init(&retry, retry_input ? retry_input : problem,
		options, propagate_bounds, feasibility_only,
		seed && seed->valid);
	retry.history = problem->execution.history;
	result->status = lp_simplex_CondUnsatisfied;
	result->iterations = 0;
	result->objective = NAN;
	result->primal_infeasibility = HUGE_VAL;
	result->dual_infeasibility = NAN;
	result->dual_residual_space = lp_simplex_RESIDUAL_UNAVAILABLE;
	if (!simplex_dual_storage(
		    problem->rows, problem->columns, (size_t)-1, &storage)) {
		result->status = lp_simplex_MemoryAllocError;
		if (problem->report)
			problem->report->stop = SIMPLEX_STOP_RESOURCE;
		return lp_simplex_EXIT_FAILURE;
	}
	if (seed && (seed->rows != problem->rows ||
			    seed->columns != problem->columns || !seed->index ||
			    !seed->status)) {
		seed->valid = seed->exported = seed->accepted = 0;
		seed->import_attempted = 0;
		seed->rejection = SIMPLEX_BASIS_LAYOUT_REJECTED;
		if (problem->report)
			problem->report->basis_rejection = seed->rejection;
		return lp_simplex_EXIT_FAILURE;
	}
	if (dual == NULL && !feasibility_only) {
		dual = (double *)lp_simplex_malloc(storage.multiplier);
		if (dual == NULL) {
			result->status = lp_simplex_MemoryAllocError;
			if (problem->report)
				problem->report->stop = SIMPLEX_STOP_RESOURCE;
			return lp_simplex_EXIT_FAILURE;
		}
		own_dual = 1;
	}
	for (j = 0; j < problem->columns; j++)
		x[j] = 0.;
	if (dual)
		for (i = 0; i < problem->rows; i++)
			dual[i] = NAN;
	if (controls->disable_scaling) {
		state = simplex_dual_backend_solve(problem, options, x, dual,
			result, propagate_bounds, problem->columns,
			problem->rows, feasibility_only, seed, &retry);
		dual_report_input_residual(problem, options, x, dual, result);
		if (own_dual)
			lp_simplex_free(dual);
		return state;
	}
	scaled = simplex_scaling_create_problem(problem, &scaling);
	if (scaled == NULL) {
		result->status = lp_simplex_MemoryAllocError;
		if (problem->report)
			problem->report->stop = SIMPLEX_STOP_RESOURCE;
		if (own_dual)
			lp_simplex_free(dual);
		return lp_simplex_EXIT_FAILURE;
	}
	/* The kernel may reject an initial factorization before producing a
	 * solution.  Initialize the recovery buffers so that failure reporting does
	 * not read indeterminate caller storage. */
	for (j = 0; j < problem->columns; j++)
		x[j] = 0.;
	if (dual != NULL)
		for (i = 0; i < problem->rows; i++)
			dual[i] = NAN;
	bridge_state = simplex_progress_bridge(&progress_bridge,
		problem->progress, NULL, &scaling, problem->rows,
		problem->columns, scaled, &problem->execution, problem);
	if (bridge_state) {
		result->status = bridge_state == SIMPLEX_PROGRESS_BRIDGE_MEMORY
					 ? lp_simplex_MemoryAllocError
					 : lp_simplex_PrecisionError;
		if (bridge_state == SIMPLEX_PROGRESS_BRIDGE_MEMORY &&
			problem->report)
			problem->report->stop = SIMPLEX_STOP_RESOURCE;
		if (bridge_state == SIMPLEX_PROGRESS_BRIDGE_INVALID &&
			problem->report) {
			problem->report->stop = SIMPLEX_STOP_NUMERICAL;
			problem->report->failure =
				SIMPLEX_FAILURE_STAGE_CONTRACT;
		}
		if (own_dual)
			lp_simplex_free(dual);
		simplex_problem_free(scaled);
		simplex_scaling_destroy(&scaling);
		return lp_simplex_EXIT_FAILURE;
	}
	scaled->progress = problem->progress ? &progress_bridge.sink : NULL;
	/* The immutable input identity precedes deterministic internal scaling.
	 * Recreated scaling buffers are not new recovery evidence. */
	retry.key.representation = scaling.active
					   ? SIMPLEX_DUAL_RETRY_EQUILIBRATED
					   : SIMPLEX_DUAL_RETRY_IDENTITY;
	state = simplex_dual_backend_solve(scaled, options, x, dual, result,
		propagate_bounds, scaled->columns, scaled->rows,
		feasibility_only, seed, &retry);
	if (scaling.active) {
		if (problem->report)
			problem->report->input_progress_space =
				SIMPLEX_SPACE_SCALED_INPUT;
		simplex_scaling_recover_primal(&scaling, x);
		if (dual != NULL)
			simplex_scaling_recover_dual(&scaling, dual);
	}
	dual_report_input_residual(problem, options, x, dual, result);
	result->objective = 0.;
	for (j = 0; j < problem->columns; j++)
		result->objective += problem->objective[j] * x[j];
	result->primal_infeasibility =
		simplex_problem_primal_infeasibility(problem, x);
	if ((result->status == lp_simplex_Success ||
		    (feasibility_only && state == lp_simplex_EXIT_SUCCESS)) &&
		result->primal_infeasibility >
			10. * options->primal_tolerance) {
		result->status = lp_simplex_PrecisionError;
		state = lp_simplex_EXIT_FAILURE;
	}
	if (own_dual)
		lp_simplex_free(dual);
	simplex_progress_bridge_destroy(&progress_bridge);
	simplex_problem_free(scaled);
	simplex_scaling_destroy(&scaling);
	return state;
}

int simplex_dual_solve_problem(const struct simplex_Problem *problem,
	const struct lp_simplex_Options *options, double *x, double *row_dual,
	struct lp_simplex_Result *result, const int propagate_bounds)
{
	return dual_solve_problem_mode(problem, options, x, row_dual, result,
		propagate_bounds, 0, NULL, NULL);
}

int simplex_dual_solve_problem_seeded(const struct simplex_Problem *problem,
	const struct lp_simplex_Options *options, double *x, double *row_dual,
	struct lp_simplex_Result *result, const int propagate_bounds,
	struct simplex_DualBasisSeed *seed)
{
	int state;

	if (seed &&
		!simplex_problem_identity_matches(&seed->identity, problem)) {
		/* Related-LP reuse requires an explicit checked projection/rebind.
		 * A foreign optional hint must not receive a cold export. */
		seed->valid = seed->exported = seed->accepted =
			seed->import_attempted = 0;
		seed->rejection = SIMPLEX_BASIS_MODEL_REJECTED;
		state = dual_solve_problem_mode(problem, options, x, row_dual,
			result, propagate_bounds, 0, NULL, NULL);
		if (problem->report)
			problem->report->basis_rejection =
				SIMPLEX_BASIS_MODEL_REJECTED;
		return state;
	}
	if (seed)
		seed->accepted = 0;
	state = dual_solve_problem_mode(problem, options, x, row_dual, result,
		propagate_bounds, 0, seed, NULL);
	if (seed && (state != lp_simplex_EXIT_SUCCESS ||
			    result->status != lp_simplex_Success)) {
		seed->exported = 0;
		/* A checked import refused only for resources did not disprove the
		 * retained input hint. It is still not this invocation's output. */
		if (result->status != lp_simplex_MemoryAllocError ||
			seed->rejection != SIMPLEX_BASIS_RESOURCE_REJECTED)
			seed->valid = 0;
	}
	if (problem->report)
		problem->report->basis_available =
			seed && seed->valid && seed->exported;
	return state;
}

int simplex_dual_solve_problem_checkpoint(const struct simplex_Problem *problem,
	const struct lp_simplex_Options *options, double *x, double *row_dual,
	struct lp_simplex_Result *result, int propagate_bounds,
	struct simplex_DualCheckpoint *checkpoint)
{
	int state, rejection = simplex_dual_checkpoint_rejection(
			   checkpoint, problem);
	if (!rejection) {
		state = dual_solve_problem_mode(problem, options, x, row_dual,
			result, propagate_bounds, 0, &checkpoint->basis, NULL);
		if (problem->report)
			problem->report->basis_available =
				checkpoint->basis.valid &&
				checkpoint->basis.exported;
		return state;
	}
	/* A rejected optional hint cannot turn a solvable LP into a global error.
	 * Do not export a cold basis into storage bound to a different input. */
	if (checkpoint) {
		checkpoint->basis.valid = checkpoint->basis.exported =
			checkpoint->basis.accepted = 0;
		checkpoint->basis.import_attempted = 0;
		checkpoint->basis.rejection = rejection;
	}
	state = dual_solve_problem_mode(problem, options, x, row_dual, result,
		propagate_bounds, 0, NULL, NULL);
	if (checkpoint && problem->report)
		problem->report->basis_rejection = rejection;
	return state;
}

int simplex_dual_find_feasible_reported(const struct simplex_Problem *input,
	const struct lp_simplex_Options *options, double *x,
	struct lp_simplex_Result *result,
	const struct simplex_ExecutionScope *execution,
	struct simplex_StageReport *report)
{
	struct simplex_Problem view = *input;
	struct simplex_ProblemIdentity binding =
		simplex_problem_identity(input);
	struct simplex_ProblemIdentity actual;
	int state;

	lp_simplex_memset(result, 0, sizeof(*result));
	result->status = lp_simplex_PrecisionError;
	result->objective = result->dual_infeasibility = NAN;
	result->primal_infeasibility = HUGE_VAL;
	result->dual_residual_space = lp_simplex_RESIDUAL_UNAVAILABLE;
	if (execution)
		view.execution = *execution;
	if (report)
		view.report = report;
	if (view.report)
		simplex_stage_reset(view.report);
	actual = simplex_problem_identity(&view);
	actual.input = binding.input;
	if (!simplex_problem_identity_equal(&actual, &binding)) {
		if (view.report) {
			view.report->stop = SIMPLEX_STOP_NUMERICAL;
			view.report->failure = SIMPLEX_FAILURE_STAGE_CONTRACT;
		}
		return lp_simplex_EXIT_FAILURE;
	}
	state = dual_solve_problem_mode(
		&view, options, x, NULL, result, 0, 1, NULL, input);
	if (!simplex_problem_identity_unchanged(&binding, input)) {
		result->status = lp_simplex_PrecisionError;
		result->objective = result->dual_infeasibility = NAN;
		result->primal_infeasibility = HUGE_VAL;
		result->dual_residual_space = lp_simplex_RESIDUAL_UNAVAILABLE;
		if (view.report) {
			simplex_stage_reset(view.report);
			view.report->stop = SIMPLEX_STOP_NUMERICAL;
			view.report->failure = SIMPLEX_FAILURE_STAGE_CONTRACT;
		}
		return lp_simplex_EXIT_FAILURE;
	}
	return state;
}

int simplex_dual_find_feasible(const struct simplex_Problem *problem,
	const struct lp_simplex_Options *options, double *x,
	struct lp_simplex_Result *result)
{
	lp_simplex_memset(result, 0, sizeof(*result));
	return dual_solve_problem_mode(
		problem, options, x, NULL, result, 0, 1, NULL, NULL);
}

int simplex_dual_solve(const struct lp_Model *model,
	const struct lp_simplex_Options *options, double *x, double *row_dual,
	struct lp_simplex_Result *result, const int propagate_bounds)
{
	struct simplex_Problem problem;
	int state;

	if (simplex_problem_from_model(&problem, model) ==
		lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	state = simplex_dual_solve_problem(
		&problem, options, x, row_dual, result, propagate_bounds);
	simplex_problem_destroy(&problem);
	return state;
}
