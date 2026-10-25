/**
 * @file simplex_original_recovery.c
 * @brief Original-model result audit, one cold retry and point selection.
 *
 * The selected algorithm owns each fresh solve; this layer does not advance
 * AUTO routes. Algebraic coordinate lifting belongs to presolve.
 */

#include "simplex_stage.h"
#include "simplex_execution.h"
#include "simplex_original_recovery.h"
#include "simplex_model_audit.h"
#include "utils.h"
#include <math.h>

int simplex_original_recovery_reason(
	int status, int work, double primal, double tolerance)
{
	int feasible = isfinite(primal) && primal >= 0. && primal <= tolerance;
	if (status == lp_simplex_Success && !feasible)
		return SIMPLEX_REASON_INVALID_POSTSOLVE;
	if (status == lp_simplex_Infeasibility && feasible)
		return SIMPLEX_REASON_INFEASIBILITY_CONFLICT;
	if (status != lp_simplex_PrecisionError)
		return SIMPLEX_REASON_NONE;
	if (work == 0)
		return SIMPLEX_REASON_STARTUP_FAILURE;
	return feasible ? SIMPLEX_REASON_FEASIBLE_CANDIDATE
			: SIMPLEX_REASON_NONE;
}

int simplex_original_recovery_run(const struct lp_Model *model,
	const struct lp_simplex_Options *options, double *x,
	struct lp_simplex_Result *result, int state,
	const struct simplex_ProgressSink *progress,
	const struct simplex_ExecutionScope *execution,
	const struct simplex_ExecutionRef *reduced_actor,
	simplex_OriginalSolve solve)
{
	double primal = simplex_primal_infeasibility(model, x);
	struct simplex_StageReport audit;
	struct lp_simplex_Result raw = *result;
	int limit = options->iteration_limit;
	int valid_work = limit >= 0 && result->iterations >= 0 &&
			 result->iterations <= limit;
	int reason = simplex_original_recovery_reason(result->status,
		result->iterations, primal, options->primal_tolerance);
	int reduced_event = reduced_actor ? simplex_execution_ref_id(
						    reduced_actor, execution, 1)
					  : 0;
	struct simplex_ExecutionRef audit_actor, retry_actor;
	struct simplex_ExecutionRef parent_actor = execution->parent_actor;
	int audit_id = simplex_execution_begin(execution,
		SIMPLEX_EXECUTION_POSTSOLVE, 0, reason, reduced_event);
	double *saved_x;
	struct lp_simplex_Result saved_result, retry_result;
	struct lp_simplex_Options retry_options;
	struct simplex_ExecutionScope retry_scope;
	int retry_state, retry_id, keep, selection, chosen_source, charge;
	int initial_feasible = isfinite(primal) && primal >= 0. &&
			       primal <= options->primal_tolerance;
	int conflict, terminal, use_retry_result;

	audit_actor = simplex_execution_ref(execution, audit_id);
	simplex_stage_reset(&audit);
	if (isfinite(primal)) {
		audit.primal.value = primal;
		audit.primal.kind = SIMPLEX_OBSERVATION_PRIMAL_MAX;
		audit.primal.space = SIMPLEX_SPACE_INPUT;
		audit.primal.quality = SIMPLEX_EVIDENCE_RECOMPUTED;
		audit.primal.work = result->iterations;
	}
	if (reason == SIMPLEX_REASON_INVALID_POSTSOLVE ||
		reason == SIMPLEX_REASON_INFEASIBILITY_CONFLICT ||
		!valid_work) {
		/* A disproved success or malformed work claim is never restorable.
		 * Unknown work conservatively exhausts this phase's whole grant. */
		result->status = lp_simplex_PrecisionError;
		state = lp_simplex_EXIT_FAILURE;
		audit.stop = SIMPLEX_STOP_NUMERICAL;
		audit.failure = valid_work ? SIMPLEX_FAILURE_CERTIFICATION
					   : SIMPLEX_FAILURE_STAGE_CONTRACT;
		if (!valid_work)
			result->iterations = limit > 0 ? limit : 0;
	}
	/* Raw producer work/status and the independent audit remain distinct. */
	simplex_execution_end(execution->history, audit_id, &raw, &audit);
	if (!valid_work || reason == SIMPLEX_REASON_NONE || !solve ||
		result->iterations >= limit)
		return state;
	saved_x =
		(double *)lp_simplex_malloc((size_t)model->n * sizeof(double));
	if (!saved_x)
		return state;
	saved_result = *result;
	retry_options = *options;
	retry_options.iteration_limit = limit - saved_result.iterations;
	retry_scope = *execution;
	retry_id = simplex_execution_begin(execution,
		SIMPLEX_EXECUTION_ORIGINAL_RETRY, options->algorithm, reason,
		audit_id);
	retry_actor = simplex_execution_ref(execution, retry_id);
	if (retry_id)
		execution->history->events[retry_id - 1].grant =
			retry_options.iteration_limit;
	simplex_execution_set_parent(&retry_scope, retry_id);
	lp_simplex_memcpy(saved_x, x, (size_t)model->n * sizeof(double));
	/* x is output storage, never an implicit inherited seed. */
	{
		int j;
		for (j = 0; j < model->n; j++)
			x[j] = NAN;
	}
	lp_simplex_memset(&retry_result, 0, sizeof(retry_result));
	retry_result.status = lp_simplex_CondUnsatisfied;
	retry_result.objective = retry_result.primal_infeasibility =
		retry_result.dual_infeasibility = NAN;
	retry_result.dual_residual_space = lp_simplex_RESIDUAL_UNAVAILABLE;
	retry_state = solve(model, &retry_options, x, &retry_result, progress,
		&retry_scope);
	raw = retry_result;
	retry_state = simplex_finalize_result(
		model, &retry_options, x, &retry_result, retry_state);
	conflict = raw.status == lp_simplex_Infeasibility &&
		   (initial_feasible ||
			   (isfinite(retry_result.primal_infeasibility) &&
				   retry_result.primal_infeasibility <=
					   options->primal_tolerance));
	if (conflict) {
		retry_result.status = lp_simplex_PrecisionError;
		retry_state = lp_simplex_EXIT_FAILURE;
	}
	valid_work = retry_result.iterations >= 0 &&
		     retry_result.iterations <= retry_options.iteration_limit;
	charge = valid_work ? retry_result.iterations
			    : retry_options.iteration_limit;
	simplex_stage_reset(&audit);
	if (isfinite(retry_result.primal_infeasibility)) {
		audit.primal.value = retry_result.primal_infeasibility;
		audit.primal.kind = SIMPLEX_OBSERVATION_PRIMAL_MAX;
		audit.primal.space = SIMPLEX_SPACE_INPUT;
		audit.primal.quality = SIMPLEX_EVIDENCE_RECOMPUTED;
		audit.primal.work = charge;
	}
	if (conflict || (raw.status == lp_simplex_Success &&
				retry_result.status != lp_simplex_Success)) {
		audit.stop = SIMPLEX_STOP_NUMERICAL;
		audit.failure = SIMPLEX_FAILURE_CERTIFICATION;
	}
	if (!valid_work) {
		retry_result.status = lp_simplex_PrecisionError;
		retry_state = lp_simplex_EXIT_FAILURE;
		audit.stop = SIMPLEX_STOP_NUMERICAL;
		audit.failure = SIMPLEX_FAILURE_STAGE_CONTRACT;
	}
	/* End the local producer before the inclusive parent is assembled.
	 * Invalid local work remains visible; its grant is charged conservatively. */
	retry_id = simplex_execution_ref_id(&retry_actor, execution, 0);
	simplex_execution_end(execution->history, retry_id, &raw, &audit);
	retry_scope = *execution;
	simplex_execution_set_parent(&retry_scope,
		simplex_execution_ref_id(&retry_actor, execution, 1));
	/* A lost retry actor must not make later audits children of a replacement
	 * slot. The original parent is usable only if its own token survived. */
	if (!retry_scope.parent)
		simplex_execution_set_parent(&retry_scope,
			simplex_execution_ref_id(&parent_actor, execution, 0));
	audit_id = simplex_execution_ref_id(&audit_actor, execution, 1);
	if (conflict) {
		/* A counterexample is audit evidence, not an inherited kernel seed.
		 * Keep this distinct from why the cold retry was launched. */
		struct simplex_StageReport counterexample;
		int counterexample_id = simplex_execution_begin(&retry_scope,
			SIMPLEX_EXECUTION_RESULT_AUDIT, 0,
			SIMPLEX_REASON_INFEASIBILITY_CONFLICT, retry_id);
		simplex_stage_reset(&counterexample);
		counterexample.stop = SIMPLEX_STOP_NUMERICAL;
		counterexample.failure = SIMPLEX_FAILURE_CERTIFICATION;
		counterexample.primal.value =
			initial_feasible ? primal
					 : retry_result.primal_infeasibility;
		counterexample.primal.kind = SIMPLEX_OBSERVATION_PRIMAL_MAX;
		counterexample.primal.space = SIMPLEX_SPACE_INPUT;
		counterexample.primal.quality = SIMPLEX_EVIDENCE_RECOMPUTED;
		counterexample.primal.work = initial_feasible
						     ? saved_result.iterations
						     : raw.iterations;
		if (counterexample_id)
			execution->history->events[counterexample_id - 1]
				.point_source =
				initial_feasible ? audit_id : retry_id;
		simplex_execution_end(execution->history, counterexample_id,
			NULL, &counterexample);
	}
	terminal = retry_result.status == lp_simplex_Success ||
		   retry_result.status == lp_simplex_Infeasibility ||
		   retry_result.status == lp_simplex_Unboundedness;
	keep = retry_result.status == lp_simplex_Success ||
	       simplex_prefer_point(model, options, x, saved_x);
	chosen_source = keep ? retry_id : audit_id;
	use_retry_result = terminal || keep;
	if (use_retry_result) {
		*result = retry_result;
		state = retry_state;
	} else
		*result = saved_result;
	if (!keep) {
		lp_simplex_memcpy(
			x, saved_x, (size_t)model->n * sizeof(double));
		if (use_retry_result) {
			/* Terminal classification and the returned point are independent.
			 * Never attach another point's residuals to restored coordinates. */
			result->objective =
				simplex_original_objective(model, x);
			result->primal_infeasibility = primal;
			result->dual_infeasibility = NAN;
			result->dual_residual_space =
				lp_simplex_RESIDUAL_UNAVAILABLE;
		}
	}
	/* Both operands are in [0, limit], with charge <= limit - prior. */
	result->iterations = saved_result.iterations + charge;
	/* Selection assembles the whole phase; its inclusive work is not a
	 * child charge of the raw retry, which reports only its own local work. */
	simplex_execution_set_parent(&retry_scope,
		simplex_execution_ref_id(&parent_actor, execution, 0));
	selection = simplex_execution_begin(&retry_scope,
		SIMPLEX_EXECUTION_POINT_SELECTION, 0,
		keep ? SIMPLEX_REASON_KEEP_RETRY_POINT
		     : SIMPLEX_REASON_RESTORE_SAVED_POINT,
		retry_id);
	if (selection)
		execution->history->events[selection - 1].point_source =
			chosen_source;
	simplex_execution_end(execution->history, selection, result, NULL);
	lp_simplex_free(saved_x);
	return state;
}
