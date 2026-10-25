/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_auto.c
 * @brief AUTO owns stage selection and the cumulative work budget.
 *
 * Kernels own their numerical safeguards. Iterations count pivots/first-order
 * steps/NNLS work, not seconds; factorization cost is deliberately not called
 * an iteration.
 */

#include "simplex_solve_context.h"
#include "simplex_execution.h"
#include "policy/controls.h"
#include "simplex_auto.h"
#include "simplex_auto_policy.h"
#include "simplex_stage.h"
#include "simplex_auto_routes.h"
#include "simplex_auto_state.h"
#include "utils.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

static void auto_reset(struct lp_simplex_Result *r)
{
	lp_simplex_memset(r, 0, sizeof(*r));
	r->status = lp_simplex_CondUnsatisfied;
	r->objective = r->primal_infeasibility = r->dual_infeasibility =
		HUGE_VAL;
}

/* Working bounds and scaled residuals never authorize passing a seed. */
static int auto_audit_point(const struct simplex_Problem *p, const double *x,
	const struct lp_simplex_Options *o, int work,
	struct simplex_StageReport *audit)
{
	simplex_stage_reset(audit);
	simplex_stage_observe(p, o, x, NULL, 1, work, audit);
	return audit->feasible_objective.quality >= SIMPLEX_EVIDENCE_RECOMPUTED;
}

/* A status code is not a substitute for a usable diagnostic. Tableau is
 * the documented exception: its internal Phase II check exports no dual. */
static int auto_usable_diagnostic(const struct lp_simplex_Result *r, int route)
{
	/* Certification belongs to the kernel's scaled/model-space checks.
	 * An exported residual in input units is not universally comparable to
	 * options->dual_tolerance; only reject malformed public diagnostics here. */
	if (isfinite(r->dual_infeasibility))
		return r->dual_infeasibility >= 0.;
	return route == SIMPLEX_AUTO_TABLEAU &&
	       r->dual_residual_space == lp_simplex_RESIDUAL_UNAVAILABLE &&
	       isnan(r->dual_infeasibility);
}

static const char *auto_reason(int status)
{
	switch (status) {
	case lp_simplex_Success:
		return "optimal";
	case lp_simplex_Infeasibility:
		return "infeasible";
	case lp_simplex_Unboundedness:
		return "unbounded";
	case lp_simplex_MemoryAllocError:
		return "memory-failure";
	case lp_simplex_ExceedIterLimit:
		return "work-budget";
	case lp_simplex_Degeneracy:
		return "stagnation";
	case lp_simplex_Singularity:
		return "singularity";
	default:
		return "numerical-or-audit-failure";
	}
}

static void auto_trace(
	const char *stage, const char *reason, int budget, int used, int seed)
{
	if (getenv(SIMPLEX_ENV_PROFILE))
		fprintf(stderr,
			"auto: stage=%s reason=%s budget=%d used=%d feasible-seed=%d\n",
			stage, reason, budget, used, seed);
}

static void auto_budget_trace(
	const char *stage, int remaining, int paths, int reserve)
{
	if (getenv(SIMPLEX_ENV_PROFILE))
		fprintf(stderr,
			"auto: allocation=%s remaining=%d recovery-paths=%d observation=%d reserve=%d grant=%d\n",
			stage, remaining, paths, SIMPLEX_AUTO_PROBE_CAP,
			reserve, remaining - reserve);
}

static const char *auto_route_name(int route)
{
	static const char *names[] = {"native-hybrid", "native-explicit",
		"classic-pan", "standard-implicit", "standard-explicit",
		"singleton-dual", "tableau"};
	return names[route];
}

/* The scheduler audits handoff contracts; mathematical certificates remain
 * kernel responsibilities. A retryable failure cannot become a terminal
 * classification merely because another route is unavailable. */
static int auto_result_terminal(struct lp_simplex_Result *r, int *state,
	int feasible, int have_best, int route, int used)
{
	if (r->status == lp_simplex_Success) {
		if (*state == lp_simplex_EXIT_SUCCESS && feasible &&
			auto_usable_diagnostic(r, route))
			return 1;
		auto_trace("result-audit",
			feasible ? "invalid-dual-diagnostic"
				 : "invalid-primal-point",
			0, used, have_best);
		r->status = lp_simplex_PrecisionError;
		*state = lp_simplex_EXIT_FAILURE;
		return 0;
	}
	if (r->status == lp_simplex_Infeasibility && (feasible || have_best)) {
		auto_trace("result-audit",
			"infeasibility-conflicts-with-feasible-point", 0, used,
			have_best);
		r->status = lp_simplex_PrecisionError;
		*state = lp_simplex_EXIT_FAILURE;
		return 0;
	}
	switch (r->status) {
	case lp_simplex_Infeasibility:
		/* Kernel-certified contradiction; no tolerance-feasible counterexample. */
	case lp_simplex_Unboundedness:
		/* Kernel Phase II owns feasibility and its improving-ray check. */
		return 1;
	default:
		return 0;
	}
}

static int auto_charge_stage(struct lp_simplex_Result *r, int *state, int grant)
{
	if (r->iterations < 0 || r->iterations > grant) {
		r->iterations = grant;
		r->status = lp_simplex_PrecisionError;
		*state = lp_simplex_EXIT_FAILURE;
		return 0;
	}
	return 1;
}

static void auto_reject_outputs(struct simplex_StageReport *report)
{
	report->stop = SIMPLEX_STOP_NUMERICAL;
	report->failure = SIMPLEX_FAILURE_STAGE_CONTRACT;
	report->lower_bound.quality = SIMPLEX_EVIDENCE_UNKNOWN;
	report->point_available = report->basis_available =
		report->factor_available = 0;
}

static const char *auto_prerequisite_name(int prerequisite)
{
	switch (prerequisite) {
	case SIMPLEX_AUTO_AUDITED_POINT:
		return "audited-primal-point";
	case SIMPLEX_AUTO_BASIS_HINT:
		return "basis-hint-needs-rebuild";
	case SIMPLEX_AUTO_IMPLICIT_ORACLE_FAILURE:
		return "diagnosed-implicit-oracle-failure";
	default:
		return "no-certified-result";
	}
}

static const char *auto_purpose_name(int purpose)
{
	switch (purpose) {
	case SIMPLEX_AUTO_SEEK_FEASIBLE:
		return "seek-feasibility";
	case SIMPLEX_AUTO_OPTIMIZE_POINT:
		return "optimize-validated-point";
	case SIMPLEX_AUTO_REBUILD_BASIS:
		return "rebuild-and-reclassify-basis";
	case SIMPLEX_AUTO_REFORMULATE:
		return "reformulate-and-audit";
	case SIMPLEX_AUTO_REPLACE_ORACLE:
		return "replace-implicit-oracle";
	case SIMPLEX_AUTO_RESTRICT_COLUMNS:
		return "restricted-column-search";
	case SIMPLEX_AUTO_ALTERNATIVE_ALGORITHM:
		return "alternative-algorithm";
	case SIMPLEX_AUTO_REPLACE_DIRECTION:
		return "replace-basis-direction-with-cone";
	default:
		return "primary-search";
	}
}

static void auto_begin_action(struct simplex_AutoAction *action, int route,
	int grant, const struct simplex_AutoPolicyEvidence *e,
	int recovery_available)
{
	struct simplex_AutoTrigger trigger =
		simplex_auto_policy_trigger(route, e);
	int uses_point = e->have_point && simplex_auto_policy_uses_point(route);
	int uses_basis = e->have_basis && route == SIMPLEX_AUTO_DUAL_FALLBACK;

	lp_simplex_memset(action, 0, sizeof(*action));
	action->failure_action = recovery_available ? SIMPLEX_AUTO_NEXT_ROUTE
						    : SIMPLEX_AUTO_NO_RECOVERY;
	action->stage = route;
	action->grant = grant;
	action->prerequisite = simplex_auto_policy_prerequisite(
		route, e->have_point, e->have_basis);
	action->purpose = simplex_auto_policy_purpose(route, e);
	action->reason = trigger.reason;
	action->cause = trigger.origin.record;
	action->cause_event = trigger.origin.event;
	action->cause_actor = trigger.origin.actor;
	action->point_source = uses_point ? e->point.record : -1;
	action->point_event = uses_point ? e->point.event : 0;
	if (uses_point)
		action->point_actor = e->point.actor;
	action->basis_source = uses_basis ? e->basis.record : -1;
	action->basis_event = uses_basis ? e->basis.event : 0;
	if (uses_basis)
		action->basis_actor = e->basis.actor;
	if (getenv(SIMPLEX_ENV_PROFILE))
		fprintf(stderr,
			"auto: transition target=%d prerequisite=%s purpose=%s "
			"reason=%d cause=%d point-source=%d basis-source=%d "
			"validation=immutable-input-audit on-retryable-failure=%s\n",
			route, auto_prerequisite_name(action->prerequisite),
			auto_purpose_name(action->purpose), action->reason,
			action->cause, action->point_source,
			action->basis_source,
			recovery_available ? "next-route-with-remaining-budget"
					   : "stop-after-audit");
}

static void auto_action_sources(struct simplex_AutoAction *action,
	const struct simplex_ExecutionScope *scope)
{
	action->cause_event =
		simplex_execution_ref_id(&action->cause_actor, scope, 1);
	action->point_event =
		simplex_execution_ref_id(&action->point_actor, scope, 1);
	action->basis_event =
		simplex_execution_ref_id(&action->basis_actor, scope, 1);
}

/* Provenance comes from runtime evidence, never from the display journal. */
static void auto_link_action(struct simplex_Problem *problem,
	const struct simplex_ExecutionScope *parent,
	const struct simplex_ExecutionRef *parent_actor,
	struct simplex_AutoAction *action)
{
	struct simplex_ExecutionEvent *event;

	problem->execution = *parent;
	simplex_execution_set_parent(&problem->execution,
		simplex_execution_ref_id(parent_actor, parent, 0));
	auto_action_sources(action, parent);
	action->execution_id = simplex_execution_begin(&problem->execution,
		SIMPLEX_EXECUTION_AUTO, action->stage, action->reason,
		action->cause_event);
	action->actor = simplex_execution_ref(parent, action->execution_id);
	simplex_execution_set_parent(&problem->execution, action->execution_id);
	if (!action->execution_id)
		return;
	event = parent->history->events + action->execution_id - 1;
	event->grant = action->grant;
	event->prerequisite = action->prerequisite;
	event->purpose = action->purpose;
	event->on_failure = action->failure_action;
	event->point_source = action->point_event;
	event->basis_source = action->basis_event;
}

static int auto_record_stage(const struct simplex_Problem *problem,
	struct simplex_AutoHistory *history,
	const struct simplex_AutoAction *action,
	const struct simplex_StageReport *report,
	const struct lp_simplex_Result *result,
	const struct lp_simplex_Result *kernel_result, int terminal,
	int feasible, int accepted_point,
	const struct simplex_StageReport *audit)
{
	struct simplex_AutoStageRecord *record;
	struct simplex_AutoAction observed = *action;
	int index = history->count;
	int claim = kernel_result->status == lp_simplex_Success ||
		    kernel_result->status == lp_simplex_Infeasibility ||
		    kernel_result->status == lp_simplex_Unboundedness;

	observed.execution_id = simplex_execution_ref_id(
		&action->actor, &problem->execution, 0);
	auto_action_sources(&observed, &problem->execution);
	simplex_execution_end(problem->execution.history, observed.execution_id,
		result, report);
	if (observed.execution_id)
		problem->execution.history->events[observed.execution_id - 1]
			.primal_audit = audit->primal;
	if (index >= SIMPLEX_AUTO_ROUTE_COUNT + 4)
		return -1;
	record = history->stages + history->count++;
	record->action = observed;
	record->work = result->iterations;
	record->kernel_status = kernel_result->status;
	record->reported_work = kernel_result->iterations;
	record->status = result->status;
	record->accepted_point = accepted_point;
	record->report = *report;
	record->primal_audit = audit->primal;
	record->objective_audit = audit->feasible_objective;
	record->validation =
		claim ? (terminal ? SIMPLEX_AUTO_TERMINAL_VALID
				  : SIMPLEX_AUTO_TERMINAL_INVALID)
		      : (feasible ? SIMPLEX_AUTO_POINT_VALID
				  : (audit->primal.quality >=
							    SIMPLEX_EVIDENCE_RECOMPUTED
						    ? SIMPLEX_AUTO_POINT_INVALID
						    : SIMPLEX_AUTO_VALIDATION_UNKNOWN));
	if (simplex_stage_observation_current(&report->lower_bound,
		    SIMPLEX_SPACE_INPUT, 0, result->iterations,
		    SIMPLEX_EVIDENCE_VERIFIED) &&
		report->lower_bound.kind == SIMPLEX_OBSERVATION_LOWER_BOUND &&
		report->lower_bound.value > history->best_bound) {
		history->best_bound = report->lower_bound.value;
		history->bound_source = index;
	}
	return index;
}

/* Snapshot before invocation. Parent action ids are deliberately excluded. */
static struct simplex_ProblemIdentity auto_model_identity(
	const struct simplex_Problem *p)
{
	return simplex_problem_identity(p);
}

static int auto_model_matches(const struct simplex_ProblemIdentity *id,
	const struct simplex_Problem *p)
{
	return simplex_problem_identity_unchanged(id, p);
}

static void auto_start_attempt(
	struct simplex_AutoAttempt *a, const struct simplex_Problem *p)
{
	auto_reset(&a->result);
	auto_reset(&a->raw_result);
	simplex_stage_reset(&a->report);
	simplex_stage_reset(&a->audit);
	a->identity = auto_model_identity(p);
	a->state = lp_simplex_EXIT_FAILURE;
	a->feasible = a->terminal = 0;
	a->work_valid = 1;
	a->objective = HUGE_VAL;
}

static int auto_audit_attempt(const struct simplex_Problem *p, const double *x,
	const struct lp_simplex_Options *o, struct simplex_AutoAttempt *a)
{
	if (!auto_model_matches(&a->identity, p)) {
		simplex_stage_reset(&a->audit);
		auto_reject_outputs(&a->report);
		a->result.status = lp_simplex_PrecisionError;
		a->state = lp_simplex_EXIT_FAILURE;
		return 0;
	}
	return auto_audit_point(p, x, o, a->result.iterations, &a->audit) &&
	       a->work_valid;
}

static int auto_candidate_better(const struct simplex_Problem *p,
	const struct simplex_AutoCandidate *best,
	const struct simplex_AutoAttempt *a)
{
	return a->work_valid && a->feasible &&
	       auto_model_matches(&a->identity, p) &&
	       p->revision != ULONG_MAX &&
	       (!best->valid || !auto_model_matches(&best->identity, p) ||
		       a->objective < best->objective);
}

/* The only candidate commit path. The journal receives a summary afterwards;
 * changing or losing that output cannot alter the retained numerical state. */
static void auto_commit_candidate(const struct simplex_Problem *p,
	struct simplex_AutoRuntime *rt, const struct simplex_AutoAttempt *a,
	const double *x, struct simplex_AutoOrigin origin,
	struct simplex_AutoHistory *history)
{
	struct simplex_AutoCandidate *best = &rt->best;

	if (best->valid && (!auto_model_matches(&best->identity, p) ||
				   p->revision == ULONG_MAX))
		best->valid = 0;
	if (auto_candidate_better(p, best, a)) {
		if (x != best->point)
			lp_simplex_memcpy(best->point, x,
				(size_t)p->columns * sizeof(double));
		best->valid = 1;
		best->identity = a->identity;
		best->objective = a->objective;
		best->result = a->result;
		best->result.objective = a->objective;
		best->result.primal_infeasibility = a->audit.primal.value;
		best->origin = origin;
	}
	rt->evidence.have_point = best->valid;
	rt->evidence.point = best->origin;
	history->have_point = best->valid;
	history->point_source = best->valid ? best->origin.record : -1;
	history->best_objective = best->valid ? best->objective : HUGE_VAL;
}

/* Enabling evidence is retained independently of the journal's capacity. */
static void auto_commit_evidence(const struct simplex_Problem *p,
	struct simplex_AutoRuntime *rt, const struct simplex_AutoAttempt *a,
	int route, struct simplex_AutoOrigin origin, int terminal)
{
	struct simplex_AutoPolicyEvidence *e = &rt->evidence;

	if (!auto_model_matches(&a->identity, p)) {
		e->have_basis = 0;
		return;
	}
	e->have_basis = rt->checkpoint.basis.valid &&
			simplex_dual_checkpoint_matches(&rt->checkpoint, p);
	if (e->have_basis)
		e->basis.event = simplex_dual_seed_source_event(
			&rt->checkpoint.basis, &p->execution);
	if (e->have_basis)
		e->basis.actor =
			simplex_execution_ref(&p->execution, e->basis.event);
	if (route == SIMPLEX_AUTO_CLASSIC && !terminal) {
		e->standard_basis_uncertified = 1;
		e->standard_basis_contract_failure =
			!a->work_valid ||
			a->report.failure == SIMPLEX_FAILURE_STAGE_CONTRACT;
		e->standard_basis = origin;
	}
	if (!a->work_valid) {
		e->retryable = origin;
		if (route == SIMPLEX_AUTO_DUAL_EXPLORE) {
			e->dual_yielded = e->dual_contract_failure = 1;
			e->dual = origin;
		}
		return;
	}
	if ((route == SIMPLEX_AUTO_DUAL_EXPLORE ||
		    route == SIMPLEX_AUTO_DUAL_FALLBACK) &&
		e->have_basis && rt->checkpoint.basis.exported) {
		e->basis.record = origin.record;
		e->basis.event = simplex_dual_seed_source_event(
			&rt->checkpoint.basis, &p->execution);
	}
	if (terminal)
		return;
	e->retryable = origin;
	if (route == SIMPLEX_AUTO_COVER) {
		e->cover_yielded = 1;
		e->cover = origin;
	}
	if (route == SIMPLEX_AUTO_DUAL_EXPLORE) {
		e->dual_yielded = 1;
		e->dual = origin;
	}
	if (simplex_auto_policy_oracle_retry(route, &a->report)) {
		if (route == SIMPLEX_AUTO_HYBRID) {
			e->native_oracle_failure = 1;
			e->native_oracle = origin;
		} else {
			e->standard_oracle_failure = 1;
			e->standard_oracle = origin;
		}
	}
}

/* One completion protocol for every action. Role-specific authority stays
 * explicit: cover needs a full certificate; an initializer can supply only a
 * point. Invocation and route selection remain in the controller below. */
static int auto_complete_attempt(const struct simplex_Problem *p,
	const struct lp_simplex_Options *o, struct simplex_AutoRuntime *rt,
	struct simplex_AutoAttempt *a, const struct simplex_AutoAction *action,
	const double *x, int proposed, struct simplex_AutoHistory *history,
	struct simplex_AutoOrigin *origin)
{
	int route = action->stage, accepted;
	int dual = route == SIMPLEX_AUTO_DUAL_EXPLORE ||
		   route == SIMPLEX_AUTO_DUAL_FALLBACK;
	int initializer = route == SIMPLEX_AUTO_INITIALIZER;

	a->raw_result = a->result;
	a->work_valid = auto_charge_stage(&a->result, &a->state, action->grant);
	if (!a->work_valid) {
		auto_reject_outputs(&a->report);
		if (dual)
			rt->checkpoint.basis.valid =
				rt->checkpoint.basis.exported =
					rt->checkpoint.basis.accepted = 0;
	}
	rt->used += a->result.iterations;
	a->feasible = auto_audit_attempt(p, x, o, a);
	a->objective = a->audit.feasible_objective.value;
	if (initializer) {
		/* No LP classification, basis, factor or dual authority. */
		a->feasible = a->feasible &&
			      a->audit.primal.quality >=
				      SIMPLEX_EVIDENCE_RECOMPUTED &&
			      a->audit.primal.value <= o->primal_tolerance;
		a->report.point_available = a->feasible;
		a->report.basis_available = a->report.factor_available = 0;
		a->report.lower_bound.quality = a->report.dual.quality =
			SIMPLEX_EVIDENCE_UNKNOWN;
	} else if (route == SIMPLEX_AUTO_COVER) {
		a->terminal = a->state == lp_simplex_EXIT_SUCCESS &&
			      a->result.status == lp_simplex_Success &&
			      a->feasible &&
			      auto_usable_diagnostic(&a->result, route) &&
			      a->result.dual_infeasibility <= o->dual_tolerance;
	} else {
		a->terminal = auto_result_terminal(&a->result, &a->state,
			a->feasible, rt->best.valid, route, rt->used);
		if (a->result.status == lp_simplex_MemoryAllocError &&
			a->report.stop == SIMPLEX_STOP_UNKNOWN)
			a->report.stop = SIMPLEX_STOP_RESOURCE;
	}
	accepted = auto_candidate_better(p, &rt->best, a);
	if (initializer) {
		if (accepted)
			a->report.stop = SIMPLEX_STOP_HANDOFF;
		else if (proposed && !a->feasible) {
			a->report.stop = SIMPLEX_STOP_NUMERICAL;
			a->report.failure = SIMPLEX_FAILURE_FEASIBILITY;
		}
		/* Contract rejection outranks a callback's claimed proposal. */
		if (!a->work_valid || !auto_model_matches(&a->identity, p))
			auto_reject_outputs(&a->report);
	}
	origin->record = auto_record_stage(p, history, action, &a->report,
		&a->result, &a->raw_result, a->terminal, a->feasible, accepted,
		&a->audit);
	origin->actor = action->actor;
	origin->event =
		simplex_execution_ref_id(&origin->actor, &p->execution, 1);
	if (dual && rt->checkpoint.basis.index && origin->record >= 0) {
		struct simplex_AutoStageRecord *record =
			history->stages + origin->record;
		record->basis_accepted = rt->checkpoint.basis.accepted;
		record->basis_exported = rt->checkpoint.basis.valid &&
					 rt->checkpoint.basis.exported;
		record->basis_producer_event =
			rt->checkpoint.basis.exported
				? simplex_dual_seed_source_event(
					  &rt->checkpoint.basis, &p->execution)
				: 0;
		record->basis_producer_actor = simplex_execution_ref(
			&p->execution, record->basis_producer_event);
	}
	auto_commit_candidate(p, rt, a, x, *origin, history);
	if (!initializer)
		auto_commit_evidence(p, rt, a, route, *origin, a->terminal);
	return accepted;
}

/* Final output is not the last kernel's report. Audit the actual returned
 * vector and describe the controller's exit, with separate point provenance.
 * Numerical phase/work/service details remain in their producing stage. */
static int auto_finish(const struct simplex_Problem *p,
	const struct simplex_ExecutionScope *root,
	const struct lp_simplex_Options *o, struct simplex_AutoRuntime *rt,
	const struct simplex_AutoAttempt *attempt,
	struct simplex_AutoOrigin cause, int exit, int restored, int state,
	const double *x, struct lp_simplex_Result *result,
	struct simplex_StageReport *output, struct simplex_AutoHistory *history)
{
	struct simplex_StageReport audit;
	struct simplex_AutoOutcome final;
	struct simplex_ExecutionScope finish_scope = *root;
	int same = auto_model_matches(&rt->input, p), j, reason, event;
	long double objective = 0.;

	lp_simplex_memset(&final, 0, sizeof(final));
	final.point.record = -1;
	final.cause = cause;
	simplex_stage_reset(&audit);
	if (!same) {
		exit = SIMPLEX_AUTO_EXIT_INPUT_CHANGED;
		result->status = lp_simplex_PrecisionError;
		result->objective = NAN;
		result->primal_infeasibility = HUGE_VAL;
		result->dual_infeasibility = NAN;
		result->dual_residual_space = lp_simplex_RESIDUAL_UNAVAILABLE;
		state = lp_simplex_EXIT_FAILURE;
		rt->best.valid = history->have_point = 0;
		history->point_source = history->bound_source = -1;
		history->best_objective = HUGE_VAL;
		history->best_bound = -HUGE_VAL;
		audit.failure = SIMPLEX_FAILURE_STAGE_CONTRACT;
	} else if (attempt) {
		simplex_stage_observe(p, o, x, NULL, 1, rt->used, &audit);
		result->primal_infeasibility = isfinite(audit.primal.value)
						       ? audit.primal.value
						       : HUGE_VAL;
		if (isfinite(audit.primal.value)) {
			for (j = 0; j < p->columns; j++)
				objective +=
					(long double)p->objective[j] * x[j];
			result->objective = (double)objective;
		} else
			result->objective = NAN;
		final.point_available =
			(restored ||
				(attempt->work_valid && attempt->feasible)) &&
			audit.feasible_objective.quality >=
				SIMPLEX_EVIDENCE_RECOMPUTED;
		if (final.point_available)
			final.point = restored ? rt->best.origin : cause;
		if (!attempt->work_valid)
			audit.failure = SIMPLEX_FAILURE_STAGE_CONTRACT;
		if (result->status == lp_simplex_Success &&
			!final.point_available) {
			result->status = lp_simplex_PrecisionError;
			state = lp_simplex_EXIT_FAILURE;
			exit = SIMPLEX_AUTO_EXIT_AUDIT;
			audit.failure = SIMPLEX_FAILURE_CERTIFICATION;
		}
	}
	result->iterations = rt->used;
	audit.stop = exit == SIMPLEX_AUTO_EXIT_CAPABILITY
			     ? SIMPLEX_STOP_NOT_APPLICABLE
			     : simplex_stage_stop_from_status(result->status);
	audit.point_available = final.point_available;
	final.exit = exit;
	final.restored = restored && final.point_available;
	final.state = state;
	final.status = result->status;
	final.work = rt->used;
	final.stop = audit.stop;
	/* Old stage indices remain local history labels. Execution references
	 * must still identify the original actors in the current journal. */
	for (j = 0; j < history->count && j < SIMPLEX_AUTO_ROUTE_COUNT + 4;
		j++) {
		struct simplex_AutoAction *action = &history->stages[j].action;
		action->execution_id =
			simplex_execution_ref_id(&action->actor, root, 1);
		auto_action_sources(action, root);
		history->stages[j].basis_producer_event =
			simplex_execution_ref_id(
				&history->stages[j].basis_producer_actor, root,
				1);
	}
	for (j = 0; j < SIMPLEX_AUTO_ROUTE_COUNT; j++)
		history->admission[j].prerequisite.event =
			simplex_execution_ref_id(
				&history->admission[j].prerequisite.actor, root,
				1);
	final.cause.event =
		simplex_execution_ref_id(&final.cause.actor, root, 1);
	final.point.event =
		simplex_execution_ref_id(&final.point.actor, root, 1);
	reason = exit == SIMPLEX_AUTO_EXIT_TERMINAL
			 ? SIMPLEX_REASON_TERMINAL_CERTIFICATE
		 : exit == SIMPLEX_AUTO_EXIT_BUDGET
			 ? SIMPLEX_REASON_TOTAL_BUDGET
		 : exit == SIMPLEX_AUTO_EXIT_INPUT_CHANGED
			 ? SIMPLEX_REASON_INPUT_CHANGED
		 : exit == SIMPLEX_AUTO_EXIT_CAPABILITY
			 ? SIMPLEX_REASON_REQUIRED_CAPABILITY
		 : exit == SIMPLEX_AUTO_EXIT_RESOURCE
			 ? SIMPLEX_REASON_REQUIRED_WORKSPACE
		 : exit == SIMPLEX_AUTO_EXIT_AUDIT
			 ? SIMPLEX_REASON_FINAL_AUDIT
			 : SIMPLEX_REASON_ROUTES_EXHAUSTED;
	simplex_execution_set_parent(&finish_scope,
		simplex_execution_ref_id(&rt->parent_actor, root, 0));
	event = simplex_execution_begin(&finish_scope,
		SIMPLEX_EXECUTION_AUTO_FINISH, exit, reason, final.cause.event);
	if (event) {
		root->history->events[event - 1].point_source =
			final.point.event;
		root->history->events[event - 1].grant = o->iteration_limit;
	}
	simplex_execution_end(root->history, event, result, &audit);
	final.execution_id = event;
	history->outcome = final;
	if (output)
		*output = audit;
	return state;
}

/* Commit diagnostic observations only after deciding from local queries.
 * A callback may corrupt output history; it cannot admit or disable a route. */
static int auto_admit_route(const struct simplex_Problem *p, int route,
	const struct simplex_AutoKernels *kernels,
	const struct simplex_AutoRoutes *routes, int condition,
	struct simplex_AutoOrigin prerequisite,
	struct simplex_AutoHistory *history)
{
	struct simplex_AutoAdmission observed;
	int binding = simplex_auto_routes_available(route, kernels, routes);
	int ready;
	prerequisite.event =
		simplex_execution_ref_id(&prerequisite.actor, &p->execution, 1);
	lp_simplex_memset(&observed, 0, sizeof(observed));
	observed.binding = binding ? SIMPLEX_AUTO_BINDING_READY
				   : SIMPLEX_AUTO_BINDING_MISSING;
	observed.model = binding ? simplex_auto_route_admission(p, route)
				 : SIMPLEX_AUTO_MODEL_UNCHECKED;
	observed.condition = condition ? SIMPLEX_AUTO_CONDITION_READY
				       : SIMPLEX_AUTO_CONDITION_ORACLE;
	observed.prerequisite.record = -1;
	observed.prerequisite.event = 0;
	if (condition && (route == SIMPLEX_AUTO_EXPLICIT ||
				 route == SIMPLEX_AUTO_STANDARD_EXPLICIT))
		observed.prerequisite = prerequisite;
	ready = binding && observed.model == SIMPLEX_AUTO_MODEL_READY &&
		condition;
	history->admission[route] = observed;
	if (getenv(SIMPLEX_ENV_PROFILE))
		fprintf(stderr,
			"auto: admission route=%s binding=%d model=%s condition=%d prerequisite=%d/%d\n",
			auto_route_name(route), observed.binding,
			simplex_auto_route_admission_name(observed.model),
			observed.condition, observed.prerequisite.record,
			observed.prerequisite.event);
	return ready;
}

static int auto_run_routes(const struct simplex_Problem *input,
	const struct lp_simplex_Options *o, double *x,
	struct lp_simplex_Result *result, int propagate_bounds,
	const struct simplex_AutoKernels *kernels,
	const struct simplex_AutoRoutes *routes,
	const struct simplex_AutoInitializer *initializer, int cover_enabled,
	struct simplex_AutoHistory *history)
{
	struct simplex_Problem problem = *input;
	const struct simplex_Problem *p = &problem;
	struct simplex_SolveContext local_context;
	struct simplex_ExecutionScope root_scope;
	struct simplex_AutoHistory local_history;
	struct simplex_AutoAction action;
	struct simplex_AutoRuntime rt;
	struct simplex_AutoAttempt attempt;
	struct simplex_AutoOrigin origin;
	struct lp_simplex_Options stage = *o;
	int j, remaining, probe, returned_best = 0, accepted;
	int owns_context = !problem.execution.context;
	int exit = SIMPLEX_AUTO_EXIT_ROUTES;
	const char *reason;
	struct simplex_AutoDecision decision;

	/* Private callers may omit the shared context. AUTO still owns one for
	 * its complete invocation; diagnostic storage is never decision state. */
	if (!problem.execution.context)
		simplex_solve_context_init(
			&local_context, NULL, &problem.execution);
	if (!problem.execution.context->configured)
		simplex_solve_context_configure(problem.execution.context, 0);
	root_scope = problem.execution;
	lp_simplex_memset(&rt, 0, sizeof(rt));
	rt.parent_actor = root_scope.parent_actor;
	rt.input = auto_model_identity(p);
	lp_simplex_memset(&origin, 0, sizeof(origin));
	origin.record = -1;
	origin.event = 0;
	rt.route = SIMPLEX_AUTO_DUAL_EXPLORE;
	rt.best.objective = HUGE_VAL;
	rt.best.origin.record = -1;
	rt.evidence.point.record = rt.evidence.basis.record = -1;
	rt.evidence.dual.record = rt.evidence.cover.record = -1;
	rt.evidence.retryable.record = -1;
	rt.evidence.native_oracle.record = rt.evidence.standard_oracle.record =
		-1;
	rt.evidence.standard_basis.record = -1;
	auto_reset(&rt.best.result);
	attempt.state = lp_simplex_EXIT_FAILURE;

	if (!history)
		history = &local_history;
	lp_simplex_memset(history, 0, sizeof(*history));
	history->best_objective = HUGE_VAL;
	history->best_bound = -HUGE_VAL;
	history->point_source = history->bound_source = -1;
	for (j = 0; j < SIMPLEX_AUTO_ROUTE_COUNT; j++)
		history->admission[j].prerequisite.record = -1;
	problem.report = &attempt.report;
	auto_reset(result);
	if (!simplex_auto_routes_available(
		    SIMPLEX_AUTO_DUAL_EXPLORE, kernels, routes)) {
		for (j = 0; j < p->columns; j++)
			x[j] = HUGE_VAL;
		result->objective = result->dual_infeasibility = NAN;
		result->dual_residual_space = lp_simplex_RESIDUAL_UNAVAILABLE;
		(void)auto_finish(p, &root_scope, o, &rt, NULL, origin,
			SIMPLEX_AUTO_EXIT_CAPABILITY, 0,
			lp_simplex_EXIT_FAILURE, x, result, input->report,
			history);
		if (owns_context)
			simplex_solve_context_destroy(&local_context);
		return lp_simplex_EXIT_FAILURE;
	}
	for (j = 0; j < SIMPLEX_AUTO_ROUTE_COUNT; j++) {
		rt.allowed[j] = auto_admit_route(p, j, kernels, routes,
			j != SIMPLEX_AUTO_EXPLICIT &&
				j != SIMPLEX_AUTO_STANDARD_EXPLICIT,
			origin, history);
		rt.recovery_count += rt.allowed[j];
	}
	rt.best.point = (double *)lp_simplex_malloc(
		((size_t)p->columns + 1) * sizeof(double));
	if (!rt.best.point) {
		result->status = lp_simplex_MemoryAllocError;
		auto_trace("finish", "memory-failure", 0, 0, 0);
		(void)auto_finish(p, &root_scope, o, &rt, NULL, origin,
			SIMPLEX_AUTO_EXIT_RESOURCE, 0, lp_simplex_EXIT_FAILURE,
			x, result, input->report, history);
		if (owns_context)
			simplex_solve_context_destroy(&local_context);
		return lp_simplex_EXIT_FAILURE;
	}
	if (simplex_auto_routes_has_checkpoint(kernels)) {
		simplex_dual_checkpoint_init(&rt.checkpoint, p, NULL, NULL);
		rt.checkpoint.basis.index = (int *)lp_simplex_malloc(
			((size_t)p->rows + 1) * sizeof(int));
		rt.checkpoint.basis.status = (unsigned char *)lp_simplex_malloc(
			(size_t)p->rows + p->columns + 1);
		if (!rt.checkpoint.basis.index || !rt.checkpoint.basis.status) {
			lp_simplex_free(rt.checkpoint.basis.index);
			lp_simplex_free(rt.checkpoint.basis.status);
			rt.checkpoint.basis.index = NULL;
			rt.checkpoint.basis.status = NULL;
		}
	}
	/* A wide positive covering LP admits restricted-column masters. Preserve
	 * half the total budget for the established routes if this optional search
	 * cannot certify the full LP; injected kernel tests keep their own schedule. */
	probe = cover_enabled ? simplex_auto_policy_cover_grant(
					p->rows, p->columns, o->iteration_limit)
			      : 0;
	if (probe > 0) {
		stage.iteration_limit = probe;
		auto_begin_action(&action, SIMPLEX_AUTO_COVER,
			stage.iteration_limit, &rt.evidence, 1);
		auto_link_action(
			&problem, &root_scope, &rt.parent_actor, &action);
		auto_start_attempt(&attempt, p);
		for (j = 0; j < p->columns; j++)
			x[j] = HUGE_VAL;
		auto_trace(
			"cover-columns", "start", stage.iteration_limit, 0, 0);
		attempt.state = simplex_auto_routes_invoke(SIMPLEX_AUTO_COVER,
			p, &stage, NULL, x, &attempt.result, propagate_bounds,
			kernels, routes, &rt.checkpoint);
		(void)auto_complete_attempt(
			p, o, &rt, &attempt, &action, x, 0, history, &origin);
		if (attempt.terminal ||
			!auto_model_matches(&attempt.identity, p)) {
			exit = attempt.terminal
				       ? SIMPLEX_AUTO_EXIT_TERMINAL
				       : SIMPLEX_AUTO_EXIT_INPUT_CHANGED;
			goto done;
		}
		auto_trace("cover-columns", "yield-without-full-certificate",
			stage.iteration_limit, rt.used, rt.best.valid);
	}
	/* Progress renews the live basis. Keep a recovery reserve, rather than
	 * imposing an absolute pivot cutoff on a progressing solve. */
	remaining = o->iteration_limit - rt.used;
	/* Each eligible recovery path and the final cold dual solve retain one
	 * observation quantum. This depends on available routes, not an opaque
	 * percentage of a possibly enormous user budget. */
	probe = simplex_auto_policy_dual_grant(remaining, rt.recovery_count);
	if (!rt.recovery_count)
		stage.algorithm = lp_simplex_ALGORITHM_DUAL_REVISED;
	stage.iteration_limit = probe;
	auto_budget_trace("dual-explore", remaining,
		rt.recovery_count ? rt.recovery_count + 1 : 0,
		remaining - probe);
	auto_start_attempt(&attempt, p);
	rt.cursor = -1;
	auto_trace("dual-explore", "start", probe, rt.used, 0);
	auto_begin_action(&action, rt.route, stage.iteration_limit,
		&rt.evidence, rt.recovery_count > 0);
	auto_link_action(&problem, &root_scope, &rt.parent_actor, &action);
	for (j = 0; j < p->columns; j++)
		x[j] = HUGE_VAL;
	attempt.state = simplex_auto_routes_invoke(rt.route, p, &stage, NULL, x,
		&attempt.result, propagate_bounds, kernels, routes,
		&rt.checkpoint);
	for (;;) {
		(void)auto_complete_attempt(
			p, o, &rt, &attempt, &action, x, 0, history, &origin);
		if (attempt.terminal) {
			exit = SIMPLEX_AUTO_EXIT_TERMINAL;
			goto done;
		}
		if (!auto_model_matches(&attempt.identity, p))
			break;
		if (getenv(SIMPLEX_ENV_PROFILE))
			fprintf(stderr,
				"auto: report route=%d stop=%s failure=%d "
				"primal=%g/%d/%d dual=%g/%d/%d "
				"bounds-epoch=%ld cost-epoch=%ld basis=%d factor=%d "
				"backend-retries=%ld factors=%ld ftran=%ld btran=%ld "
				"pdhg-stop=%s repairs=%d/%d cpu=%.6f\n",
				rt.route,
				simplex_stage_stop_name(attempt.report.stop),
				attempt.report.failure,
				attempt.report.primal.value,
				attempt.report.primal.quality,
				attempt.report.primal.space,
				attempt.report.dual.value,
				attempt.report.dual.quality,
				attempt.report.dual.space,
				attempt.report.bounds_epoch,
				attempt.report.cost_epoch,
				attempt.report.basis_available,
				attempt.report.factor_available,
				attempt.report.backend_retries,
				attempt.report.work.factorizations,
				attempt.report.work.ftran,
				attempt.report.work.btran,
				simplex_stage_stop_name(
					attempt.report.first_order_stop),
				attempt.report.work.repair_successes,
				attempt.report.work.repair_attempts,
				attempt.report.work.cpu_seconds);
		remaining = o->iteration_limit - rt.used;
		if (remaining <= 0 || rt.route == SIMPLEX_AUTO_DUAL_FALLBACK ||
			!rt.recovery_count)
			break;
		if (rt.route == SIMPLEX_AUTO_DUAL_EXPLORE && !rt.best.valid &&
			initializer && initializer->run) {
			struct lp_simplex_Options initialization = *o;
			struct simplex_AutoAttempt initialized;
			struct simplex_AutoAction initialization_action;
			int work = 0, proposed;

			initialization.iteration_limit =
				simplex_auto_policy_initializer_grant(
					remaining, rt.recovery_count + 1);
			auto_budget_trace("feasibility-initializer", remaining,
				rt.recovery_count + 1,
				remaining - initialization.iteration_limit);
			auto_begin_action(&initialization_action,
				SIMPLEX_AUTO_INITIALIZER,
				initialization.iteration_limit, &rt.evidence,
				1);
			auto_link_action(&problem, &root_scope,
				&rt.parent_actor, &initialization_action);
			/* A failed initializer may not write a point at all. Never audit
			 * uninitialized storage, and never interpret its return code as
			 * either original-model feasibility or infeasibility. */
			for (j = 0; j < p->columns; j++)
				rt.best.point[j] = HUGE_VAL;
			auto_start_attempt(&initialized, p);
			problem.report = &initialized.report;
			proposed = simplex_auto_routes_initialize(initializer,
				p, &initialization, rt.best.point, &work);
			problem.report = &attempt.report;
			initialized.result.iterations = work;
			initialized.result.dual_infeasibility = NAN;
			initialized.result.dual_residual_space =
				lp_simplex_RESIDUAL_UNAVAILABLE;
			accepted = auto_complete_attempt(p, o, &rt,
				&initialized, &initialization_action,
				rt.best.point, proposed, history, &origin);
			work = initialized.result.iterations;
			remaining = o->iteration_limit - rt.used;
			if (accepted)
				/* A validated point supersedes singleton feasibility work. */
				rt.cursor = simplex_auto_policy_seed_cursor(
					rt.cursor, accepted);
			auto_trace("feasibility-initializer",
				accepted ? "validated-seed"
					 : "no-validated-seed",
				work, rt.used, rt.best.valid);
			if (!auto_model_matches(&initialized.identity, p)) {
				attempt.result.status =
					lp_simplex_PrecisionError;
				attempt.state = lp_simplex_EXIT_FAILURE;
				break;
			}
			if (remaining <= 0)
				break;
		}
		/* The kernel supplies a stop cause. Short work alone is not evidence
		 * of a voluntary yield, nor does a generic precision error identify
		 * an implicit-oracle defect. Unknown reports retain that distinction. */
		reason = simplex_stage_stop_name(attempt.report.stop);
		if (rt.route == SIMPLEX_AUTO_HYBRID)
			rt.allowed[SIMPLEX_AUTO_EXPLICIT] = auto_admit_route(p,
				SIMPLEX_AUTO_EXPLICIT, kernels, routes,
				simplex_auto_policy_oracle_retry(
					rt.route, &attempt.report),
				rt.evidence.native_oracle, history);
		if (rt.route == SIMPLEX_AUTO_STANDARD_IMPLICIT)
			rt.allowed[SIMPLEX_AUTO_STANDARD_EXPLICIT] =
				auto_admit_route(p,
					SIMPLEX_AUTO_STANDARD_EXPLICIT, kernels,
					routes,
					simplex_auto_policy_oracle_retry(
						rt.route, &attempt.report),
					rt.evidence.standard_oracle, history);
		/* Policy returns a proposal; only the controller advances the cursor
		 * and installs the next lease. No history or candidates enter policy. */
		decision = simplex_auto_policy_next(
			rt.allowed, rt.cursor, remaining);
		rt.cursor = decision.index;
		rt.route = decision.route;
		stage = *o;
		stage.iteration_limit = decision.grant;
		if (rt.route == SIMPLEX_AUTO_DUAL_FALLBACK) {
			stage.algorithm = lp_simplex_ALGORITHM_DUAL_REVISED;
			auto_trace(
				"dual-fallback", reason, remaining, rt.used, 0);
		} else {
			if (rt.route == SIMPLEX_AUTO_HYBRID)
				auto_budget_trace("native-hybrid", remaining,
					decision.recovery_paths,
					decision.reserve);
			auto_trace(auto_route_name(rt.route), reason,
				stage.iteration_limit, rt.used, rt.best.valid);
		}
		auto_begin_action(&action, rt.route, stage.iteration_limit,
			&rt.evidence, rt.route != SIMPLEX_AUTO_DUAL_FALLBACK);
		auto_link_action(
			&problem, &root_scope, &rt.parent_actor, &action);
		auto_start_attempt(&attempt, p);
		for (j = 0; j < p->columns; j++)
			x[j] = HUGE_VAL;
		attempt.state = simplex_auto_routes_invoke(rt.route, p, &stage,
			rt.best.valid ? rt.best.point : NULL, x,
			&attempt.result, propagate_bounds, kernels, routes,
			&rt.checkpoint);
	}
	if (rt.best.valid) {
		lp_simplex_memcpy(
			x, rt.best.point, (size_t)p->columns * sizeof(double));
		attempt.result = rt.best.result;
		returned_best = 1;
	}
	auto_trace("recovery",
		rt.used >= o->iteration_limit ? "total-budget-exhausted"
					      : "admissible-routes-exhausted",
		0, rt.used, rt.best.valid);
	attempt.result.status = rt.used >= o->iteration_limit
					? lp_simplex_ExceedIterLimit
					: lp_simplex_PrecisionError;
	exit = rt.used >= o->iteration_limit ? SIMPLEX_AUTO_EXIT_BUDGET
					     : SIMPLEX_AUTO_EXIT_ROUTES;
	attempt.state = lp_simplex_EXIT_FAILURE;
done:
	if (rt.best.valid && !auto_model_matches(&rt.best.identity, p))
		rt.best.valid = 0;
	if (attempt.result.status != lp_simplex_Success && rt.best.valid &&
		(!attempt.feasible || attempt.objective > rt.best.objective ||
			attempt.result.status == lp_simplex_MemoryAllocError)) {
		lp_simplex_memcpy(
			x, rt.best.point, (size_t)p->columns * sizeof(double));
		returned_best = 1;
		attempt.result.objective = rt.best.objective;
		attempt.result.primal_infeasibility =
			rt.best.result.primal_infeasibility;
		attempt.result.dual_infeasibility =
			rt.best.result.dual_infeasibility;
		attempt.result.dual_residual_space =
			rt.best.result.dual_residual_space;
	}
	*result = attempt.result;
	attempt.state = auto_finish(p, &root_scope, o, &rt, &attempt, origin,
		exit, returned_best, attempt.state, x, result, input->report,
		history);
	auto_trace("finish", auto_reason(result->status), 0, rt.used, 0);
	lp_simplex_free(rt.checkpoint.basis.index);
	lp_simplex_free(rt.checkpoint.basis.status);
	lp_simplex_free(rt.best.point);
	if (owns_context)
		simplex_solve_context_destroy(&local_context);
	return attempt.state;
}

int simplex_auto_solve(const struct simplex_Problem *p,
	const struct lp_simplex_Options *o, double *x,
	struct lp_simplex_Result *result, int propagate_bounds,
	const struct simplex_AutoBindings *bindings,
	const struct simplex_AutoInitializer *initializer,
	struct simplex_AutoHistory *history)
{
	struct simplex_AutoKernels kernels = {NULL, NULL};
	struct simplex_AutoRoutes routes = {NULL};
	int cover_enabled = 1;

	if (!bindings)
		simplex_auto_routes_defaults(&kernels, &routes);
	else {
		/* Freeze the callable set for this invocation before callbacks run.
		 * Missing injected bindings remain absent capabilities. */
		if (bindings->kernels)
			kernels = *bindings->kernels;
		if (bindings->routes)
			routes = *bindings->routes;
		cover_enabled = bindings->cover_enabled;
	}
	return auto_run_routes(p, o, x, result, propagate_bounds, &kernels,
		&routes, initializer, cover_enabled, history);
}
