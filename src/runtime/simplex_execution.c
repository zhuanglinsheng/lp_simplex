/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_execution.c
 * @brief Bounded execution journal and validated causal references.
 */

#include "simplex_stage.h"
#include "simplex_solve_context.h"
#include "simplex_execution.h"
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

struct simplex_ExecutionRef simplex_execution_ref(
	const struct simplex_ExecutionScope *scope, int id)
{
	struct simplex_ExecutionRef reference;
	const struct simplex_ExecutionEvent *event;
	memset(&reference, 0, sizeof(reference));
	if (!scope || !scope->history || !scope->context || id <= 0 ||
		id > scope->history->count || id > SIMPLEX_EXECUTION_CAPACITY)
		return reference;
	event = scope->history->events + id - 1;
	if (!event->serial || event->context != scope->context)
		return reference;
	reference.history = scope->history;
	reference.context = scope->context;
	reference.serial = event->serial;
	reference.id = id;
	reference.kind = event->kind;
	reference.detail = event->detail;
	reference.model = event->model;
	return reference;
}

int simplex_execution_ref_id(const struct simplex_ExecutionRef *reference,
	const struct simplex_ExecutionScope *scope, int completed)
{
	const struct simplex_ExecutionEvent *event;
	if (!scope || !scope->history || !reference->serial ||
		reference->history != scope->history ||
		reference->context != scope->context || reference->id <= 0 ||
		reference->id > scope->history->count ||
		reference->id > SIMPLEX_EXECUTION_CAPACITY)
		return 0;
	event = scope->history->events + reference->id - 1;
	return event->context == reference->context &&
			       event->serial == reference->serial &&
			       event->kind == reference->kind &&
			       event->detail == reference->detail &&
			       event->model == reference->model &&
			       (!completed || event->complete)
		       ? reference->id
		       : 0;
}

void simplex_execution_set_parent(struct simplex_ExecutionScope *scope, int id)
{
	scope->parent = id;
	scope->parent_actor = simplex_execution_ref(scope, id);
}

int simplex_execution_parent(const struct simplex_ExecutionScope *scope)
{
	return scope && scope->parent > 0 &&
			       scope->parent == scope->parent_actor.id
		       ? simplex_execution_ref_id(
				 &scope->parent_actor, scope, 0)
		       : 0;
}

void simplex_execution_inherit_parent(struct simplex_ExecutionScope *target,
	const struct simplex_ExecutionScope *source)
{
	target->parent = source->parent;
	target->parent_actor = source->parent_actor;
}

int simplex_execution_begin(const struct simplex_ExecutionScope *scope,
	int kind, int detail, int reason, int cause)
{
	struct simplex_Execution *h = scope->history;
	struct simplex_ExecutionEvent *event;
	if (!h)
		return 0;
	if (h->count == SIMPLEX_EXECUTION_CAPACITY) {
		if (h->dropped < INT_MAX)
			h->dropped++;
		return 0;
	}
	event = h->events + h->count++;
	memset(event, 0, sizeof(*event));
	event->context = scope->context;
	if (scope->context && scope->context->execution_serial < ULONG_MAX)
		event->serial = ++scope->context->execution_serial;
	event->kind = kind;
	event->detail = detail;
	event->reason = reason;
	event->model = event->source_model = scope->model;
	event->parent = simplex_execution_parent(scope);
	event->cause = cause;
	event->primal.value = event->primal_audit.value = NAN;
	return h->count;
}

void simplex_execution_end(struct simplex_Execution *h, int id,
	const struct lp_simplex_Result *result,
	const struct simplex_StageReport *r)
{
	struct simplex_ExecutionEvent *event;
	if (!h || id <= 0 || id > h->count)
		return;
	event = h->events + id - 1;
	event->complete = 1;
	if (result) {
		event->result_available = 1;
		event->status = result->status;
		event->work = result->iterations;
	}
	if (r) {
		event->tableau = r->tableau;
		event->singleton = r->singleton;
		event->stop = r->stop;
		event->failure = r->failure;
		event->basis_import_attempted = r->basis_import_attempted;
		event->basis_rejection = r->basis_rejection;
		event->basis_exported = r->basis_available;
		event->recovery_skipped = r->recovery_skipped;
		event->bounds_epoch = r->bounds_epoch;
		event->cost_epoch = r->cost_epoch;
		event->primal = r->primal;
	}
}

static const char *execution_kind_name(int kind)
{
	static const char *const names[] = {"unknown", "solve",
		"presolve-model", "reduced-solve", "scaling-model", "auto",
		"backend", "stable-initialize", "postsolve-audit",
		"original-retry", "point-selection", "recovery-skip",
		"pan-service", "pan-model", "column-model", "cover-master",
		"column-handoff", "auto-finish", "point-archive",
		"result-audit", "tableau", "tableau-model", "tableau-standard",
		"singleton-dual", "singleton-dual-model", "pan-trigger",
		"pan-face-opportunity"};
	return kind >= 0 && (size_t)kind < sizeof(names) / sizeof(names[0])
		       ? names[kind]
		       : names[0];
}

static const char *execution_reason_name(int reason)
{
	static const char *const names[] = {"none", "numerical-failure",
		"invalid-postsolve", "startup-failure",
		"original-feasible-candidate", "keep-retry-point",
		"restore-saved-point", "restore-archived-point",
		"repeated-failure", "dual-handoff", "implicit-oracle-failure",
		"basis-rebuild", "alternative-representation",
		"recovery-exhausted", "cover-uncertified",
		"alternative-algorithm", "stage-contract-failure",
		"column-admission", "dual-refinement", "terminal-certificate",
		"total-work-budget", "routes-exhausted", "input-changed",
		"required-workspace", "final-output-audit",
		"infeasibility-conflict", "required-capability",
		"alternative-direction", "pdhg-stagnation", "face-initial",
		"face-low-progress", "face-after-multiplier-search",
		"primal-margin", "pdhg-primal-stagnation",
		"face-pdhg-improvement"};
	return reason >= 0 && (size_t)reason < sizeof(names) / sizeof(names[0])
		       ? names[reason]
		       : "unknown";
}

static const char *execution_face_action_name(int action)
{
	static const char *const names[] = {"not-run", "no-change", "adjusted",
		"exact-request-skip", "optional-resource-refused"};
	return action >= 0 && (size_t)action < sizeof(names) / sizeof(names[0])
		       ? names[action]
		       : names[0];
}

void simplex_execution_print(const struct simplex_Execution *h)
{
	int i;
	for (i = 0; i < h->count; i++) {
		const struct simplex_ExecutionEvent *e = h->events + i;
		fprintf(stderr,
			"execution: id=%d kind=%s detail=%d model=%d source-model=%d "
			"parent=%d cause=%d point-source=%d basis-source=%d grant=%d reason=%s complete=%d result=%d status=%d work=%d "
			"stop=%s failure=%d basis-import=%d basis-rejection=%d basis-export=%d recovery-skipped=%d bounds-epoch=%ld cost-epoch=%ld primal=%g/%d audit=%g/%d serial=%lu\n",
			i + 1, execution_kind_name(e->kind), e->detail,
			e->model, e->source_model, e->parent, e->cause,
			e->point_source, e->basis_source, e->grant,
			execution_reason_name(e->reason), e->complete,
			e->result_available, e->status, e->work,
			simplex_stage_stop_name(e->stop), e->failure,
			e->basis_import_attempted, e->basis_rejection,
			e->basis_exported, e->recovery_skipped, e->bounds_epoch,
			e->cost_epoch, e->primal.value, e->primal.quality,
			e->primal_audit.value, e->primal_audit.quality,
			e->serial);
		if (e->singleton.available)
			fprintf(stderr,
				"execution-singleton: id=%d input-model=%d auxiliary-model=%d auxiliary-result=%d auxiliary-status=%d mapping=%d structure=%d point-checked=%d\n",
				i + 1, e->singleton.input_model,
				e->singleton.auxiliary_model,
				e->singleton.auxiliary_available,
				e->singleton.auxiliary_status,
				e->singleton.mapping,
				e->singleton.structure_verified,
				e->singleton.input_point_checked);
		if (e->tableau.available)
			fprintf(stderr,
				"execution-tableau: id=%d phase=%d phase1-work=%d cleanup-work=%d phase2-work=%d point=%d input-model=%d standard-model=%d rows=%d columns=%d converted=%d recovered=%d\n",
				i + 1, e->tableau.phase,
				e->tableau.phase_one_work,
				e->tableau.cleanup_work,
				e->tableau.phase_two_work,
				e->tableau.point_available,
				e->tableau.input_model,
				e->tableau.standard_model,
				e->tableau.standard_rows,
				e->tableau.standard_columns,
				e->tableau.converted,
				e->tableau.recovered_point);
		if (e->kind == SIMPLEX_EXECUTION_SOLVE)
			fprintf(stderr, "execution-config: mode=%s\n",
				e->experimental_environment
					? "experimental-environment"
					: "production");
		if (e->kind == SIMPLEX_EXECUTION_POINT_ARCHIVE)
			fprintf(stderr,
				"execution-archive: id=%d stage=%s producer-work=%d transforms=%d rows=%d columns=%d\n",
				i + 1, e->observation_stage,
				e->observation_work, e->detail, e->model_rows,
				e->model_columns);
		if (e->kind == SIMPLEX_EXECUTION_PAN_SERVICE)
			fprintf(stderr,
				"execution-service: id=%d role=%s code=%d output=%s valid=%d raw-stop=%s reported-end=%d reported-work=%d relative=%d implementation=%d rows=%d columns=%d\n",
				i + 1, simplex_pan_service_name(e->detail),
				e->service_code,
				simplex_pan_service_output_name(
					e->service_output),
				e->service_valid,
				simplex_stage_stop_name(e->service_raw_stop),
				e->service_reported_end,
				e->service_reported_work,
				e->service_relative_work,
				e->service_implementation, e->model_rows,
				e->model_columns);
		if (e->kind == SIMPLEX_EXECUTION_PAN_TRIGGER &&
			(e->reason == SIMPLEX_REASON_PDHG_STAGNATION ||
				e->reason ==
					SIMPLEX_REASON_PDHG_PRIMAL_STAGNATION ||
				e->reason == SIMPLEX_REASON_PRIMAL_MARGIN))
			fprintf(stderr,
				"execution-trigger: id=%d producer-work=%d metric-epoch=%ld working-primal=%g quality=%d\n",
				i + 1, e->observation_work, e->primal.epoch,
				e->primal.value, e->primal.quality);
		if ((e->kind == SIMPLEX_EXECUTION_PAN_FACE_OPPORTUNITY ||
			    (e->kind == SIMPLEX_EXECUTION_PAN_SERVICE &&
				    e->detail ==
					    SIMPLEX_PAN_SERVICE_FACE_RECOVERY)) &&
			e->face.opportunity)
			fprintf(stderr,
				"execution-face: id=%d opportunity=%d producer-work=%d steps=%d low-progress=%d action=%s changed=%d remembered=%d\n",
				i + 1, e->face.opportunity, e->observation_work,
				e->face.steps, e->face.low_progress,
				execution_face_action_name(
					e->face.result.action),
				e->face.result.changed,
				e->face.result.remembered);
		if (e->kind == SIMPLEX_EXECUTION_COLUMN_HANDOFF)
			fprintf(stderr,
				"execution-handoff: id=%d role=%d rejection=%d rows=%d columns=%d\n",
				i + 1, e->detail, e->handoff_rejection,
				e->model_rows, e->model_columns);
	}
	if (h->dropped)
		fprintf(stderr, "execution: dropped=%d\n", h->dropped);
}
