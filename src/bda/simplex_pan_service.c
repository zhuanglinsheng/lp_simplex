/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_pan_service.c
 * @brief BDA subprocess roles, model bindings and bounded work handoffs.
 */

#include "simplex_solve_context.h"
#include "simplex_execution.h"
#include "simplex_pan_service.h"
#include "simplex_pan_config.h"
#include "simplex_pan_basis.h"
#include "simplex_pan_bounded.h"
#include "simplex_pan_face_recovery.h"
#include "simplex_pan_cone.h"
#include "simplex_pan_recovery.h"
#include "simplex_pan_phase1.h"
#include "simplex_dual.h"
#include "simplex_stage.h"
#include "utils.h"
#include <limits.h>
#include <math.h>
#include <string.h>

const char *simplex_pan_service_name(int service)
{
	static const char *names[] = {"nnls-initialize", "dual-feasibility",
		"pdhg-advance", "cone-direction", "feasibility-restore",
		"dual-refine", "direction-predict", "native-face-recovery"};
	return service >= 0 && service < SIMPLEX_PAN_SERVICE_COUNT
		       ? names[service]
		       : "unknown";
}

const char *simplex_pan_service_output_name(int output)
{
	static const char *names[] = {"none", "nonnegative-equality-point",
		"input-feasible-point", "native-feasible-point",
		"direction-proposal", "stationary-residual",
		"checked-direction", "checked-dual", "native-kkt",
		"native-face-adjustment"};
	return output >= 0 && output <=
				       SIMPLEX_PAN_OUTPUT_NATIVE_FACE_ADJUSTMENT
		       ? names[output]
		       : "unknown";
}

static void service_bind_model(struct simplex_PanStandard *a,
	struct simplex_PanServiceContext *c, int source_model)
{
	struct simplex_ExecutionScope scope;
	int count;
	a->services = c;
	a->service_revision = 0;
	a->service_model = 0;
	if (!c)
		return;
	if (c->representations == ULONG_MAX) {
		c->failed = 1;
		return;
	}
	a->service_revision = ++c->representations;
	scope = c->scope;
	scope.model = source_model;
	scope.parent = c->parent;
	scope.parent_actor = c->parent_actor;
	count = scope.history ? scope.history->count : 0;
	scope = simplex_solve_context_derive(
		&scope, SIMPLEX_EXECUTION_PAN_MODEL);
	a->service_model = scope.model;
	if (scope.history && scope.history->count > count) {
		scope.history->events[count].model_rows = a->rows;
		scope.history->events[count].model_columns = a->columns;
	}
}

void simplex_pan_service_bind(
	struct simplex_PanStandard *a, struct simplex_PanServiceContext *c)
{
	service_bind_model(a, c, c ? c->scope.model : 0);
}

void simplex_pan_service_derive(
	struct simplex_PanStandard *a, const struct simplex_PanStandard *source)
{
	service_bind_model(a, source->services, source->service_model);
}

void simplex_pan_service_context_init(
	struct simplex_PanServiceContext *c, const struct simplex_Problem *p)
{
	memset(c, 0, sizeof(*c));
	c->input = p;
	if (p) {
		c->scope = p->execution;
		c->identity = simplex_problem_identity(p);
		c->parent = p->execution.parent;
		c->parent_actor = p->execution.parent_actor;
	}
}

int simplex_pan_service_context_valid(const struct simplex_PanServiceContext *c)
{
	const struct simplex_Problem *p;
	if (!c)
		return 1;
	if (c->failed)
		return 0;
	if (!c->input)
		return !c->identity.input;
	p = c->input;
	return simplex_problem_identity_matches(&c->identity, p) &&
	       c->scope.context == c->identity.context &&
	       c->scope.model == c->identity.model;
}

static int service_model_matches(const struct simplex_PanServiceFrame *f)
{
	const struct simplex_PanStandard *a = f->outcome.representation;
	const struct simplex_PanStandard *b = &f->identity;
	return a && a->rows == b->rows && a->columns == b->columns &&
	       a->original_columns == b->original_columns &&
	       a->service_revision == b->service_revision &&
	       a->service_model == b->service_model &&
	       a->nonzeros == b->nonzeros && a->value == b->value &&
	       a->normalized_value == b->normalized_value &&
	       a->column_norm == b->column_norm && a->rhs == b->rhs &&
	       a->objective == b->objective &&
	       a->column_start == b->column_start &&
	       a->row_index == b->row_index &&
	       a->original_column == b->original_column &&
	       a->original_scale == b->original_scale &&
	       a->original_shift == b->original_shift &&
	       a->objective_offset == b->objective_offset &&
	       a->audit_objective == b->audit_objective &&
	       a->audit_objective_offset == b->audit_objective_offset &&
	       a->services == f->context &&
	       simplex_pan_service_context_valid(f->context) &&
	       (!f->context || simplex_problem_identity_equal(
				       &f->input, &f->context->identity)) &&
	       (!f->bounded ||
		       (&f->bounded->a == a &&
			       f->bounded->problem == f->bounded_input &&
			       f->bounded->lower == f->lower &&
			       f->bounded->upper == f->upper)) &&
	       (!f->basis || (f->basis->standard == a && f->basis->count >= 0 &&
				     f->basis->count <= a->rows &&
				     f->basis->count <= f->basis->capacity));
}

int simplex_pan_service_begin(struct simplex_PanServiceFrame *f,
	const struct simplex_PanStandard *a, int service, int begin, int limit,
	int requested)
{
	struct simplex_ExecutionScope scope;
	int counter_valid =
		limit >= 0 && begin >= 0 && begin <= limit && requested >= 0;

	memset(f, 0, sizeof(*f));
	f->outcome.service = service;
	f->outcome.representation = a;
	f->outcome.before =
		begin >= 0 && limit >= 0 && begin <= limit ? begin : 0;
	f->outcome.stop = SIMPLEX_STOP_UNKNOWN;
	if (!a)
		return 0;
	f->context = a->services;
	f->sink = a->report;
	f->outcome.rows = a->rows;
	f->outcome.working_model = a->service_model;
	f->outcome.representation_revision = a->service_revision;
	f->outcome.columns = a->columns;
	f->identity = *a;
	if (f->context)
		f->input = f->context->identity;
	f->outcome.valid = counter_valid && service >= 0 &&
			   service < SIMPLEX_PAN_SERVICE_COUNT &&
			   simplex_pan_service_context_valid(f->context);
	if (f->outcome.valid)
		f->outcome.grant =
			requested < limit - begin ? requested : limit - begin;
	if (f->context) {
		f->outcome.input = f->context->input;
		f->outcome.revision = f->context->identity.revision;
		f->outcome.input_model = f->context->scope.model;
		f->saved_parent = f->context->parent;
		f->saved_parent_actor = f->context->parent_actor;
		scope = f->context->scope;
		scope.parent = f->saved_parent;
		scope.parent_actor = f->saved_parent_actor;
		scope.model = service == SIMPLEX_PAN_SERVICE_DUAL_FEASIBILITY
				      ? f->outcome.input_model
				      : a->service_model;
		f->outcome.event = simplex_execution_begin(&scope,
			SIMPLEX_EXECUTION_PAN_SERVICE, service,
			SIMPLEX_REASON_NONE, 0);
		f->actor = simplex_execution_ref(&scope, f->outcome.event);
		if (f->outcome.event) {
			struct simplex_ExecutionEvent *e =
				scope.history->events + f->outcome.event - 1;
			e->source_model = f->outcome.input_model;
			e->grant = f->outcome.grant;
			e->model_rows = a->rows;
			e->model_columns = a->columns;
		}
		f->context->parent = f->outcome.event;
		f->context->parent_actor = f->actor;
	}
	return f->outcome.valid;
}

static int service_output_allowed(int service, int output)
{
	if (output == SIMPLEX_PAN_OUTPUT_NONE)
		return 1;
	switch (service) {
	case SIMPLEX_PAN_SERVICE_NNLS:
		return output == SIMPLEX_PAN_OUTPUT_NONNEGATIVE_EQUALITY_POINT;
	case SIMPLEX_PAN_SERVICE_DUAL_FEASIBILITY:
		return output == SIMPLEX_PAN_OUTPUT_INPUT_FEASIBLE_POINT;
	case SIMPLEX_PAN_SERVICE_PDHG:
		return output == SIMPLEX_PAN_OUTPUT_NATIVE_KKT ||
		       output == SIMPLEX_PAN_OUTPUT_NATIVE_FEASIBLE_POINT;
	case SIMPLEX_PAN_SERVICE_RESTORE:
		return output == SIMPLEX_PAN_OUTPUT_NATIVE_FEASIBLE_POINT ||
		       output == SIMPLEX_PAN_OUTPUT_NONNEGATIVE_EQUALITY_POINT;
	case SIMPLEX_PAN_SERVICE_FACE_RECOVERY:
		return output == SIMPLEX_PAN_OUTPUT_NATIVE_FACE_ADJUSTMENT;
	case SIMPLEX_PAN_SERVICE_CONE:
		return output == SIMPLEX_PAN_OUTPUT_DIRECTION_PROPOSAL ||
		       output == SIMPLEX_PAN_OUTPUT_STATIONARY_RESIDUAL ||
		       output == SIMPLEX_PAN_OUTPUT_CHECKED_DIRECTION;
	case SIMPLEX_PAN_SERVICE_DUAL_REFINE:
		return output == SIMPLEX_PAN_OUTPUT_CHECKED_DUAL;
	case SIMPLEX_PAN_SERVICE_PREDICT:
		return output == SIMPLEX_PAN_OUTPUT_CHECKED_DUAL ||
		       output == SIMPLEX_PAN_OUTPUT_CHECKED_DIRECTION;
	default:
		return 0;
	}
}

static int service_trigger_reason(int trigger)
{
	switch (trigger) {
	case SIMPLEX_PAN_TRIGGER_PDHG_STAGNATION:
		return SIMPLEX_REASON_PDHG_STAGNATION;
	case SIMPLEX_PAN_TRIGGER_PRIMAL_MARGIN:
		return SIMPLEX_REASON_PRIMAL_MARGIN;
	case SIMPLEX_PAN_TRIGGER_PDHG_PRIMAL_STAGNATION:
		return SIMPLEX_REASON_PDHG_PRIMAL_STAGNATION;
	default:
		return SIMPLEX_REASON_NONE;
	}
}

static int service_face_reason(int opportunity)
{
	switch (opportunity) {
	case SIMPLEX_PAN_FACE_INITIAL:
		return SIMPLEX_REASON_FACE_INITIAL;
	case SIMPLEX_PAN_FACE_LOW_PROGRESS:
		return SIMPLEX_REASON_FACE_LOW_PROGRESS;
	case SIMPLEX_PAN_FACE_AFTER_MULTIPLIER_SEARCH:
		return SIMPLEX_REASON_FACE_MULTIPLIER_SEARCH;
	case SIMPLEX_PAN_FACE_PDHG_IMPROVEMENT:
		return SIMPLEX_REASON_FACE_PDHG_IMPROVEMENT;
	default:
		return SIMPLEX_REASON_NONE;
	}
}

int simplex_pan_service_finish(struct simplex_PanServiceFrame *f,
	int *iterations, int code, int status, int stop, int output,
	struct simplex_PanServiceOutcome *out)
{
	struct simplex_PanServiceOutcome *r = &f->outcome;
	int entered = r->valid, end = r->before + r->grant;
	int work_valid = *iterations >= r->before && *iterations <= end;

	r->reported_end = *iterations;
	if (!r->relative_work)
		r->reported_work =
			work_valid ? *iterations - r->before : *iterations;
	r->raw_code = code;
	r->raw_status = status;
	r->raw_stop = stop;
	r->raw_output = output;
	r->valid = entered && work_valid && service_model_matches(f) &&
		   service_output_allowed(r->service, output) &&
		   stop >= SIMPLEX_STOP_UNKNOWN &&
		   stop <= SIMPLEX_STOP_NOT_APPLICABLE;
	r->end = !entered ? r->before : work_valid ? *iterations : end;
	*iterations = r->end;
	r->stop = r->valid ? stop : SIMPLEX_STOP_NUMERICAL;
	r->output = r->valid ? output : SIMPLEX_PAN_OUTPUT_NONE;
	if (f->context) {
		struct simplex_Execution *h = f->context->scope.history;
		f->context->parent = f->saved_parent;
		f->context->parent_actor = f->saved_parent_actor;
		r->event = simplex_execution_ref_id(
			&f->actor, &f->context->scope, 0);
		if (!r->valid)
			f->context->failed = 1;
		if (r->event) {
			struct simplex_ExecutionEvent *e =
				h->events + r->event - 1;
			e->complete = e->result_available = 1;
			e->status =
				r->valid ? status : lp_simplex_PrecisionError;
			e->work = r->end - r->before;
			e->stop = r->stop;
			e->failure = r->valid ? SIMPLEX_FAILURE_UNKNOWN
					      : SIMPLEX_FAILURE_STAGE_CONTRACT;
			e->service_output = r->output;
			e->service_code = code;
			e->service_valid = r->valid;
			e->service_raw_stop = stop;
			e->service_reported_end = r->reported_end;
			e->service_reported_work = r->reported_work;
			e->service_relative_work = r->relative_work;
			e->service_implementation = r->implementation;
			if (r->service == SIMPLEX_PAN_SERVICE_RESTORE &&
				service_trigger_reason(r->trigger.kind) !=
					SIMPLEX_REASON_NONE) {
				e->reason =
					service_trigger_reason(r->trigger.kind);
				e->cause = simplex_execution_ref_id(
					&f->trigger_actor, &f->context->scope,
					1);
			}
			if (r->service == SIMPLEX_PAN_SERVICE_FACE_RECOVERY) {
				e->face = r->face;
				e->observation_work = r->before;
				e->reason = service_face_reason(
					r->face.opportunity);
				e->cause = simplex_execution_ref_id(
					&f->trigger_actor, &f->context->scope,
					1);
				e->recovery_skipped =
					r->face.result.action ==
					SIMPLEX_PAN_FACE_EXACT_SKIP;
			}
		}
	}
	if (f->sink && r->service >= 0 &&
		r->service < SIMPLEX_PAN_SERVICE_COUNT) {
		struct simplex_PanServiceSummary *s = &f->sink->services;
		int k = r->service, work = r->end - r->before;
		if (s->calls[k] < LONG_MAX)
			s->calls[k]++;
		s->work[k] = s->work[k] <= LONG_MAX - work ? s->work[k] + work
							   : LONG_MAX;
		s->stop[k] = r->stop;
		s->output[k] = r->output;
		s->implementation[k] = r->implementation;
		if (!r->valid && s->invalid_contracts < LONG_MAX)
			s->invalid_contracts++;
	}
	r->actor = f->actor;
	if (out)
		*out = *r;
	return r->valid;
}

static int service_remaining(int limit, int work)
{
	return work >= 0 && work <= limit ? limit - work : 0;
}

int simplex_pan_service_nnls(struct simplex_PanBasis *basis, double *point,
	double primal_tolerance, double dual_tolerance, int limit,
	int *iterations, int *status, struct simplex_PanServiceOutcome *out)
{
	struct simplex_PanServiceFrame f;
	int code = lp_simplex_EXIT_FAILURE, output = SIMPLEX_PAN_OUTPUT_NONE;
	if (simplex_pan_service_begin(&f, basis->standard,
		    SIMPLEX_PAN_SERVICE_NNLS, *iterations, limit,
		    service_remaining(limit, *iterations)))
		code = simplex_pan_phase1_run(basis, point, primal_tolerance,
			dual_tolerance, f.outcome.before + f.outcome.grant,
			iterations, status);
	else
		*status = lp_simplex_PrecisionError;
	f.basis = basis;
	if (code == lp_simplex_EXIT_SUCCESS)
		output = SIMPLEX_PAN_OUTPUT_NONNEGATIVE_EQUALITY_POINT;
	if (!simplex_pan_service_finish(&f, iterations, code, *status,
		    code == 0 ? SIMPLEX_STOP_HANDOFF
			      : simplex_stage_stop_from_status(*status),
		    output, out)) {
		*status = lp_simplex_PrecisionError;
		return lp_simplex_EXIT_FAILURE;
	}
	return code;
}

int simplex_pan_service_dual_feasible(struct simplex_PanServiceContext *context,
	const struct simplex_Problem *input, const struct lp_simplex_Options *o,
	double *point, int *iterations, struct lp_simplex_Result *result,
	struct simplex_PanServiceOutcome *out)
{
	struct simplex_PanServiceFrame f;
	struct simplex_PanStandard view;
	struct simplex_ExecutionScope scope = input->execution;
	struct simplex_StageReport raw_report;
	struct lp_simplex_Options lease = *o;
	int code = lp_simplex_EXIT_FAILURE, ready;

	lp_simplex_memset(&view, 0, sizeof(view));
	view.services = context;
	view.service_model = input->execution.model;
	view.rows = input->rows;
	view.columns = input->columns;
	view.nonzeros = input->matrix.nonzeros;
	view.value = input->matrix.value;
	view.rhs = input->rhs;
	view.objective = input->objective;
	view.report = input->report;
	view.column_start = input->matrix.column_start;
	view.row_index = input->matrix.row_index;
	ready = simplex_pan_service_begin(&f, &view,
		SIMPLEX_PAN_SERVICE_DUAL_FEASIBILITY, *iterations,
		o->iteration_limit,
		service_remaining(o->iteration_limit, *iterations));
	lp_simplex_memset(result, 0, sizeof(*result));
	result->status = lp_simplex_PrecisionError;
	if (ready && context && context->input == input) {
		lease.iteration_limit = f.outcome.grant;
		lease.algorithm = lp_simplex_ALGORITHM_DUAL_REVISED;
		lease.pricing = lp_simplex_PRICING_DUAL_STEEPEST_EDGE;
		simplex_stage_reset(&raw_report);
		simplex_execution_set_parent(&scope, f.outcome.event);
		code = simplex_dual_find_feasible_reported(
			input, &lease, point, result, &scope, &raw_report);
		if (result->iterations >= 0 &&
			result->iterations <= f.outcome.grant)
			*iterations = f.outcome.before + result->iterations;
		else
			*iterations = -1;
	} else
		f.outcome.valid = 0;
	f.outcome.relative_work = 1;
	f.outcome.reported_work = result->iterations;
	ready = simplex_pan_service_finish(&f, iterations, code, result->status,
		code == 0 ? SIMPLEX_STOP_HANDOFF
			  : simplex_stage_stop_from_status(result->status),
		code == 0 ? SIMPLEX_PAN_OUTPUT_INPUT_FEASIBLE_POINT
			  : SIMPLEX_PAN_OUTPUT_NONE,
		NULL);
	/* The public-input-space output must not expose the stack metadata view. */
	f.outcome.representation = NULL;
	f.outcome.input = input;
	if (out)
		*out = f.outcome;
	result->iterations = f.outcome.end - f.outcome.before;
	if (!ready) {
		result->status = lp_simplex_PrecisionError;
		return lp_simplex_EXIT_FAILURE;
	}
	return code;
}

static int service_bind_bounded(
	struct simplex_PanServiceFrame *f, const struct simplex_PanBounded *b)
{
	f->bounded = b;
	f->bounded_input = b->problem;
	f->lower = b->lower;
	f->upper = b->upper;
	if (f->context && f->context->input && f->context->input != b->problem)
		f->outcome.valid = 0;
	return f->outcome.valid;
}

/* Opportunities describe the native loop's policy. The causal source for
 * a multiplier opportunity is its typed, completed producer, not whichever
 * diagnostic event happens to precede this call. */
static struct simplex_ExecutionRef service_observe_face(
	const struct simplex_PanStandard *a, int work,
	const struct simplex_PanFaceReport *face,
	const struct simplex_PanServiceOutcome *search)
{
	struct simplex_ExecutionRef actor = {0};
	struct simplex_PanServiceContext *c = a->services;
	struct simplex_ExecutionScope scope;
	int id, cause = 0;
	if (!c || !c->scope.history || !simplex_pan_service_context_valid(c))
		return actor;
	scope = c->scope;
	scope.model = a->service_model;
	scope.parent = c->parent;
	scope.parent_actor = c->parent_actor;
	if (face->opportunity == SIMPLEX_PAN_FACE_AFTER_MULTIPLIER_SEARCH &&
		search && search->valid && search->input == c->input &&
		search->input_model == c->scope.model && search->end <= work &&
		(search->service == SIMPLEX_PAN_SERVICE_DUAL_REFINE ||
			search->service == SIMPLEX_PAN_SERVICE_PREDICT) &&
		search->actor.kind == SIMPLEX_EXECUTION_PAN_SERVICE &&
		search->actor.detail == search->service &&
		search->actor.model == search->working_model)
		cause = simplex_execution_ref_id(&search->actor, &scope, 1);
	id = simplex_execution_begin(&scope,
		SIMPLEX_EXECUTION_PAN_FACE_OPPORTUNITY, face->opportunity,
		service_face_reason(face->opportunity), cause);
	if (id) {
		struct simplex_ExecutionEvent *e =
			scope.history->events + id - 1;
		e->source_model = c->scope.model;
		e->observation_work = work;
		e->face = *face;
	}
	simplex_execution_end(scope.history, id, NULL, NULL);
	return simplex_execution_ref(&scope, id);
}

int simplex_pan_service_face(struct simplex_PanFaceRecovery *owner,
	const struct simplex_PanBounded *b, double *point,
	const double *multiplier, const struct lp_simplex_Options *o,
	int *iterations, double *saved, double *residual, double *work,
	int opportunity, int steps, int low_progress,
	const struct simplex_PanServiceOutcome *search,
	struct simplex_PanServiceOutcome *out)
{
	struct simplex_PanServiceFrame f;
	struct simplex_PanFaceReport face;
	struct simplex_ExecutionRef cause;
	int code = -1, entered, stop = SIMPLEX_STOP_UNKNOWN;
	memset(&face, 0, sizeof(face));
	face.opportunity = opportunity;
	face.steps = steps;
	face.low_progress = low_progress;
	cause = service_observe_face(&b->a, *iterations, &face, search);
	entered = simplex_pan_service_begin(&f, &b->a,
		SIMPLEX_PAN_SERVICE_FACE_RECOVERY, *iterations,
		o->iteration_limit, 0);
	f.trigger_actor = cause;
	f.outcome.face = face;
	if (entered && service_bind_bounded(&f, b)) {
		code = simplex_pan_face_recover_reported(owner, b, point,
			multiplier, o, saved, residual, work,
			&f.outcome.face.result);
		stop = code < 0	  ? SIMPLEX_STOP_RESOURCE
		       : code > 0 ? SIMPLEX_STOP_HANDOFF
		       : f.outcome.face.result.action ==
				       SIMPLEX_PAN_FACE_EXACT_SKIP
			       ? SIMPLEX_STOP_NOT_APPLICABLE
			       : SIMPLEX_STOP_REASSESS;
	}
	if (!simplex_pan_service_finish(&f, iterations, code,
		    lp_simplex_CondUnsatisfied, stop,
		    code > 0 ? SIMPLEX_PAN_OUTPUT_NATIVE_FACE_ADJUSTMENT
			     : SIMPLEX_PAN_OUTPUT_NONE,
		    out))
		return -1;
	return code;
}

int simplex_pan_service_pdhg(const struct simplex_PanBounded *b, double *point,
	double *multiplier, const struct lp_simplex_Options *o, int *iterations,
	int budget, const struct simplex_PanConfig *cfg,
	struct simplex_PanServiceOutcome *out)
{
	struct simplex_PanServiceFrame f;
	struct lp_simplex_Options lease = *o;
	int code = 0, stop = SIMPLEX_STOP_UNKNOWN,
	    output = SIMPLEX_PAN_OUTPUT_NONE;
	if (simplex_pan_service_begin(&f, &b->a, SIMPLEX_PAN_SERVICE_PDHG,
		    *iterations, o->iteration_limit, budget) &&
		service_bind_bounded(&f, b)) {
		lease.iteration_limit = f.outcome.before + f.outcome.grant;
		code = simplex_pan_bounded_pdhg_configured_reported(b, point,
			multiplier, &lease, iterations, f.outcome.grant, cfg,
			&stop, &output);
	}
	if (!simplex_pan_service_finish(&f, iterations, code,
		    code == 1 ? lp_simplex_Success : lp_simplex_CondUnsatisfied,
		    stop, output, out))
		return 0;
	return code;
}

/* A completed producer observation is a cause, not another numerical action.
 * Its journal may be absent or full without changing the repair decision. */
static struct simplex_ExecutionRef service_observe_trigger(
	const struct simplex_PanStandard *a,
	const struct simplex_PanServiceTrigger *trigger)
{
	struct simplex_ExecutionRef actor = {0};
	struct simplex_PanServiceContext *c = a->services;
	struct simplex_ExecutionScope scope;
	int id, reason = trigger ? service_trigger_reason(trigger->kind)
				 : SIMPLEX_REASON_NONE;

	if (!trigger || reason == SIMPLEX_REASON_NONE || !c ||
		!simplex_pan_service_context_valid(c))
		return actor;
	scope = c->scope;
	scope.model = a->service_model;
	scope.parent = c->parent;
	scope.parent_actor = c->parent_actor;
	id = simplex_execution_begin(&scope, SIMPLEX_EXECUTION_PAN_TRIGGER,
		trigger->kind, reason, 0);
	if (id) {
		struct simplex_ExecutionEvent *e =
			scope.history->events + id - 1;
		e->source_model = c->scope.model;
		e->observation_work = trigger->work;
		e->primal.value = trigger->primal;
		e->primal.quality =
			isfinite(trigger->primal) && trigger->primal >= 0.
				? SIMPLEX_EVIDENCE_RECOMPUTED
				: SIMPLEX_EVIDENCE_UNKNOWN;
		e->primal.space = SIMPLEX_SPACE_WORKING;
		e->primal.kind = SIMPLEX_OBSERVATION_PRIMAL_MAX;
		e->primal.work = trigger->work;
		e->primal.epoch = trigger->epoch;
	}
	simplex_execution_end(scope.history, id, NULL, NULL);
	return simplex_execution_ref(&scope, id);
}

int simplex_pan_service_restore_standard(const struct simplex_PanStandard *a,
	double *point, const struct lp_simplex_Options *o, long double ceiling,
	int *iterations, const struct simplex_PanConfig *cfg,
	const struct simplex_PanServiceTrigger *trigger,
	struct simplex_PanServiceOutcome *out)
{
	struct simplex_PanServiceFrame f;
	struct lp_simplex_Options lease = *o;
	int code = -1, stop = SIMPLEX_STOP_UNKNOWN;
	struct simplex_ExecutionRef cause = service_observe_trigger(a, trigger);
	if (simplex_pan_service_begin(&f, a, SIMPLEX_PAN_SERVICE_RESTORE,
		    *iterations, o->iteration_limit,
		    service_remaining(o->iteration_limit, *iterations))) {
		lease.iteration_limit = f.outcome.before + f.outcome.grant;
		code = simplex_pan_restore_configured_reported(
			a, point, &lease, ceiling, iterations, cfg, &stop);
	}
	f.trigger_actor = cause;
	if (trigger && trigger->kind == SIMPLEX_PAN_TRIGGER_PRIMAL_MARGIN)
		f.outcome.trigger = *trigger;
	if (!simplex_pan_service_finish(&f, iterations, code,
		    lp_simplex_CondUnsatisfied, stop,
		    code == 1 ? SIMPLEX_PAN_OUTPUT_NONNEGATIVE_EQUALITY_POINT
			      : SIMPLEX_PAN_OUTPUT_NONE,
		    out))
		return -1;
	return code;
}

int simplex_pan_service_restore(const struct simplex_PanBounded *b,
	double *point, const struct lp_simplex_Options *o, int *iterations,
	long double ceiling, const struct simplex_PanConfig *cfg,
	const struct simplex_PanServiceTrigger *trigger,
	struct simplex_PanServiceOutcome *out)
{
	struct simplex_PanServiceFrame f;
	struct lp_simplex_Options lease = *o;
	int code = -1, stop = SIMPLEX_STOP_UNKNOWN, entered;
	struct simplex_ExecutionRef cause =
		service_observe_trigger(&b->a, trigger);
	entered = simplex_pan_service_begin(&f, &b->a,
		SIMPLEX_PAN_SERVICE_RESTORE, *iterations, o->iteration_limit,
		service_remaining(o->iteration_limit, *iterations));
	f.trigger_actor = cause;
	if (trigger &&
		(trigger->kind == SIMPLEX_PAN_TRIGGER_PDHG_STAGNATION ||
			trigger->kind ==
				SIMPLEX_PAN_TRIGGER_PDHG_PRIMAL_STAGNATION))
		f.outcome.trigger = *trigger;
	if (entered && service_bind_bounded(&f, b)) {
		lease.iteration_limit = f.outcome.before + f.outcome.grant;
		code = simplex_pan_bounded_restore_configured_reported(b, point,
			&lease, iterations, ceiling, cfg,
			trigger &&
				trigger->kind ==
					SIMPLEX_PAN_TRIGGER_PDHG_PRIMAL_STAGNATION,
			&stop);
	}
	if (!simplex_pan_service_finish(&f, iterations, code,
		    lp_simplex_CondUnsatisfied, stop,
		    code == 1 ? SIMPLEX_PAN_OUTPUT_NATIVE_FEASIBLE_POINT
			      : SIMPLEX_PAN_OUTPUT_NONE,
		    out))
		return -1;
	return code;
}

int simplex_pan_service_first_order(const struct simplex_PanStandard *a,
	const double *point, double *multiplier,
	const struct lp_simplex_Options *o, int *iterations, double *violation,
	double *direction, const unsigned char *free_mask, int warm,
	const struct simplex_PanConfig *cfg,
	struct simplex_PanServiceOutcome *out)
{
	struct simplex_PanServiceFrame f;
	struct lp_simplex_Options lease = *o;
	int code = 0, stop = SIMPLEX_STOP_UNKNOWN;
	int budget = direction ? cfg->first_order_iterations
			       : cfg->dual_polish_iterations;
	if (simplex_pan_service_begin(&f, a,
		    direction ? SIMPLEX_PAN_SERVICE_PREDICT
			      : SIMPLEX_PAN_SERVICE_DUAL_REFINE,
		    *iterations, o->iteration_limit, budget)) {
		lease.iteration_limit = f.outcome.before + f.outcome.grant;
		code = simplex_pan_cone_first_order_configured_reported(a,
			point, multiplier, &lease, iterations, violation,
			direction, free_mask, warm, cfg, &stop);
	}
	if (!simplex_pan_service_finish(&f, iterations, code,
		    lp_simplex_CondUnsatisfied, stop,
		    code == 1	? SIMPLEX_PAN_OUTPUT_CHECKED_DUAL
		    : code == 2 ? SIMPLEX_PAN_OUTPUT_CHECKED_DIRECTION
				: SIMPLEX_PAN_OUTPUT_NONE,
		    out))
		return 0;
	return code;
}

int simplex_pan_service_predict(const struct simplex_PanStandard *a,
	const double *point, double *multiplier,
	const struct lp_simplex_Options *o, int *iterations, double *violation,
	double *direction, const unsigned char *free_mask, int budget,
	int (*accept)(void *, const double *), void *accept_context,
	const struct simplex_PanConfig *cfg,
	struct simplex_PanServiceOutcome *out)
{
	struct simplex_PanServiceFrame f;
	struct lp_simplex_Options lease = *o;
	int code = 0, stop = SIMPLEX_STOP_UNKNOWN;
	if (simplex_pan_service_begin(&f, a, SIMPLEX_PAN_SERVICE_PREDICT,
		    *iterations, o->iteration_limit, budget)) {
		lease.iteration_limit = f.outcome.before + f.outcome.grant;
		code = simplex_pan_cone_predict_audited_configured_reported(a,
			point, multiplier, &lease, iterations, violation,
			direction, free_mask, budget, accept, accept_context,
			cfg, &stop);
	}
	if (!simplex_pan_service_finish(&f, iterations, code,
		    lp_simplex_CondUnsatisfied, stop,
		    code == 1	? SIMPLEX_PAN_OUTPUT_CHECKED_DUAL
		    : code == 2 ? SIMPLEX_PAN_OUTPUT_CHECKED_DIRECTION
				: SIMPLEX_PAN_OUTPUT_NONE,
		    out))
		return 0;
	return code;
}

int simplex_pan_service_cone(struct simplex_PanBasis *basis, int implicit,
	double *direction, const unsigned char *free_mask, double *residual,
	double tolerance, int limit, int *iterations, int *status,
	int (*accept)(void *, const double *), void *accept_context,
	struct simplex_PanServiceOutcome *out)
{
	struct simplex_PanServiceFrame f;
	int code = -1, output = SIMPLEX_PAN_OUTPUT_NONE;
	if (simplex_pan_service_begin(&f, basis->standard,
		    SIMPLEX_PAN_SERVICE_CONE, *iterations, limit,
		    service_remaining(limit, *iterations)))
		code = implicit ? simplex_pan_cone_implicit_audited(basis,
					  direction, free_mask, residual,
					  tolerance,
					  f.outcome.before + f.outcome.grant,
					  iterations, status, accept,
					  accept_context)
				: simplex_pan_phase1_mixed_audited(basis,
					  direction, free_mask, residual,
					  tolerance, tolerance,
					  f.outcome.before + f.outcome.grant,
					  iterations, status, accept,
					  accept_context);
	else
		*status = lp_simplex_PrecisionError;
	f.basis = basis;
	f.outcome.implementation = implicit ? SIMPLEX_PAN_IMPLICIT_ORACLE
					    : SIMPLEX_PAN_EXPLICIT_ORACLE;
	if (code == 0)
		output = SIMPLEX_PAN_OUTPUT_DIRECTION_PROPOSAL;
	if (code == 1)
		output = SIMPLEX_PAN_OUTPUT_STATIONARY_RESIDUAL;
	if (code == 2)
		output = SIMPLEX_PAN_OUTPUT_CHECKED_DIRECTION;
	if (!simplex_pan_service_finish(&f, iterations, code, *status,
		    code >= 0 ? SIMPLEX_STOP_HANDOFF
			      : simplex_stage_stop_from_status(*status),
		    output, out)) {
		*status = lp_simplex_PrecisionError;
		return -1;
	}
	return code;
}
