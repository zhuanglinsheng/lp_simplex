/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_progress.c
 * @brief Audited input-model progress and producer provenance.
 */

#include "simplex_stage.h"
#include "simplex_execution.h"
#include "simplex_progress.h"
#include "policy/defaults.h"
#include "utils.h"
#include <math.h>
#include <float.h>
#include <string.h>

/* Outward long-double intervals also work where long double == double.
 * Exact zero operations are kept exact, notably for unbounded/free columns. */
static long double down_add(long double a, long double b)
{
	long double sum, v, error;
	if (b == 0.)
		return a;
	if (a == 0.)
		return b;
	sum = a + b;
	if (!isfinite(sum))
		return -HUGE_VALL;
	v = sum - a;
	error = (a - (sum - v)) + (b - v);
	return error >= 0. ? sum : nextafterl(sum, -HUGE_VALL);
}

static long double down_product(long double a, long double b)
{
	long double product;
	if (a == 0. || b == 0.)
		return 0.;
	product = a * b;
	if (!isfinite(product))
		return -HUGE_VALL;
	/* Close to underflow, even the FMA residual can round to zero. */
	if (fabsl(product) >= 4.L * (LDBL_MIN / LDBL_EPSILON) &&
		fmal(a, b, -product) >= 0.)
		return product;
	return nextafterl(product, -HUGE_VALL);
}

static long double multiplier(
	const struct simplex_Problem *p, const double *y, double scale, int i)
{
	long double v = y ? (long double)y[i] * scale : 0.;
	if (p->row_type[i] == optm_CONS_T_LE && v < 0.)
		v = 0.;
	if (p->row_type[i] == optm_CONS_T_GE && v > 0.)
		v = 0.;
	return v;
}

double simplex_progress_dual_bound(
	const struct simplex_Problem *p, const double *y, double scale)
{
	long double total = 0.;
	int i, j, k;
	if (!isfinite(scale))
		return -HUGE_VAL;
	for (i = 0; i < p->rows; i++) {
		long double v = multiplier(p, y, scale, i);
		if (!isfinite(v))
			return -HUGE_VAL;
		total = down_add(total, down_product(-p->rhs[i], v));
	}
	for (j = 0; j < p->columns; j++) {
		const struct optm_VariableBound *b = p->bounds + j;
		long double lo = p->objective[j], hi = lo, term;
		double lower = (b->b_type == optm_BOUND_T_LO ||
				       b->b_type == optm_BOUND_T_BS)
				       ? b->lb
				       : -HUGE_VAL;
		double upper = (b->b_type == optm_BOUND_T_UP ||
				       b->b_type == optm_BOUND_T_BS)
				       ? b->ub
				       : HUGE_VAL;
		for (k = p->matrix.column_start[j];
			k < p->matrix.column_start[j + 1]; k++) {
			long double v =
				multiplier(p, y, scale, p->matrix.row_index[k]);
			long double a = p->matrix.value[k];
			lo = down_add(lo, down_product(a, v));
			hi = -down_add(-hi, down_product(-a, v));
		}
		if (!isfinite(lo) || !isfinite(hi))
			return -HUGE_VAL;
		/* Minimise g*x over both the box and the enclosing g interval.
		 * Crossing zero with an infinite endpoint cannot yield a bound. */
		if ((!isfinite(lower) && hi > 0.) ||
			(!isfinite(upper) && lo < 0.))
			return -HUGE_VAL;
		term = 0.;
		if (isfinite(lower)) {
			term = down_product(lower, lower >= 0. ? lo : hi);
			if (!isfinite(upper) && lo <= 0. && hi >= 0.)
				term = fminl(term, 0.);
		}
		if (isfinite(upper)) {
			long double other =
				down_product(upper, upper >= 0. ? lo : hi);
			term = isfinite(lower) ? fminl(term, other) : other;
			if (!isfinite(lower) && lo <= 0. && hi >= 0.)
				term = fminl(term, 0.);
		}
		total = down_add(total, term);
	}
	if (!isfinite(total) || !isfinite((double)total))
		return -HUGE_VAL;
	{
		double rounded = (double)total;
		return (long double)rounded > total
			       ? nextafter(rounded, -HUGE_VAL)
			       : rounded;
	}
}

static int tracker_due(void *context)
{
	struct simplex_ProgressTracker *t = context;
	/* Observation is deterministic and independent of UI refresh time. */
	t->polls++;
	return t->polls % SIMPLEX_PROGRESS_OBSERVATION_PERIOD == 0;
}

static int tracker_archive(struct simplex_ProgressTracker *t,
	const struct simplex_ProgressSource *source, const char *stage,
	int work, double primal)
{
	struct simplex_StageReport audit;
	struct simplex_ExecutionEvent *event;
	int id, producer = 0, model = 0;
	simplex_execution_set_parent(&t->execution,
		simplex_execution_ref_id(&t->parent_actor, &t->execution, 0));
	t->point_serial = 0;
	if (!source || !t->execution.context ||
		source->coordinates.context != t->execution.context ||
		source->coordinates.model != t->execution.model ||
		source->work != work)
		return 0;
	if (source->producer.context == t->execution.context) {
		model = source->producer.model;
		if (source->producer.history == t->execution.history)
			producer = simplex_execution_parent(&source->producer);
		if (producer &&
			t->execution.history->events[producer - 1].model !=
				model)
			producer = 0;
	}
	id = simplex_execution_begin(&t->execution,
		SIMPLEX_EXECUTION_POINT_ARCHIVE, (int)source->transforms,
		SIMPLEX_REASON_FEASIBLE_CANDIDATE, producer);
	if (!id)
		return 0;
	event = t->execution.history->events + id - 1;
	event->point_source = producer;
	event->source_model = model;
	event->observation_work = work;
	event->model_rows = t->original.rows;
	event->model_columns = t->original.columns;
	strncpy(event->observation_stage, stage ? stage : "unknown",
		sizeof(event->observation_stage) - 1);
	event->observation_stage[sizeof(event->observation_stage) - 1] = 0;
	simplex_stage_reset(&audit);
	audit.stop = SIMPLEX_STOP_HANDOFF;
	audit.point_available = 1;
	audit.primal.value = primal;
	audit.primal.kind = SIMPLEX_OBSERVATION_PRIMAL_MAX;
	audit.primal.space = SIMPLEX_SPACE_INPUT;
	audit.primal.quality = SIMPLEX_EVIDENCE_RECOMPUTED;
	audit.primal.work = work;
	simplex_execution_end(t->execution.history, id, NULL, &audit);
	t->point_serial = event->serial;
	return id;
}

static void tracker_report(void *context, const double *x, const double *y,
	double scale, const char *stage, int iterations, int force,
	const struct simplex_ProgressSource *source)
{
	struct simplex_ProgressTracker *t = context;
	struct lp_simplex_Progress event;
	long double objective = t->model->objective_offset;
	double primal = HUGE_VAL, bound, candidate;
	int j, improved = 0;
	if (source &&
		(source->rows != t->original.rows ||
			source->columns != t->original.columns ||
			(t->execution.context && source->coordinates.context &&
				(source->coordinates.context !=
						t->execution.context ||
					source->coordinates.model !=
						t->execution.model))))
		return;

	if (x) {
		primal = simplex_problem_primal_infeasibility(&t->original, x);
		if (isfinite(primal) &&
			primal <= t->options->primal_tolerance) {
			for (j = 0; j < t->original.columns; j++)
				objective = -down_add(-objective,
					down_product(-t->original.objective[j],
						x[j]));
			candidate = (double)objective;
			if ((long double)candidate < objective)
				candidate = nextafter(candidate, HUGE_VAL);
			if (isfinite(objective) && isfinite(candidate) &&
				(!t->have_primal ||
					candidate < t->best_primal)) {
				t->best_primal = candidate;
				t->primal_error = primal;
				t->have_primal = improved = 1;
				memcpy(t->point, x,
					(size_t)t->model->n * sizeof(double));
				t->point_source = tracker_archive(
					t, source, stage, iterations, primal);
			}
		}
	}
	bound = simplex_progress_dual_bound(&t->original, y, scale);
	if (isfinite(bound)) {
		long double v = down_add(bound, t->model->objective_offset);
		bound = (double)v;
		if ((long double)bound > v)
			bound = nextafter(bound, -HUGE_VAL);
		if (isfinite(bound) && bound > t->best_dual)
			t->best_dual = bound;
	}
	event.stage = stage;
	event.stage_iterations = iterations;
	event.cpu_seconds = (double)(clock() - t->begun) / CLOCKS_PER_SEC;
	if (t->model->objective_sense == 1.) {
		event.lower_bound = t->best_dual;
		event.upper_bound = t->best_primal;
	} else {
		event.lower_bound = -t->best_primal;
		event.upper_bound = -t->best_dual;
	}
	event.absolute_gap = event.relative_gap = HUGE_VAL;
	/* A crossed interval is unavailable rather than a false zero gap.
	 * Finite-tolerance primal feasibility can cause a tiny crossing. */
	if (isfinite(event.lower_bound) && isfinite(event.upper_bound) &&
		event.lower_bound <= event.upper_bound) {
		event.absolute_gap = event.upper_bound - event.lower_bound;
		event.relative_gap = event.absolute_gap /
				     fmax(1., fmax(fabs(event.lower_bound),
						      fabs(event.upper_bound)));
	}
	event.point = improved ? t->point : NULL;
	event.primal_infeasibility = t->primal_error;
	t->point_pending |= improved;
	if (t->options->progress_callback &&
		(force || (double)(clock() - t->last) / CLOCKS_PER_SEC >=
				  t->options->progress_interval)) {
		event.point = t->point_pending ? t->point : NULL;
		t->options->progress_callback(
			t->options->progress_context, &event);
		t->point_pending = 0;
		t->last = clock();
	}
}

int simplex_progress_init(struct simplex_ProgressTracker *t,
	const struct lp_Model *model, const struct lp_simplex_Options *o)
{
	return simplex_progress_init_reported(t, model, o, NULL);
}

int simplex_progress_init_reported(struct simplex_ProgressTracker *t,
	const struct lp_Model *model, const struct lp_simplex_Options *o,
	const struct simplex_ExecutionScope *execution)
{
	memset(t, 0, sizeof(*t));
	if (execution)
		t->execution = *execution;
	t->parent_actor = t->execution.parent_actor;

	t->model = model;
	t->options = o;
	t->best_primal = t->primal_error = HUGE_VAL;
	t->best_dual = -HUGE_VAL;
	t->begun = t->last = clock();
	if (simplex_problem_from_model(&t->original, model))
		return -1;
	t->storage = lp_simplex_malloc((size_t)model->n * sizeof(double));
	if (!t->storage) {
		simplex_problem_destroy(&t->original);
		return -1;
	}
	t->point = t->storage;
	t->sink.context = t;
	t->sink.due = tracker_due;
	t->sink.report = tracker_report;
	tracker_report(t, NULL, NULL, 1., "start", 0, 1, NULL);
	return 0;
}

void simplex_progress_finish(
	struct simplex_ProgressTracker *t, const double *x, int work)
{
	struct simplex_ProgressSource source;
	memset(&source, 0, sizeof(source));
	simplex_execution_set_parent(&t->execution,
		simplex_execution_ref_id(&t->parent_actor, &t->execution, 0));
	source.producer = source.coordinates = t->execution;
	source.rows = t->original.rows;
	source.columns = t->original.columns;
	source.work = work;
	tracker_report(t, x, NULL, 1., "finish", work, 1, &source);
}

int simplex_progress_point_source(const struct simplex_ProgressTracker *t,
	const struct simplex_ExecutionScope *execution)
{
	const struct simplex_ExecutionEvent *event;
	int id;
	if (!t || !execution || !t->have_primal || !t->point_serial ||
		!t->execution.context ||
		t->execution.context != execution->context ||
		t->execution.model != execution->model ||
		t->execution.history != execution->history ||
		!execution->history || execution->history->count < 0 ||
		execution->history->count > SIMPLEX_EXECUTION_CAPACITY)
		return 0;
	id = t->point_source;
	if (id <= 0 || id > execution->history->count)
		return 0;
	event = execution->history->events + id - 1;
	return event->context == execution->context &&
			       event->serial == t->point_serial &&
			       event->complete &&
			       event->kind == SIMPLEX_EXECUTION_POINT_ARCHIVE &&
			       event->model == execution->model &&
			       event->stop == SIMPLEX_STOP_HANDOFF &&
			       event->primal.quality ==
				       SIMPLEX_EVIDENCE_RECOMPUTED
		       ? id
		       : 0;
}

int simplex_progress_record_restore(const struct simplex_ProgressTracker *t,
	const struct simplex_ExecutionScope *execution,
	const struct lp_simplex_Result *result)
{
	int source = simplex_progress_point_source(t, execution);
	int id = simplex_execution_begin(execution,
		SIMPLEX_EXECUTION_POINT_SELECTION, 0,
		SIMPLEX_REASON_RESTORE_ARCHIVED_POINT, source);
	if (id)
		execution->history->events[id - 1].point_source = source;
	simplex_execution_end(execution->history, id, result, NULL);
	return id;
}

void simplex_progress_destroy(struct simplex_ProgressTracker *t)
{
	if (!t->sink.report)
		return;
	simplex_problem_destroy(&t->original);
	lp_simplex_free(t->storage);
}

static int bridge_due(void *context)
{
	struct simplex_ProgressBridge *b = context;
	return b->parent->due(b->parent->context);
}

static void bridge_report(void *context, const double *x, const double *y,
	double scale, const char *stage, int iterations, int force,
	const struct simplex_ProgressSource *source)
{
	struct simplex_ProgressBridge *b = context;
	struct simplex_ProgressSource mapped;
	struct simplex_ProblemIdentity incoming;
	int i;
	/* Check before reading vectors at the bound representation's dimensions.
	 * An absent legacy source remains unknown; a known wrong model is rejected. */
	if (!simplex_problem_identity_unchanged(&b->input, b->input.input) ||
		(source &&
			(source->rows != b->input.rows ||
				source->columns != b->input.columns ||
				source->coordinates.context !=
					b->input.context ||
				source->coordinates.model != b->input.model)))
		return;
	if (source) {
		incoming = source->input;
		/* Read-only views may change the owner/parent but must retain every
		 * storage pointer, model identity, dimension and revision. This is
		 * observation compatibility, not checkpoint reuse permission. */
		incoming.input = b->input.input;
		if (!simplex_problem_identity_equal(&incoming, &b->input))
			return;
	}
	if (b->output.input && !simplex_problem_identity_unchanged(
				       &b->output, b->output.input))
		return;
	if (b->presolve) {
		if (x)
			simplex_presolve_postsolve(b->presolve, x, b->x);
		memset(b->y, 0, (size_t)b->rows * sizeof(double));
		if (y)
			for (i = 0; i < b->presolve->reduced->rows; i++)
				b->y[b->presolve->row_map[i]] = y[i];
	} else {
		if (x)
			for (i = 0; i < b->columns; i++)
				b->x[i] = x[i] * b->scaling->column[i];
		if (y)
			for (i = 0; i < b->rows; i++)
				b->y[i] = y[i] * b->scaling->row[i];
	}
	if (source) {
		mapped = *source;
		mapped.coordinates = b->target;
		mapped.input = b->output;
		mapped.rows = b->rows;
		mapped.columns = b->columns;
		mapped.transforms |= b->presolve ? SIMPLEX_PROGRESS_PRESOLVE
						 : SIMPLEX_PROGRESS_SCALING;
	}
	b->parent->report(b->parent->context, x ? b->x : NULL, y ? b->y : NULL,
		scale, stage, iterations, force, source ? &mapped : NULL);
}

int simplex_progress_bridge(struct simplex_ProgressBridge *b,
	const struct simplex_ProgressSink *parent,
	const struct simplex_Presolve *presolve,
	const struct simplex_Scaling *scaling, int rows, int columns,
	const struct simplex_Problem *input,
	const struct simplex_ExecutionScope *target,
	const struct simplex_Problem *output)
{
	memset(b, 0, sizeof(*b));
	if (!parent)
		return 0;
	if (!input || !target || (presolve != NULL) == (scaling != NULL) ||
		input->execution.context != target->context || rows < 0 ||
		columns < 0 ||
		(output &&
			(output->rows != rows || output->columns != columns ||
				output->execution.context != target->context ||
				output->execution.model != target->model)) ||
		(scaling && !output) ||
		(presolve &&
			(presolve->reduced != input || !presolve->original ||
				presolve->original->m != rows ||
				presolve->original->n != columns)) ||
		(scaling &&
			(scaling->rows != rows || scaling->columns != columns ||
				input->rows != rows ||
				input->columns != columns)))
		return SIMPLEX_PROGRESS_BRIDGE_INVALID;
	b->input = simplex_problem_identity(input);
	b->output = simplex_problem_identity(output);
	b->target = *target;
	b->storage =
		lp_simplex_malloc(((size_t)rows + columns) * sizeof(double));
	if (!b->storage)
		return SIMPLEX_PROGRESS_BRIDGE_MEMORY;
	b->parent = parent;
	b->presolve = presolve;
	b->scaling = scaling;
	b->rows = rows;
	b->columns = columns;
	b->x = b->storage;
	b->y = b->x + columns;
	b->sink.context = b;
	b->sink.due = bridge_due;
	b->sink.report = bridge_report;
	return 0;
}

void simplex_progress_bridge_destroy(struct simplex_ProgressBridge *b)
{
	lp_simplex_free(b->storage);
}

void simplex_progress_emit(const struct simplex_Problem *p, const double *x,
	const double *y, double scale, const char *stage, int iterations,
	int force)
{
	const struct simplex_ProgressSink *s = p->progress;
	if (s && (force || s->due(s->context))) {
		struct simplex_ProgressSource source;
		source.producer = source.coordinates = p->execution;
		source.input = simplex_problem_identity(p);
		source.rows = p->rows;
		source.columns = p->columns;
		source.work = iterations;
		source.transforms = 0;
		s->report(s->context, x, y, scale, stage, iterations, force,
			&source);
	}
}

void simplex_progress_observe(const struct simplex_Problem *p, const double *x,
	const double *y, double scale, const char *stage, int work)
{
	const struct simplex_ProgressSink *s = p->progress;
	if (s) {
		struct simplex_ProgressSource source;
		source.producer = source.coordinates = p->execution;
		source.input = simplex_problem_identity(p);
		source.rows = p->rows;
		source.columns = p->columns;
		source.work = work;
		source.transforms = 0;
		s->report(s->context, x, y, scale, stage, work, 0, &source);
	}
}
