/**
 * @file test_progress_origin.c
 * @brief Provenance crosses reversible observation bridges without becoming
 * state.
 */

#include "simplex_stage.h"
#include "simplex_solve_context.h"
#include "simplex_execution.h"
#include "simplex_progress.h"
#include "simplex_dual.h"
#include <lp_simplex/lp_simplex.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#define CHECK(c) \
	do { \
		if (!(c)) { \
			fprintf(stderr, "progress origin line %d: %s\n", \
				__LINE__, #c); \
			return 1; \
		} \
	} while (0)

static struct lp_Model *fixture(void)
{
	struct lp_Model *m = lp_model_create(3, 3);
	int i;
	if (!m)
		return NULL;
	for (i = 0; i < 3; i++) {
		m->objective[i] = 1.;
		m->bounds[i].b_type = optm_BOUND_T_BS;
		m->bounds[i].ub = 1.;
		m->constraints[i].type = optm_CONS_T_GE;
		m->constraints[i].rhs = i == 0 ? 1048576. : 1.;
	}
	m->constraints[0].type = optm_CONS_T_EQ;
	m->constraints[0].coef[0] = m->constraints[0].coef[1] = 1048576.;
	m->constraints[1].coef[1] = m->constraints[1].coef[2] = 1.;
	m->constraints[2].coef[0] = m->constraints[2].coef[2] = 1.;
	return m;
}

static int actual_chain(int unscaled, int tracing)
{
	struct lp_Model *m = fixture();
	struct simplex_SolveContext context;
	struct simplex_Execution history;
	struct simplex_ExecutionScope root, wrong;
	struct simplex_ProgressTracker tracker;
	struct simplex_ProgressBridge bridge;
	struct simplex_Presolve pre;
	struct lp_simplex_Options o;
	struct lp_simplex_Result result;
	double x[3], worse[3] = {0., 1., 1.};
	int archive, selection, i;
	CHECK(m);
	simplex_solve_context_init(&context, tracing ? &history : NULL, &root);
	simplex_solve_context_configure(&context, 0);
	context.controls.disable_scaling = unscaled;
	simplex_execution_set_parent(
		&root, simplex_execution_begin(
			       &root, SIMPLEX_EXECUTION_SOLVE, 0, 0, 0));
	if (tracing == 2)
		for (i = history.count; i < SIMPLEX_EXECUTION_CAPACITY; i++)
			(void)simplex_execution_begin(
				&root, SIMPLEX_EXECUTION_SOLVE, 0, 0, 0);
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_DUAL_REVISED);
	CHECK(!simplex_progress_init_reported(&tracker, m, &o, &root));
	CHECK(!simplex_presolve_run(&pre, m, o.primal_tolerance));
	CHECK(pre.reduced && pre.reduced->columns < m->n);
	pre.reduced->execution =
		simplex_solve_context_derive(&root, SIMPLEX_EXECUTION_PRESOLVE);
	CHECK(!simplex_progress_bridge(&bridge, &tracker.sink, &pre, NULL, m->m,
		m->n, pre.reduced, &root, NULL));
	pre.reduced->progress = &bridge.sink;
	CHECK(!simplex_dual_solve_problem(
		pre.reduced, &o, x, NULL, &result, 0));
	CHECK(tracker.have_primal && fabs(tracker.best_primal - 1.5) < 1e-8);
	archive = simplex_progress_point_source(&tracker, &root);
	if (tracing == 1) {
		const struct simplex_ExecutionEvent *e;
		CHECK(archive > 0);
		e = history.events + archive - 1;
		CHECK(e->kind == SIMPLEX_EXECUTION_POINT_ARCHIVE &&
			e->model == root.model);
		CHECK(!e->result_available && e->complete &&
			e->stop == SIMPLEX_STOP_HANDOFF);
		CHECK(e->point_source > 0 && e->point_source < archive &&
			e->cause == e->point_source);
		CHECK(history.events[e->point_source - 1].kind ==
			SIMPLEX_EXECUTION_BACKEND);
		CHECK(e->source_model ==
			history.events[e->point_source - 1].model);
		CHECK(e->model_rows == m->m && e->model_columns == m->n);
		CHECK(e->observation_work <= result.iterations &&
			!strcmp(e->observation_stage, "dual"));
		CHECK(e->detail ==
			(SIMPLEX_PROGRESS_PRESOLVE |
				(unscaled ? 0 : SIMPLEX_PROGRESS_SCALING)));
		CHECK(e->primal.quality == SIMPLEX_EVIDENCE_RECOMPUTED &&
			e->primal.value <= o.primal_tolerance);
	} else
		CHECK(!archive);
	/* Destroy all transformed model buffers: provenance never borrows those. */
	simplex_progress_bridge_destroy(&bridge);
	simplex_presolve_destroy(&pre);
	simplex_progress_finish(&tracker, worse, 50);
	CHECK(simplex_progress_point_source(&tracker, &root) == archive);
	CHECK(fabs(tracker.best_primal - 1.5) < 1e-8);
	wrong = root;
	wrong.context = NULL;
	CHECK(!simplex_progress_point_source(&tracker, &wrong));
	memcpy(worse, tracker.point, sizeof(worse));
	CHECK(simplex_problem_primal_infeasibility(&tracker.original, worse) <=
		o.primal_tolerance);
	result.status = lp_simplex_PrecisionError;
	result.iterations = 50;
	result.objective = tracker.best_primal;
	result.primal_infeasibility = tracker.primal_error;
	result.dual_infeasibility = NAN;
	result.dual_residual_space = lp_simplex_RESIDUAL_UNAVAILABLE;
	selection = simplex_progress_record_restore(&tracker, &root, &result);
	if (tracing == 1) {
		CHECK(selection > archive &&
			history.events[selection - 1].point_source == archive);
		CHECK(history.events[selection - 1].cause == archive &&
			history.events[selection - 1].work == 50);
	} else
		CHECK(!selection);
	simplex_progress_destroy(&tracker);
	simplex_solve_context_destroy(&context);
	lp_model_free(m);
	return 0;
}

static int rejected_bridges(void)
{
	struct lp_Model *m = fixture();
	struct simplex_SolveContext context;
	struct simplex_Execution h;
	struct simplex_ExecutionScope root;
	struct simplex_Presolve pre;
	struct simplex_ProgressTracker t;
	struct simplex_ProgressBridge b;
	struct simplex_ProgressSource source;
	struct lp_simplex_Options o;
	double reduced[3], original[3] = {.5, .5, .5};
	int count, j;
	CHECK(m);
	simplex_solve_context_init(&context, &h, &root);
	simplex_execution_set_parent(
		&root, simplex_execution_begin(
			       &root, SIMPLEX_EXECUTION_SOLVE, 0, 0, 0));
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_AUTO);
	CHECK(!simplex_presolve_run(&pre, m, o.primal_tolerance));
	CHECK(pre.reduced);
	pre.reduced->execution =
		simplex_solve_context_derive(&root, SIMPLEX_EXECUTION_PRESOLVE);
	CHECK(!simplex_progress_init_reported(&t, m, &o, &root));
	{
		struct simplex_ExecutionScope foreign = root;
		foreign.context = NULL;
		CHECK(simplex_progress_bridge(&b, &t.sink, &pre, NULL, 3, 3,
			      pre.reduced, &foreign,
			      NULL) == SIMPLEX_PROGRESS_BRIDGE_INVALID);
		CHECK(!b.storage);
		CHECK(simplex_progress_bridge(&b, &t.sink, &pre, NULL, INT_MAX,
			      3, pre.reduced, &root,
			      NULL) == SIMPLEX_PROGRESS_BRIDGE_INVALID);
		CHECK(!b.storage);
	}
	CHECK(!simplex_progress_bridge(
		&b, &t.sink, &pre, NULL, 3, 3, pre.reduced, &root, NULL));
	for (j = 0; j < pre.reduced->columns; j++)
		reduced[j] = original[pre.column_map[j]];
	memset(&source, 0, sizeof(source));
	source.producer = source.coordinates = pre.reduced->execution;
	source.input = simplex_problem_identity(pre.reduced);
	source.rows = pre.reduced->rows;
	source.columns = pre.reduced->columns;
	source.work = 7;
	count = h.count;
	source.columns = INT_MAX;
	b.sink.report(b.sink.context, reduced, NULL, 1., "bad-dimensions", 7, 1,
		&source);
	CHECK(!t.have_primal && h.count == count);
	source.columns = pre.reduced->columns;
	source.coordinates.model++;
	b.sink.report(b.sink.context, reduced, NULL, 1., "foreign-model", 7, 1,
		&source);
	CHECK(!t.have_primal && h.count == count);
	source.coordinates = pre.reduced->execution;
	{
		struct simplex_Problem foreign = *pre.reduced;
		double rhs[3];
		for (j = 0; j < pre.reduced->rows; j++)
			rhs[j] = pre.reduced->rhs[j];
		foreign.rhs = rhs;
		foreign.progress = &b.sink;
		simplex_progress_emit(
			&foreign, reduced, NULL, 1., "foreign-storage", 7, 1);
		CHECK(!t.have_primal && h.count == count);
	}
	simplex_problem_mark_changed(pre.reduced);
	b.sink.report(b.sink.context, reduced, NULL, 1., "stale-bridge", 7, 1,
		&source);
	CHECK(!t.have_primal && h.count == count);
	simplex_progress_bridge_destroy(&b);
	CHECK(!simplex_progress_bridge(
		&b, &t.sink, &pre, NULL, 3, 3, pre.reduced, &root, NULL));
	b.sink.report(b.sink.context, reduced, NULL, 1., "stale-source", 7, 1,
		&source);
	CHECK(!t.have_primal && h.count == count);
	source.input = simplex_problem_identity(pre.reduced);
	b.sink.report(
		b.sink.context, reduced, NULL, 1., "rebound", 7, 1, &source);
	CHECK(t.have_primal && fabs(t.best_primal - 1.5) < 1e-8);
	simplex_progress_bridge_destroy(&b);
	simplex_progress_destroy(&t);
	simplex_presolve_destroy(&pre);
	simplex_solve_context_destroy(&context);
	lp_model_free(m);
	return 0;
}

static int unknown_source(int full)
{
	struct lp_Model *m = fixture();
	struct simplex_SolveContext c;
	struct simplex_Execution h, foreign;
	struct simplex_ExecutionScope root, other;
	struct simplex_ProgressTracker t;
	struct simplex_ProgressSource source;
	struct lp_simplex_Options o;
	double worse[3] = {0., 1., 1.}, better[3] = {.5, .5, .5};
	int archive, i;
	char stage[100];
	CHECK(m);
	simplex_solve_context_init(&c, &h, &root);
	simplex_execution_set_parent(
		&root, simplex_execution_begin(
			       &root, SIMPLEX_EXECUTION_SOLVE, 0, 0, 0));
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_AUTO);
	CHECK(!simplex_progress_init_reported(&t, m, &o, &root));
	memset(&source, 0, sizeof(source));
	source.coordinates = source.producer = root;
	source.rows = source.columns = 3;
	source.work = 1;
	memset(stage, 'A', sizeof(stage) - 1);
	stage[sizeof(stage) - 1] = 0;
	t.sink.report(t.sink.context, worse, NULL, 1., stage, 1, 0, &source);
	archive = simplex_progress_point_source(&t, &root);
	CHECK(archive > 0);
	stage[0] = 'Z';
	CHECK(h.events[archive - 1].observation_stage[0] == 'A');
	CHECK(strlen(h.events[archive - 1].observation_stage) ==
		sizeof(h.events[archive - 1].observation_stage) - 1);
	other = root;
	memset(&foreign, 0, sizeof(foreign));
	other.history = &foreign;
	CHECK(!simplex_progress_point_source(&t, &other));
	if (full == 1)
		for (i = h.count; i < SIMPLEX_EXECUTION_CAPACITY; i++)
			(void)simplex_execution_begin(
				&root, SIMPLEX_EXECUTION_SOLVE, 0, 0, 0);
	if (full == 2)
		c.execution_serial = ULONG_MAX;
	source.work = 2;
	t.sink.report(t.sink.context, better, NULL, 1., "legacy", 2, 0,
		full ? &source : NULL);
	CHECK(t.have_primal && fabs(t.best_primal - 1.5) < 1e-8 &&
		!simplex_progress_point_source(&t, &root));
	CHECK(h.count == (full == 1	     ? SIMPLEX_EXECUTION_CAPACITY
				 : full == 2 ? archive + 1
					     : archive));
	/* Never reuse an older point's provenance after losing a new source. */
	simplex_progress_destroy(&t);
	simplex_solve_context_destroy(&c);
	lp_model_free(m);
	return 0;
}

/* Archive lifetime is independent of the real solver point's lifetime. */
static int archive_lifetime(int mode)
{
	struct lp_Model *m = fixture();
	struct simplex_Problem p;
	struct simplex_SolveContext c, other;
	struct simplex_Execution h;
	struct simplex_ExecutionScope root, foreign;
	struct simplex_ProgressTracker t;
	struct simplex_ExecutionEvent saved, *e;
	struct lp_simplex_Options o;
	struct lp_simplex_Result r;
	double x[3], point[3], objective;
	int id, i, selection;
	CHECK(m && !simplex_problem_from_model(&p, m));
	simplex_solve_context_init(&c, &h, &root);
	simplex_solve_context_init(&other, NULL, &foreign);
	simplex_solve_context_configure(&c, 0);
	p.execution = root;
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_DUAL_REVISED);
	CHECK(!simplex_progress_init_reported(&t, m, &o, &root));
	p.progress = &t.sink;
	CHECK(!simplex_dual_solve_problem(&p, &o, x, NULL, &r, 0));
	id = simplex_progress_point_source(&t, &root);
	CHECK(id > 0 && t.have_primal);
	saved = h.events[id - 1];
	memcpy(point, t.point, sizeof(point));
	objective = t.best_primal;
	if (mode == 1 || mode == 2) {
		memset(&h, 0, sizeof(h));
		foreign.history = &h;
		for (i = 0; i < id; i++)
			CHECK(simplex_execution_begin(
				      mode == 1 ? &root : &foreign,
				      SIMPLEX_EXECUTION_POINT_ARCHIVE, 0, 0,
				      0) == i + 1);
		/* Same complete payload, newly minted owner and serial. */
		e = h.events + id - 1;
		saved.serial = e->serial;
		saved.context = e->context;
		*e = saved;
	}
	e = h.events + id - 1;
	switch (mode) {
	case 3:
		e->complete = 0;
		break;
	case 4:
		e->serial = 0;
		break;
	case 5:
		e->context = &other;
		break;
	case 6:
		e->model++;
		break;
	case 7:
		e->kind = SIMPLEX_EXECUTION_POINT_SELECTION;
		break;
	case 8:
		h.count = INT_MAX;
		break;
	case 9:
		e->primal.quality = SIMPLEX_EVIDENCE_UNKNOWN;
		break;
	}
	CHECK(simplex_progress_point_source(&t, &root) == (mode ? 0 : id));
	CHECK(t.have_primal && t.best_primal == objective &&
		!memcmp(point, t.point, sizeof(point)) &&
		simplex_problem_primal_infeasibility(&t.original, t.point) <=
			o.primal_tolerance);
	/* Restore a valid append boundary for the diagnostic selection event. */
	if (mode == 8)
		h.count = id;
	selection = simplex_progress_record_restore(&t, &root, &r);
	CHECK(selection > id && h.events[selection - 1].point_source ==
					(mode && mode != 8 ? 0 : id));
	simplex_progress_destroy(&t);
	simplex_solve_context_destroy(&other);
	simplex_solve_context_destroy(&c);
	simplex_problem_destroy(&p);
	lp_model_free(m);
	return 0;
}

static int stale_producer(void)
{
	struct lp_Model *m = fixture();
	struct simplex_SolveContext context;
	struct simplex_Execution history, saved;
	struct simplex_ExecutionScope root, producer;
	struct simplex_ProgressTracker tracker;
	struct simplex_ProgressSource source;
	struct lp_simplex_Options options;
	double worse[3] = {0., 1., 1.}, better[3] = {.5, .5, .5};
	int i, archive;
	CHECK(m);
	simplex_solve_context_init(&context, &history, &root);
	simplex_execution_set_parent(
		&root, simplex_execution_begin(
			       &root, SIMPLEX_EXECUTION_SOLVE, 0, 0, 0));
	producer = root;
	simplex_execution_set_parent(
		&producer, simplex_execution_begin(
				   &root, SIMPLEX_EXECUTION_BACKEND, 0, 0, 0));
	lp_simplex_default_options(&options, lp_simplex_ALGORITHM_AUTO);
	CHECK(!simplex_progress_init_reported(&tracker, m, &options, &root));
	memset(&source, 0, sizeof(source));
	source.coordinates = root;
	source.producer = producer;
	source.rows = source.columns = 3;
	source.work = 1;
	tracker.sink.report(
		tracker.sink.context, worse, NULL, 1., "before", 1, 1, &source);
	saved = history;
	memset(&history, 0, sizeof(history));
	for (i = 0; i < saved.count; i++)
		(void)simplex_execution_begin(&root, saved.events[i].kind,
			saved.events[i].detail, 0, 0);
	source.work = 2;
	tracker.sink.report(
		tracker.sink.context, better, NULL, 1., "after", 2, 1, &source);
	archive = simplex_progress_point_source(&tracker, &root);
	CHECK(archive > 0 && tracker.have_primal &&
		fabs(tracker.best_primal - 1.5) < 1e-8);
	CHECK(!history.events[archive - 1].cause &&
		!history.events[archive - 1].point_source);
	CHECK(!history.events[archive - 1].parent);
	simplex_progress_destroy(&tracker);
	simplex_solve_context_destroy(&context);
	lp_model_free(m);
	return 0;
}

int main(void)
{
	int unscaled, tracing, mode;
	for (unscaled = 0; unscaled < 2; unscaled++)
		for (tracing = 0; tracing < 3; tracing++)
			CHECK(!actual_chain(unscaled, tracing));
	for (mode = 0; mode < 10; mode++)
		CHECK(!archive_lifetime(mode));
	CHECK(!rejected_bridges() && !unknown_source(0) && !unknown_source(1) &&
		!unknown_source(2) && !stale_producer());
	return 0;
}
