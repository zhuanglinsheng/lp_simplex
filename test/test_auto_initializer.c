/**
 * @file test_auto_initializer.c
 * @brief Exact initialization lease, report isolation and real NNLS model
 * binding.
 */

#include "simplex_stage.h"
#include "simplex_solve_context.h"
#include "simplex_execution.h"
#include "simplex_auto_routes.h"
#include "simplex_pan_config.h"
#include "simplex_presolve.h"
#include <lp_simplex/lp_simplex.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#define CHECK(c) \
	do { \
		if (!(c)) { \
			fprintf(stderr, "initializer line %d mode=%d: %s\n", \
				__LINE__, mode, #c); \
			return 1; \
		} \
	} while (0)
static int mode, dual_calls, initialized, seed_ok, seen_grant, report_ok;
static struct simplex_StageReport *dual_report;

static int dual(const struct simplex_Problem *p,
	const struct lp_simplex_Options *o, double *x, double *y,
	struct lp_simplex_Result *r, int propagate)
{
	(void)o;
	(void)y;
	(void)propagate;
	dual_calls++;
	dual_report = p->report;
	p->report->phase = 11;
	x[0] = x[1] = 0.;
	r->status = lp_simplex_PrecisionError;
	r->iterations = 1;
	return -1;
}

static int initialize(void *context, const struct simplex_Problem *p,
	const struct lp_simplex_Options *o, double *x, int *work)
{
	(void)context;
	initialized++;
	seen_grant = o->iteration_limit;
	report_ok =
		p->report && p->report != dual_report && p->report->phase == -1;
	p->report->phase = 22;
	p->report->lower_bound.quality = SIMPLEX_EVIDENCE_VERIFIED;
	p->report->lower_bound.kind = SIMPLEX_OBSERVATION_LOWER_BOUND;
	p->report->lower_bound.space = SIMPLEX_SPACE_INPUT;
	p->report->lower_bound.work = 1;
	p->report->lower_bound.value = 999.;
	p->report->basis_available = p->report->factor_available = 1;
	p->report->dual.quality = SIMPLEX_EVIDENCE_VERIFIED;
	x[0] = x[1] = .5;
	*work = mode == 1   ? o->iteration_limit
		: mode == 2 ? o->iteration_limit + 1
			    : 1;
	if (mode == 3)
		((struct simplex_Problem *)p)->columns = INT_MAX;
	return 0; /* Candidate validity is independent of this proposal flag. */
}

static int hybrid(const struct simplex_Problem *p,
	const struct lp_simplex_Options *o, const double *seed, double *x,
	struct lp_simplex_Result *r)
{
	(void)o;
	seed_ok = mode == 2 || mode == 4
			  ? !seed
			  : seed && seed[0] == .5 && seed[1] == .5;
	p->report->phase = 33;
	x[0] = x[1] = .5;
	r->iterations = 1;
	r->status = lp_simplex_PrecisionError;
	return -1;
}

static int controller(int tracing)
{
	struct lp_Model *m = lp_model_create(1, 2);
	struct simplex_Problem p;
	struct simplex_SolveContext context;
	struct simplex_Execution trace;
	struct simplex_AutoHistory history;
	struct simplex_AutoInitializer initializer;
	struct simplex_AutoKernels kernels;
	struct simplex_AutoBindings bindings = {NULL, NULL, 0};
	struct lp_simplex_Options o;
	struct lp_simplex_Result r;
	double x[2];
	int i;
	bindings.kernels = &kernels;
	CHECK(m);
	m->constraints[0].type = optm_CONS_T_EQ;
	m->constraints[0].rhs = 1.;
	m->constraints[0].coef[0] = m->constraints[0].coef[1] = 1.;
	CHECK(!simplex_problem_from_model(&p, m));
	simplex_solve_context_init(
		&context, tracing ? &trace : NULL, &p.execution);
	if (tracing == 2)
		for (i = 0; i < SIMPLEX_EXECUTION_CAPACITY; i++)
			(void)simplex_execution_begin(
				&p.execution, SIMPLEX_EXECUTION_SOLVE, 0, 0, 0);
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_AUTO);
	o.iteration_limit = 100000;
	kernels.dual = dual;
	kernels.hybrid = hybrid;
	initializer.context = NULL;
	initializer.run = mode == 4 ? NULL : initialize;
	dual_calls = initialized = seen_grant = seed_ok = report_ok = 0;
	CHECK(simplex_auto_solve(&p, &o, x, &r, 0, &bindings, &initializer,
		      &history) != 0);
	CHECK(r.status != lp_simplex_Success &&
		r.iterations <= o.iteration_limit);
	CHECK(history.stages[0].report.phase == 11);
	if (mode == 3) {
		CHECK(history.count == 2 && dual_calls == 1);
		CHECK(history.outcome.exit == SIMPLEX_AUTO_EXIT_INPUT_CHANGED);
		CHECK(history.stages[1].report.failure ==
			SIMPLEX_FAILURE_STAGE_CONTRACT);
	} else {
		CHECK(dual_calls == 2 && seed_ok);
		CHECK(history.count == (mode == 4 ? 3 : 4));
		CHECK(history.stages[mode == 4 ? 1 : 2].report.phase == 33);
	}
	if (mode == 4)
		CHECK(!initialized);
	else {
		const struct simplex_AutoStageRecord *record =
			history.stages + 1;
		CHECK(initialized == 1 && report_ok);
		CHECK(seen_grant == SIMPLEX_AUTO_INITIALIZER_OBSERVATIONS *
					    SIMPLEX_AUTO_PROBE_CAP);
		CHECK(record->action.grant == seen_grant &&
			record->report.phase == 22);
		CHECK(record->report.lower_bound.quality ==
			SIMPLEX_EVIDENCE_UNKNOWN);
		CHECK(record->report.dual.quality == SIMPLEX_EVIDENCE_UNKNOWN);
		CHECK(!record->report.basis_available &&
			!record->report.factor_available);
		CHECK(history.best_bound == -HUGE_VAL);
		if (mode == 1 || mode == 2)
			CHECK(record->work == seen_grant);
		if (mode == 2) {
			CHECK(record->reported_work == seen_grant + 1 &&
				!record->accepted_point);
			CHECK(record->report.failure ==
				SIMPLEX_FAILURE_STAGE_CONTRACT);
		}
		if (tracing == 1)
			CHECK(trace.events[record->action.execution_id - 1]
					.grant == seen_grant);
	}
	simplex_solve_context_destroy(&context);
	simplex_problem_destroy(&p);
	lp_model_free(m);
	return 0;
}

static int real_adapter(int tracing, int fault)
{
	struct lp_Model *m = lp_model_create(1, 6);
	struct simplex_Problem p, view;
	struct simplex_Presolve pre;
	struct simplex_AutoOriginalInitializer c;
	struct simplex_PanConfig controls;
	struct simplex_SolveContext context, foreign;
	struct simplex_Execution trace;
	struct simplex_ExecutionScope root;
	struct simplex_StageReport report;
	struct lp_simplex_Options o;
	double original[6], point[6], other_rhs = 1., other_cost[6] = {0.};
	int map[6], j, work = 999, valid, action;
	unsigned long revision;
	CHECK(m);
	m->constraints[0].type = optm_CONS_T_EQ;
	m->constraints[0].rhs = 1.;
	for (j = 0; j < 6; j++)
		m->constraints[0].coef[j] = 1.;
	CHECK(!simplex_problem_from_model(&p, m));
	simplex_solve_context_init(&context, tracing ? &trace : NULL, &root);
	p.execution =
		simplex_solve_context_derive(&root, SIMPLEX_EXECUTION_PRESOLVE);
	for (j = 0; j < 6; j++) {
		p.bounds[j].b_type = optm_BOUND_T_BS;
		p.bounds[j].ub = 1.;
		map[j] = j;
		original[j] = point[j] = NAN;
	}
	simplex_problem_mark_changed(&p);
	memset(&pre, 0, sizeof(pre));
	pre.original = m;
	pre.reduced = &p;
	pre.column_map = map;
	CHECK(simplex_pan_config_default(&controls, SIMPLEX_PAN_STANDARD,
		SIMPLEX_PAN_BASIS_DIRECTION, SIMPLEX_PAN_NO_ORACLE));
	c.model = m;
	c.binding = *m;
	c.presolve = &pre;
	c.controls = &controls;
	c.original = original;
	c.scope = root;
	if (fault == 9)
		p.revision =
			ULONG_MAX; /* Current-call binding, not checkpoint reuse. */
	revision = p.revision;
	c.target = simplex_problem_identity(&p);
	view = p;
	simplex_stage_reset(&report);
	view.report = &report;
	action = simplex_execution_begin(&p.execution, SIMPLEX_EXECUTION_AUTO,
		SIMPLEX_AUTO_INITIALIZER, 0, 0);
	simplex_execution_set_parent(&view.execution, action);
	if (tracing == 2)
		while (trace.count < SIMPLEX_EXECUTION_CAPACITY)
			(void)simplex_execution_begin(
				&root, SIMPLEX_EXECUTION_SOLVE, 0, 0, 0);
	switch (fault) {
	case 1:
		view.columns = INT_MAX;
		break;
	case 2:
		view.rhs = &other_rhs;
		break;
	case 3:
		p.revision++;
		break;
	case 4:
		m->n = INT_MAX;
		break;
	case 5:
		map[2] = INT_MAX;
		break;
	case 6:
		m->objective = other_cost;
		break;
	case 7:
		c.scope.context = &foreign;
		break;
	case 8:
		view.matrix.columns = INT_MAX;
		break;
	case 10:
		c.original = point;
		break;
	}
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_AUTO);
	o.iteration_limit = (fault == 11 ? 2 : 1) *
			    SIMPLEX_AUTO_INITIALIZER_OBSERVATIONS *
			    SIMPLEX_AUTO_PROBE_CAP;
	valid = simplex_auto_routes_initialize_original(
		&c, &view, &o, point, &work);
	if (fault == 0 || fault == 9 || fault == 11) {
		CHECK(valid && work >= 0 && work <= o.iteration_limit);
		CHECK(simplex_problem_primal_infeasibility(&p, point) <=
			o.primal_tolerance);
		CHECK(report.services.calls[SIMPLEX_PAN_SERVICE_NNLS] == 1);
		CHECK(report.services.work[SIMPLEX_PAN_SERVICE_NNLS] == work);
		CHECK(report.stop == SIMPLEX_STOP_HANDOFF);
		if (tracing == 1) {
			int service = 0, representation = 0;
			for (j = 0; j < trace.count; j++) {
				const struct simplex_ExecutionEvent *e =
					trace.events + j;
				if (e->kind == SIMPLEX_EXECUTION_PAN_MODEL) {
					representation = e->model;
					CHECK(e->source_model == root.model &&
						e->parent == action);
				}
				if (e->kind == SIMPLEX_EXECUTION_PAN_SERVICE) {
					service++;
					CHECK(e->detail ==
							SIMPLEX_PAN_SERVICE_NNLS &&
						e->parent == action);
					CHECK(e->model == representation &&
						e->source_model == root.model);
					CHECK(e->grant == o.iteration_limit &&
						e->work == work &&
						e->service_valid);
					CHECK(e->service_output ==
						SIMPLEX_PAN_OUTPUT_NONNEGATIVE_EQUALITY_POINT);
				}
			}
			CHECK(service == 1 && representation != root.model &&
				representation != p.execution.model);
		}
	} else {
		CHECK(!valid && work == 0 &&
			report.failure == SIMPLEX_FAILURE_STAGE_CONTRACT);
		for (j = 0; j < 6; j++)
			CHECK(isnan(point[j]) && isnan(original[j]));
		CHECK(report.services.calls[SIMPLEX_PAN_SERVICE_NNLS] == 0);
	}
	p.revision = revision;
	m->n = c.binding.n;
	m->objective = c.binding.objective;
	simplex_solve_context_destroy(&context);
	simplex_problem_destroy(&p);
	lp_model_free(m);
	return 0;
}

int main(void)
{
	int tracing, fault;
	for (mode = 0; mode < 5; mode++)
		for (tracing = 0; tracing < 3; tracing++)
			if (controller(tracing))
				return 1;
	for (tracing = 0; tracing < 3; tracing++)
		for (fault = 0; fault < 12; fault++)
			if (real_adapter(tracing, fault))
				return 1;
	return 0;
}
