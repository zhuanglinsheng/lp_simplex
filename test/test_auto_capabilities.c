/**
 * @file test_auto_capabilities.c
 * @brief Route objects are bindings, not proof of executable capabilities.
 */

#include "simplex_stage.h"
#include "simplex_solve_context.h"
#include "simplex_execution.h"
#include "simplex_auto.h"
#include "simplex_auto_routes.h"
#include "simplex_auto_policy.h"
#include <lp_simplex/lp_simplex.h>
#include <math.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>

#define CHECK(c) \
	do { \
		if (!(c)) { \
			fprintf(stderr, "AUTO capabilities %d: %s\n", \
				__LINE__, #c); \
			return 1; \
		} \
	} while (0)
static int first_grant, calls, hybrid_calls, route_calls, oracle_failure;
static struct simplex_AutoHistory *poison_output;
static struct simplex_AutoKernels *changed_kernels;
static struct simplex_AutoRoutes *changed_routes;

static int dual(const struct simplex_Problem *p,
	const struct lp_simplex_Options *o, double *x, double *y,
	struct lp_simplex_Result *r, int propagate)
{
	(void)x;
	(void)y;
	(void)propagate;
	if (!calls++) {
		first_grant = o->iteration_limit;
		if (poison_output)
			memset(poison_output->admission, 0,
				sizeof(poison_output->admission));
		if (changed_kernels) {
			changed_kernels->dual = NULL;
			changed_kernels->hybrid = NULL;
			changed_routes->solve = NULL;
		}
	}
	r->iterations = 0;
	r->status = lp_simplex_PrecisionError;
	p->report->stop = SIMPLEX_STOP_NUMERICAL;
	return lp_simplex_EXIT_FAILURE;
}

static int hybrid(const struct simplex_Problem *p,
	const struct lp_simplex_Options *o, const double *seed, double *x,
	struct lp_simplex_Result *r)
{
	hybrid_calls++;
	(void)o;
	(void)seed;
	(void)x;
	r->iterations = 0;
	r->status = lp_simplex_PrecisionError;
	p->report->stop = SIMPLEX_STOP_NUMERICAL;
	if (oracle_failure)
		p->report->failure = SIMPLEX_FAILURE_IMPLICIT_ORACLE;
	return lp_simplex_EXIT_FAILURE;
}

static int route(int kind, const struct simplex_Problem *p,
	const struct lp_simplex_Options *o, const double *seed, double *x,
	struct lp_simplex_Result *r)
{
	(void)kind;
	(void)o;
	(void)seed;
	(void)x;
	route_calls++;
	r->iterations = 0;
	r->status = lp_simplex_PrecisionError;
	p->report->stop = SIMPLEX_STOP_NUMERICAL;
	p->report->failure = SIMPLEX_FAILURE_IMPLICIT_ORACLE;
	return lp_simplex_EXIT_FAILURE;
}

static int empty_binding(void)
{
	struct lp_Model *m = lp_model_create(1, 2);
	struct simplex_Problem p;
	struct simplex_AutoKernels kernels;
	struct simplex_AutoBindings route_bindings = {NULL, NULL, 0};
	struct simplex_AutoRoutes empty;
	struct simplex_AutoBindings bindings = {NULL, NULL, 0};
	struct simplex_AutoHistory history[2];
	struct lp_simplex_Options o;
	struct lp_simplex_Result r;
	double x[2];
	int i, grant[2];
	route_bindings.kernels = &kernels;
	bindings.kernels = &kernels;
	bindings.routes = &empty;
	CHECK(m);
	m->constraints[0].type = optm_CONS_T_EQ;
	m->constraints[0].rhs = 1.;
	m->constraints[0].coef[0] = m->constraints[0].coef[1] = 1.;
	CHECK(!simplex_problem_from_model(&p, m));
	kernels.dual = dual;
	kernels.hybrid = hybrid;
	memset(&empty, 0, sizeof(empty));
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_AUTO);
	o.iteration_limit = 10000;
	for (i = 0; i < 2; i++) {
		calls = first_grant = 0;
		CHECK(simplex_auto_solve(&p, &o, x, &r, 0,
			      (i ? &bindings : &route_bindings), NULL,
			      history + i) != 0);
		grant[i] = first_grant;
		printf("empty=%d first-grant=%d actions=%d\n", i, grant[i],
			history[i].count);
		CHECK(calls == 2 && !isfinite(x[0]) && !history[i].have_point);
	}
	CHECK(grant[0] == grant[1] && history[0].count == history[1].count);
	simplex_problem_destroy(&p);
	lp_model_free(m);
	return 0;
}

static int bindings(void)
{
	const struct simplex_AutoBindings absent_bindings = {NULL, NULL, 0};
	struct lp_Model *m = lp_model_create(1, 2);
	struct simplex_Problem p;
	struct simplex_SolveContext context;
	struct simplex_Execution events;
	struct simplex_AutoKernels kernels;
	struct simplex_AutoRoutes routes;
	struct simplex_AutoBindings bindings = {NULL, NULL, 0};
	struct simplex_AutoHistory history;
	struct simplex_StageReport report;
	struct lp_simplex_Options o;
	struct lp_simplex_Result r;
	double x[2];
	int mask, tracing, i, available, expected;
	bindings.kernels = &kernels;
	bindings.routes = &routes;
	CHECK(m);
	m->constraints[0].type = optm_CONS_T_EQ;
	m->constraints[0].rhs = 1.;
	m->constraints[0].coef[0] = m->constraints[0].coef[1] = 1.;
	CHECK(!simplex_problem_from_model(&p, m));
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_AUTO);
	o.iteration_limit = 10000;
	oracle_failure = 1;
	for (tracing = 0; tracing < 3; tracing++) {
		for (mask = 0; mask < 8; mask++) {
			simplex_solve_context_init(&context,
				tracing ? &events : NULL, &p.execution);
			if (tracing == 2)
				for (i = 0; i < SIMPLEX_EXECUTION_CAPACITY; i++)
					(void)simplex_execution_begin(
						&p.execution,
						SIMPLEX_EXECUTION_SOLVE, 0, 0,
						0);
			kernels.dual = mask & 1 ? dual : NULL;
			kernels.hybrid = mask & 2 ? hybrid : NULL;
			routes.solve = mask & 4 ? route : NULL;
			CHECK(simplex_auto_routes_available(
				      SIMPLEX_AUTO_DUAL_EXPLORE, &kernels,
				      &routes) == !!(mask & 1));
			CHECK(simplex_auto_routes_available(
				      SIMPLEX_AUTO_DUAL_FALLBACK, &kernels,
				      &routes) == !!(mask & 1));
			for (i = 0; i < SIMPLEX_AUTO_ROUTE_COUNT; i++)
				CHECK(simplex_auto_routes_available(
					      i, &kernels, &routes) ==
					(i == SIMPLEX_AUTO_HYBRID
							? !!(mask & 2)
							: !!(mask & 4)));
			calls = hybrid_calls = route_calls = first_grant = 0;
			p.report = &report;
			CHECK(simplex_auto_solve(&p, &o, x, &r, 0, &bindings,
				      NULL, &history) != 0);
			CHECK(!isfinite(x[0]) && !history.have_point &&
				!r.iterations);
			if (!(mask & 1)) {
				for (i = 0; i < SIMPLEX_AUTO_ROUTE_COUNT; i++)
					CHECK(history.admission[i].binding ==
							SIMPLEX_AUTO_BINDING_UNCHECKED &&
						history.admission[i].model ==
							SIMPLEX_AUTO_MODEL_UNCHECKED &&
						history.admission[i]
								.condition ==
							SIMPLEX_AUTO_CONDITION_UNCHECKED &&
						history.admission[i]
								.prerequisite
								.record == -1);
				CHECK(!calls && !hybrid_calls && !route_calls &&
					!history.count);
				CHECK(r.status == lp_simplex_CondUnsatisfied &&
					history.outcome.exit ==
						SIMPLEX_AUTO_EXIT_CAPABILITY &&
					report.stop ==
						SIMPLEX_STOP_NOT_APPLICABLE);
				if (tracing == 1)
					CHECK(events.count == 1 &&
						events.events[0].reason ==
							SIMPLEX_REASON_REQUIRED_CAPABILITY);
			} else {
				available = !!(mask & 2) + (mask & 4 ? 3 : 0);
				CHECK(first_grant ==
					o.iteration_limit -
						simplex_auto_recovery_reserve(
							o.iteration_limit,
							available
								? available + 1
								: 0));
				CHECK(calls == (available ? 2 : 1));
				CHECK(hybrid_calls == !!(mask & 2));
				expected = mask & 4 ? 4 + !!(mask & 2) : 0;
				CHECK(route_calls == expected);
				CHECK(history.count ==
					calls + hybrid_calls + route_calls);
				for (i = 0; i < SIMPLEX_AUTO_ROUTE_COUNT; i++) {
					const struct simplex_AutoAdmission *a =
						history.admission + i;
					int present = i == SIMPLEX_AUTO_HYBRID
							      ? !!(mask & 2)
							      : !!(mask & 4);
					int conditional =
						i == SIMPLEX_AUTO_EXPLICIT ||
						i == SIMPLEX_AUTO_STANDARD_EXPLICIT;
					int opened = i == SIMPLEX_AUTO_EXPLICIT
							     ? !!(mask & 2)
							     : !!(mask & 4);
					int k, producer = -1;
					CHECK(a->binding ==
						(present ? SIMPLEX_AUTO_BINDING_READY
							 : SIMPLEX_AUTO_BINDING_MISSING));
					CHECK(a->model ==
						(!present ? SIMPLEX_AUTO_MODEL_UNCHECKED
							: i == SIMPLEX_AUTO_SINGLETON
								? SIMPLEX_AUTO_MODEL_SINGLETON_PAIRS
								: SIMPLEX_AUTO_MODEL_READY));
					CHECK(a->condition ==
						(!conditional || opened
								? SIMPLEX_AUTO_CONDITION_READY
								: SIMPLEX_AUTO_CONDITION_ORACLE));
					if (conditional && opened) {
						for (k = 0; k < history.count;
							k++)
							if (history.stages[k]
									.action
									.stage ==
								(i == SIMPLEX_AUTO_EXPLICIT
										? SIMPLEX_AUTO_HYBRID
										: SIMPLEX_AUTO_STANDARD_IMPLICIT))
								producer = k;
						CHECK(producer >= 0 &&
							a->prerequisite.record ==
								producer &&
							a->prerequisite.event ==
								history.stages[producer]
									.action
									.execution_id);
					} else
						CHECK(a->prerequisite.record ==
								-1 &&
							!a->prerequisite.event);
				}
			}
			simplex_solve_context_destroy(&context);
		}
	}
	/* Poisoning output admission cannot change the committed route schedule. */
	memset(&p.execution, 0, sizeof(p.execution));
	kernels.dual = dual;
	kernels.hybrid = hybrid;
	routes.solve = route;
	calls = hybrid_calls = route_calls = 0;
	poison_output = &history;
	CHECK(simplex_auto_solve(&p, &o, x, &r, 0, &bindings, NULL, &history) !=
		0);
	poison_output = NULL;
	CHECK(calls == 2 && hybrid_calls == 1 && route_calls == 5 &&
		history.count == 8);
	/* Entry copies callback capabilities before the first callback. Mutating
	 * the caller's tables cannot replace already admitted implementations. */
	for (tracing = 0; tracing < 3; tracing++) {
		simplex_solve_context_init(
			&context, tracing ? &events : NULL, &p.execution);
		if (tracing == 2)
			for (i = 0; i < SIMPLEX_EXECUTION_CAPACITY; i++)
				(void)simplex_execution_begin(&p.execution,
					SIMPLEX_EXECUTION_SOLVE, 0, 0, 0);
		kernels.dual = dual;
		kernels.hybrid = hybrid;
		routes.solve = route;
		calls = hybrid_calls = route_calls = 0;
		changed_kernels = &kernels;
		changed_routes = &routes;
		CHECK(simplex_auto_solve(&p, &o, x, &r, 0, &bindings, NULL,
			      &history) != 0);
		changed_kernels = NULL;
		changed_routes = NULL;
		CHECK(!kernels.dual && !kernels.hybrid && !routes.solve);
		CHECK(calls == 2 && hybrid_calls == 1 && route_calls == 5 &&
			history.count == 8 && !r.iterations);
		simplex_solve_context_destroy(&context);
	}
	/* A wholly absent required kernel object has the same explicit exit. */
	memset(&p.execution, 0, sizeof(p.execution));
	calls = 0;
	CHECK(simplex_auto_solve(
		      &p, &o, x, &r, 0, &absent_bindings, NULL, &history) != 0);
	CHECK(!calls && history.outcome.exit == SIMPLEX_AUTO_EXIT_CAPABILITY);
	CHECK(!simplex_auto_routes_has_checkpoint(NULL));
	CHECK(!simplex_auto_routes_available(
		SIMPLEX_AUTO_INITIALIZER, NULL, NULL));
	CHECK(!simplex_auto_routes_available(100, NULL, NULL));
	/* The invoke guard cannot call absent functions or invent an output. */
	simplex_stage_reset(&report);
	report.point_available = report.basis_available = 1;
	x[0] = .25;
	CHECK(simplex_auto_routes_invoke(SIMPLEX_AUTO_HYBRID, &p, &o, NULL, x,
		      &r, 0, NULL, NULL, NULL) != 0);
	CHECK(!r.iterations && r.status == lp_simplex_CondUnsatisfied &&
		isnan(r.objective) && x[0] == .25 && !report.point_available &&
		!report.basis_available &&
		report.stop == SIMPLEX_STOP_NOT_APPLICABLE);
	/* A production Dual binding remains usable without recovery callbacks. */
	memset(&p.execution, 0, sizeof(p.execution));
	kernels.dual = simplex_dual_solve_problem;
	kernels.hybrid = NULL;
	routes.solve = NULL;
	CHECK(simplex_auto_routes_has_checkpoint(&kernels));
	CHECK(!simplex_auto_solve(&p, &o, x, &r, 0, &bindings, NULL, &history));
	CHECK(r.status == lp_simplex_Success && history.count == 1 &&
		history.stages[0].action.grant == o.iteration_limit &&
		history.stages[0].action.failure_action ==
			SIMPLEX_AUTO_NO_RECOVERY &&
		simplex_problem_primal_infeasibility(&p, x) <=
			o.primal_tolerance);
	simplex_problem_destroy(&p);
	lp_model_free(m);
	return 0;
}

static int model_reasons(void)
{
	struct lp_Model *m = lp_model_create(1, 2);
	struct simplex_Problem p, changed;
	int reason;
	CHECK(m);
	m->constraints[0].type = optm_CONS_T_EQ;
	m->constraints[0].coef[0] = 1.;
	m->constraints[0].coef[1] = -1.;
	m->objective[0] = m->objective[1] = 1.;
	CHECK(!simplex_problem_from_model(&p, m));
	CHECK(simplex_auto_route_admission(&p, SIMPLEX_AUTO_SINGLETON) ==
		SIMPLEX_AUTO_MODEL_READY);
	CHECK(simplex_auto_route_admission(&p, -100) ==
		SIMPLEX_AUTO_MODEL_ROUTE);
	CHECK(simplex_auto_route_admission(NULL, SIMPLEX_AUTO_HYBRID) ==
		SIMPLEX_AUTO_MODEL_SHAPE);
	changed = p;
	changed.rows = INT_MAX / 4;
	CHECK(simplex_auto_route_admission(&changed, SIMPLEX_AUTO_HYBRID) ==
		SIMPLEX_AUTO_MODEL_INDEX);
	changed = p;
	changed.matrix.nonzeros = INT_MAX;
	CHECK(simplex_auto_route_admission(&changed, SIMPLEX_AUTO_CLASSIC) ==
		SIMPLEX_AUTO_MODEL_INDEX);
	changed = p;
	changed.rows = 100000000;
	CHECK(simplex_auto_route_admission(&changed, SIMPLEX_AUTO_CLASSIC) ==
		SIMPLEX_AUTO_MODEL_CONVERSION);
	CHECK(simplex_auto_route_admission(&changed, SIMPLEX_AUTO_TABLEAU) ==
		SIMPLEX_AUTO_MODEL_DENSE);
	changed = p;
	changed.rows = 100000;
	reason = simplex_auto_route_admission(&changed, SIMPLEX_AUTO_HYBRID);
	CHECK(reason == SIMPLEX_AUTO_MODEL_READY ||
		reason == SIMPLEX_AUTO_MODEL_DENSE_QR);
	CHECK(simplex_auto_route_allowed(&changed, SIMPLEX_AUTO_HYBRID) ==
		(reason == SIMPLEX_AUTO_MODEL_READY));
	printf("native-resource=%s\n",
		simplex_auto_route_admission_name(reason));
	p.bounds[0].lb = 1.;
	CHECK(simplex_auto_route_admission(&p, SIMPLEX_AUTO_SINGLETON) ==
		SIMPLEX_AUTO_MODEL_SINGLETON_BOUNDS);
	p.bounds[0].lb = 0.;
	p.row_type[0] = optm_CONS_T_GE;
	CHECK(simplex_auto_route_admission(&p, SIMPLEX_AUTO_SINGLETON) ==
		SIMPLEX_AUTO_MODEL_SINGLETON_ROWS);
	p.row_type[0] = optm_CONS_T_EQ;
	p.objective[0] = 0.;
	CHECK(simplex_auto_route_admission(&p, SIMPLEX_AUTO_SINGLETON) ==
		SIMPLEX_AUTO_MODEL_SINGLETON_PAIRS);
	CHECK(!strcmp(simplex_auto_route_admission_name(INT_MAX), "unknown"));
	simplex_problem_destroy(&p);
	lp_model_free(m);
	return 0;
}

int main(void)
{
	return empty_binding() || bindings() || model_reasons();
}
