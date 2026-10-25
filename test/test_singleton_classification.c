/**
 * @file test_singleton_classification.c
 * @brief Dualized singleton recovery must classify its input model, not the
 * dual LP.
 */

#include "simplex_stage.h"
#include "simplex_solve_context.h"
#include "simplex_singleton_dual.h"
#include "simplex_auto.h"
#include "simplex_auto_routes.h"
#include "simplex_execution.h"
#include <lp_simplex/lp_simplex.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define CHECK(c) \
	do { \
		if (!(c)) { \
			fprintf(stderr, "singleton classification %d: %s\n", \
				__LINE__, #c); \
			return 1; \
		} \
	} while (0)

static struct lp_Model *fixture(int core, double rhs, double cost)
{
	struct lp_Model *m = lp_model_create(1, core ? 3 : 2);
	if (!m)
		return NULL;
	m->constraints[0].type = optm_CONS_T_EQ;
	m->constraints[0].rhs = rhs;
	m->constraints[0].coef[0] = 1.;
	m->constraints[0].coef[1] = -1.;
	m->objective[0] = m->objective[1] = 1.;
	if (core) {
		m->constraints[0].coef[2] = -1.;
		m->objective[2] = cost;
	}
	return m;
}

static int actual(void)
{
	double costs[4] = {-3., -1., -.5, 0.};
	struct lp_simplex_Options o;
	struct lp_simplex_Result r;
	struct simplex_StageReport report;
	struct simplex_SolveContext context;
	struct simplex_Execution h;
	struct simplex_ExecutionScope scope;
	struct lp_Model *m;
	double x[3], rhs, objective;
	int core, c, b, d, i, applicable, state, outer, model, backend;
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_DUAL_REVISED);
	o.iteration_limit = 100;
	for (core = 0; core < 2; core++)
		for (c = 0; c < (core ? 4 : 1); c++)
			for (b = -1; b <= 1; b++)
				for (d = 0; d < 3; d++) {
					rhs = 2. * b;
					m = fixture(core, rhs, costs[c]);
					CHECK(m);
					simplex_solve_context_init(&context,
						d ? &h : NULL, &scope);
					if (d == 2)
						while (h.count <
							SIMPLEX_EXECUTION_CAPACITY)
							(void)simplex_execution_begin(
								&scope,
								SIMPLEX_EXECUTION_SOLVE,
								0, 0, 0);
					memset(&report, 0xff, sizeof(report));
					memset(&r, 0xff, sizeof(r));
					x[0] = x[1] = x[2] = NAN;
					state = simplex_singleton_dual_solve_reported(
						m, &o, x, &r, &applicable,
						&report, &scope);
					CHECK(applicable &&
						report.singleton
							.structure_verified &&
						report.singleton.input_model ==
							1);
					CHECK(report.singleton
							.input_point_checked &&
						report.point_available &&
						r.primal_infeasibility <=
							o.primal_tolerance);
					CHECK(fabs(x[0] - x[1] -
						      (core ? x[2] : 0.) -
						      rhs) <=
						o.primal_tolerance);
					CHECK(!report.basis_available &&
						!report.factor_available &&
						!report.lower_bound.quality);
					CHECK(r.iterations >= 0 &&
						r.iterations <=
							o.iteration_limit &&
						report.phase == -1);
					if (core && c == 0) {
						CHECK(state != 0 &&
							r.status ==
								lp_simplex_Unboundedness &&
							report.stop ==
								SIMPLEX_STOP_UNBOUNDED);
						CHECK(report.singleton
								.auxiliary_available &&
							report.singleton.auxiliary_status ==
								lp_simplex_Infeasibility);
						CHECK(report.singleton
								.mapping ==
							SIMPLEX_SINGLETON_UNBOUNDED);
						CHECK(isnan(r.dual_infeasibility) &&
							r.dual_residual_space ==
								lp_simplex_RESIDUAL_UNAVAILABLE);
					} else {
						objective =
							core && rhs < 0.
								? -rhs * costs[c]
								: fabs(rhs);
						CHECK(!state &&
							r.status ==
								lp_simplex_Success &&
							fabs(r.objective -
								objective) <
								1e-10);
						CHECK(report.singleton.mapping ==
								SIMPLEX_SINGLETON_OPTIMUM &&
							report.stop ==
								SIMPLEX_STOP_OPTIMAL);
						CHECK(isfinite(
							      r.dual_infeasibility) &&
							r.dual_infeasibility <=
								o.dual_tolerance &&
							r.dual_residual_space ==
								lp_simplex_RESIDUAL_ORIGINAL);
					}
					CHECK(report.singleton
							.auxiliary_available ==
						core);
					CHECK(report.singleton
							.auxiliary_model ==
						(core ? 2 : 0));
					if (!core)
						CHECK(r.iterations == 0 &&
							context.models == 1);
					if (d == 1) {
						outer = model = backend = 0;
						for (i = 0; i < h.count; i++) {
							const struct simplex_ExecutionEvent
								*e = h.events +
								     i;
							CHECK(e->complete &&
								e->context ==
									&context);
							if (e->kind ==
								SIMPLEX_EXECUTION_SINGLETON) {
								outer = i + 1;
								CHECK(e->model ==
										1 &&
									e->status ==
										r.status &&
									e->singleton.mapping ==
										report.singleton
											.mapping);
							}
							if (e->kind ==
								SIMPLEX_EXECUTION_SINGLETON_MODEL) {
								model = i + 1;
								CHECK(e->model ==
										2 &&
									e->source_model ==
										1 &&
									e->parent ==
										outer &&
									!e->result_available);
							}
							if (e->kind == SIMPLEX_EXECUTION_BACKEND &&
								e->result_available &&
								e->status ==
									report.singleton
										.auxiliary_status)
								backend++;
						}
						CHECK(outer &&
							(!!model == core));
						if (core)
							CHECK(backend);
					}
					simplex_solve_context_destroy(&context);
					lp_model_free(m);
				}
	return 0;
}

static int failed_dual(const struct simplex_Problem *p,
	const struct lp_simplex_Options *o, double *x, double *y,
	struct lp_simplex_Result *r, int propagate)
{
	(void)o;
	(void)x;
	(void)y;
	(void)propagate;
	r->status = lp_simplex_PrecisionError;
	r->iterations = 1;
	p->report->stop = SIMPLEX_STOP_NUMERICAL;
	return -1;
}

static int failed_hybrid(const struct simplex_Problem *p,
	const struct lp_simplex_Options *o, const double *seed, double *x,
	struct lp_simplex_Result *r)
{
	(void)seed;
	return failed_dual(p, o, x, NULL, r, 0);
}

static int controller(void)
{
	struct lp_Model *m = fixture(1, -2., -3.);
	struct simplex_Problem p;
	struct simplex_StageReport final;
	struct simplex_SolveContext context;
	struct simplex_Execution h;
	struct simplex_AutoHistory history;
	struct simplex_AutoKernels k;
	struct simplex_AutoRoutes routes;
	struct simplex_AutoBindings bindings = {NULL, NULL, 0};
	struct lp_simplex_Options o;
	struct lp_simplex_Result r;
	double x[3];
	int d, i, found;
	bindings.kernels = &k;
	bindings.routes = &routes;
	CHECK(m && !simplex_problem_from_model(&p, m));
	p.report = &final;
	k.dual = failed_dual;
	k.hybrid = failed_hybrid;
	routes.solve = simplex_auto_route_solve;
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_AUTO);
	o.iteration_limit = 100;
	for (d = 0; d < 3; d++) {
		simplex_solve_context_init(
			&context, d ? &h : NULL, &p.execution);
		if (d == 2)
			while (h.count < SIMPLEX_EXECUTION_CAPACITY)
				(void)simplex_execution_begin(&p.execution,
					SIMPLEX_EXECUTION_SOLVE, 0, 0, 0);
		CHECK(simplex_auto_solve(&p, &o, x, &r, 0, &bindings, NULL,
			      &history) != 0);
		CHECK(r.status == lp_simplex_Unboundedness &&
			r.iterations == 2 &&
			r.primal_infeasibility <= o.primal_tolerance);
		CHECK(!final.singleton.available &&
			final.stop == SIMPLEX_STOP_UNBOUNDED &&
			final.point_available);
		found = 0;
		for (i = 0; i < history.count; i++)
			if (history.stages[i].action.stage ==
				SIMPLEX_AUTO_SINGLETON) {
				CHECK(history.stages[i].kernel_status ==
					lp_simplex_Unboundedness);
				CHECK(history.stages[i]
						.report.singleton
						.auxiliary_status ==
					lp_simplex_Infeasibility);
				CHECK(history.stages[i]
						.report.singleton.mapping ==
					SIMPLEX_SINGLETON_UNBOUNDED);
				CHECK(history.stages[i].validation ==
					SIMPLEX_AUTO_TERMINAL_VALID);
				found++;
			}
		CHECK(found == 1);
		simplex_solve_context_destroy(&context);
	}
	simplex_problem_destroy(&p);
	lp_model_free(m);
	return 0;
}

static int stale_rejection(void)
{
	struct lp_Model *m = fixture(1, 0., -3.);
	struct simplex_Problem p;
	struct simplex_StageReport report;
	struct lp_simplex_Options o;
	struct lp_simplex_Result r;
	double x[3] = {17., 17., 17.};
	int state;
	CHECK(m && !simplex_problem_from_model(&p, m));
	p.report = &report;
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_AUTO);
	memset(&r, 0xff, sizeof(r));
	memset(&report, 0xff, sizeof(report));
	state = simplex_auto_route_solve(-1, &p, &o, NULL, x, &r);
	CHECK(state && r.status == lp_simplex_CondUnsatisfied &&
		!r.iterations && isnan(r.objective) &&
		isnan(r.dual_infeasibility));
	CHECK(r.dual_residual_space == lp_simplex_RESIDUAL_UNAVAILABLE &&
		report.stop == SIMPLEX_STOP_NOT_APPLICABLE);
	CHECK(!report.point_available && !report.basis_available &&
		!report.singleton.available && x[0] == 17.);
	simplex_problem_destroy(&p);
	lp_model_free(m);
	return 0;
}

/* Test only: compile the production adapter under a namespaced subprocess.
 * These checks validate rejection/mapping contracts, not stubbed mathematics. */
static int raw_status, raw_work, raw_state, raw_calls;
static double seen_lower, seen_upper;

static int singleton_test_kernel(const struct simplex_Problem *p,
	const struct lp_simplex_Options *o, double *x, double *y,
	struct lp_simplex_Result *r, int propagate)
{
	int i;
	raw_calls++;
	seen_lower = p->bounds[0].lb;
	seen_upper = p->bounds[0].ub;
	(void)o;
	(void)propagate;
	for (i = 0; i < p->columns; i++)
		x[i] = 0.;
	for (i = 0; i < p->rows; i++)
		y[i] = 0.;
	r->status = raw_status;
	r->iterations = raw_work;
	r->objective = r->dual_infeasibility = r->primal_infeasibility = 1234.;
	r->dual_residual_space = lp_simplex_RESIDUAL_AUXILIARY;
	return raw_state;
}

/* Inject only adapter-owned vectors; model/CSC and kernel allocations retain
 * their real library ownership. Poison storage to expose premature reads. */
static void *owned[5];
static int allocation_requests, allocation_failure, allocation_live, bad_free;

static void *singleton_test_malloc(size_t bytes)
{
	void *p;
	int i;
	allocation_requests++;
	if (allocation_requests == allocation_failure)
		return NULL;
	p = malloc(bytes);
	if (!p)
		return NULL;
	memset(p, 0xa5, bytes);
	for (i = 0; i < 5; i++)
		if (!owned[i]) {
			owned[i] = p;
			allocation_live++;
			return p;
		}
	free(p);
	return NULL;
}

static void singleton_test_free(void *p)
{
	int i;
	if (!p)
		return;
	for (i = 0; i < 5; i++)
		if (owned[i] == p) {
			owned[i] = NULL;
			allocation_live--;
			free(p);
			return;
		}
	bad_free++;
}

#define lp_simplex_malloc singleton_test_malloc
#define lp_simplex_free singleton_test_free
#define simplex_dual_solve_problem singleton_test_kernel
#define simplex_singleton_dual_solve_reported singleton_test_reported
#define simplex_singleton_dual_solve singleton_test_solve
#include "../src/dual/simplex_singleton_dual.c"

#undef lp_simplex_malloc
#undef lp_simplex_free
#undef simplex_dual_solve_problem
#undef simplex_singleton_dual_solve_reported
#undef simplex_singleton_dual_solve

static int subprocess_contract(void)
{
	struct lp_Model *m;
	struct lp_simplex_Options o;
	struct lp_simplex_Result r;
	struct simplex_StageReport report;
	double x[3];
	int kind, applicable, state;
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_DUAL_REVISED);
	o.iteration_limit = 10;
	for (kind = 0; kind < 9; kind++) {
		m = fixture(1, 0., -3.);
		CHECK(m);
		raw_state = -1;
		raw_work = 1;
		raw_status = lp_simplex_Unboundedness;
		if (kind == 1)
			raw_status = lp_simplex_ExceedIterLimit;
		if (kind == 2)
			raw_status = lp_simplex_MemoryAllocError;
		if (kind == 3) {
			raw_status = lp_simplex_Infeasibility;
			m->constraints[0].rhs = DBL_MAX;
			m->constraints[0].coef[0] = DBL_MIN;
			m->objective[2] = -DBL_MAX;
		}
		if (kind == 4)
			raw_status = lp_simplex_Success;
		if (kind == 5) {
			raw_status = lp_simplex_Infeasibility;
			raw_work = 11;
		}
		if (kind == 6)
			raw_status = 99;
		if (kind == 7) {
			raw_status = lp_simplex_Infeasibility;
			raw_state = 2;
		}
		if (kind == 8) {
			raw_status = lp_simplex_Success;
			raw_state = 0;
		}
		x[0] = x[1] = x[2] = NAN;
		state = singleton_test_reported(
			m, &o, x, &r, &applicable, &report, NULL);
		CHECK(state && applicable &&
			report.singleton.auxiliary_available &&
			report.singleton.auxiliary_status == raw_status);
		CHECK(!report.point_available &&
			!report.singleton.input_point_checked &&
			!report.basis_available && !report.factor_available &&
			!report.lower_bound.quality);
		CHECK(isnan(r.objective) && isnan(r.dual_infeasibility) &&
			r.dual_residual_space ==
				lp_simplex_RESIDUAL_UNAVAILABLE);
		CHECK(r.iterations == raw_work);
		if (kind == 1 || kind == 2)
			CHECK(r.status == raw_status &&
				!report.singleton.mapping);
		else
			CHECK(r.status == lp_simplex_PrecisionError &&
				report.singleton.mapping ==
					SIMPLEX_SINGLETON_REJECTED);
		if (kind >= 4 && kind < 8)
			CHECK(report.failure == SIMPLEX_FAILURE_STAGE_CONTRACT);
		lp_model_free(m);
	}
	return 0;
}

static int ratio_enclosure(void)
{
	struct lp_Model *m;
	struct lp_simplex_Options o;
	struct lp_simplex_Result r;
	struct simplex_StageReport report;
	double x[3];
	int kind, applicable;
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_DUAL_REVISED);
	raw_status = lp_simplex_ExceedIterLimit;
	raw_state = -1;
	raw_work = 1;
	for (kind = 0; kind < 3; kind++) {
		m = fixture(1, 0., -3.);
		CHECK(m);
		if (kind == 0) {
			m->constraints[0].coef[0] = 49.;
			m->constraints[0].coef[1] = -49.;
		}
		if (kind == 1) {
			m->objective[0] = m->objective[1] = DBL_MIN;
			m->constraints[0].coef[0] = DBL_MAX;
			m->constraints[0].coef[1] = -DBL_MAX;
		}
		if (kind == 2) {
			m->objective[0] = DBL_MAX;
			m->constraints[0].coef[0] = DBL_MIN;
		}
		raw_calls = 0;
		CHECK(singleton_test_reported(
			      m, &o, x, &r, &applicable, &report, NULL) != 0);
		if (kind == 2)
			CHECK(!raw_calls &&
				r.status == lp_simplex_PrecisionError &&
				!report.singleton.auxiliary_available &&
				!r.iterations);
		else {
			CHECK(raw_calls == 1 &&
				r.status == lp_simplex_ExceedIterLimit);
			if (kind == 0)
				CHECK(seen_lower < -1. / 49. &&
					seen_upper > 1. / 49.);
			else
				CHECK(seen_lower < 0. && seen_upper > 0.);
		}
		CHECK(!report.point_available && isnan(r.objective) &&
			isnan(r.dual_infeasibility));
		lp_model_free(m);
	}
	return 0;
}

static int resource_exits(void)
{
	struct lp_Model *m;
	struct lp_simplex_Options o;
	struct lp_simplex_Result r;
	struct simplex_StageReport report;
	struct simplex_SolveContext context;
	struct simplex_Execution h;
	struct simplex_ExecutionScope scope;
	double x[3];
	int core, diagnostics, failure, state, applicable, i;
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_DUAL_REVISED);
	for (core = 0; core < 2; core++)
		for (diagnostics = 0; diagnostics < 3; diagnostics++)
			for (failure = 1; failure <= 4 + core; failure++) {
				m = fixture(core, 2., -.5);
				CHECK(m && !allocation_live && !bad_free);
				simplex_solve_context_init(&context,
					diagnostics ? &h : NULL, &scope);
				if (diagnostics == 2)
					while (h.count <
						SIMPLEX_EXECUTION_CAPACITY)
						(void)simplex_execution_begin(
							&scope,
							SIMPLEX_EXECUTION_SOLVE,
							0, 0, 0);
				allocation_requests = raw_calls = 0;
				allocation_failure = failure;
				x[0] = x[1] = x[2] = 17.;
				state = singleton_test_reported(m, &o, x, &r,
					&applicable, &report, &scope);
				CHECK(state &&
					r.status ==
						lp_simplex_MemoryAllocError &&
					!r.iterations &&
					report.stop == SIMPLEX_STOP_RESOURCE &&
					!raw_calls);
				CHECK(applicable == (failure > 3) &&
					report.singleton.structure_verified ==
						applicable);
				CHECK(allocation_requests ==
						(failure < 3 ? 3 : failure) &&
					!allocation_live && !bad_free);
				CHECK(!report.point_available &&
					!report.singleton.input_point_checked &&
					!report.singleton.auxiliary_available &&
					!report.singleton.mapping &&
					!report.basis_available &&
					!report.factor_available &&
					!report.lower_bound.quality);
				CHECK(isnan(r.objective) &&
					isnan(r.dual_infeasibility));
				for (i = 0; i < 3; i++)
					CHECK(x[i] == 17.);
				if (diagnostics == 1)
					CHECK(h.count == 1 &&
						h.events[0].complete &&
						h.events[0].stop ==
							SIMPLEX_STOP_RESOURCE &&
						!h.events[0].point_source);
				simplex_solve_context_destroy(&context);
				lp_model_free(m);
			}
	allocation_failure = 0;
	return 0;
}

int main(void)
{
	return actual() || controller() || stale_rejection() ||
	       subprocess_contract() || ratio_enclosure() || resource_exits();
}
