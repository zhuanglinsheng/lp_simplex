/* Mechanism tests: no benchmark names, reference trajectories or time targets. */
#include "simplex_auto.h"
#include "simplex_degeneracy.h"
#include <lp_simplex/lp_simplex.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#define CHECK(c) \
	do { \
		if (!(c)) { \
			printf("AUTO policy line %d: %s\n", __LINE__, #c); \
			return 0; \
		} \
	} while (0)

static int leases(void)
{
	struct simplex_DualExplore h;
	int i;
	memset(&h, 0, sizeof(h));
	/* A decreasing feasibility residual renews the current live solve. */
	for (i = 0; i <= 1024; i += 128)
		CHECK(!simplex_dual_explore_yield(
			&h, i, 16384, 0, 0, 100. - .05 * i, 1e-7, 0));
	CHECK(h.limit == 2048);
	/* The same stage on a sustained plateau must not keep renewing. */
	for (i = 1152; i < 2048; i += 128)
		CHECK(!simplex_dual_explore_yield(
			&h, i, 16384, 0, 0, 48.8, 1e-7, 0));
	CHECK(simplex_dual_explore_yield(&h, 2048, 16384, 0, 0, 48.8, 1e-7, 0));
	/* A new measurement space gets a bounded baseline, not a comparison
	 * of objective residuals against the old feasibility residual. */
	CHECK(!simplex_dual_explore_yield(
		&h, 2048, 16384, 10, 0, 500., 1e-7, 0));
	CHECK(h.limit == 4096);
	CHECK(simplex_dual_explore_yield(
		&h, 4096, 16384, 10, 0, 400., 1e-7, 2));
	memset(&h, 0, sizeof(h));
	for (i = 0; i <= 1024; i += 128)
		CHECK(!simplex_dual_explore_yield(
			&h, i, 1024, 0, 0, 1., 1e-7, 0));
	CHECK(h.limit <= 1024);
	CHECK(simplex_dual_explore_yield(&h, 1, 1024, 0, 0, HUGE_VAL, 1e-7, 0));
	/* Reduction inside aggregate tolerance noise cannot buy another lease. */
	memset(&h, 0, sizeof(h));
	for (i = 0; i < 1024; i += 128)
		CHECK(!simplex_dual_explore_yield(
			&h, i, 16384, 0, 0, 100., 3., 0));
	CHECK(simplex_dual_explore_yield(&h, 1024, 16384, 0, 0, 98., 3., 0));
	/* Oscillations may supply local recovery evidence but never exceed the
	 * hard cap. Every grant is finite and the cap need not be a power of 2. */
	memset(&h, 0, sizeof(h));
	for (i = 0; i <= 5000; i += 8) {
		(void)simplex_dual_explore_yield(
			&h, i, 5000, i % 3, i / 128, 1., 1e-7, 0);
		CHECK(h.limit <= 5000);
	}
	return 1;
}

static struct lp_Model *problem(void)
{
	struct lp_Model *m = lp_model_create(2, 5);
	int i;
	if (!m)
		return NULL;
	for (i = 0; i < 2; i++) {
		m->constraints[i].type = optm_CONS_T_EQ;
		m->constraints[i].rhs = 1.;
		m->constraints[i].coef[2 * i] = 1.;
		m->constraints[i].coef[2 * i + 1] = -1.;
		m->objective[2 * i] = m->objective[2 * i + 1] = 2.;
		m->constraints[i].coef[4] = 1.;
	}
	m->objective[4] = 1.;
	return m;
}

static int adapters(void)
{
	struct lp_Model *m = problem();
	struct simplex_Problem p, huge;
	struct lp_simplex_Options o;
	struct lp_simplex_Result r;
	double x[5], seed[5] = {0., 0., 0., 0., 1.};
	int route, warm;
	CHECK(m && simplex_problem_from_model(&p, m) == 0);
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_AUTO);
	o.iteration_limit = 1000;
	for (route = 0; route < SIMPLEX_AUTO_ROUTE_COUNT; route++) {
		CHECK(simplex_auto_route_allowed(&p, route));
		for (warm = 0; warm <= 1; warm++) {
			memset(&r, 0, sizeof(r));
			CHECK(simplex_auto_route_solve(route, &p, &o,
				      warm ? seed : NULL, x, &r) == 0);
			CHECK(r.status == lp_simplex_Success &&
				r.iterations <= o.iteration_limit);
			CHECK(simplex_problem_primal_infeasibility(&p, x) <=
				o.primal_tolerance);
			CHECK(fabs(x[4] - 1.) < 1e-7);
			if (route == SIMPLEX_AUTO_TABLEAU)
				CHECK(isnan(r.dual_infeasibility) &&
					r.dual_residual_space ==
						lp_simplex_RESIDUAL_UNAVAILABLE);
		}
	}
	huge = p;
	huge.rows = INT_MAX - 1;
	for (route = 0; route < SIMPLEX_AUTO_ROUTE_COUNT; route++)
		CHECK(!simplex_auto_route_allowed(&huge, route));
	huge = p;
	huge.rows = 100000;
	CHECK(!simplex_auto_route_allowed(&huge, SIMPLEX_AUTO_TABLEAU));
	CHECK(!simplex_auto_route_allowed(&huge, SIMPLEX_AUTO_SINGLETON));
	simplex_problem_destroy(&p);
	lp_model_free(m);
	return 1;
}

static int target, calls[SIMPLEX_AUTO_ROUTE_COUNT], dual_calls, spent, numeric;
static int real_target, bad_diagnostic;
static int routed_mock(int route, const struct simplex_Problem *p,
	const struct lp_simplex_Options *o, const double *seed, double *x,
	struct lp_simplex_Result *r)
{
	int j, state;
	calls[route]++;
	if (real_target && route == target) {
		state = simplex_auto_route_solve(route, p, o, seed, x, r);
		spent += r->iterations;
		return state;
	}
	r->iterations = o->iteration_limit > 2 ? 2 : o->iteration_limit;
	spent += r->iterations;
	for (j = 0; j < p->columns; j++)
		x[j] = j == 4 ? 1. : 0.;
	r->status = route == target ? lp_simplex_Success
				    : (numeric ? lp_simplex_PrecisionError
					       : lp_simplex_ExceedIterLimit);
	r->dual_infeasibility = 0.;
	if (bad_diagnostic && route == SIMPLEX_AUTO_TABLEAU) {
		r->status = lp_simplex_Success;
		r->dual_infeasibility = NAN;
		r->dual_residual_space = lp_simplex_RESIDUAL_ORIGINAL;
	}
	return r->status == lp_simplex_Success ? 0 : -1;
}
static int dual_mock(const struct simplex_Problem *p,
	const struct lp_simplex_Options *o, double *x, double *y,
	struct lp_simplex_Result *r, int propagate)
{
	int j;
	(void)y;
	(void)propagate;
	dual_calls++;
	for (j = 0; j < p->columns; j++)
		x[j] = HUGE_VAL;
	r->iterations = o->iteration_limit > 2 ? 2 : o->iteration_limit;
	spent += r->iterations;
	r->status = lp_simplex_PrecisionError;
	return -1;
}
static int hybrid_mock(const struct simplex_Problem *p,
	const struct lp_simplex_Options *o, const double *seed, double *x,
	struct lp_simplex_Result *r)
{
	return routed_mock(SIMPLEX_AUTO_HYBRID, p, o, seed, x, r);
}
static int routing(void)
{
	struct lp_Model *m = problem();
	struct simplex_Problem p;
	struct lp_simplex_Options o;
	struct lp_simplex_Result r;
	struct simplex_AutoKernels k;
	struct simplex_AutoRoutes routes;
	double x[5];
	int state, i, budget;
	CHECK(m && simplex_problem_from_model(&p, m) == 0);
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_AUTO);
	k.dual = dual_mock;
	k.hybrid = hybrid_mock;
	routes.solve = routed_mock;
	numeric = 1;
	for (target = 0; target < SIMPLEX_AUTO_ROUTE_COUNT; target++) {
		memset(calls, 0, sizeof(calls));
		dual_calls = spent = 0;
		o.iteration_limit = 160;
		state = simplex_auto_run_routes(&p, &o, x, &r, 0, &k, &routes);
		if (state != 0 || r.status != lp_simplex_Success ||
			calls[target] != 1)
			printf("target=%d real=%d status=%d state=%d calls=%d work=%d\n",
				target, real_target, r.status, state,
				calls[target], r.iterations);
		CHECK(state == 0 && r.status == lp_simplex_Success &&
			calls[target] == 1);
		CHECK(r.iterations == spent && spent <= o.iteration_limit);
		for (i = 0; i < SIMPLEX_AUTO_ROUTE_COUNT; i++)
			CHECK(calls[i] <= 1);
	}
	/* Drive every actual adapter through the dispatcher as well: this
	 * includes successful tableau with an honestly unavailable dual. */
	real_target = 1;
	for (target = 0; target < SIMPLEX_AUTO_ROUTE_COUNT; target++) {
		memset(calls, 0, sizeof(calls));
		dual_calls = spent = 0;
		o.iteration_limit = 10000;
		state = simplex_auto_run_routes(&p, &o, x, &r, 0, &k, &routes);
		if (state != 0 || r.status != lp_simplex_Success ||
			calls[target] != 1)
			printf("target=%d real=%d status=%d state=%d calls=%d work=%d\n",
				target, real_target, r.status, state,
				calls[target], r.iterations);
		CHECK(state == 0 && r.status == lp_simplex_Success &&
			calls[target] == 1);
		CHECK(r.iterations == spent && spent <= o.iteration_limit);
		if (target == SIMPLEX_AUTO_TABLEAU)
			CHECK(isnan(r.dual_infeasibility) &&
				r.dual_residual_space ==
					lp_simplex_RESIDUAL_UNAVAILABLE);
	}
	real_target = 0;
	bad_diagnostic = 1;
	target = -1;
	memset(calls, 0, sizeof(calls));
	dual_calls = spent = 0;
	state = simplex_auto_run_routes(&p, &o, x, &r, 0, &k, &routes);
	CHECK(state != 0 && r.status != lp_simplex_Success);
	CHECK(calls[SIMPLEX_AUTO_TABLEAU] == 1 && r.iterations == spent);
	bad_diagnostic = 0;
	/* Exhaustion alone is not a representation failure: skip explicit QR.
	 * Small budgets may end anywhere but must always conserve work. */
	target = -1;
	numeric = 0;
	for (budget = 1; budget <= 64; budget++) {
		memset(calls, 0, sizeof(calls));
		dual_calls = spent = 0;
		o.iteration_limit = budget;
		state = simplex_auto_run_routes(&p, &o, x, &r, 0, &k, &routes);
		CHECK(state != 0 && r.iterations == spent && spent <= budget);
		CHECK(calls[SIMPLEX_AUTO_EXPLICIT] == 0 &&
			calls[SIMPLEX_AUTO_STANDARD_EXPLICIT] == 0);
		CHECK(dual_calls <= 2);
		for (i = 0; i < SIMPLEX_AUTO_ROUTE_COUNT; i++)
			CHECK(calls[i] <= 1);
	}
	simplex_problem_destroy(&p);
	lp_model_free(m);
	return 1;
}
int main(void)
{
	if (!leases() || !adapters() || !routing())
		return 1;
	puts("AUTO leases, real adapters, route reachability, resource admission and shared budgets passed");
	return 0;
}
