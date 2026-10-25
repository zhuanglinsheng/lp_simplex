/**
 * @file test_auto_policy.c
 * @brief Mechanism tests: no benchmark names, reference trajectories or time
 * targets.
 */

#include "simplex_auto.h"
#include "simplex_auto_routes.h"
#include "simplex_auto_policy.h"
#include "simplex_stage.h"
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

static int reserves(void)
{
	int work, paths, reserve, previous;
	for (paths = 0; paths <= SIMPLEX_AUTO_ROUTE_COUNT + 1; paths++) {
		previous = 0;
		for (work = 0; work < 10000; work++) {
			reserve = simplex_auto_recovery_reserve(work, paths);
			CHECK(reserve >= 0 && reserve <= work / 2);
			CHECK(reserve >= previous);
			if (work / 2 >= paths * SIMPLEX_AUTO_PROBE_CAP)
				CHECK(reserve ==
					paths * SIMPLEX_AUTO_PROBE_CAP);
			previous = reserve;
		}
	}
	for (paths = 0; paths <= SIMPLEX_AUTO_ROUTE_COUNT; paths++)
		for (work = 0; work < 20000; work++) {
			reserve = simplex_auto_policy_dual_grant(work, paths);
			CHECK(reserve >= 0 && reserve <= work);
			CHECK(reserve ==
				work - simplex_auto_recovery_reserve(
					       work, paths ? paths + 1 : 0));
			CHECK(!paths || reserve >= work - work / 2);
		}
	CHECK(simplex_auto_policy_dual_grant(INT_MIN, INT_MAX) == 0);
	CHECK(simplex_auto_policy_dual_grant(INT_MAX, INT_MAX) ==
		INT_MAX - INT_MAX / 2);
	CHECK(simplex_auto_policy_dual_grant(INT_MAX, 0) == INT_MAX);
	CHECK(simplex_auto_policy_dual_grant(1, 1) == 1);
	CHECK(simplex_auto_recovery_reserve(INT_MAX, INT_MAX) == INT_MAX / 2);
	CHECK(simplex_auto_recovery_reserve(-1, 2) == 0);
	CHECK(simplex_auto_recovery_reserve(100, -1) == 0);
	CHECK(simplex_auto_policy_initializer_grant(INT_MIN, INT_MAX) == 0);
	CHECK(simplex_auto_policy_initializer_grant(1, INT_MAX) == 1);
	CHECK(simplex_auto_policy_initializer_grant(100, 2) == 50);
	CHECK(simplex_auto_policy_initializer_grant(INT_MAX, INT_MAX) ==
		SIMPLEX_AUTO_INITIALIZER_OBSERVATIONS * SIMPLEX_AUTO_PROBE_CAP);
	CHECK(simplex_auto_initializer_limit(-1) == 0);
	CHECK(simplex_auto_initializer_limit(0) == 0);
	CHECK(simplex_auto_initializer_limit(1) == 1);
	CHECK(simplex_auto_initializer_limit(INT_MAX) ==
		SIMPLEX_AUTO_INITIALIZER_OBSERVATIONS * SIMPLEX_AUTO_PROBE_CAP);
	for (work = 0; work < 20000; work++) {
		reserve = simplex_auto_initializer_limit(work);
		CHECK(reserve >= 0 && reserve <= work);
		CHECK(reserve <= SIMPLEX_AUTO_INITIALIZER_OBSERVATIONS *
					 SIMPLEX_AUTO_PROBE_CAP);
	}
	return 1;
}

static int cover_grants(void)
{
	static const int values[] = {INT_MIN, -1, 0, 1, 2, 63, 1024, INT_MAX};
	int i, j, k, grant, expected;
	for (i = 0; i < 8; i++)
		for (j = 0; j < 8; j++)
			for (k = 0; k < 8; k++) {
				expected =
					values[i] > 0 &&
							values[j] / (2 * SIMPLEX_AUTO_COVER_COLUMN_MARGIN) >=
								values[i] &&
							values[k] / (2 * SIMPLEX_AUTO_COVER_MIN_SWEEPS) >=
								values[i]
						? values[k] / 2
						: 0;
				grant = simplex_auto_policy_cover_grant(
					values[i], values[j], values[k]);
				CHECK(grant == expected && grant >= 0);
				CHECK(!grant || grant <= values[k] / 2);
			}
	return 1;
}

/* Querying policy must not consume routes or mutate the caller's runtime.
 * Exhaust every admission subset, including zero-work and one-work leases. */
static int decisions(void)
{
	int allowed[SIMPLEX_AUTO_ROUTE_COUNT], saved[SIMPLEX_AUTO_ROUTE_COUNT];
	int seen[SIMPLEX_AUTO_ROUTE_COUNT], mask, i, work, cursor, steps;
	struct simplex_AutoDecision d, repeated;
	struct simplex_StageReport report;

	for (mask = 0; mask < (1 << SIMPLEX_AUTO_ROUTE_COUNT); mask++) {
		for (i = 0; i < SIMPLEX_AUTO_ROUTE_COUNT; i++)
			allowed[i] = (mask >> i) & 1;
		memcpy(saved, allowed, sizeof(saved));
		for (work = 0; work <= 16; work++) {
			memset(seen, 0, sizeof(seen));
			cursor = -1;
			for (steps = 0; steps <= SIMPLEX_AUTO_ROUTE_COUNT;
				steps++) {
				d = simplex_auto_policy_next(
					allowed, cursor, work);
				repeated = simplex_auto_policy_next(
					allowed, cursor, work);
				CHECK(d.index == repeated.index &&
					d.route == repeated.route &&
					d.grant == repeated.grant &&
					d.reserve == repeated.reserve);
				CHECK(memcmp(saved, allowed, sizeof(saved)) ==
					0);
				CHECK(d.grant >= 0 && d.grant <= work);
				if (d.route == SIMPLEX_AUTO_DUAL_FALLBACK) {
					CHECK(d.grant == work);
					break;
				}
				CHECK(d.index > cursor &&
					d.index < SIMPLEX_AUTO_ROUTE_COUNT);
				CHECK(allowed[d.route] && !seen[d.route]);
				CHECK(d.grant > 0 && d.reserve <= work / 2);
				seen[d.route] = 1;
				cursor = d.index; /* Explicit caller commit. */
			}
			CHECK(steps <= SIMPLEX_AUTO_ROUTE_COUNT);
			if (work > 0)
				for (i = 0; i < SIMPLEX_AUTO_ROUTE_COUNT; i++)
					CHECK(seen[i] == allowed[i]);
		}
	}
	for (cursor = -1; cursor <= SIMPLEX_AUTO_ROUTE_COUNT; cursor++) {
		CHECK(simplex_auto_policy_seed_cursor(cursor, 0) == cursor);
		CHECK(simplex_auto_policy_seed_cursor(cursor, 1) ==
			(cursor < 0 ? 0 : cursor));
	}
	CHECK(simplex_auto_policy_seed_cursor(INT_MAX, 1) == INT_MAX);
	/* An initializer's accepted point bypasses singleton feasibility work. */
	for (i = 0; i < SIMPLEX_AUTO_ROUTE_COUNT; i++)
		allowed[i] = 1;
	d = simplex_auto_policy_next(allowed, 0, INT_MAX);
	CHECK(d.route == SIMPLEX_AUTO_HYBRID && d.grant > 0);
	d = simplex_auto_policy_next(
		allowed, SIMPLEX_AUTO_ROUTE_COUNT, INT_MAX);
	CHECK(d.route == SIMPLEX_AUTO_DUAL_FALLBACK && d.grant == INT_MAX);
	simplex_stage_reset(&report);
	report.stop = SIMPLEX_STOP_NUMERICAL;
	CHECK(!simplex_auto_policy_oracle_retry(SIMPLEX_AUTO_HYBRID, &report));
	report.failure = SIMPLEX_FAILURE_IMPLICIT_ORACLE;
	CHECK(simplex_auto_policy_oracle_retry(SIMPLEX_AUTO_HYBRID, &report));
	CHECK(simplex_auto_policy_oracle_retry(
		SIMPLEX_AUTO_STANDARD_IMPLICIT, &report));
	CHECK(!simplex_auto_policy_oracle_retry(
		SIMPLEX_AUTO_DUAL_EXPLORE, &report));
	return 1;
}

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
	/* Progress survives the former 16,384-pivot cutoff in one live solve. */
	memset(&h, 0, sizeof(h));
	for (i = 0; i <= 32768; i += 128)
		CHECK(!simplex_dual_explore_yield(
			&h, i, 60000, 0, 0, 100000. - i, 1e-7, 0));
	CHECK(h.limit > 32768);
	/* An upward excursion followed by return to the starting value is
	 * not improvement over the beginning of the observation lease. */
	memset(&h, 0, sizeof(h));
	for (i = 0; i < 1024; i += 128)
		CHECK(!simplex_dual_explore_yield(
			&h, i, 10000, 0, 0, i == 512 ? 200. : 100., 1e-7, 0));
	CHECK(simplex_dual_explore_yield(&h, 1024, 10000, 0, 0, 100., 1e-7, 0));
	/* Changing working-bound epochs must not erase a plateau. */
	memset(&h, 0, sizeof(h));
	for (i = 0; i < 1024; i += 128)
		CHECK(!simplex_dual_explore_yield(
			&h, i, 10000, 0, i / 128, 100., 1e-7, 0));
	CHECK(simplex_dual_explore_yield(&h, 1024, 10000, 0, 9, 100., 1e-7, 0));
	CHECK(h.reason == SIMPLEX_DUAL_EXPLORE_MEASUREMENT);
	CHECK(!simplex_dual_explore_reassess(&h, 10000, 0U, 0));
	CHECK(h.limit == 1024);
	CHECK(simplex_dual_explore_reassess(
		&h, 10000, SIMPLEX_INPUT_PRIMAL_PROGRESS, 0));
	CHECK(h.limit == 2048 && h.reason == SIMPLEX_DUAL_EXPLORE_CONTINUE);
	/* A smaller value in a changed working problem is not a reduction of
	 * the old problem. Without independent evidence it must be reassessed. */
	memset(&h, 0, sizeof(h));
	CHECK(!simplex_dual_explore_yield(&h, 0, 10000, 0, 0, 100., 1e-7, 0));
	CHECK(simplex_dual_explore_yield(&h, 1024, 10000, 0, 1, 1., 1e-7, 0));
	CHECK(h.reason == SIMPLEX_DUAL_EXPLORE_MEASUREMENT);
	CHECK(!simplex_dual_explore_reassess(
		&h, 1024, SIMPLEX_INPUT_BOUND_PROGRESS, 0));

	/* Fewer violated rows is useful feasibility progress even when their
	 * aggregate magnitude grows. A later plateau still has to yield. */
	memset(&h, 0, sizeof(h));
	for (i = 0; i <= 1024; i += 128)
		CHECK(!simplex_dual_explore_yield_pair(&h, i, 10000, 0, 0,
			100. + i, 1e-7, 0, 100. - .02 * i, .5));
	CHECK(h.limit == 2048);
	for (i = 1152; i < 2048; i += 128)
		CHECK(!simplex_dual_explore_yield_pair(
			&h, i, 10000, 0, i, 2000., 1e-7, 0, 79.52, .5));
	CHECK(simplex_dual_explore_yield_pair(
		&h, 2048, 10000, 0, 2048, 2000., 1e-7, 0, 79.52, .5));
	/* Objective progress works with negative objectives, too. */
	memset(&h, 0, sizeof(h));
	for (i = 0; i <= 1024; i += 128)
		CHECK(!simplex_dual_explore_yield_pair(&h, i, 10000, 10, 0,
			100., 1e-7, 0, -100. - .02 * i, 1e-5));
	CHECK(h.limit == 2048);
	/* A numerical Phase-I repair may get finite grace, never a fresh
	 * allowance for every epoch. Genuine progress remains uncapped. */
	memset(&h, 0, sizeof(h));
	h.recovery_until = 4096;
	for (i = 0; i < 4096; i += 128)
		CHECK(!simplex_dual_explore_yield(
			&h, i, 10000, 0, i, 100., 1e-7, 0));
	CHECK(simplex_dual_explore_yield(
		&h, 4096, 10000, 0, 4096, 100., 1e-7, 0));
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
	/* CSR admission must count only positive-cost singleton columns. Column
	 * 4 touches both rows and must not replace a missing singleton. */
	CHECK(simplex_auto_route_allowed(&p, SIMPLEX_AUTO_SINGLETON));
	p.objective[0] = 0.;
	CHECK(!simplex_auto_route_allowed(&p, SIMPLEX_AUTO_SINGLETON));
	p.objective[0] = -1.;
	CHECK(!simplex_auto_route_allowed(&p, SIMPLEX_AUTO_SINGLETON));
	p.objective[0] = 2.;
	p.matrix.value[p.matrix.column_start[0]] = -1.;
	CHECK(!simplex_auto_route_allowed(&p, SIMPLEX_AUTO_SINGLETON));
	p.matrix.value[p.matrix.column_start[0]] = 1.;
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
	if (p->report) {
		p->report->stop = simplex_stage_stop_from_status(r->status);
		p->report->failure = numeric ? SIMPLEX_FAILURE_IMPLICIT_ORACLE
					     : SIMPLEX_FAILURE_UNKNOWN;
	}
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
	struct simplex_AutoBindings bindings = {NULL, NULL, 0};
	struct simplex_AutoHistory history;
	double x[5];
	int state, i, budget;
	bindings.kernels = &k;
	bindings.routes = &routes;
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
		state = simplex_auto_solve(
			&p, &o, x, &r, 0, &bindings, NULL, &history);
		for (i = 0; i < history.count; i++) {
			const struct simplex_AutoAction *action =
				&history.stages[i].action;
			if (action->stage == SIMPLEX_AUTO_TABLEAU ||
				action->stage == SIMPLEX_AUTO_SINGLETON) {
				CHECK(action->point_source == -1);
				CHECK(action->prerequisite ==
					SIMPLEX_AUTO_UNCERTIFIED);
			}
			if (action->prerequisite == SIMPLEX_AUTO_AUDITED_POINT)
				CHECK(action->point_source >= 0 &&
					action->point_source < i);
			if (action->prerequisite ==
				SIMPLEX_AUTO_IMPLICIT_ORACLE_FAILURE) {
				CHECK(action->cause >= 0 && action->cause < i);
				CHECK(history.stages[action->cause]
						.report.failure ==
					SIMPLEX_FAILURE_IMPLICIT_ORACLE);
			}
		}
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
		state = simplex_auto_solve(
			&p, &o, x, &r, 0, &bindings, NULL, NULL);
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
	state = simplex_auto_solve(&p, &o, x, &r, 0, &bindings, NULL, NULL);
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
		state = simplex_auto_solve(
			&p, &o, x, &r, 0, &bindings, NULL, NULL);
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
	if (!reserves() || !cover_grants() || !decisions() || !leases() ||
		!adapters() || !routing())
		return 1;
	puts("AUTO leases, real adapters, route reachability, resource admission and shared budgets passed");
	return 0;
}
