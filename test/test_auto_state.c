/**
 * @file test_auto_state.c
 * @brief Runtime, candidate identity and diagnostic history are independent
 * contracts.
 */

#include "simplex_stage.h"
#include "simplex_solve_context.h"
#include "simplex_execution.h"
#include "simplex_auto.h"
#include "simplex_auto_policy.h"
#include <lp_simplex/lp_simplex.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#define CHECK(c) \
	do { \
		if (!(c)) { \
			fprintf(stderr, "AUTO state line %d: %s\n", __LINE__, \
				#c); \
			return 1; \
		} \
	} while (0)

static int triggers(void)
{
	struct simplex_AutoPolicyEvidence e;
	struct simplex_AutoTrigger t;

	memset(&e, 0, sizeof(e));
	e.dual_yielded = e.have_point = e.have_basis = 1;
	e.native_oracle_failure = e.standard_oracle_failure = 1;
	e.dual.record = 0;
	e.dual.event = 10;
	e.point.record = 1;
	e.point.event = 20;
	e.basis.record = 2;
	e.basis.event = 30;
	e.native_oracle.record = 3;
	e.native_oracle.event = 40;
	e.standard_oracle.record = 5;
	e.standard_oracle.event = 60;
	e.retryable.record = 8;
	e.retryable.event = 90;
	t = simplex_auto_policy_trigger(SIMPLEX_AUTO_HYBRID, &e);
	CHECK(t.reason == SIMPLEX_REASON_DUAL_HANDOFF && t.origin.event == 10);
	t = simplex_auto_policy_trigger(SIMPLEX_AUTO_EXPLICIT, &e);
	CHECK(t.reason == SIMPLEX_REASON_IMPLICIT_ORACLE_FAILURE &&
		t.origin.event == 40);
	t = simplex_auto_policy_trigger(SIMPLEX_AUTO_STANDARD_EXPLICIT, &e);
	CHECK(t.reason == SIMPLEX_REASON_IMPLICIT_ORACLE_FAILURE &&
		t.origin.event == 60);
	t = simplex_auto_policy_trigger(SIMPLEX_AUTO_DUAL_FALLBACK, &e);
	CHECK(t.reason == SIMPLEX_REASON_BASIS_REBUILD && t.origin.record == 2);
	e.have_basis = 0;
	t = simplex_auto_policy_trigger(SIMPLEX_AUTO_DUAL_FALLBACK, &e);
	CHECK(t.reason == SIMPLEX_REASON_RECOVERY_EXHAUSTED &&
		t.origin.record == 8);
	t = simplex_auto_policy_trigger(SIMPLEX_AUTO_TABLEAU, &e);
	CHECK(t.reason == SIMPLEX_REASON_ALTERNATIVE_ALGORITHM);
	CHECK(simplex_auto_policy_purpose(SIMPLEX_AUTO_TABLEAU, &e) ==
		SIMPLEX_AUTO_ALTERNATIVE_ALGORITHM);
	/* Changing direction within the same representation is not another
	 * representation change. A later unrelated output is not its cause. */
	CHECK(simplex_auto_policy_purpose(SIMPLEX_AUTO_STANDARD_IMPLICIT, &e) ==
		SIMPLEX_AUTO_REFORMULATE);
	e.standard_basis_uncertified = 1;
	e.standard_basis.record = 4;
	e.standard_basis.event = 50;
	t = simplex_auto_policy_trigger(SIMPLEX_AUTO_STANDARD_IMPLICIT, &e);
	CHECK(t.reason == SIMPLEX_REASON_ALTERNATIVE_DIRECTION &&
		t.origin.record == 4 && t.origin.event == 50);
	CHECK(simplex_auto_policy_purpose(SIMPLEX_AUTO_STANDARD_IMPLICIT, &e) ==
		SIMPLEX_AUTO_REPLACE_DIRECTION);
	e.standard_basis_contract_failure = 1;
	t = simplex_auto_policy_trigger(SIMPLEX_AUTO_STANDARD_IMPLICIT, &e);
	CHECK(t.reason == SIMPLEX_REASON_STAGE_CONTRACT &&
		t.origin.event == 50);
	e.standard_basis_contract_failure = 0;
	e.standard_basis.record = -1;
	e.standard_basis.event = 0;
	t = simplex_auto_policy_trigger(SIMPLEX_AUTO_STANDARD_IMPLICIT, &e);
	CHECK(t.reason == SIMPLEX_REASON_ALTERNATIVE_DIRECTION &&
		!t.origin.event);
	e.native_oracle_failure = 0;
	t = simplex_auto_policy_trigger(SIMPLEX_AUTO_EXPLICIT, &e);
	CHECK(t.reason == SIMPLEX_REASON_NONE && t.origin.record == -1);
	/* Truncated diagnostic sources do not erase real enabling evidence. */
	e.have_basis = 1;
	e.basis.record = -1;
	e.basis.event = 0;
	t = simplex_auto_policy_trigger(SIMPLEX_AUTO_DUAL_FALLBACK, &e);
	CHECK(t.reason == SIMPLEX_REASON_BASIS_REBUILD && t.origin.event == 0);
	return 0;
}

static int mode, dual_calls, seed_ok;
static struct simplex_AutoHistory *output;

static int dual(const struct simplex_Problem *p,
	const struct lp_simplex_Options *o, double *x, double *y,
	struct lp_simplex_Result *r, int propagate)
{
	(void)o;
	(void)y;
	(void)propagate;
	dual_calls++;
	x[0] = x[1] = 0.;
	r->iterations = 1;
	r->status = lp_simplex_ExceedIterLimit;
	p->report->stop = SIMPLEX_STOP_BUDGET;
	return -1;
}

static int initialize(void *context, const struct simplex_Problem *input,
	const struct lp_simplex_Options *o, double *x, int *work)
{
	(void)input;
	(void)context;
	(void)o;
	x[0] = x[1] = .5;
	*work = 1;
	/* Failure injection: poison diagnostic summaries while the controller
	 * is producing a new candidate. These fields must not govern handoff. */
	output->have_point = 1;
	output->point_source = INT_MAX;
	output->best_objective = -HUGE_VAL;
	return 0;
}

static int hybrid(const struct simplex_Problem *p,
	const struct lp_simplex_Options *o, const double *seed, double *x,
	struct lp_simplex_Result *r)
{
	(void)o;
	seed_ok = mode == 3 ? !seed : seed && seed[0] == .5 && seed[1] == .5;
	x[0] = .75;
	x[1] = .25;
	r->iterations = 1;
	r->status = lp_simplex_PrecisionError;
	if (mode == 1 || mode == 2) {
		/* Deliberately violate const input to exercise the handoff guard.
		 * Shape and feasible coordinates are unchanged; identity is not. */
		struct simplex_Problem *changed = (struct simplex_Problem *)p;
		if (mode == 1)
			changed->execution.model++;
		else
			simplex_problem_mark_changed(changed);
		r->status = lp_simplex_Success;
		r->dual_infeasibility = 0.;
		return 0;
	}
	if (mode == 3) {
		r->status = lp_simplex_Success;
		r->dual_infeasibility = 0.;
		return 0;
	}
	p->report->stop = SIMPLEX_STOP_NUMERICAL;
	p->report->failure = SIMPLEX_FAILURE_DIRECTION;
	return -1;
}

static int runtime(void)
{
	struct lp_Model *m = lp_model_create(1, 2);
	struct simplex_Problem p;
	struct simplex_SolveContext context;
	struct simplex_Execution execution;
	struct simplex_AutoHistory journal;
	struct simplex_AutoKernels kernels;
	struct simplex_AutoBindings bindings = {NULL, NULL, 0};
	struct simplex_AutoInitializer initializer;
	struct lp_simplex_Options o;
	struct lp_simplex_Result r;
	double x[2];
	int diagnostics, i, state;

	bindings.kernels = &kernels;
	CHECK(m);
	m->constraints[0].type = optm_CONS_T_EQ;
	m->constraints[0].rhs = 1.;
	m->constraints[0].coef[0] = m->constraints[0].coef[1] = 1.;
	m->objective[1] = 1.;
	CHECK(!simplex_problem_from_model(&p, m));
	kernels.dual = dual;
	kernels.hybrid = hybrid;
	initializer.context = NULL;
	initializer.run = initialize;
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_AUTO);
	o.iteration_limit = 160;
	output = &journal;
	for (diagnostics = 0; diagnostics < 3; diagnostics++) {
		for (mode = 0; mode < 4; mode++) {
			simplex_solve_context_init(
				&context, &execution, &p.execution);
			if (diagnostics == 0)
				p.execution.history = NULL;
			if (diagnostics == 2)
				for (i = 0; i < SIMPLEX_EXECUTION_CAPACITY; i++)
					(void)simplex_execution_begin(
						&p.execution,
						SIMPLEX_EXECUTION_SOLVE, 0, 0,
						0);
			p.revision = mode == 3 ? ULONG_MAX : 1;
			dual_calls = seed_ok = 0;
			state = simplex_auto_solve(&p, &o, x, &r, 0, &bindings,
				&initializer, &journal);
			CHECK(seed_ok);
			CHECK(journal.stages[2].action.reason ==
				SIMPLEX_REASON_DUAL_HANDOFF);
			CHECK(journal.stages[2].action.cause == 0);
			CHECK(journal.stages[2].action.point_source ==
				(mode == 3 ? -1 : 1));
			if (mode == 0) {
				CHECK(state != 0 && dual_calls == 2 &&
					r.iterations == 4);
				CHECK(x[0] == .75 && x[1] == .25 &&
					r.objective == .25);
				CHECK(journal.have_point &&
					journal.point_source == 2);
				CHECK(journal.stages[3].action.reason ==
					SIMPLEX_REASON_RECOVERY_EXHAUSTED);
			} else if (mode != 3) {
				CHECK(state != 0 && dual_calls == 1 &&
					r.iterations == 3);
				CHECK(r.status == lp_simplex_PrecisionError &&
					!journal.have_point);
				CHECK(journal.stages[2].validation ==
					SIMPLEX_AUTO_TERMINAL_INVALID);
				CHECK(!journal.stages[2].accepted_point);
			} else {
				CHECK(state == 0 &&
					r.status == lp_simplex_Success &&
					r.iterations == 3);
				CHECK(!journal.have_point &&
					!journal.stages[1].accepted_point);
				CHECK(journal.stages[2].validation ==
					SIMPLEX_AUTO_TERMINAL_VALID);
			}
			if (diagnostics == 1) {
				CHECK(execution.events[2].reason ==
					SIMPLEX_REASON_DUAL_HANDOFF);
				CHECK(execution.events[2].cause == 1);
				CHECK(execution.events[2].point_source ==
					(mode == 3 ? 0 : 2));
			} else
				CHECK(!journal.stages[2].action.cause_event);
			if (diagnostics == 2)
				CHECK(execution.dropped == journal.count + 1);
			simplex_solve_context_destroy(&context);
		}
	}
	simplex_problem_destroy(&p);
	lp_model_free(m);
	return 0;
}

static int basis_contract;
static int trace_mutation, fake_id;
static struct simplex_Execution replacement;
static struct simplex_Execution *fake_history;

static int direction_route(int route, const struct simplex_Problem *p,
	const struct lp_simplex_Options *o, const double *seed, double *x,
	struct lp_simplex_Result *r)
{
	(void)o;
	if (route == SIMPLEX_AUTO_CLASSIC && trace_mutation &&
		p->execution.history) {
		struct simplex_ExecutionScope scope = p->execution;
		int i, count = scope.parent;
		if (trace_mutation == 2) {
			scope.history = &replacement;
			((struct simplex_Problem *)p)->execution.history =
				&replacement;
		}
		memset(scope.history, 0, sizeof(*scope.history));
		for (i = 0; i < count; i++)
			fake_id = simplex_execution_begin(&scope,
				SIMPLEX_EXECUTION_AUTO, SIMPLEX_AUTO_CLASSIC, 0,
				0);
		fake_history = scope.history;
		fake_history->events[fake_id - 1].grant = 777;
	}
	if (route == SIMPLEX_AUTO_STANDARD_IMPLICIT &&
		(!seed || seed[0] != .75 || seed[1] != .25))
		seed_ok = 0;
	x[0] = .75;
	x[1] = .25;
	r->status = lp_simplex_PrecisionError;
	r->iterations =
		route == SIMPLEX_AUTO_CLASSIC && basis_contract ? -1 : 1;
	p->report->stop = SIMPLEX_STOP_NUMERICAL;
	return -1;
}

static int direction_hybrid(const struct simplex_Problem *p,
	const struct lp_simplex_Options *o, const double *seed, double *x,
	struct lp_simplex_Result *r)
{
	return direction_route(SIMPLEX_AUTO_HYBRID, p, o, seed, x, r);
}

static int direction_handoff(void)
{
	struct lp_Model *m = lp_model_create(1, 2);
	struct simplex_Problem p;
	struct simplex_SolveContext context;
	struct simplex_Execution execution;
	struct simplex_AutoHistory history;
	struct simplex_AutoKernels kernels;
	struct simplex_AutoRoutes routes;
	struct simplex_AutoBindings bindings = {NULL, NULL, 0};
	struct lp_simplex_Options o;
	struct lp_simplex_Result r;
	double x[2];
	int diagnostics, i, classic, cone, point_source;
	int expected_work[2] = {0, 0};

	bindings.kernels = &kernels;
	bindings.routes = &routes;
	CHECK(m);
	m->constraints[0].type = optm_CONS_T_EQ;
	m->constraints[0].rhs = 1.;
	m->constraints[0].coef[0] = m->constraints[0].coef[1] = 1.;
	m->objective[1] = 1.;
	CHECK(!simplex_problem_from_model(&p, m));
	kernels.dual = dual;
	kernels.hybrid = direction_hybrid;
	routes.solve = direction_route;
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_AUTO);
	o.iteration_limit = 1000;
	for (diagnostics = 0; diagnostics < 5; diagnostics++)
		for (basis_contract = 0; basis_contract < 2; basis_contract++) {
			simplex_solve_context_init(
				&context, &execution, &p.execution);
			if (diagnostics != 0 && diagnostics != 2)
				simplex_execution_set_parent(&p.execution,
					simplex_execution_begin(&p.execution,
						SIMPLEX_EXECUTION_SOLVE, 0, 0,
						0));
			if (diagnostics == 0)
				p.execution.history = NULL;
			if (diagnostics == 2)
				for (i = 0; i < SIMPLEX_EXECUTION_CAPACITY; i++)
					(void)simplex_execution_begin(
						&p.execution,
						SIMPLEX_EXECUTION_SOLVE, 0, 0,
						0);
			seed_ok = 1;
			trace_mutation = diagnostics == 3   ? 1
					 : diagnostics == 4 ? 2
							    : 0;
			(void)simplex_auto_solve(
				&p, &o, x, &r, 0, &bindings, NULL, &history);
			classic = cone = point_source = -1;
			for (i = 0; i < history.count; i++) {
				if (history.stages[i].action.stage ==
					SIMPLEX_AUTO_CLASSIC)
					classic = i;
				if (history.stages[i].action.stage ==
					SIMPLEX_AUTO_HYBRID)
					point_source = i;
				if (history.stages[i].action.stage ==
					SIMPLEX_AUTO_STANDARD_IMPLICIT)
					cone = i;
			}
			CHECK(classic >= 0 && cone > classic && seed_ok);
			CHECK(history.stages[cone].action.purpose ==
				SIMPLEX_AUTO_REPLACE_DIRECTION);
			CHECK(history.stages[cone].action.reason ==
				(basis_contract ? SIMPLEX_REASON_STAGE_CONTRACT
						: SIMPLEX_REASON_ALTERNATIVE_DIRECTION));
			CHECK(history.stages[cone].action.cause == classic);
			CHECK(point_source >= 0 && point_source != classic);
			CHECK(history.stages[cone].action.point_source ==
				point_source);
			CHECK(r.iterations > 0 &&
				r.iterations <= o.iteration_limit);
			if (diagnostics == 0)
				expected_work[basis_contract] = r.iterations;
			CHECK(r.iterations == expected_work[basis_contract]);
			CHECK(x[0] == .75 && x[1] == .25 && r.objective == .25);
			if (trace_mutation) {
				CHECK(fake_id > 0 &&
					!fake_history->events[fake_id - 1]
						.complete);
				CHECK(fake_history->events[fake_id - 1].grant ==
					777);
				CHECK(!history.stages[classic]
						.action.execution_id);
				CHECK(!history.stages[cone].action.cause_event);
				CHECK(history.stages[cone].action.point_event ==
					(trace_mutation == 1
							? 0
							: history.stages[point_source]
								  .action
								  .execution_id));
				CHECK(execution.events[history.stages[cone]
							       .action
							       .execution_id -
						       1]
						.parent ==
					(trace_mutation == 1 ? 0 : 1));
			}
			if (diagnostics == 1) {
				int id = history.stages[cone]
						 .action.execution_id;
				CHECK(id > 0);
				CHECK(execution.events[id - 1].cause ==
					history.stages[classic]
						.action.execution_id);
				CHECK(execution.events[id - 1].purpose ==
					SIMPLEX_AUTO_REPLACE_DIRECTION);
			} else
				CHECK(!history.stages[cone].action.cause_event);
			simplex_solve_context_destroy(&context);
		}
	simplex_problem_destroy(&p);
	lp_model_free(m);
	return 0;
}

int main(void)
{
	if (triggers() || runtime() || direction_handoff())
		return 1;
	puts("AUTO runtime, candidate identity, journal independence and semantic causes passed");
	return 0;
}
