/**
 * @file test_postsolve_classification.c
 * @brief Mathematical classification, candidate retention and counterexample
 * origin.
 */

#include "simplex_stage.h"
#include "simplex_solve_context.h"
#include "simplex_execution.h"
#include "simplex_original_recovery.h"
#include "simplex_model_audit.h"
#include <lp_simplex/lp_simplex.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#define CHECK(c) \
	do { \
		if (!(c)) { \
			fprintf(stderr, \
				"classification line %d case=%d: %s\n", \
				__LINE__, which, #c); \
			return 1; \
		} \
	} while (0)

struct Case {
	const char *name;
	int initial_status, used, limit;
	double saved[2];
	int raw_status, raw_work;
	double produced[2];
	int expected_status, expected_work;
	double returned[2];
	int counterexample, point_source, overflow;
};
static const struct Case cases[] = {
	{"old infeasible claim conflicts with saved point",
		lp_simplex_Infeasibility, 3, 10, {.5, .5},
		lp_simplex_Infeasibility, 1, {NAN, NAN},
		lp_simplex_PrecisionError, 4, {.5, .5}, 1, 1, 0},
	{"new infeasible claim with better feasible output",
		lp_simplex_PrecisionError, 3, 10, {.5, .5},
		lp_simplex_Infeasibility, 1, {1., 0.},
		lp_simplex_PrecisionError, 4, {1., 0.}, 1, 2, 0},
	{"no points does not erase kernel infeasible status",
		lp_simplex_PrecisionError, 0, 10, {NAN, NAN},
		lp_simplex_Infeasibility, 1, {NAN, NAN},
		lp_simplex_Infeasibility, 1, {NAN, NAN}, 0, 1, 0},
	{"new point alone contradicts infeasibility", lp_simplex_PrecisionError,
		0, 10, {NAN, NAN}, lp_simplex_Infeasibility, 1, {1., 0.},
		lp_simplex_PrecisionError, 1, {1., 0.}, 2, 2, 0},
	{"unbounded status survives restoring better saved point",
		lp_simplex_PrecisionError, 3, 10, {.5, .5},
		lp_simplex_Unboundedness, 1, {0., 1.}, lp_simplex_Unboundedness,
		4, {.5, .5}, 0, 1, 0},
	{"unbounded status retains new feasible point",
		lp_simplex_PrecisionError, 0, 10, {NAN, NAN},
		lp_simplex_Unboundedness, 1, {1., 0.}, lp_simplex_Unboundedness,
		1, {1., 0.}, 0, 2, 0},
	{"raw zero residual cannot make infeasible output a counterexample",
		lp_simplex_PrecisionError, 0, 10, {NAN, NAN},
		lp_simplex_Infeasibility, 1, {2., 2.}, lp_simplex_Infeasibility,
		1, {NAN, NAN}, 0, 1, 0},
	{"saved feasible point contradicts infeasibility despite bad output",
		lp_simplex_PrecisionError, 3, 10, {.5, .5},
		lp_simplex_Infeasibility, 1, {2., 2.},
		lp_simplex_PrecisionError, 4, {.5, .5}, 1, 1, 0},
	{"objective overflow does not erase feasibility counterexample",
		lp_simplex_PrecisionError, 0, 10, {NAN, NAN},
		lp_simplex_Infeasibility, 1, {1., 1.},
		lp_simplex_PrecisionError, 1, {NAN, NAN}, 2, 1, 1},
	{"exhausted budget still rejects contradictory old claim",
		lp_simplex_Infeasibility, 3, 3, {.5, .5}, 0, 0, {NAN, NAN},
		lp_simplex_PrecisionError, 3, {.5, .5}, 0, 0, 0},
	{"negative work rejects terminal status", lp_simplex_PrecisionError, 3,
		10, {.5, .5}, lp_simplex_Unboundedness, INT_MIN, {NAN, NAN},
		lp_simplex_PrecisionError, 10, {.5, .5}, 0, 1, 0},
	{"excessive work rejects terminal status", lp_simplex_PrecisionError, 3,
		10, {.5, .5}, lp_simplex_Unboundedness, INT_MAX, {NAN, NAN},
		lp_simplex_PrecisionError, 10, {.5, .5}, 0, 1, 0},
	{"failed retry cannot resurrect rejected old claim",
		lp_simplex_Infeasibility, 3, 10, {.5, .5},
		lp_simplex_MemoryAllocError, 1, {NAN, NAN},
		lp_simplex_PrecisionError, 4, {.5, .5}, 0, 1, 0},
	{"no feasible counterexample preserves old kernel status",
		lp_simplex_Infeasibility, 3, 3, {2., 2.}, 0, 0, {NAN, NAN},
		lp_simplex_Infeasibility, 3, {2., 2.}, 0, 0, 0}};
static int which, calls, grant;
static int reset_history, replaced;

static int invoke(const struct lp_Model *m, const struct lp_simplex_Options *o,
	double *x, struct lp_simplex_Result *r,
	const struct simplex_ProgressSink *p,
	const struct simplex_ExecutionScope *s)
{
	const struct Case *c = cases + which;
	if (reset_history && s->history) {
		int j, count = s->parent;
		memset(s->history, 0, sizeof(*s->history));
		for (j = 0; j < count; j++)
			replaced = simplex_execution_begin(s,
				SIMPLEX_EXECUTION_ORIGINAL_RETRY, o->algorithm,
				0, 0);
		s->history->events[replaced - 1].grant = 777;
	}
	(void)m;
	(void)p;
	(void)s;
	calls++;
	grant = o->iteration_limit;
	r->iterations = c->raw_work;
	r->status = c->raw_status;
	r->primal_infeasibility = 0.;
	r->objective = -999.;
	r->dual_infeasibility = .25;
	r->dual_residual_space = lp_simplex_RESIDUAL_ORIGINAL;
	x[0] = c->produced[0];
	x[1] = c->produced[1];
	return lp_simplex_EXIT_FAILURE;
}

static int same_point(const double *actual, const double *expected)
{
	int j;
	for (j = 0; j < 2; j++)
		if (!(isnan(expected[j]) ? isnan(actual[j])
					 : actual[j] == expected[j]))
			return 0;
	return 1;
}

static int run(int tracing, int algorithm)
{
	const struct Case *c = cases + which;
	struct lp_Model *m = lp_model_create(1, 2);
	struct lp_simplex_Options o;
	struct lp_simplex_Result r;
	struct simplex_SolveContext context;
	struct simplex_Execution history;
	struct simplex_ExecutionScope root;
	double x[2];
	int state, i, retried = c->point_source != 0;
	CHECK(m);
	m->constraints[0].type = optm_CONS_T_EQ;
	m->constraints[0].rhs = c->overflow ? 2. : 1.;
	m->constraints[0].coef[0] = m->constraints[0].coef[1] = 1.;
	m->objective[0] = c->overflow ? 1e308 : 1.;
	m->objective[1] = c->overflow ? 1e308 : 2.;
	lp_simplex_default_options(&o, algorithm);
	o.iteration_limit = c->limit;
	memset(&r, 0, sizeof(r));
	r.status = c->initial_status;
	r.iterations = c->used;
	x[0] = c->saved[0];
	x[1] = c->saved[1];
	simplex_solve_context_init(&context, tracing ? &history : NULL, &root);
	if (tracing == 2)
		for (i = 0; i < SIMPLEX_EXECUTION_CAPACITY; i++)
			(void)simplex_execution_begin(
				&root, SIMPLEX_EXECUTION_SOLVE, 0, 0, 0);
	calls = grant = 0;
	reset_history = tracing == 3;
	replaced = 0;
	state = simplex_original_recovery_run(
		m, &o, x, &r, lp_simplex_EXIT_FAILURE, NULL, &root, 0, invoke);
	CHECK(calls == retried && state == lp_simplex_EXIT_FAILURE);
	CHECK(r.status == c->expected_status &&
		r.iterations == c->expected_work);
	CHECK(same_point(x, c->returned));
	if (tracing == 3 && retried) {
		CHECK(replaced > 0 && !history.events[replaced - 1].complete);
		CHECK(history.events[replaced - 1].grant == 777);
		for (i = replaced; i < history.count; i++) {
			CHECK(!history.events[i].parent &&
				!history.events[i].cause);
			CHECK(!history.events[i].point_source);
		}
	}
	if (retried)
		CHECK(grant == c->limit - c->used);
	if (c->point_source == 2)
		CHECK(r.objective == 1.);
	if (c->expected_status == lp_simplex_Unboundedness &&
		c->point_source == 1) {
		CHECK(r.objective == 1.5 && r.primal_infeasibility == 0.);
		CHECK(isnan(r.dual_infeasibility) &&
			r.dual_residual_space ==
				lp_simplex_RESIDUAL_UNAVAILABLE);
	}
	if (tracing == 1) {
		CHECK(history.count ==
			1 + (retried ? 2 : 0) + (c->counterexample ? 1 : 0));
		if (c->initial_status == lp_simplex_Infeasibility &&
			c->expected_status == lp_simplex_PrecisionError) {
			CHECK(history.events[0].reason ==
				SIMPLEX_REASON_INFEASIBILITY_CONFLICT);
			CHECK(history.events[0].status ==
				lp_simplex_Infeasibility);
			CHECK(history.events[0].failure ==
				SIMPLEX_FAILURE_CERTIFICATION);
		}
		if (retried) {
			CHECK(history.events[1].point_source == 0 &&
				history.events[1].basis_source == 0);
			CHECK(history.events[history.count - 1].point_source ==
				c->point_source);
			CHECK(history.events[history.count - 1].status ==
				r.status);
			CHECK(history.events[history.count - 1].parent ==
				simplex_execution_parent(&root));
		}
		if (c->counterexample) {
			const struct simplex_ExecutionEvent *e =
				history.events + 2;
			CHECK(e->kind == SIMPLEX_EXECUTION_RESULT_AUDIT &&
				e->reason ==
					SIMPLEX_REASON_INFEASIBILITY_CONFLICT);
			CHECK(e->parent == 2 && e->cause == 2 &&
				e->point_source == c->counterexample);
			CHECK(e->model == root.model && !e->result_available &&
				e->work == 0);
			CHECK(e->primal.quality ==
					SIMPLEX_EVIDENCE_RECOMPUTED &&
				e->primal.value == 0.);
			CHECK(history.events[1].status ==
				lp_simplex_Infeasibility);
			CHECK(history.events[1].failure ==
				SIMPLEX_FAILURE_CERTIFICATION);
		}
	}
	if (tracing == 2)
		CHECK(history.dropped ==
			1 + (retried ? 2 : 0) + (c->counterexample ? 1 : 0));
	simplex_solve_context_destroy(&context);
	lp_model_free(m);
	return 0;
}

static int final_gate(void)
{
	struct lp_Model *m = lp_model_create(1, 2);
	struct lp_simplex_Options o;
	struct lp_simplex_Result r;
	double x[2] = {1., 1.};
	CHECK(m);
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_TABLEAU);
	m->constraints[0].type = optm_CONS_T_LE;
	m->constraints[0].rhs = 3.;
	m->constraints[0].coef[0] = m->constraints[0].coef[1] = 1.;
	m->objective[0] = m->objective[1] = 1e308;
	memset(&r, 0, sizeof(r));
	r.status = lp_simplex_Infeasibility;
	CHECK(simplex_finalize_result(m, &o, x, &r, lp_simplex_EXIT_FAILURE) ==
		lp_simplex_EXIT_FAILURE);
	CHECK(r.status == lp_simplex_PrecisionError &&
		r.primal_infeasibility == 0.);
	CHECK(!isfinite(
		r.objective)); /* Objective is irrelevant to feasibility. */
	x[0] = NAN;
	r.status = lp_simplex_Infeasibility;
	(void)simplex_finalize_result(m, &o, x, &r, lp_simplex_EXIT_FAILURE);
	CHECK(r.status == lp_simplex_Infeasibility);
	x[0] = 1.;
	m->constraints[0].rhs = 1.;
	r.status = lp_simplex_Infeasibility;
	(void)simplex_finalize_result(m, &o, x, &r, lp_simplex_EXIT_FAILURE);
	CHECK(r.status == lp_simplex_Infeasibility);
	lp_model_free(m);
	return 0;
}

int main(void)
{
	int tracing, algorithm;
	for (which = 0; which < (int)(sizeof(cases) / sizeof(cases[0]));
		which++)
		for (tracing = 0; tracing < 4; tracing++)
			for (algorithm = lp_simplex_ALGORITHM_DUAL_REVISED;
				algorithm <= lp_simplex_ALGORITHM_AUTO;
				algorithm++)
				if (run(tracing, algorithm))
					return 1;
	return final_gate();
}
