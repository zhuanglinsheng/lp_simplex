/**
 * @file test_postsolve.c
 * @brief Original-coordinate recovery contract with adversarial solve outputs.
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
			fprintf(stderr, "postsolve line %d mode=%d: %s\n", \
				__LINE__, mode, #c); \
			return 1; \
		} \
	} while (0)
static int mode, calls, expected_algorithm, expected_grant, adapter_error;
static const struct lp_Model *expected_model;
static struct simplex_SolveContext *expected_context;
static int original_model;
static int reset_history, replaced;

static int attempt(const struct lp_Model *model,
	const struct lp_simplex_Options *o, double *x,
	struct lp_simplex_Result *r,
	const struct simplex_ProgressSink *progress,
	const struct simplex_ExecutionScope *scope)
{
	(void)progress;
	calls++;
	if (model != expected_model || o->algorithm != expected_algorithm ||
		o->iteration_limit != expected_grant || !o->presolve ||
		scope->context != expected_context ||
		scope->model != original_model || !isnan(x[0]) || !isnan(x[1]))
		adapter_error = 1;
	if (scope->history && scope->parent) {
		const struct simplex_ExecutionEvent *e =
			scope->history->events + scope->parent - 1;
		if (e->kind != SIMPLEX_EXECUTION_ORIGINAL_RETRY ||
			e->grant != expected_grant || e->point_source ||
			e->basis_source)
			adapter_error = 1;
	}
	if (reset_history && scope->history) {
		int j, count = scope->parent;
		memset(scope->history, 0, sizeof(*scope->history));
		for (j = 0; j < count; j++)
			replaced = simplex_execution_begin(scope,
				SIMPLEX_EXECUTION_ORIGINAL_RETRY, o->algorithm,
				0, 0);
		scope->history->events[replaced - 1].grant = 777;
	}
	r->iterations = mode == 7 || mode == 11 ? 0
			: mode == 4		? INT_MIN
			: mode == 5		? INT_MAX
			: mode == 8		? 3
						: 2;
	r->status = mode == 1 || mode == 2    ? lp_simplex_PrecisionError
		    : mode == 7 || mode == 11 ? lp_simplex_MemoryAllocError
					      : lp_simplex_Success;
	r->objective = -999.; /* Must be recomputed before selection. */
	r->primal_infeasibility = 0.;
	r->dual_infeasibility = NAN;
	r->dual_residual_space = lp_simplex_RESIDUAL_UNAVAILABLE;
	if (mode != 6 && mode != 7 && mode != 11) {
		x[0] = mode == 1 ? 0. : mode == 3 ? 2. : 1.;
		x[1] = mode == 1 ? 1. : mode == 3 ? 2. : 0.;
	}
	return r->status == lp_simplex_Success ? lp_simplex_EXIT_SUCCESS
					       : lp_simplex_EXIT_FAILURE;
}

static int scenario(int tracing, int algorithm)
{
	struct lp_Model *m = lp_model_create(1, 2);
	struct lp_simplex_Options o;
	struct lp_simplex_Result r;
	struct simplex_SolveContext context;
	struct simplex_Execution history;
	struct simplex_ExecutionScope root;
	double x[2] = {.5, .5};
	int i, state, retried, kept,
		initial_work = mode == 15   ? 0
			       : mode == 8  ? INT_MAX - 2
			       : mode == 9  ? -1
			       : mode == 10 ? 11
					    : 3;
	CHECK(m);
	m->constraints[0].type = optm_CONS_T_EQ;
	m->constraints[0].rhs = 1.;
	m->constraints[0].coef[0] = m->constraints[0].coef[1] = 1.;
	m->objective[0] = 1.;
	m->objective[1] = 2.;
	lp_simplex_default_options(&o, algorithm);
	o.iteration_limit = mode == 8 ? INT_MAX : mode == 12 ? 3 : 10;
	o.presolve = 1;
	memset(&r, 0, sizeof(r));
	r.status = mode == 11 || mode == 13 ? lp_simplex_Success
					    : lp_simplex_PrecisionError;
	r.iterations = initial_work;
	r.objective = 1.5;
	r.dual_residual_space = lp_simplex_RESIDUAL_UNAVAILABLE;
	r.dual_infeasibility = NAN;
	if (mode == 11 || mode == 14)
		x[0] = x[1] = 2.;
	if (mode == 15)
		x[0] = x[1] = NAN;
	simplex_solve_context_init(&context, tracing ? &history : NULL, &root);
	if (tracing == 2)
		for (i = 0; i < SIMPLEX_EXECUTION_CAPACITY; i++)
			(void)simplex_execution_begin(
				&root, SIMPLEX_EXECUTION_SOLVE, 0, 0, 0);
	calls = adapter_error = 0;
	reset_history = tracing == 3;
	replaced = 0;
	expected_model = m;
	expected_algorithm = algorithm;
	expected_context = &context;
	original_model = root.model;
	expected_grant =
		o.iteration_limit - (initial_work >= 0 ? initial_work : 0);
	state = simplex_original_recovery_run(m, &o, x, &r,
		r.status == lp_simplex_Success ? lp_simplex_EXIT_SUCCESS
					       : lp_simplex_EXIT_FAILURE,
		NULL, &root, 0, attempt);
	retried = mode < 9 || mode == 11 || mode == 15;
	kept = mode == 0 || mode == 2 || mode == 4 || mode == 5 || mode == 8 ||
	       mode == 15;
	CHECK(calls == retried && !adapter_error);
	CHECK(r.iterations >= 0 && r.iterations <= o.iteration_limit);
	if (mode == 4 || mode == 5 || mode == 8 || mode == 9 || mode == 10)
		CHECK(r.iterations == o.iteration_limit);
	else
		CHECK(r.iterations ==
			initial_work +
				(retried ? (mode == 7 || mode == 11 ? 0 : 2)
					 : 0));
	if (kept) {
		CHECK(x[0] == 1. && x[1] == 0. && r.objective == 1.);
	} else if (mode == 11 || mode == 14) {
		CHECK(x[0] == 2. && x[1] == 2.);
	} else {
		CHECK(x[0] == .5 && x[1] == .5);
	}
	CHECK((state == lp_simplex_EXIT_SUCCESS) ==
		(mode == 0 || mode == 13 || mode == 15));
	CHECK(r.status == (state == lp_simplex_EXIT_SUCCESS
					  ? lp_simplex_Success
					  : lp_simplex_PrecisionError));
	if (tracing == 3 && retried) {
		CHECK(replaced > 0 && !history.events[replaced - 1].complete);
		CHECK(history.events[replaced - 1].grant == 777);
		CHECK(history.events[history.count - 1].cause == 0);
		CHECK(history.events[history.count - 1].point_source == 0);
		CHECK(history.events[history.count - 1].parent == 0);
	}
	if (tracing == 1) {
		CHECK(history.count == (retried ? 3 : 1));
		if (retried) {
			CHECK(history.events[1].grant == expected_grant);
			CHECK(history.events[1].cause == 1 &&
				history.events[1].model == root.model);
			CHECK(history.events[2].point_source == (kept ? 2 : 1));
			CHECK(history.events[2].work == r.iterations);
			CHECK(history.events[2].parent ==
				simplex_execution_parent(&root));
			if (mode == 4 || mode == 5 || mode == 8) {
				CHECK(history.events[1].status ==
					lp_simplex_Success);
				CHECK(history.events[1].work ==
					(mode == 4	    ? INT_MIN
						: mode == 5 ? INT_MAX
							    : 3));
				CHECK(history.events[1].failure ==
					SIMPLEX_FAILURE_STAGE_CONTRACT);
			}
			if (mode == 3 || mode == 6)
				CHECK(history.events[1].failure ==
					SIMPLEX_FAILURE_CERTIFICATION);
		}
		if (mode == 9 || mode == 10)
			CHECK(history.events[0].failure ==
				SIMPLEX_FAILURE_STAGE_CONTRACT);
		if (mode == 11)
			CHECK(history.events[0].status == lp_simplex_Success &&
				history.events[0].failure ==
					SIMPLEX_FAILURE_CERTIFICATION);
	}
	if (tracing == 2)
		CHECK(history.dropped == (retried ? 3 : 1));
	simplex_solve_context_destroy(&context);
	lp_model_free(m);
	return 0;
}

static int refused_requests;

static void *refuse_saved_point(size_t bytes)
{
	(void)bytes;
	refused_requests++;
	return NULL;
}

#define lp_simplex_malloc refuse_saved_point
#define simplex_original_recovery_reason refused_postsolve_reason
#define simplex_original_recovery_run refused_postsolve_recover
#include "../src/core/simplex_original_recovery.c"
#undef lp_simplex_malloc
#undef simplex_original_recovery_reason
#undef simplex_original_recovery_run

static int resource_refusal(void)
{
	struct lp_Model *m = lp_model_create(1, 2);
	struct lp_simplex_Options o;
	struct lp_simplex_Result r;
	struct simplex_SolveContext context;
	struct simplex_Execution h;
	struct simplex_ExecutionScope scope;
	double x[2], before[2];
	int tracing, kind, state;
	CHECK(m);
	m->constraints[0].type = optm_CONS_T_EQ;
	m->constraints[0].coef[0] = m->constraints[0].coef[1] = 1.;
	m->constraints[0].rhs = 1.;
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_AUTO);
	o.iteration_limit = 10;
	for (tracing = 0; tracing < 3; tracing++)
		for (kind = 0; kind < 4; kind++) {
			simplex_solve_context_init(
				&context, tracing ? &h : NULL, &scope);
			if (tracing == 2)
				while (h.count < SIMPLEX_EXECUTION_CAPACITY)
					(void)simplex_execution_begin(&scope,
						SIMPLEX_EXECUTION_SOLVE, 0, 0,
						0);
			memset(&r, 0, sizeof(r));
			r.status = kind == 1   ? lp_simplex_Success
				   : kind == 2 ? lp_simplex_Infeasibility
					       : lp_simplex_PrecisionError;
			r.iterations = kind == 3 ? 10 : 3;
			x[0] = x[1] = kind == 1 ? 2. : .5;
			memcpy(before, x, sizeof(x));
			calls = refused_requests = 0;
			state = refused_postsolve_recover(m, &o, x, &r,
				kind == 1 ? lp_simplex_EXIT_SUCCESS
					  : lp_simplex_EXIT_FAILURE,
				NULL, &scope, NULL, attempt);
			CHECK(state && r.status == lp_simplex_PrecisionError &&
				r.iterations == (kind == 3 ? 10 : 3));
			CHECK(!calls &&
				refused_requests == (kind == 3 ? 0 : 1) &&
				!memcmp(before, x, sizeof(x)));
			if (tracing == 1)
				CHECK(h.count == 1 && h.events[0].complete &&
					h.events[0].kind ==
						SIMPLEX_EXECUTION_POSTSOLVE &&
					h.events[0].status ==
						(kind == 1 ? lp_simplex_Success
							: kind == 2
								? lp_simplex_Infeasibility
								: lp_simplex_PrecisionError));
			simplex_solve_context_destroy(&context);
		}
	lp_model_free(m);
	return 0;
}

int main(void)
{
	int tracing, algorithm;
	for (mode = 0; mode < 16; mode++)
		for (tracing = 0; tracing < 4; tracing++)
			for (algorithm = lp_simplex_ALGORITHM_DUAL_REVISED;
				algorithm <= lp_simplex_ALGORITHM_AUTO;
				algorithm++)
				if (scenario(tracing, algorithm))
					return 1;
	return resource_refusal();
}
