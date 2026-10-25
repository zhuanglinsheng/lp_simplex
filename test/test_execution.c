/**
 * @file test_execution.c
 * @brief Solve-wide causal history and original-model recovery admission.
 */

#include "simplex_stage.h"
#include "simplex_solve_context.h"
#include "simplex_execution.h"
#include "simplex_original_recovery.h"
#include "simplex_auto.h"
#include "simplex_scaling.h"
#include <lp_simplex/lp_simplex.h>
#include <math.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>

#define CHECK(c) \
	do { \
		if (!(c)) { \
			fprintf(stderr, "execution line %d: %s\n", __LINE__, \
				#c); \
			return 1; \
		} \
	} while (0)

static int admission(void)
{
	CHECK(simplex_original_recovery_reason(lp_simplex_Success, 10, 1.,
		      1e-7) == SIMPLEX_REASON_INVALID_POSTSOLVE);
	CHECK(simplex_original_recovery_reason(lp_simplex_Success, 10, NAN,
		      1e-7) == SIMPLEX_REASON_INVALID_POSTSOLVE);
	CHECK(simplex_original_recovery_reason(lp_simplex_PrecisionError, 0,
		      NAN, 1e-7) == SIMPLEX_REASON_STARTUP_FAILURE);
	CHECK(simplex_original_recovery_reason(lp_simplex_PrecisionError, 10,
		      0., 1e-7) == SIMPLEX_REASON_FEASIBLE_CANDIDATE);
	/* Missing, stale (not supplied), negative and nonfinite reduced residuals
	 * cannot authorize original-feasibility recovery. */
	CHECK(!simplex_original_recovery_reason(
		lp_simplex_PrecisionError, 10, NAN, 1e-7));
	CHECK(!simplex_original_recovery_reason(
		lp_simplex_PrecisionError, 10, HUGE_VAL, 1e-7));
	CHECK(!simplex_original_recovery_reason(
		lp_simplex_PrecisionError, 10, -1., 1e-7));
	CHECK(!simplex_original_recovery_reason(
		lp_simplex_ExceedIterLimit, 10, 0., 1e-7));
	CHECK(simplex_original_recovery_reason(lp_simplex_Infeasibility, 0, 0.,
		      1e-7) == SIMPLEX_REASON_INFEASIBILITY_CONFLICT);
	CHECK(!simplex_original_recovery_reason(
		lp_simplex_Infeasibility, 0, NAN, 1e-7));
	CHECK(!simplex_original_recovery_reason(
		lp_simplex_Infeasibility, 0, -1., 1e-7));
	return 0;
}

static int identities(void)
{
	struct simplex_SolveContext context;
	struct simplex_Execution h;
	struct simplex_ExecutionScope original, reduced, scaled;
	struct lp_simplex_Result r;
	int root, action, i, count;
	simplex_solve_context_init(&context, &h, &original);
	root = simplex_execution_begin(
		&original, SIMPLEX_EXECUTION_SOLVE, 0, 0, 0);
	simplex_execution_set_parent(&original, root);
	reduced = simplex_solve_context_derive(
		&original, SIMPLEX_EXECUTION_PRESOLVE);
	scaled = simplex_solve_context_derive(
		&reduced, SIMPLEX_EXECUTION_SCALING);
	CHECK(original.model != reduced.model && reduced.model != scaled.model);
	CHECK(h.events[1].source_model == original.model);
	CHECK(h.events[2].source_model == reduced.model);
	CHECK(!h.events[1].result_available && isnan(h.events[1].primal.value));
	action = simplex_execution_begin(&scaled, SIMPLEX_EXECUTION_BACKEND, 0,
		SIMPLEX_REASON_NUMERICAL_FAILURE, root);
	CHECK(h.events[action - 1].model == scaled.model);
	CHECK(h.events[action - 1].parent == root);
	memset(&r, 0, sizeof(r));
	r.status = lp_simplex_PrecisionError;
	r.iterations = 7;
	simplex_execution_end(&h, action, &r, NULL);
	CHECK(h.events[action - 1].result_available &&
		h.events[action - 1].work == 7);
	count = h.count;
	for (i = count; i < SIMPLEX_EXECUTION_CAPACITY + 3; i++)
		(void)simplex_execution_begin(
			&scaled, SIMPLEX_EXECUTION_AUTO, 0, 0, 0);
	CHECK(h.count == SIMPLEX_EXECUTION_CAPACITY && h.dropped == 3);
	/* Exhaustion cannot prevent completing an existing parent. */
	simplex_execution_end(&h, root, &r, NULL);
	CHECK(h.events[0].complete && h.events[0].work == 7);
	/* Model namespaces are solve state, even when no trace was ever created. */
	simplex_solve_context_destroy(&context);
	simplex_solve_context_init(&context, NULL, &original);
	reduced = simplex_solve_context_derive(
		&original, SIMPLEX_EXECUTION_PRESOLVE);
	scaled = simplex_solve_context_derive(
		&reduced, SIMPLEX_EXECUTION_SCALING);
	CHECK(!scaled.history && original.model == 1 && reduced.model == 2 &&
		scaled.model == 3);
	context.models = INT_MAX;
	scaled = simplex_solve_context_derive(
		&reduced, SIMPLEX_EXECUTION_SCALING);
	CHECK(scaled.model == 0 && context.models == INT_MAX);
	simplex_solve_context_destroy(&context);
	return 0;
}

static int journal_dual(const struct simplex_Problem *p,
	const struct lp_simplex_Options *o, double *x, double *y,
	struct lp_simplex_Result *r, int propagate)
{
	(void)o;
	(void)y;
	(void)propagate;
	CHECK(p->execution.parent > 0);
	x[0] = x[1] = 0.;
	r->status = lp_simplex_ExceedIterLimit;
	r->iterations = 1;
	return -1;
}

static int journal_hybrid(const struct simplex_Problem *p,
	const struct lp_simplex_Options *o, const double *seed, double *x,
	struct lp_simplex_Result *r)
{
	(void)o;
	CHECK(p->execution.parent > 0 && seed && seed[0] == .5);
	x[0] = x[1] = .5;
	r->status = lp_simplex_MemoryAllocError;
	r->iterations = 1;
	return -1;
}

static int journal_initializer(void *context,
	const struct simplex_Problem *input, const struct lp_simplex_Options *o,
	double *x, int *work)
{
	(void)input;
	(void)context;
	(void)o;
	x[0] = x[1] = .5;
	*work = 1;
	return 0; /* Useful point despite a failed initializer. */
}

static int provenance(void)
{
	struct lp_Model *m = lp_model_create(1, 2);
	struct simplex_Problem p;
	struct simplex_SolveContext context;
	struct simplex_Execution h;
	struct simplex_AutoHistory local;
	struct simplex_AutoKernels kernels;
	struct simplex_AutoBindings bindings = {NULL, NULL, 0};
	struct simplex_AutoInitializer initializer;
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
	simplex_solve_context_init(&context, &h, &p.execution);
	kernels.dual = journal_dual;
	kernels.hybrid = journal_hybrid;
	initializer.context = NULL;
	initializer.run = journal_initializer;
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_AUTO);
	CHECK(simplex_auto_solve(
		      &p, &o, x, &r, 0, &bindings, &initializer, &local) != 0);
	CHECK(local.count == 4 && h.count == 5);
	for (i = 0; i < local.count; i++) {
		CHECK(local.stages[i].action.execution_id == i + 1);
		CHECK(h.events[i].complete && h.events[i].work == 1);
		CHECK(h.events[i].grant == local.stages[i].action.grant);
	}
	/* Failure and point have different sources; an initialized feasible point
	 * does not erase the failed dual action that caused recovery. */
	CHECK(h.events[2].cause == 1 && h.events[2].point_source == 2);
	CHECK(h.events[2].reason == SIMPLEX_REASON_DUAL_HANDOFF);
	CHECK(!h.events[3].point_source); /* Cold dual ignores primal seed. */
	CHECK(h.events[3].cause == 3);
	CHECK(h.events[3].reason == SIMPLEX_REASON_RECOVERY_EXHAUSTED);
	CHECK(r.iterations == 4);
	simplex_problem_destroy(&p);
	lp_model_free(m);
	simplex_solve_context_destroy(&context);
	return 0;
}

static int integration(void)
{
	struct simplex_Execution h;
	struct lp_simplex_Options o;
	struct lp_simplex_Result r;
	struct lp_Model *m = lp_model_create(3, 3);
	double x[3];
	int i, presolve, autos, backends, audits, reduced, transformed;
	CHECK(m);
	for (i = 0; i < 3; i++) {
		double scale = i == 0 ? 1048576. : i == 1 ? .0009765625 : 32.;
		m->constraints[i].type = optm_CONS_T_GE;
		m->constraints[i].rhs = scale;
		m->constraints[i].coef[i] = scale;
		m->constraints[i].coef[(i + 1) % 3] = scale;
		m->bounds[i].b_type = optm_BOUND_T_BS;
		m->bounds[i].ub = 1.;
		m->objective[i] = 1.;
	}
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_AUTO);
	CHECK(!o.progress_callback);
	for (presolve = 0; presolve < 2; presolve++) {
		o.presolve = presolve;
		CHECK(!simplex_solve_reported(m, &o, x, &r, &h));
		CHECK(r.status == lp_simplex_Success &&
			fabs(r.objective - 1.5) < 1e-7);
		CHECK(!h.dropped && h.events[0].work == r.iterations);
		autos = backends = audits = reduced = transformed = 0;
		for (i = 0; i < h.count; i++) {
			const struct simplex_ExecutionEvent *e = h.events + i;
			CHECK(e->complete && e->model > 0);
			CHECK(e->parent >= 0 && e->parent <= i);
			CHECK(e->cause >= 0 && e->cause <= i);
			if (e->kind == SIMPLEX_EXECUTION_AUTO) {
				autos++;
				CHECK(e->parent > 0);
			}
			if (e->kind == SIMPLEX_EXECUTION_BACKEND) {
				backends++;
				CHECK(e->parent > 0);
				CHECK(h.events[e->parent - 1].kind ==
					SIMPLEX_EXECUTION_AUTO);
			}
			if (e->kind == SIMPLEX_EXECUTION_SCALING) {
				transformed++;
				CHECK(e->model != e->source_model);
			}
			if (e->kind == SIMPLEX_EXECUTION_REDUCED_SOLVE)
				reduced++;
			if (e->kind == SIMPLEX_EXECUTION_POSTSOLVE) {
				audits++;
				CHECK(e->model == 1 && e->cause > 0);
				CHECK(h.events[e->cause - 1].kind ==
					SIMPLEX_EXECUTION_REDUCED_SOLVE);
				CHECK(e->primal.quality ==
					SIMPLEX_EVIDENCE_RECOMPUTED);
				CHECK(e->primal.value <= o.primal_tolerance);
			}
		}
		CHECK(autos > 0 && backends > 0 && transformed > 0);
		CHECK(audits == presolve && reduced == presolve);
	}
	lp_model_free(m);
	return 0;
}

struct root_reset {
	struct simplex_Execution *history;
	const char *stage;
	int done, id, count;
};

static void reset_root(
	void *context, const struct lp_simplex_Progress *progress)
{
	struct root_reset *reset = context;
	struct simplex_ExecutionScope scope;
	struct simplex_Execution previous;
	int i;
	if (reset->done || !reset->stage ||
		strcmp(progress->stage, reset->stage))
		return;
	memset(&scope, 0, sizeof(scope));
	scope.history = reset->history;
	scope.context = (struct simplex_SolveContext *)reset->history->events[0]
				.context;
	scope.model = reset->history->events[0].model;
	previous = *reset->history;
	reset->count = previous.count;
	memset(reset->history, 0, sizeof(*reset->history));
	for (i = 0; i < previous.count; i++) {
		int id;
		scope.model = previous.events[i].model;
		id = simplex_execution_begin(&scope, previous.events[i].kind,
			previous.events[i].detail, 0, 0);
		reset->history->events[id - 1].grant = 777;
	}
	reset->id = 1;
	reset->done = 1;
}

static int root_lifetime(void)
{
	struct lp_Model *m = lp_model_create(3, 4);
	struct simplex_Execution history;
	struct root_reset reset;
	struct lp_simplex_Options options;
	struct lp_simplex_Result result;
	double x[4];
	int algorithm, pre, mode, i, expected;
	CHECK(m);
	for (i = 0; i < 3; i++) {
		m->constraints[i].type = optm_CONS_T_GE;
		m->constraints[i].rhs = 1.;
		m->constraints[i].coef[i] =
			m->constraints[i].coef[(i + 1) % 3] = 1.;
		m->bounds[i].b_type = optm_BOUND_T_BS;
		m->bounds[i].ub = 1.;
		m->objective[i] = 1.;
	}
	m->bounds[3].b_type = optm_BOUND_T_BS;
	m->bounds[3].lb = m->bounds[3].ub = 0.;
	m->objective[3] = 2.;
	for (algorithm = 0; algorithm <= lp_simplex_ALGORITHM_AUTO; algorithm++)
		for (pre = 0; pre < 2; pre++) {
			expected = -1;
			for (mode = 0; mode < 5; mode++) {
				if (mode == 3 &&
					algorithm !=
						lp_simplex_ALGORITHM_PAN_BDA)
					continue;
				if (mode == 4 &&
					algorithm !=
						lp_simplex_ALGORITHM_DUAL_REVISED &&
					algorithm != lp_simplex_ALGORITHM_AUTO)
					continue;
				lp_simplex_default_options(&options, algorithm);
				options.presolve = pre;
				memset(&reset, 0, sizeof(reset));
				reset.history = &history;
				reset.stage = mode == 1	  ? "start"
					      : mode == 2 ? "finish"
					      : mode == 3 ? "pan-finish"
					      : mode == 4 ? "dual"
							  : NULL;
				options.progress_callback = reset_root;
				options.progress_context = &reset;
				CHECK(!simplex_solve_reported(
					m, &options, x, &result, &history));
				CHECK(result.status == lp_simplex_Success &&
					fabs(result.objective - 1.5) < 1e-7);
				if (mode == 0)
					expected = result.iterations;
				CHECK(result.iterations == expected);
				if (mode) {
					CHECK(reset.done && reset.id > 0);
					CHECK(!history.events[reset.id - 1]
							.complete);
					CHECK(history.events[reset.id - 1]
							.grant == 777);
					for (i = 0; i < reset.count; i++)
						CHECK(!history.events[i]
								.complete);
					for (i = reset.count; i < history.count;
						i++) {
						CHECK(!history.events[i]
								.parent ||
							history.events[i]
									.parent >
								reset.count);
						CHECK(!history.events[i]
								.cause ||
							history.events[i]
									.cause >
								reset.count);
						CHECK(!history.events[i]
								.point_source ||
							history.events[i]
									.point_source >
								reset.count);
					}
				}
			}
		}
	lp_model_free(m);
	return 0;
}

static int references(void)
{
	struct simplex_SolveContext context, other;
	struct simplex_Execution history, replacement;
	struct simplex_ExecutionScope scope, alternate;
	struct simplex_ExecutionRef reference;
	int id;

	simplex_solve_context_init(&context, &history, &scope);
	id = simplex_execution_begin(&scope, SIMPLEX_EXECUTION_AUTO, 3, 0, 0);
	reference = simplex_execution_ref(&scope, id);
	CHECK(simplex_execution_ref_id(&reference, &scope, 0) == id);
	CHECK(!simplex_execution_ref_id(&reference, &scope, 1));
	simplex_execution_end(&history, id, NULL, NULL);
	CHECK(simplex_execution_ref_id(&reference, &scope, 1) == id);
	alternate = scope;
	replacement = history;
	alternate.history = &replacement;
	CHECK(!simplex_execution_ref_id(&reference, &alternate, 1));
	simplex_solve_context_init(&other, NULL, &alternate);
	alternate.history = &history;
	CHECK(!simplex_execution_ref_id(&reference, &alternate, 1));
	/* Same slot, context, model, kind and detail; only the actor token
	 * distinguishes this replacement after an in-place journal reset. */
	memset(&history, 0, sizeof(history));
	CHECK(simplex_execution_begin(
		      &scope, SIMPLEX_EXECUTION_AUTO, 3, 0, 0) == id);
	simplex_execution_end(&history, id, NULL, NULL);
	CHECK(!simplex_execution_ref_id(&reference, &scope, 0));
	context.execution_serial = ULONG_MAX;
	id = simplex_execution_begin(&scope, SIMPLEX_EXECUTION_AUTO, 3, 0, 0);
	reference = simplex_execution_ref(&scope, id);
	CHECK(!simplex_execution_ref_id(&reference, &scope, 0));
	simplex_solve_context_destroy(&other);
	simplex_solve_context_destroy(&context);
	return 0;
}

int main(void)
{
	return admission() || identities() || provenance() || integration() ||
	       references() || root_lifetime();
}
