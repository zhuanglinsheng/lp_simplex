/**
 * @file test_retry.c
 * @brief Recovery admission depends on consumed state, not caller depth or
 * clocks.
 */

#include "simplex_stage.h"
#include "simplex_solve_context.h"
#include "simplex_execution.h"
#include "simplex_retry.h"
#include "simplex_dual.h"
#include "simplex_basis.h"
#include "policy/controls.h"
#include <lp_simplex/lp_simplex.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define CHECK(c) \
	do { \
		if (!(c)) { \
			fprintf(stderr, "retry line %d: %s\n", __LINE__, #c); \
			return 1; \
		} \
	} while (0)

static int ledger_contract(void)
{
	struct simplex_RetryLedger ledger;
	struct simplex_RetryKey key, changed;
	struct simplex_Problem input;
	struct simplex_SolveContext context;
	struct lp_simplex_Result r;
	int i;
	memset(&ledger, 0, sizeof(ledger));
	memset(&key, 0, sizeof(key));
	memset(&input, 0, sizeof(input));
	simplex_solve_context_init(&context, NULL, &input.execution);
	key.input = simplex_problem_identity(&input);
	memset(&r, 0, sizeof(r));
	key.known = 1;
	key.model = 1;
	key.action = SIMPLEX_EXECUTION_BACKEND;
	key.backend = SIMPLEX_BASIS_BACKEND_INTERNAL;
	key.algorithm = lp_simplex_ALGORITHM_AUTO;
	key.pricing = lp_simplex_PRICING_DUAL_STEEPEST_EDGE;
	key.primal_tolerance = key.dual_tolerance = 1e-7;
	key.pivot_tolerance = 1e-9;
	r.status = lp_simplex_PrecisionError;
	r.iterations = 3;
	simplex_retry_record(
		&ledger, &key, NULL, 7, lp_simplex_EXIT_FAILURE, &r, 100, 0);
	CHECK(ledger.count == 1 && simplex_retry_failed(&ledger, &key) == 1);
	CHECK(ledger.records[0].event == 7 && ledger.records[0].work == 3);
	/* More unused budget or a different caller is not a new cold start. */
	simplex_retry_record(
		&ledger, &key, NULL, 20, lp_simplex_EXIT_FAILURE, &r, 1000, 0);
	CHECK(ledger.count == 1);
#define DIFFERENT(field, value) \
	do { \
		changed = key; \
		changed.field = value; \
		CHECK(!simplex_retry_failed(&ledger, &changed)); \
	} while (0)
	DIFFERENT(model, 2);
	DIFFERENT(revision, 1);
	DIFFERENT(action, SIMPLEX_EXECUTION_STABLE_INITIALIZE);
	DIFFERENT(backend, SIMPLEX_BASIS_BACKEND_KLU);
	DIFFERENT(algorithm, lp_simplex_ALGORITHM_DUAL_REVISED);
	DIFFERENT(pricing, lp_simplex_PRICING_BLAND);
	DIFFERENT(representation, SIMPLEX_DUAL_RETRY_EQUILIBRATED);
	DIFFERENT(propagate, 1);
	DIFFERENT(feasibility_only, 1);
	DIFFERENT(controls, 1);
	DIFFERENT(primal_tolerance, 1e-6);
	DIFFERENT(dual_tolerance, 1e-6);
	DIFFERENT(pivot_tolerance, 1e-8);
	DIFFERENT(input.input, NULL);
	DIFFERENT(input.context, NULL);
	DIFFERENT(input.model, 2);
	DIFFERENT(input.revision, 1);
	DIFFERENT(input.rows, 1);
	DIFFERENT(input.columns, 1);
	DIFFERENT(input.nonzeros, 1);
	DIFFERENT(input.values, &r);
	DIFFERENT(input.starts, &r);
	DIFFERENT(input.indices, &r);
	DIFFERENT(input.row_starts, &r);
	DIFFERENT(input.column_indices, &r);
	DIFFERENT(input.row_values, &r);
	DIFFERENT(input.rhs, &r);
	DIFFERENT(input.cost, &r);
	DIFFERENT(input.bounds, &r);
	DIFFERENT(input.row_type, &r);
	DIFFERENT(known, 0);
	DIFFERENT(model, 0);
	DIFFERENT(revision, ULONG_MAX);
	DIFFERENT(primal_tolerance, NAN);
#undef DIFFERENT
	/* Budget, allocation, mathematical conclusions and malformed work are
	 * not completed numerical failures of an otherwise identical start. */
	changed = key;
	changed.model = 2;
	changed.input.model = 2;
	r.status = lp_simplex_ExceedIterLimit;
	simplex_retry_record(&ledger, &changed, NULL, 0, -1, &r, 100, 0);
	r.status = lp_simplex_MemoryAllocError;
	simplex_retry_record(&ledger, &changed, NULL, 0, -1, &r, 100, 0);
	r.status = lp_simplex_Infeasibility;
	simplex_retry_record(&ledger, &changed, NULL, 0, -1, &r, 100, 0);
	r.status = lp_simplex_PrecisionError;
	simplex_retry_record(&ledger, &changed, NULL, 0, 0, &r, 100, 0);
	r.iterations = 100;
	simplex_retry_record(&ledger, &changed, NULL, 0, -1, &r, 100, 0);
	r.iterations = -1;
	simplex_retry_record(&ledger, &changed, NULL, 0, -1, &r, 100, 0);
	CHECK(!simplex_retry_failed(&ledger, &changed) && ledger.count == 1);
	/* Recording without a trace id still protects against a known repeat. */
	r.iterations = 0;
	simplex_retry_record(&ledger, &changed, NULL, 0, -1, &r, 100, 0);
	CHECK(simplex_retry_failed(&ledger, &changed) == 2);
	for (i = 3; i <= SIMPLEX_RETRY_CAPACITY + 1; i++) {
		changed.model = i;
		changed.input.model = i;
		simplex_retry_record(
			&ledger, &changed, NULL, 0, -1, &r, 100, 0);
	}
	CHECK(ledger.count == SIMPLEX_RETRY_CAPACITY && ledger.unrecorded == 1);
	CHECK(!simplex_retry_failed(&ledger, &changed));
	CHECK(simplex_retry_failed(&ledger, &key) == 1);
	simplex_retry_destroy(&ledger);
	simplex_solve_context_destroy(&context);
	return 0;
}

static int shared_evidence(int tracing)
{
	struct simplex_SolveContext context;
	struct simplex_Execution execution;
	struct simplex_Problem p;
	struct simplex_Problem other;
	struct simplex_DualRetry first, later;
	struct simplex_StageReport report;
	struct lp_simplex_Options o;
	struct lp_simplex_Result r;
	double changed_rhs = 1.;
	int i;
	memset(&p, 0, sizeof(p));
	simplex_solve_context_init(
		&context, tracing ? &execution : NULL, &p.execution);
	p.report = &report;
	simplex_stage_reset(&report);
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_AUTO);
	simplex_dual_retry_init(&first, &p, &o, 0, 0, 0);
	memset(&r, 0, sizeof(r));
	r.status = lp_simplex_PrecisionError;
	r.iterations = 3;
	simplex_dual_retry_record(
		&first, SIMPLEX_BASIS_BACKEND_INTERNAL, 0, -1, &r, 100, 0);
	CHECK(first.key.known);
	CHECK(context.retries.count == 1);
	/* A short action can stop before the previous failure and export a useful
	 * checkpoint. Equal and shorter grants are not classified by that failure. */
	CHECK(!simplex_dual_retry_blocked(
		&first, &p, SIMPLEX_BASIS_BACKEND_INTERNAL, 2));
	CHECK(!simplex_dual_retry_blocked(
		&first, &p, SIMPLEX_BASIS_BACKEND_INTERNAL, 3));
	CHECK(simplex_dual_retry_blocked(
		&first, &p, SIMPLEX_BASIS_BACKEND_INTERNAL, 4));
	/* Exhausting optional trace storage cannot erase recovery decision state. */
	if (tracing)
		for (i = execution.count; i < SIMPLEX_EXECUTION_CAPACITY; i++)
			(void)simplex_execution_begin(
				&p.execution, SIMPLEX_EXECUTION_AUTO, 0, 0, 0);
	simplex_stage_reset(&report);
	o.iteration_limit *= 2;
	o.progress_interval *= 2.;
	simplex_execution_set_parent(&p.execution, 1);
	simplex_dual_retry_init(&later, &p, &o, 0, 0, 0);
	CHECK(simplex_dual_retry_blocked(
		&later, &p, SIMPLEX_BASIS_BACKEND_INTERNAL, 100));
	CHECK(report.recovery_skipped == 1);
	if (tracing)
		CHECK(execution.dropped == 1);
	CHECK(report.stop ==
		SIMPLEX_STOP_UNKNOWN); /* No invented kernel stop/result. */
	/* Removing the entire output journal cannot erase an earlier failure. */
	p.execution.history = NULL;
	simplex_dual_retry_init(&later, &p, &o, 0, 0, 0);
	CHECK(simplex_dual_retry_blocked(
		&later, &p, SIMPLEX_BASIS_BACKEND_INTERNAL, 100));
	CHECK(report.recovery_skipped == 2);
	/* Input storage and the live owner are part of the exact start. A copied
	 * view or replacement pointer cannot inherit another start's failure. */
	other = p;
	simplex_dual_retry_init(&later, &other, &o, 0, 0, 0);
	CHECK(!simplex_dual_retry_blocked(
		&later, &other, SIMPLEX_BASIS_BACKEND_INTERNAL, 100));
	p.rhs = &changed_rhs;
	CHECK(!simplex_dual_retry_blocked(
		&first, &p, SIMPLEX_BASIS_BACKEND_INTERNAL, 100));
	simplex_dual_retry_init(&later, &p, &o, 0, 0, 0);
	CHECK(!simplex_dual_retry_blocked(
		&later, &p, SIMPLEX_BASIS_BACKEND_INTERNAL, 100));
	/* A callback/input mutation cannot record a failure under the old input. */
	r.iterations = 1;
	simplex_dual_retry_record(
		&first, SIMPLEX_BASIS_BACKEND_KLU, 0, -1, &r, 100, 0);
	CHECK(context.retries.count == 1);
	p.rhs = NULL;
	simplex_dual_retry_init(&later, &p, &o, 0, 0, 1);
	CHECK(!simplex_dual_retry_blocked(
		&later, &p, SIMPLEX_BASIS_BACKEND_INTERNAL, 100));
	simplex_problem_mark_changed(&p);
	simplex_dual_retry_init(&later, &p, &o, 0, 0, 0);
	CHECK(!simplex_dual_retry_blocked(
		&later, &p, SIMPLEX_BASIS_BACKEND_INTERNAL, 100));
	simplex_solve_context_destroy(&context);
	return 0;
}

/* Inject one known default-backend failure into the numerical admission seam.
 * The actual second backend must still solve, with no finite placeholder point
 * inherited from the skipped primary. No backend means no fabricated success. */
static int remaining_backend(int tracing)
{
	struct lp_Model *m = lp_model_create(1, 1);
	struct simplex_Problem p;
	struct simplex_SolveContext context;
	struct simplex_Execution execution;
	struct simplex_StageReport report;
	struct simplex_DualRetry retry;
	struct simplex_Basis factor;
	struct lp_simplex_Options o;
	struct lp_simplex_Result failed, result;
	int index[1], backend, source, state, i, backends = 0, skips = 0;
	double x[1], dual[1];
	CHECK(m);
	m->constraints[0].type = optm_CONS_T_EQ;
	m->constraints[0].rhs = m->constraints[0].coef[0] = 1.;
	m->objective[0] = 1.;
	CHECK(!simplex_problem_from_model(&p, m));
	simplex_solve_context_init(
		&context, tracing ? &execution : NULL, &p.execution);
	simplex_solve_context_configure(&context, 0);
	context.controls.disable_scaling = 1;
	p.report = &report;
	index[0] = 0;
	CHECK(!simplex_basis_create_configured(
		&factor, &p.matrix, p.columns, index, 0));
	backend = simplex_basis_resolved_backend(&factor);
	simplex_basis_destroy(&factor);
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_DUAL_REVISED);
	o.iteration_limit = 100;
	simplex_dual_retry_init(&retry, &p, &o, 0, 0, 0);
	memset(&failed, 0, sizeof(failed));
	failed.status = lp_simplex_Singularity;
	failed.iterations = 1;
	source =
		simplex_execution_begin(&p.execution, SIMPLEX_EXECUTION_BACKEND,
			SIMPLEX_BASIS_BACKEND_DEFAULT, SIMPLEX_REASON_NONE, 0);
	simplex_execution_end(p.execution.history, source, &failed, NULL);
	simplex_dual_retry_record(&retry, SIMPLEX_BASIS_BACKEND_DEFAULT, source,
		-1, &failed, 100, 0);
	CHECK(context.retries.count == 1);
	x[0] = dual[0] = 123.;
	state = simplex_dual_solve_problem(&p, &o, x, dual, &result, 0);
	CHECK(report.recovery_skipped == 1);
	if (backend == SIMPLEX_BASIS_BACKEND_KLU) {
		CHECK(state == 0 && result.status == lp_simplex_Success);
		CHECK(fabs(x[0] - 1.) <= o.primal_tolerance &&
			isfinite(dual[0]));
		CHECK(report.backend_retries == 1 && report.point_available);
	} else {
		CHECK(state == -1 && result.status == lp_simplex_Singularity);
		CHECK(result.iterations == 0 && !isfinite(x[0]) &&
			!isfinite(dual[0]));
		CHECK(!report.point_available && !report.basis_available);
	}
	for (i = source; tracing && i < execution.count; i++) {
		const struct simplex_ExecutionEvent *e = execution.events + i;
		if (e->kind == SIMPLEX_EXECUTION_BACKEND) {
			backends++;
			CHECK(e->detail == SIMPLEX_BASIS_BACKEND_INTERNAL &&
				e->cause == source);
		}
		if (e->kind == SIMPLEX_EXECUTION_RECOVERY_SKIP) {
			skips++;
			CHECK(e->cause == source && !e->result_available);
		}
	}
	if (tracing)
		CHECK(skips == 1 &&
			backends == (backend == SIMPLEX_BASIS_BACKEND_KLU));
	/* A shorter/equal prefix is useful new work: it may yield before the
	 * recorded breakdown. Admission must not suppress that actual backend. */
	source = tracing ? execution.count : 0;
	o.iteration_limit = failed.iterations;
	(void)simplex_dual_solve_problem(&p, &o, x, dual, &result, 0);
	CHECK(result.iterations <= o.iteration_limit &&
		!report.recovery_skipped);
	backends = 0;
	for (i = source; tracing && i < execution.count; i++)
		backends +=
			execution.events[i].kind == SIMPLEX_EXECUTION_BACKEND;
	if (tracing)
		CHECK(backends == 1);
	simplex_problem_destroy(&p);
	lp_model_free(m);
	simplex_solve_context_destroy(&context);
	return 0;
}

static int backend_events(const struct simplex_Execution *execution, int begin)
{
	int i, count = 0;
	for (i = begin; i < execution->count; i++)
		count += execution->events[i].kind == SIMPLEX_EXECUTION_BACKEND;
	return count;
}

static int kernel_recovery(int tracing)
{
	static struct lp_simplex_Result no_trace_first, no_trace_second;
	struct lp_Model *m;
	struct simplex_Problem p;
	struct simplex_SolveContext context;
	struct simplex_Execution execution, alternate;
	struct simplex_Execution *journal = &execution;
	struct simplex_StageReport report;
	struct lp_simplex_Options o;
	struct lp_simplex_Result first, second;
	double x[1], replacement = 1., *original_rhs;
	int calls, begin, i, skips = 0;
	if (!getenv(SIMPLEX_ENV_DISABLE_SCALING))
		return 0;
	m = lp_model_create(1, 1);
	CHECK(m);
	/* Feasible x=1e300. With scaling disabled, the usable pivot lies below
	 * kernel resolution; rejecting its free-column infeasibility proof must
	 * yield a numerical failure rather than a mathematical conclusion. */
	m->constraints[0].type = optm_CONS_T_EQ;
	m->constraints[0].rhs = 1.;
	m->constraints[0].coef[0] = 1e-300;
	m->bounds[0].b_type = optm_BOUND_T_FR;
	CHECK(!simplex_problem_from_model(&p, m));
	simplex_solve_context_init(
		&context, tracing ? &execution : NULL, &p.execution);
	p.report = &report;
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_DUAL_REVISED);
	o.iteration_limit = 100;
	CHECK(simplex_dual_solve_problem(&p, &o, x, NULL, &first, 0));
	CHECK(first.status == lp_simplex_PrecisionError ||
		first.status == lp_simplex_Singularity);
	calls = context.retries.count;
	CHECK(calls >= 1 && context.retries.count == calls);
	if (tracing)
		CHECK(calls == backend_events(&execution, 0));
	if (tracing == 2)
		for (i = execution.count; i < SIMPLEX_EXECUTION_CAPACITY; i++)
			(void)simplex_execution_begin(
				&p.execution, SIMPLEX_EXECUTION_AUTO, 0, 0, 0);
	if (tracing == 3) {
		memset(&alternate, 0, sizeof(alternate));
		p.execution.history = journal = &alternate;
	}
	begin = tracing ? journal->count : 0;
	CHECK(simplex_dual_solve_problem(&p, &o, x, NULL, &second, 0));
	CHECK(second.status == first.status &&
		second.iterations <= first.iterations);
	if (!tracing) {
		no_trace_first = first;
		no_trace_second = second;
	} else {
		CHECK(first.status == no_trace_first.status &&
			first.iterations == no_trace_first.iterations);
		CHECK(second.status == no_trace_second.status &&
			second.iterations == no_trace_second.iterations);
	}
	/* Every known identical cold backend is skipped, including the primary.
	 * A skipped call cannot manufacture a new point or mathematical result. */
	if (tracing)
		CHECK(backend_events(journal, begin) == 0);
	CHECK(second.iterations == 0 && !isfinite(x[0]));
	CHECK(report.recovery_skipped == calls);
	for (i = begin; tracing && i < journal->count; i++) {
		const struct simplex_ExecutionEvent *e = journal->events + i;
		if (e->kind == SIMPLEX_EXECUTION_RECOVERY_SKIP) {
			skips++;
			if (tracing == 3)
				CHECK(e->cause ==
					0); /* Old journal ids are not portable. */
			else
				CHECK(e->cause > 0 && e->cause <= begin);
			CHECK(e->reason == SIMPLEX_REASON_REPEATED_FAILURE &&
				!e->result_available);
		}
	}
	if (tracing == 1 || tracing == 3)
		CHECK(skips == calls);
	if (tracing == 2)
		CHECK(skips == 0 && execution.dropped == calls);
	/* Replacing RHS storage is a new immutable input identity, even if its
	 * values happen to be equal. Both numerical backends must run again. */
	original_rhs = p.rhs;
	p.rhs = &replacement;
	begin = tracing ? journal->count : 0;
	CHECK(simplex_dual_solve_problem(&p, &o, x, NULL, &second, 0));
	CHECK(second.status == first.status && report.recovery_skipped == 0);
	CHECK(context.retries.count == 2 * calls);
	if (tracing == 1 || tracing == 3)
		CHECK(backend_events(journal, begin) == calls);
	p.rhs = original_rhs;
	simplex_problem_destroy(&p);
	lp_model_free(m);
	simplex_solve_context_destroy(&context);
	return 0;
}

int main(void)
{
	return ledger_contract() || shared_evidence(0) || shared_evidence(1) ||
	       kernel_recovery(0) || kernel_recovery(1) || kernel_recovery(2) ||
	       kernel_recovery(3) || remaining_backend(0) ||
	       remaining_backend(1);
}
