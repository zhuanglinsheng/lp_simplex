/**
 * @file test_basis_retry.c
 * @brief Exact basis requests, backend adaptation and snapshot ownership.
 */

#include "simplex_stage.h"
#include "simplex_solve_context.h"
#include "simplex_execution.h"
#include "simplex_dual.h"
#include "simplex_dual_internal.h"
#include "simplex_basis.h"
#include "policy/resources.h"
#include <lp_simplex/lp_simplex.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define CHECK(c) \
	do { \
		if (!(c)) { \
			fprintf(stderr, "basis retry line %d: %s\n", __LINE__, \
				#c); \
			return 1; \
		} \
	} while (0)

static int snapshots(void)
{
	struct lp_Model *m = lp_model_create(1, 2);
	struct simplex_Problem p;
	struct simplex_SolveContext c;
	struct simplex_DualBasisSeed seed;
	struct simplex_DualRetry a, b;
	struct lp_simplex_Options o;
	struct lp_simplex_Result r;
	int index[1] = {0}, copy[1] = {0};
	unsigned char status[3] = {0, 1, 3}, copied[3] = {0, 1, 3};
	CHECK(m && !simplex_problem_from_model(&p, m));
	simplex_solve_context_init(&c, NULL, &p.execution);
	simplex_solve_context_configure(&c, 0);
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_DUAL_REVISED);
	simplex_dual_seed_init(&seed, &p, index, status);
	seed.valid = 1;
	simplex_dual_retry_init(&a, &p, &o, 0, 0, 1);
	simplex_dual_retry_bind_basis(&a, &seed);
	memset(&r, 0, sizeof(r));
	r.status = lp_simplex_PrecisionError;
	r.iterations = 3;
	simplex_dual_retry_record(
		&a, SIMPLEX_BASIS_BACKEND_INTERNAL, 0, -1, &r, 100, 1);
	CHECK(c.retries.count == 1 &&
		c.retries.basis_bytes == sizeof(index) + sizeof(status));
	CHECK(c.retries.records[0].key.basis != index &&
		c.retries.records[0].key.status != status);
	CHECK(simplex_dual_retry_blocked(
		&a, &p, SIMPLEX_BASIS_BACKEND_INTERNAL, 100));
	index[0] = 1;
	CHECK(!simplex_dual_retry_blocked(
		&a, &p, SIMPLEX_BASIS_BACKEND_INTERNAL, 100));
	CHECK(c.retries.records[0].key.basis[0] == 0);
	index[0] = 0;
	status[1] = 2;
	CHECK(!simplex_dual_retry_blocked(
		&a, &p, SIMPLEX_BASIS_BACKEND_INTERNAL, 100));
	status[1] = 1;
	seed.index = copy;
	seed.status = copied;
	simplex_dual_retry_init(&b, &p, &o, 0, 0, 1);
	simplex_dual_retry_bind_basis(&b, &seed);
	CHECK(simplex_dual_retry_blocked(
		&b, &p, SIMPLEX_BASIS_BACKEND_INTERNAL, 100));
	seed.prefer_dual = 0;
	simplex_dual_retry_bind_basis(&b, &seed);
	CHECK(!simplex_dual_retry_blocked(
		&b, &p, SIMPLEX_BASIS_BACKEND_INTERNAL, 100));
	seed.prefer_dual = 1;
	seed.valid = 0;
	simplex_dual_retry_bind_basis(&b, &seed);
	CHECK(!simplex_dual_retry_blocked(
		&b, &p, SIMPLEX_BASIS_BACKEND_INTERNAL, 100));
	simplex_solve_context_destroy(&c);
	CHECK(!c.retries.count && !c.retries.basis_bytes);
	simplex_solve_context_destroy(&c);
	simplex_problem_destroy(&p);
	lp_model_free(m);
	return 0;
}

static int storage_limit(void)
{
	struct simplex_RetryLedger ledger;
	struct simplex_RetryKey key;
	struct simplex_Problem p;
	struct simplex_SolveContext c;
	struct lp_simplex_Result r;
	int basis = 0;
	unsigned char *status =
		(unsigned char *)calloc(SIMPLEX_RETRY_BASIS_BYTES, 1);
	CHECK(status);
	memset(&ledger, 0, sizeof(ledger));
	memset(&p, 0, sizeof(p));
	memset(&key, 0, sizeof(key));
	simplex_solve_context_init(&c, NULL, &p.execution);
	p.rows = 1;
	p.columns = (int)(SIMPLEX_RETRY_BASIS_BYTES - sizeof(int) - 1);
	key.input = simplex_problem_identity(&p);
	key.known = key.model = key.action = key.backend = 1;
	key.primal_tolerance = key.dual_tolerance = key.pivot_tolerance = 1e-7;
	key.seeded = 1;
	key.basis = &basis;
	key.status = status;
	memset(&r, 0, sizeof(r));
	r.status = lp_simplex_Singularity;
	simplex_retry_record(&ledger, &key, NULL, 0, -1, &r, 100, 1);
	CHECK(ledger.count == 1 &&
		ledger.basis_bytes == SIMPLEX_RETRY_BASIS_BYTES);
	key.prefer_dual = 1;
	simplex_retry_record(&ledger, &key, NULL, 0, -1, &r, 100, 1);
	CHECK(ledger.count == 1 && ledger.unrecorded == 1);
	key.input.rows = key.input.columns = INT_MAX;
	CHECK(!simplex_retry_failed(&ledger, &key));
	simplex_retry_record(&ledger, &key, NULL, 0, -1, &r, 100, 1);
	CHECK(ledger.count == 1);
	simplex_retry_destroy(&ledger);
	CHECK(!ledger.basis_bytes && !ledger.count);
	simplex_solve_context_destroy(&c);
	free(status);
	return 0;
}

static int backend_count(const struct simplex_Execution *h, int begin)
{
	int i, n = 0;
	for (i = begin; i < h->count; i++)
		n += h->events[i].kind == SIMPLEX_EXECUTION_BACKEND;
	return n;
}

static int numerical(int invalid, int tracing)
{
	struct lp_Model *m = lp_model_create(1, 1);
	struct simplex_Problem p;
	struct simplex_SolveContext c;
	struct simplex_Execution h, alternate;
	struct simplex_StageReport report;
	struct simplex_DualBasisSeed seed;
	struct lp_simplex_Options o;
	struct lp_simplex_Result r;
	int index[1], begin, first_calls, i;
	unsigned char status[2] = {DUAL_STATUS_BASIC, DUAL_STATUS_FIXED};
	double x[1];
	CHECK(m);
	m->constraints[0].type = optm_CONS_T_EQ;
	m->constraints[0].rhs = 1.;
	m->constraints[0].coef[0] = 1e-320;
	m->bounds[0].b_type = optm_BOUND_T_FR;
	CHECK(!simplex_problem_from_model(&p, m));
	simplex_solve_context_init(&c, tracing ? &h : NULL, &p.execution);
	simplex_solve_context_configure(&c, 0);
	c.controls.disable_scaling = 1;
	p.report = &report;
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_DUAL_REVISED);
	o.iteration_limit = 100;
	index[0] = invalid ? -1 : 0;
	simplex_dual_seed_init(&seed, &p, index, status);
	seed.valid = 1;
	CHECK(simplex_dual_solve_problem_seeded(&p, &o, x, NULL, &r, 0, &seed));
	CHECK((r.status == lp_simplex_PrecisionError ||
		      r.status == lp_simplex_Singularity) &&
		c.retries.count >= 1);
	first_calls = c.retries.count;
	CHECK(c.retries.records[0].key.seeded &&
		c.retries.records[0].basis_retry == !invalid);
	if (first_calls > 1)
		CHECK(c.retries.records[1].key.seeded == !invalid);
	CHECK(c.retries.basis_bytes > 0);
	/* Retain only the real primary failure: the alternate has no known
	 * failure now and must receive the original hint, not a cold substitute. */
	if (!invalid && tracing == 1 && first_calls == 2) {
		free(c.retries.records[1].basis_storage);
		c.retries.basis_bytes -= sizeof(index) + sizeof(status);
		c.retries.count = 1;
		seed.valid = 1;
		begin = h.count;
		CHECK(simplex_dual_solve_problem_seeded(
			&p, &o, x, NULL, &r, 0, &seed));
		CHECK(backend_count(&h, begin) == 1 &&
			report.recovery_skipped == 1);
		CHECK(c.retries.count == 2);
		for (i = begin; i < h.count; i++)
			if (h.events[i].kind == SIMPLEX_EXECUTION_BACKEND)
				CHECK(h.events[i].basis_import_attempted);
		CHECK(c.retries.records[1].key.seeded);
	}
	if (tracing == 2)
		for (i = h.count; i < SIMPLEX_EXECUTION_CAPACITY; i++)
			(void)simplex_execution_begin(
				&p.execution, SIMPLEX_EXECUTION_AUTO, 0, 0, 0);
	if (tracing == 3) {
		memset(&alternate, 0, sizeof(alternate));
		p.execution.history = &alternate;
	}
	begin = tracing ? p.execution.history->count : 0;
	seed.valid = 1;
	CHECK(simplex_dual_solve_problem_seeded(&p, &o, x, NULL, &r, 0, &seed));
	CHECK(r.iterations == 0 && !isfinite(x[0]) && !seed.valid &&
		!report.basis_available);
	CHECK(report.recovery_skipped == first_calls &&
		c.retries.count == first_calls);
	if (tracing)
		CHECK(!backend_count(p.execution.history, begin));
	for (i = begin; tracing && i < p.execution.history->count; i++) {
		const struct simplex_ExecutionEvent *e =
			p.execution.history->events + i;
		if (e->kind == SIMPLEX_EXECUTION_RECOVERY_SKIP) {
			CHECK(!e->result_available);
			CHECK(tracing == 3 ? e->cause == 0
					   : e->cause > 0 && e->cause <= begin);
		}
	}
	/* Same storage, genuinely changed basis content: the primary must run. */
	index[0] = 1;
	status[0] = DUAL_STATUS_FREE;
	status[1] = DUAL_STATUS_BASIC;
	seed.valid = 1;
	begin = tracing ? p.execution.history->count : 0;
	CHECK(simplex_dual_solve_problem_seeded(&p, &o, x, NULL, &r, 0, &seed));
	CHECK(c.retries.count > first_calls);
	if (tracing && tracing != 2)
		CHECK(backend_count(p.execution.history, begin) >= 1);
	simplex_solve_context_destroy(&c);
	CHECK(!c.retries.basis_bytes);
	simplex_problem_destroy(&p);
	lp_model_free(m);
	return 0;
}

/* A finite bounded LP with zero as a feasible witness. Its ill-conditioned
 * powers-of-two matrix exercises nonzero-work numerical fallback in current
 * KLU builds, without an externally downloaded benchmark or fixed pivot trace. */
static int backend_selection(void)
{
	static const int exponent[5][8] = {{-7, 50, -3, -17, -4, -22, -6, -27},
		{41, 19, -15, 42, -41, 43, -48, 40},
		{36, 46, -9, 22, 29, 42, -28, 11},
		{-22, 33, 46, 0, -13, -1, -45, 39},
		{-35, -9, -24, 11, -32, -16, -5, -39}};
	static const unsigned int negative[5] = {0, 170, 0, 85, 0};
	static const int rhs_exponent[5] = {-12, -16, 4, 18, -10};
	static const double cost[8] = {
		-4., -15., -10., -17., -20., -11., -6., -13.};
	struct lp_Model *m = lp_model_create(5, 8);
	struct simplex_Problem p;
	struct simplex_SolveContext context;
	struct simplex_Execution history;
	struct simplex_ExecutionRef owner;
	struct simplex_StageReport report;
	struct lp_simplex_Options o;
	struct lp_simplex_Result result, reference;
	double x[8], saved[8], zero[8] = {0.};
	int mode, i, j, root, backends, work, selected, retry, point;

	CHECK(m);
	for (i = 0; i < 5; i++) {
		m->constraints[i].type = optm_CONS_T_LE;
		m->constraints[i].rhs = ldexp(1., rhs_exponent[i]);
		for (j = 0; j < 8; j++)
			m->constraints[i].coef[j] =
				ldexp(negative[i] & (1U << j) ? -1. : 1.,
					exponent[i][j]);
	}
	for (j = 0; j < 8; j++) {
		m->objective[j] = cost[j];
		m->bounds[j].b_type = optm_BOUND_T_BS;
		m->bounds[j].lb = 0.;
		m->bounds[j].ub = 2.;
	}
	CHECK(!simplex_problem_from_model(&p, m));
	CHECK(simplex_problem_primal_infeasibility(&p, zero) == 0.);
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_DUAL_REVISED);
	o.iteration_limit = 80;
	p.report = &report;
	for (mode = 0; mode < 4; mode++) {
		simplex_solve_context_init(
			&context, mode == 1 ? NULL : &history, &p.execution);
		simplex_solve_context_configure(&context, 0);
		context.controls.disable_scaling = 1;
		root = simplex_execution_begin(&p.execution,
			SIMPLEX_EXECUTION_SOLVE, o.algorithm, 0, 0);
		owner = simplex_execution_ref(&p.execution, root);
		simplex_execution_set_parent(&p.execution, root);
		if (root)
			history.events[root - 1].grant = o.iteration_limit;
		if (mode == 2)
			for (i = history.count; i < SIMPLEX_EXECUTION_CAPACITY;
				i++)
				(void)simplex_execution_begin(&p.execution,
					SIMPLEX_EXECUTION_SOLVE, 0, 0, 0);
		if (mode == 3)
			context.execution_serial = ULONG_MAX;
		(void)simplex_dual_solve_problem(&p, &o, x, NULL, &result, 0);
		CHECK(result.status != lp_simplex_Infeasibility &&
			result.status != lp_simplex_Unboundedness &&
			result.iterations >= 0 &&
			result.iterations <= o.iteration_limit);
		if (!mode) {
			backends = work = selected = retry = 0;
			for (i = root; i < history.count; i++) {
				const struct simplex_ExecutionEvent *e =
					history.events + i;
				if (e->kind == SIMPLEX_EXECUTION_BACKEND) {
					backends++;
					work += e->work;
					retry = i + 1;
					CHECK(e->parent == root &&
						e->work >= 0);
				}
				if (e->kind ==
					SIMPLEX_EXECUTION_POINT_SELECTION)
					selected = i + 1;
			}
			CHECK(backends >= 1 && work == result.iterations);
			if (selected) {
				CHECK(report.backend_retries && backends == 2 &&
					work > 0);
				CHECK(history.events[selected - 1].work ==
					work);
				CHECK(history.events[selected - 1].parent ==
						root &&
					history.events[selected - 1].cause ==
						retry);
				point = history.events[selected - 1]
						.point_source;
				CHECK(point > 0 &&
					history.events[point - 1].kind ==
						SIMPLEX_EXECUTION_BACKEND &&
					history.events[point - 1].complete);
			} else
				/* A future backend may certify this bounded fixture;
				 * no failed-attempt point selection is then needed. */
				CHECK(backends == 1 ||
					result.status == lp_simplex_Success);
			reference = result;
			memcpy(saved, x, sizeof(saved));
		} else
			CHECK(result.status == reference.status &&
				result.iterations == reference.iterations &&
				!memcmp(saved, x, sizeof(saved)));
		root = simplex_execution_ref_id(&owner, &p.execution, 0);
		simplex_execution_end(p.execution.history, root, &result, NULL);
		simplex_solve_context_destroy(&context);
	}
	simplex_problem_destroy(&p);
	lp_model_free(m);
	return 0;
}

int main(void)
{
	int invalid, tracing;
	CHECK(!snapshots() && !storage_limit() && !backend_selection());
	for (invalid = 0; invalid < 2; invalid++)
		for (tracing = 0; tracing < 4; tracing++)
			CHECK(!numerical(invalid, tracing));
	return 0;
}
