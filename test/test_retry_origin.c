/**
 * @file test_retry_origin.c
 * @brief Numerical failure memory survives lost diagnostics; actor references
 * do not.
 */

#include "simplex_stage.h"
#include "simplex_solve_context.h"
#include "simplex_execution.h"
#include "simplex_retry.h"
#include "simplex_dual.h"
#include "simplex_basis.h"
#include <lp_simplex/lp_simplex.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#define CHECK(c) \
	do { \
		if (!(c)) { \
			fprintf(stderr, "retry origin line %d case=%d: %s\n", \
				__LINE__, scenario, #c); \
			return 1; \
		} \
	} while (0)
static int scenario;

enum Case {
	LIVE,
	SCALED,
	RESET_SAME,
	RESET_FOREIGN,
	OTHER_JOURNAL,
	NO_JOURNAL,
	SERIAL_SATURATED,
	UNFINISHED,
	WRONG_STATUS,
	WRONG_WORK,
	WRONG_KIND,
	FOREIGN_AT_CAPTURE,
	CHANGED_MODEL,
	CHANGED_DETAIL,
	CHANGED_STATUS,
	CHANGED_WORK,
	CHANGED_SERIAL,
	CHANGED_CONTEXT,
	CLEARED_COMPLETE,
	CLEARED_RESULT,
	BAD_COUNT,
	BAD_REFERENCE,
	CASE_COUNT
};

static int finish_failure(const struct simplex_ExecutionScope *scope, int kind,
	const struct lp_simplex_Result *r, int finish)
{
	int id = simplex_execution_begin(
		scope, kind, SIMPLEX_BASIS_BACKEND_INTERNAL, 0, 0);
	if (finish)
		simplex_execution_end(scope->history, id, r, NULL);
	return id;
}

static int contract(int seeded)
{
	struct lp_Model *m = lp_model_create(1, 1);
	struct simplex_Problem p;
	struct simplex_SolveContext c, foreign;
	struct simplex_Execution h, other;
	struct simplex_ExecutionScope actor, foreign_scope;
	struct simplex_DualRetry retry;
	struct simplex_DualBasisSeed seed;
	struct simplex_RetryKey key;
	struct lp_simplex_Options o;
	struct lp_simplex_Result r, raw;
	struct simplex_ExecutionEvent *e;
	int index = 0, id, expected, reference;
	unsigned char status[2] = {0, 3};
	CHECK(m);
	m->constraints[0].rhs = m->constraints[0].coef[0] = m->objective[0] =
		1.;
	CHECK(!simplex_problem_from_model(&p, m));
	simplex_solve_context_init(
		&c, scenario == NO_JOURNAL ? NULL : &h, &p.execution);
	simplex_solve_context_configure(&c, 0);
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_DUAL_REVISED);
	simplex_dual_retry_init(&retry, &p, &o, 0, 0, seeded);
	if (seeded) {
		simplex_dual_seed_init(&seed, &p, &index, status);
		seed.valid = 1;
		simplex_dual_retry_bind_basis(&retry, &seed);
	}
	key = retry.key;
	key.backend = SIMPLEX_BASIS_BACKEND_INTERNAL;
	memset(&r, 0, sizeof(r));
	r.status = lp_simplex_PrecisionError;
	r.iterations = 3;
	raw = r;
	actor = p.execution;
	simplex_solve_context_init(&foreign, NULL, &foreign_scope);
	memset(&other, 0, sizeof(other));
	if (scenario == SCALED)
		actor = simplex_solve_context_derive(
			&p.execution, SIMPLEX_EXECUTION_SCALING);
	if (scenario == SERIAL_SATURATED)
		c.execution_serial = ULONG_MAX;
	if (scenario == FOREIGN_AT_CAPTURE) {
		actor = foreign_scope;
		actor.history = &h;
	}
	if (scenario == WRONG_STATUS)
		raw.status = lp_simplex_Singularity;
	if (scenario == WRONG_WORK)
		raw.iterations++;
	id = finish_failure(&actor,
		scenario == WRONG_KIND ? SIMPLEX_EXECUTION_AUTO
				       : SIMPLEX_EXECUTION_BACKEND,
		&raw, scenario != UNFINISHED);
	simplex_retry_record(
		&c.retries, &key, p.execution.history, id, -1, &r, 100, seeded);
	CHECK(c.retries.count == 1 &&
		simplex_retry_failed(&c.retries, &key) == 1);
	if (scenario == UNFINISHED)
		simplex_execution_end(&h, id, &r, NULL);
	if (scenario == RESET_SAME) {
		memset(&h, 0, sizeof(h));
		(void)finish_failure(&actor, SIMPLEX_EXECUTION_BACKEND, &r, 1);
	}
	if (scenario == RESET_FOREIGN) {
		simplex_solve_context_init(&foreign, &h, &foreign_scope);
		(void)finish_failure(
			&foreign_scope, SIMPLEX_EXECUTION_BACKEND, &r, 1);
	}
	if (id > 0) {
		e = h.events + id - 1;
		switch (scenario) {
		case CHANGED_MODEL:
			e->model++;
			break;
		case CHANGED_DETAIL:
			e->detail++;
			break;
		case CHANGED_STATUS:
			e->status = lp_simplex_Singularity;
			break;
		case CHANGED_WORK:
			e->work++;
			break;
		case CHANGED_SERIAL:
			e->serial++;
			break;
		case CHANGED_CONTEXT:
			e->context = &foreign;
			break;
		case CLEARED_COMPLETE:
			e->complete = 0;
			break;
		case CLEARED_RESULT:
			e->result_available = 0;
			break;
		case BAD_COUNT:
			h.count = INT_MAX;
			break;
		case BAD_REFERENCE:
			c.retries.records[0].event = INT_MAX;
			break;
		}
	}
	expected = scenario == LIVE || scenario == SCALED;
	reference = simplex_retry_cause(c.retries.records,
		scenario == OTHER_JOURNAL ? &other : p.execution.history);
	CHECK(reference == (expected ? id : 0));
	CHECK(simplex_retry_failed(&c.retries, &key) == 1 &&
		c.retries.records[0].work == 3);
	CHECK(c.retries.basis_bytes == (seeded ? sizeof(int) + 2 : 0));
	simplex_retry_record(&c.retries, &key, NULL, 0, -1, &r, 1000, seeded);
	CHECK(c.retries.count ==
		1); /* No trace state enters numerical equality. */
	simplex_solve_context_destroy(&c);
	simplex_solve_context_destroy(&foreign);
	simplex_problem_destroy(&p);
	lp_model_free(m);
	return 0;
}

static void replay_prefix(struct simplex_Execution *journal,
	const struct simplex_Execution *old,
	struct simplex_ExecutionScope *owner)
{
	int i;
	memset(journal, 0, sizeof(*journal));
	for (i = 0; i < old->count; i++) {
		int id = simplex_execution_begin(owner, old->events[i].kind,
			old->events[i].detail, 0, 0);
		unsigned long serial = journal->events[id - 1].serial;
		const struct simplex_SolveContext *context =
			journal->events[id - 1].context;
		journal->events[id - 1] = old->events[i];
		journal->events[id - 1].serial = serial;
		journal->events[id - 1].context = context;
	}
}

static int numerical(int mode)
{
	struct lp_Model *m = lp_model_create(1, 1);
	struct simplex_Problem p;
	struct simplex_SolveContext c, foreign;
	struct simplex_Execution h, old;
	struct simplex_ExecutionScope foreign_scope;
	struct simplex_StageReport report;
	struct lp_simplex_Options o;
	struct lp_simplex_Result first, second;
	double x[1];
	int calls, begin, i, skips = 0;
	CHECK(m);
	m->constraints[0].type = optm_CONS_T_EQ;
	m->constraints[0].rhs = 1.;
	m->constraints[0].coef[0] = 1e-300;
	m->bounds[0].b_type = optm_BOUND_T_FR;
	CHECK(!simplex_problem_from_model(&p, m));
	simplex_solve_context_init(&c, mode ? &h : NULL, &p.execution);
	simplex_solve_context_configure(&c, 0);
	c.controls.disable_scaling = 1;
	p.report = &report;
	simplex_solve_context_init(&foreign, NULL, &foreign_scope);
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_DUAL_REVISED);
	o.iteration_limit = 100;
	CHECK(simplex_dual_solve_problem(&p, &o, x, NULL, &first, 0) != 0);
	CHECK(first.status == lp_simplex_PrecisionError ||
		first.status == lp_simplex_Singularity);
	calls = c.retries.count;
	CHECK(calls > 0);
	if (mode == 2) {
		old = h;
		replay_prefix(&h, &old, &p.execution);
	}
	if (mode == 3) {
		old = h;
		foreign_scope.history = &h;
		replay_prefix(&h, &old, &foreign_scope);
	}
	begin = mode ? h.count : 0;
	CHECK(simplex_dual_solve_problem(&p, &o, x, NULL, &second, 0) != 0);
	CHECK(second.status == first.status && second.iterations == 0 &&
		report.recovery_skipped == calls);
	CHECK(c.retries.count == calls && !isfinite(x[0]) &&
		!report.point_available);
	for (i = begin; mode && i < h.count; i++) {
		const struct simplex_ExecutionEvent *e = h.events + i;
		CHECK(e->kind != SIMPLEX_EXECUTION_BACKEND);
		if (e->kind == SIMPLEX_EXECUTION_RECOVERY_SKIP) {
			skips++;
			CHECK(!e->result_available && !e->work);
			if (mode == 1)
				CHECK(e->cause > 0);
			else
				CHECK(e->cause == 0);
		}
	}
	if (mode)
		CHECK(skips == calls);
	simplex_solve_context_destroy(&c);
	simplex_solve_context_destroy(&foreign);
	simplex_problem_destroy(&p);
	lp_model_free(m);
	return 0;
}

static int basis_owner(void)
{
	struct simplex_SolveContext a, b;
	struct simplex_Execution h;
	struct simplex_ExecutionScope one, two;
	struct simplex_DualBasisSeed seed;
	struct lp_simplex_Result r;
	int id;
	simplex_solve_context_init(&a, &h, &one);
	memset(&seed, 0, sizeof(seed));
	memset(&r, 0, sizeof(r));
	r.status = lp_simplex_ExceedIterLimit;
	id = finish_failure(&one, SIMPLEX_EXECUTION_BACKEND, &r, 1);
	simplex_execution_set_parent(&one, id);
	h.events[id - 1].basis_exported = 1;
	simplex_dual_seed_export(&seed, &one);
	CHECK(simplex_dual_seed_source_event(&seed, &one) == id);
	simplex_solve_context_init(&b, &h, &two);
	id = finish_failure(&two, SIMPLEX_EXECUTION_BACKEND, &r, 1);
	h.events[id - 1].basis_exported = 1;
	CHECK(h.events[id - 1].serial ==
		seed.export_serial); /* Natural per-context collision. */
	CHECK(!simplex_dual_seed_source_event(&seed, &one) && seed.valid &&
		seed.exported);
	simplex_solve_context_destroy(&a);
	simplex_solve_context_destroy(&b);
	return 0;
}

int main(void)
{
	int seeded, mode;
	for (scenario = 0; scenario < CASE_COUNT; scenario++)
		for (seeded = 0; seeded < 2; seeded++)
			if (contract(seeded))
				return 1;
	for (mode = 0; mode < 4; mode++)
		if (numerical(mode))
			return 1;
	return basis_owner();
}
