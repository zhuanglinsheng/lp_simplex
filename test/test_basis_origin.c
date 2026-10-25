/**
 * @file test_basis_origin.c
 * @brief Basis output and point provenance are independent, with actual Dual
 * exports.
 */

#include "simplex_solve_context.h"
#include "simplex_execution.h"
#include "simplex_auto.h"
#include "simplex_dual.h"
#include <lp_simplex/lp_simplex.h>
#include <math.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>

#define CHECK(c) \
	do { \
		if (!(c)) { \
			fprintf(stderr, "basis origin line %d: %s\n", \
				__LINE__, #c); \
			return 1; \
		} \
	} while (0)
static int offer_point;
static int reset_history;

static int hybrid(const struct simplex_Problem *p,
	const struct lp_simplex_Options *o, const double *seed, double *x,
	struct lp_simplex_Result *r)
{
	int j;
	(void)o;
	(void)seed;
	if (reset_history && p->execution.history) {
		int count = p->execution.parent;
		memset(p->execution.history, 0, sizeof(*p->execution.history));
		for (j = 0; j < count; j++)
			(void)simplex_execution_begin(&p->execution,
				SIMPLEX_EXECUTION_AUTO, SIMPLEX_AUTO_HYBRID, 0,
				0);
	}
	for (j = 0; j < p->columns; j++)
		x[j] = offer_point ? 1. : 0.;
	r->iterations = 1;
	r->status = lp_simplex_PrecisionError;
	return -1;
}

static int controller(int tracing, int scaled, int point)
{
	struct lp_Model *m = lp_model_create(20, 20);
	struct simplex_Problem p;
	struct simplex_SolveContext context;
	struct simplex_Execution trace;
	struct simplex_AutoHistory h;
	struct simplex_AutoKernels k;
	struct simplex_AutoBindings bindings = {NULL, NULL, 0};
	struct lp_simplex_Options o;
	struct lp_simplex_Result r;
	double x[20];
	int i, producer;
	bindings.kernels = &k;
	CHECK(m);
	offer_point = point;
	reset_history = tracing == 3;
	for (i = 0; i < 20; i++) {
		double scale = scaled ? (i % 2 ? 1048576. : .0009765625) : 1.;
		m->constraints[i].type = optm_CONS_T_GE;
		m->constraints[i].rhs = scale;
		m->constraints[i].coef[i] = scale;
		m->objective[i] = 1.;
	}
	CHECK(!simplex_problem_from_model(&p, m));
	simplex_solve_context_init(
		&context, tracing ? &trace : NULL, &p.execution);
	if (tracing == 2)
		for (i = 0; i < SIMPLEX_EXECUTION_CAPACITY; i++)
			(void)simplex_execution_begin(
				&p.execution, SIMPLEX_EXECUTION_SOLVE, 0, 0, 0);
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_AUTO);
	o.iteration_limit = 16;
	k.dual = simplex_dual_solve_problem;
	k.hybrid = hybrid;
	CHECK(simplex_auto_solve(&p, &o, x, &r, 0, &bindings, NULL, &h) != 0);
	CHECK(h.count == 3 && r.iterations == 16 &&
		r.status == lp_simplex_ExceedIterLimit);
	CHECK(h.stages[0].basis_exported && !h.stages[0].basis_accepted);
	CHECK(h.stages[2].basis_exported && h.stages[2].basis_accepted);
	CHECK(h.stages[2].action.basis_source == 0 &&
		h.stages[2].action.cause == 0);
	CHECK(h.stages[2].action.reason == SIMPLEX_REASON_BASIS_REBUILD);
	CHECK(h.stages[2].action.point_source ==
		-1); /* No primal seed to cold/warm Dual. */
	CHECK(h.stages[1].action.stage == SIMPLEX_AUTO_HYBRID);
	if (point) {
		CHECK(h.have_point && h.point_source == 1 &&
			h.outcome.point.record == 1);
		CHECK(r.objective == 20. && h.stages[1].accepted_point &&
			!h.stages[2].accepted_point);
		for (i = 0; i < 20; i++)
			CHECK(x[i] == 1.);
	} else
		CHECK(!h.have_point);
	producer = h.stages[0].basis_producer_event;
	if (tracing == 1) {
		const struct simplex_ExecutionEvent *e;
		int next = h.stages[2].basis_producer_event;
		CHECK(producer > 0 && next > producer &&
			producer != h.stages[0].action.execution_id);
		e = trace.events + producer - 1;
		CHECK(e->kind == SIMPLEX_EXECUTION_BACKEND &&
			e->basis_exported && e->complete);
		CHECK(h.stages[2].action.basis_event == producer &&
			h.stages[2].action.cause_event == producer);
		CHECK(trace.events[next - 1].basis_source == producer);
		CHECK(trace.events[next - 1].basis_import_attempted &&
			trace.events[next - 1].basis_exported);
		CHECK(!trace.events[next - 1].basis_rejection);
	} else if (tracing == 3) {
		CHECK(!producer && !h.stages[2].action.basis_event);
		CHECK(h.stages[2].basis_producer_event > 0);
	} else
		CHECK(!producer && !h.stages[2].action.basis_event &&
			!h.stages[2].basis_producer_event);
	simplex_solve_context_destroy(&context);
	simplex_problem_destroy(&p);
	lp_model_free(m);
	return 0;
}

static int source_contract(void)
{
	struct simplex_DualBasisSeed s;
	struct simplex_SolveContext c;
	struct simplex_Execution h, other;
	struct simplex_ExecutionScope scope, observer;
	struct lp_simplex_Result r;
	int id;
	simplex_solve_context_init(&c, &h, &scope);
	memset(&s, 0, sizeof(s));
	memset(&r, 0, sizeof(r));
	r.status = lp_simplex_ExceedIterLimit;
	id = simplex_execution_begin(
		&scope, SIMPLEX_EXECUTION_BACKEND, 0, 0, 0);
	simplex_execution_set_parent(&scope, id);
	simplex_dual_seed_export(&s, &scope);
	CHECK(s.valid && s.exported &&
		!simplex_dual_seed_source_event(&s, &scope));
	simplex_execution_end(&h, id, &r, NULL);
	h.events[id - 1].basis_exported = 1;
	CHECK(simplex_dual_seed_source_event(&s, &scope) == id);
	s.accepted = s.import_attempted = 1;
	s.rejection = SIMPLEX_BASIS_LAYOUT_REJECTED;
	simplex_dual_seed_begin(&s);
	CHECK(s.valid && !s.exported && !s.accepted && !s.import_attempted &&
		!s.rejection);
	CHECK(simplex_dual_seed_source_event(&s, &scope) ==
		id); /* Retained producer, no new export. */
	observer = scope;
	observer.history = &other;
	CHECK(!simplex_dual_seed_source_event(&s, &observer));
	observer = scope;
	observer.context = NULL;
	CHECK(!simplex_dual_seed_source_event(&s, &observer));
	h.events[id - 1].basis_exported = 0;
	CHECK(!simplex_dual_seed_source_event(&s, &scope) && s.valid);
	h.events[id - 1].basis_exported = 1;
	h.events[id - 1].model++;
	CHECK(!simplex_dual_seed_source_event(&s, &scope) && s.valid);
	/* Reinitializing the same journal buffer cannot make an old actor id
	 * refer to a new, otherwise identical backend record. */
	memset(&h, 0, sizeof(h));
	id = simplex_execution_begin(
		&scope, SIMPLEX_EXECUTION_BACKEND, 0, 0, 0);
	simplex_execution_end(&h, id, &r, NULL);
	h.events[id - 1].basis_exported = 1;
	CHECK(!simplex_dual_seed_source_event(&s, &scope) && s.valid);
	/* A late export cannot relabel the retained numerical basis as output
	 * of the replacement actor occupying the old parent's integer slot. */
	simplex_dual_seed_export(&s, &scope);
	CHECK(s.valid && s.exported && !s.export_serial &&
		!simplex_dual_seed_source_event(&s, &scope));
	c.execution_serial = ULONG_MAX;
	id = simplex_execution_begin(
		&scope, SIMPLEX_EXECUTION_BACKEND, 0, 0, 0);
	simplex_execution_set_parent(&scope, id);
	simplex_dual_seed_export(&s, &scope);
	simplex_execution_end(&h, id, &r, NULL);
	h.events[id - 1].basis_exported = 1;
	CHECK(!simplex_dual_seed_source_event(&s, &scope) && s.valid &&
		s.exported);
	s.valid = 0;
	simplex_dual_seed_begin(&s);
	CHECK(!s.export_source.parent && !s.export_source.history);
	simplex_solve_context_destroy(&c);
	return 0;
}

int main(void)
{
	int tracing, scaled, point;
	if (source_contract())
		return 1;
	for (tracing = 0; tracing < 4; tracing++)
		for (scaled = 0; scaled < 2; scaled++)
			for (point = 0; point < 2; point++)
				if (controller(tracing, scaled, point))
					return 1;
	return 0;
}
