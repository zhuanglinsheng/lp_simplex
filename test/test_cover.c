/**
 * @file test_cover.c
 * @brief Regression tests for cover.
 */

#include "simplex_stage.h"
#include "simplex_solve_context.h"
#include "simplex_execution.h"
#include "simplex_cover.h"
#include "simplex_dual.h"
#include "simplex_column_handoff.h"
#include "simplex_auto.h"
#include "simplex_progress.h"
#include <lp_simplex/lp_simplex.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define CHECK(c) \
	do { \
		if (!(c)) { \
			printf("cover test line %d\n", __LINE__); \
			return 0; \
		} \
	} while (0)

static int certified_cover(void)
{
	struct lp_Model *m = lp_model_create(2, 3);
	struct simplex_Problem p, huge;
	struct lp_simplex_Options o;
	struct lp_simplex_Result r;
	double x[3], y[2];
	int i;
	CHECK(m);
	for (i = 0; i < 2; i++) {
		m->constraints[i].type = optm_CONS_T_GE;
		m->constraints[i].rhs = 1.;
		m->constraints[i].coef[i] = 1.;
		m->constraints[i].coef[2] = 1.;
	}
	for (i = 0; i < 3; i++) {
		m->bounds[i].b_type = optm_BOUND_T_BS;
		m->bounds[i].lb = 0.;
		m->bounds[i].ub = 1.;
		m->objective[i] = i < 2 ? 1. : 1.5;
	}
	CHECK(simplex_problem_from_model(&p, m) == 0);
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_AUTO);
	/* Initial per-row columns cost 2. Full pricing must discover the shared
	 * column and certify 1.5 on the complete LP, not stop at the first master. */
	CHECK(simplex_cover_solve(&p, &o, x, y, &r) == 0);
	CHECK(r.status == lp_simplex_Success && fabs(r.objective - 1.5) < 1e-9);
	CHECK(x[2] > .999 && r.iterations <= o.iteration_limit);
	CHECK(simplex_problem_primal_infeasibility(&p, x) <=
		o.primal_tolerance);
	CHECK(simplex_problem_dual_infeasibility(
		      &p, x, y, o.primal_tolerance, -1) <= o.dual_tolerance);
	o.iteration_limit = 1;
	CHECK(simplex_cover_solve(&p, &o, x, y, &r) != 0);
	CHECK(r.status == lp_simplex_ExceedIterLimit && r.iterations <= 1);
	o.iteration_limit = 1000;
	/* All-zero RHS: no selected columns are needed. */
	p.rhs[0] = p.rhs[1] = 0.;
	CHECK(simplex_cover_solve(&p, &o, x, y, &r) == 0);
	CHECK(r.objective == 0. && x[0] == 0. && x[1] == 0. && x[2] == 0.);
	/* Ineligible models must not receive a mathematical terminal status. */
	p.objective[0] = -1.;
	CHECK(simplex_cover_solve(&p, &o, x, y, &r) != 0);
	CHECK(r.status == lp_simplex_CondUnsatisfied && r.iterations == 0);
	p.objective[0] = 1.;
	p.row_type[0] = optm_CONS_T_EQ;
	CHECK(simplex_cover_solve(&p, &o, x, y, &r) != 0);
	CHECK(r.status == lp_simplex_CondUnsatisfied);
	p.row_type[0] = optm_CONS_T_GE;
	huge = p;
	huge.rows = INT_MAX;
	CHECK(simplex_cover_solve(&huge, &o, x, y, &r) != 0);
	CHECK(r.status == lp_simplex_CondUnsatisfied);
	simplex_problem_destroy(&p);
	lp_model_free(m);
	return 1;
}

static int insufficient_single_column(void)
{
	struct lp_Model *m = lp_model_create(1, 2);
	struct simplex_Problem p;
	struct lp_simplex_Options o;
	struct lp_simplex_Result r;
	double x[2], y[1];
	int j;
	CHECK(m);
	m->constraints[0].type = optm_CONS_T_GE;
	m->constraints[0].rhs = 1.;
	for (j = 0; j < 2; j++) {
		m->constraints[0].coef[j] = m->objective[j] = 1.;
		m->bounds[j].b_type = optm_BOUND_T_BS;
		m->bounds[j].lb = 0.;
		m->bounds[j].ub = .6;
	}
	CHECK(simplex_problem_from_model(&p, m) == 0);
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_DUAL_REVISED);
	CHECK(simplex_cover_solve(&p, &o, x, y, &r) != 0);
	CHECK(r.status == lp_simplex_CondUnsatisfied && r.iterations == 0);
	/* Combined capacity is sufficient, so rejecting the initializer must not
	 * incorrectly mark the original LP infeasible. */
	CHECK(simplex_dual_solve_problem(&p, &o, x, y, &r, 0) == 0);
	CHECK(fabs(r.objective - 1.) < 1e-9);
	simplex_problem_destroy(&p);
	lp_model_free(m);
	return 1;
}

static int bounded_heap_and_checkpoint(void)
{
	struct lp_Model *m = lp_model_create(1, 6);
	struct simplex_Problem p;
	struct lp_simplex_Options o;
	struct lp_simplex_Result r;
	double costs[] = {10., 1., .5, .2, .1, 0.};
	double x[6], y[1];
	int j;
	CHECK(m);
	m->constraints[0].type = optm_CONS_T_GE;
	m->constraints[0].rhs = 1.;
	for (j = 0; j < 6; j++) {
		m->constraints[0].coef[j] = 1.;
		m->objective[j] = costs[j];
		m->bounds[j].b_type = optm_BOUND_T_BS;
		m->bounds[j].lb = 0.;
		m->bounds[j].ub = j == 0 ? 1. : (j == 5 ? 0. : .25);
	}
	CHECK(simplex_problem_from_model(&p, m) == 0);
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_AUTO);
	o.iteration_limit = 3;
	/* Capacity-one heap must retain the best of several negative columns.
	 * Exhausting the next pricing budget retains a feasible incumbent, never
	 * certifies that incomplete master or the fixed zero-capacity column. */
	CHECK(simplex_cover_solve(&p, &o, x, y, &r) != 0);
	CHECK(r.status == lp_simplex_ExceedIterLimit && r.iterations <= 3);
	CHECK(fabs(r.objective - 7.525) < 1e-9 && x[4] == .25 && x[5] == 0.);
	CHECK(simplex_problem_primal_infeasibility(&p, x) <=
		o.primal_tolerance);
	o.iteration_limit = 1000;
	CHECK(simplex_cover_solve(&p, &o, x, y, &r) == 0);
	CHECK(fabs(r.objective - .45) < 1e-9 && x[5] == 0.);
	CHECK(simplex_problem_dual_infeasibility(
		      &p, x, y, o.primal_tolerance, -1) <= o.dual_tolerance);
	simplex_problem_destroy(&p);
	lp_model_free(m);
	return 1;
}

static int always_due(void *context)
{
	(void)context;
	return 1;
}

static void count_cover(void *context, const double *x, const double *y,
	double scale, const char *stage, int iterations, int force,
	const struct simplex_ProgressSource *source)
{
	(void)x;
	(void)y;
	(void)scale;
	(void)iterations;
	(void)force;
	(void)source;
	if (strcmp(stage, "cover-columns") == 0)
		(*(int *)context)++;
}

static int auto_route_and_fallback(void)
{
	struct lp_Model *m = lp_model_create(2, 64);
	struct simplex_Problem p;
	struct simplex_ProgressSink sink;
	struct lp_simplex_Options o;
	struct lp_simplex_Result r;
	double x[64];
	int i, j, reports = 0;
	CHECK(m);
	for (i = 0; i < 2; i++) {
		m->constraints[i].type = optm_CONS_T_GE;
		m->constraints[i].rhs = 1.;
		m->constraints[i].coef[i] = 1.;
		m->constraints[i].coef[63] = 1.;
	}
	for (j = 0; j < 64; j++) {
		m->objective[j] = j < 2 ? 1. : (j == 63 ? 1.5 : 10.);
		m->bounds[j].b_type = optm_BOUND_T_BS;
		m->bounds[j].lb = 0.;
		m->bounds[j].ub = 1.;
		if (j >= 2 && j < 63)
			m->constraints[j % 2].coef[j] = 1.;
	}
	CHECK(simplex_problem_from_model(&p, m) == 0);
	sink.context = &reports;
	sink.due = always_due;
	sink.report = count_cover;
	p.progress = &sink;
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_AUTO);
	o.iteration_limit = 100;
	CHECK(simplex_auto_solve(&p, &o, x, &r, 0, NULL, NULL, NULL) == 0);
	CHECK(reports > 0 && r.status == lp_simplex_Success &&
		fabs(r.objective - 1.5) < 1e-9);
	CHECK(r.iterations <= 100);
	reports = 0;
	o.iteration_limit = 1;
	(void)simplex_auto_solve(&p, &o, x, &r, 0, NULL, NULL, NULL);
	CHECK(reports == 0 && r.iterations <= 1);
	o.iteration_limit = 6;
	CHECK(simplex_auto_solve(&p, &o, x, &r, 0, NULL, NULL, NULL) == 0);
	CHECK(reports == 0 && r.iterations <= 6 &&
		fabs(r.objective - 1.5) < 1e-9);
	simplex_problem_destroy(&p);
	lp_model_free(m);
	/* This wide LP is feasible only by combining smaller capacities. The
	 * optional initializer refuses it, then the established dual route solves. */
	m = lp_model_create(1, 64);
	CHECK(m);
	m->constraints[0].type = optm_CONS_T_GE;
	m->constraints[0].rhs = 1.;
	for (j = 0; j < 64; j++) {
		m->objective[j] = m->constraints[0].coef[j] = 1.;
		m->bounds[j].b_type = optm_BOUND_T_BS;
		m->bounds[j].lb = 0.;
		m->bounds[j].ub = .05;
	}
	CHECK(simplex_problem_from_model(&p, m) == 0);
	p.progress = &sink;
	o.iteration_limit = 100;
	CHECK(simplex_auto_solve(&p, &o, x, &r, 0, NULL, NULL, NULL) == 0);
	CHECK(reports == 0 && fabs(r.objective - 1.) < 1e-9 &&
		r.iterations <= 100);
	simplex_problem_destroy(&p);
	lp_model_free(m);
	return 1;
}

/* An individually sub-tolerance reduced cost can carry a large objective
 * improvement over a wide domain. Check both explicit and covering-implied
 * capacities; a restricted master is not a full objective certificate. */
static int wide_domain_objective_accuracy(int bounded)
{
	struct lp_Model *m = lp_model_create(2, 3);
	struct simplex_Problem p;
	struct lp_simplex_Options o;
	struct lp_simplex_Result r;
	double x[3], y[2];
	int i;
	CHECK(m);
	for (i = 0; i < 2; i++) {
		m->constraints[i].type = optm_CONS_T_GE;
		m->constraints[i].rhs = 1e10;
		m->constraints[i].coef[i] = 1.;
		m->constraints[i].coef[2] = 1.;
	}
	for (i = 0; i < 3; i++) {
		m->bounds[i].b_type =
			bounded ? optm_BOUND_T_BS : optm_BOUND_T_LO;
		m->bounds[i].lb = 0.;
		m->bounds[i].ub = 1e10;
		m->objective[i] = i < 2 ? 1e-8 : 1.5e-8;
	}
	CHECK(simplex_problem_from_model(&p, m) == 0);
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_AUTO);
	CHECK(simplex_cover_solve(&p, &o, x, y, &r) == 0);
	CHECK(r.status == lp_simplex_Success &&
		fabs(r.objective - 150.) < 1e-7);
	CHECK(simplex_problem_primal_infeasibility(&p, x) <=
		o.primal_tolerance);
	CHECK(simplex_problem_dual_infeasibility(
		      &p, x, y, o.primal_tolerance, -1) <= o.dual_tolerance);
	/* The first master costs 200. Exhausting its pricing budget must retain
	 * that feasible point as a limit result, never call it optimal. */
	o.iteration_limit = 3;
	CHECK(simplex_cover_solve(&p, &o, x, y, &r) != 0);
	CHECK(r.status == lp_simplex_ExceedIterLimit && r.iterations <= 3);
	CHECK(fabs(r.objective - 200.) < 1e-7);
	simplex_problem_destroy(&p);
	lp_model_free(m);
	return 1;
}

static int many_small_reduced_costs(void)
{
	struct lp_Model *m = lp_model_create(2, 1026);
	struct simplex_Problem p;
	struct lp_simplex_Options o;
	struct lp_simplex_Result result;
	double x[1026], y[2];
	int i, j;
	CHECK(m);
	for (i = 0; i < 2; i++) {
		m->constraints[i].type = optm_CONS_T_GE;
		m->constraints[i].rhs = 1e10;
		m->constraints[i].coef[i] = 1.;
		for (j = 2; j < 1026; j++)
			m->constraints[i].coef[j] = 1.;
	}
	for (j = 0; j < 1026; j++) {
		m->bounds[j].b_type = optm_BOUND_T_BS;
		m->bounds[j].lb = 0.;
		m->bounds[j].ub = 1e10;
		m->objective[j] = j < 2 ? 1e-8 : 2e-8 - 1e-12;
	}
	CHECK(simplex_problem_from_model(&p, m) == 0);
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_AUTO);
	/* Tiny reduced costs still have material objective impact. Admission
	 * filtering must find them; success requires the complete objective gap. */
	CHECK(simplex_cover_solve(&p, &o, x, y, &result) == 0);
	CHECK(fabs(result.objective - 199.99) < 1e-6);
	CHECK(result.primal_infeasibility <= o.primal_tolerance &&
		result.dual_infeasibility <= o.dual_tolerance);
	simplex_problem_destroy(&p);
	lp_model_free(m);
	return 1;
}

static int traced_cover(void)
{
	struct lp_Model *m = lp_model_create(2, 3);
	struct simplex_Problem p;
	struct simplex_SolveContext context;
	struct simplex_Execution history;
	struct simplex_StageReport report;
	struct lp_simplex_Options o;
	struct lp_simplex_Result r;
	double x[3], y[2];
	int i, j, diagnostics, masters, restrictions, models, first_work = -1;
	CHECK(m);
	for (i = 0; i < 2; i++) {
		m->constraints[i].type = optm_CONS_T_GE;
		m->constraints[i].rhs = 1.;
		m->constraints[i].coef[i] = m->constraints[i].coef[2] = 1.;
	}
	for (j = 0; j < 3; j++) {
		m->objective[j] = j < 2 ? 1. : 1.5;
		m->bounds[j].b_type = optm_BOUND_T_BS;
		m->bounds[j].lb = 0.;
		m->bounds[j].ub = 1.;
	}
	CHECK(!simplex_problem_from_model(&p, m));
	p.report = &report;
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_AUTO);
	for (diagnostics = 0; diagnostics < 3; diagnostics++) {
		simplex_solve_context_init(
			&context, diagnostics ? &history : NULL, &p.execution);
		if (diagnostics == 2)
			for (i = 0; i < SIMPLEX_EXECUTION_CAPACITY; i++)
				(void)simplex_execution_begin(&p.execution,
					SIMPLEX_EXECUTION_SOLVE, 0, 0, 0);
		CHECK(!simplex_cover_solve(&p, &o, x, y, &r));
		CHECK(r.objective == 1.5 && report.point_available &&
			report.stop == SIMPLEX_STOP_OPTIMAL);
		if (!diagnostics)
			first_work = r.iterations;
		else
			CHECK(r.iterations == first_work);
		if (diagnostics != 1)
			continue;
		masters = restrictions = models = 0;
		for (i = 0; i < history.count; i++) {
			const struct simplex_ExecutionEvent *e =
				history.events + i;
			CHECK(e->complete && e->model > 0 && e->parent <= i);
			if (e->kind == SIMPLEX_EXECUTION_COLUMN_MODEL) {
				models++;
				CHECK(e->source_model == p.execution.model &&
					e->model != p.execution.model);
				CHECK(e->model_rows == 2 &&
					e->model_columns >= 2);
			}
			if (e->kind == SIMPLEX_EXECUTION_COVER_MASTER) {
				masters++;
				CHECK(e->work <= e->grant &&
					e->result_available &&
					e->status == lp_simplex_Success);
				CHECK(e->reason ==
						SIMPLEX_REASON_COLUMN_ADMISSION ||
					e->reason ==
						SIMPLEX_REASON_DUAL_REFINEMENT);
			}
			if (e->kind == SIMPLEX_EXECUTION_COLUMN_HANDOFF) {
				CHECK(!e->handoff_rejection &&
					!e->result_available);
				if (e->detail ==
					SIMPLEX_COLUMN_BASIS_RESTRICT) {
					restrictions++;
					CHECK(e->cause > 0 && e->cause <= i);
					CHECK(history.events[e->cause - 1]
							.kind ==
						SIMPLEX_EXECUTION_COLUMN_HANDOFF);
					CHECK(history.events[e->cause - 1]
							.detail ==
						SIMPLEX_COLUMN_BASIS_LIFT);
					CHECK(e->parent > 0 &&
						history.events[e->parent - 1]
								.kind ==
							SIMPLEX_EXECUTION_COVER_MASTER);
				}
			}
		}
		CHECK(masters >= 2 && models == masters &&
			restrictions == masters - 1);
		simplex_solve_context_destroy(&context);
	}
	simplex_problem_destroy(&p);
	lp_model_free(m);
	return 1;
}

/* Only cover-owned roots/cursors are injected. Real master models, Dual and
 * numerical backends retain their normal library allocation paths. */
static void *cover_live[5];
static int cover_requests, cover_failure, cover_active, cover_bad_free;

static void *cover_fault_malloc(size_t bytes)
{
	void *p;
	int i;
	cover_requests++;
	if (cover_requests == cover_failure)
		return NULL;
	p = malloc(bytes ? bytes : 1);
	if (!p)
		return NULL;
	memset(p, 0xa5, bytes);
	for (i = 0; i < 5; i++)
		if (!cover_live[i]) {
			cover_live[i] = p;
			cover_active++;
			return p;
		}
	free(p);
	return NULL;
}

static void cover_fault_free(void *p)
{
	int i;
	if (!p)
		return;
	for (i = 0; i < 5; i++)
		if (cover_live[i] == p) {
			cover_live[i] = NULL;
			cover_active--;
			free(p);
			return;
		}
	cover_bad_free++;
}

#define lp_simplex_malloc cover_fault_malloc
#define lp_simplex_free cover_fault_free
#define simplex_cover_solve cover_fault_solve
#include "../src/dual/simplex_cover.c"
#undef lp_simplex_malloc
#undef lp_simplex_free
#undef simplex_cover_solve

static int resource_exits(void)
{
	struct lp_Model *m = lp_model_create(2, 3);
	struct simplex_Problem p;
	struct simplex_SolveContext context;
	struct simplex_Execution h;
	struct simplex_StageReport report;
	struct lp_simplex_Options o;
	struct lp_simplex_Result r;
	double x[3], y[2];
	int tracing, failure, total = 0, i, code;
	CHECK(m);
	for (i = 0; i < 2; i++) {
		m->constraints[i].type = optm_CONS_T_GE;
		m->constraints[i].rhs = 1.;
		m->constraints[i].coef[i] = m->constraints[i].coef[2] = 1.;
	}
	for (i = 0; i < 3; i++) {
		m->objective[i] = i < 2 ? 1. : 1.5;
		m->bounds[i].b_type = optm_BOUND_T_BS;
		m->bounds[i].ub = 1.;
	}
	CHECK(!simplex_problem_from_model(&p, m));
	p.report = &report;
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_AUTO);
	for (tracing = 0; tracing < 3; tracing++)
		for (failure = 0; failure <= total; failure++) {
			simplex_solve_context_init(
				&context, tracing ? &h : NULL, &p.execution);
			simplex_solve_context_configure(&context, 0);
			if (tracing == 2)
				while (h.count < SIMPLEX_EXECUTION_CAPACITY)
					(void)simplex_execution_begin(
						&p.execution,
						SIMPLEX_EXECUTION_SOLVE, 0, 0,
						0);
			cover_requests = 0;
			cover_failure = failure;
			memset(&report, 0xff, sizeof(report));
			x[0] = x[1] = x[2] = y[0] = y[1] = 17.;
			code = cover_fault_solve(&p, &o, x, y, &r);
			CHECK(!cover_active && !cover_bad_free);
			if (!failure) {
				if (!total)
					total = cover_requests;
				CHECK(total >= 6 && total <= 16 &&
					cover_requests == total);
				CHECK(!code && r.status == lp_simplex_Success &&
					fabs(r.objective - 1.5) < 1e-9 &&
					report.point_available);
			} else {
				CHECK(code &&
					r.status ==
						lp_simplex_MemoryAllocError &&
					report.stop == SIMPLEX_STOP_RESOURCE);
				CHECK(cover_requests ==
					(failure < 4 ? 4 : failure));
				CHECK(!report.basis_available &&
					!report.factor_available &&
					y[0] == 17. && y[1] == 17.);
				if (failure <= 5)
					CHECK(!r.iterations &&
						!report.point_available);
				else
					CHECK(r.iterations > 0 &&
						report.point_available &&
						simplex_problem_primal_infeasibility(
							&p, x) <=
							o.primal_tolerance);
			}
			simplex_solve_context_destroy(&context);
		}
	cover_failure = 0;
	simplex_problem_destroy(&p);
	lp_model_free(m);
	return 1;
}

int main(void)
{
	if (!resource_exits() || !certified_cover() ||
		!insufficient_single_column() ||
		!bounded_heap_and_checkpoint() || !auto_route_and_fallback() ||
		!wide_domain_objective_accuracy(1) ||
		!wide_domain_objective_accuracy(0) ||
		!many_small_reduced_costs() || !traced_cover())
		return 1;
	puts("cover pricing certifies the full LP and preserves fallback/budget semantics");
	return 0;
}
