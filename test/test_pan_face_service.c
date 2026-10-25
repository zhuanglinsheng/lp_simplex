/**
 * @file test_pan_face_service.c
 * @brief Native-face roles, producer causes and journal-independent solving.
 */

#include "simplex_pan_face_recovery.h"
#include "simplex_pan_service.h"
#include "simplex_pan_config.h"
#include "simplex_stage.h"
#include <lp_simplex/lp_simplex.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>

#define CHECK(c) \
	do { \
		if (!(c)) { \
			fprintf(stderr, "face service line %d: %s\n", \
				__LINE__, #c); \
			return 1; \
		} \
	} while (0)

static int lifecycle(void)
{
	struct lp_Model *m = lp_model_create(1, 2);
	struct simplex_Problem p;
	struct simplex_PanBounded b;
	struct simplex_PanFaceRecovery owner;
	struct simplex_PanServiceContext services;
	struct simplex_PanServiceOutcome out, search, stale;
	struct simplex_PanServiceFrame frame;
	struct simplex_SolveContext context;
	struct simplex_Execution h;
	struct simplex_StageReport report;
	struct lp_simplex_Options o;
	double x[3], y[1] = {0.}, saved[3], r[1], scratch[1];
	int mode, work, root, unrelated, mutation;
	CHECK(m);
	m->constraints[0].type = optm_CONS_T_GE;
	m->constraints[0].rhs = 1.;
	m->constraints[0].coef[0] = m->constraints[0].coef[1] = 1.;
	m->objective[1] = 1.;
	CHECK(!simplex_problem_from_model(&p, m));
	p.report = &report;
	CHECK(!simplex_pan_bounded_create(&b, &p));
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_PAN_BDA);
	o.iteration_limit = 7;
	for (mode = 0; mode < 6; mode++) {
		memset(&owner, 0, sizeof(owner));
		simplex_stage_reset(&report);
		simplex_solve_context_init(&context, &h, &p.execution);
		if (mode == 0)
			p.execution.history = NULL;
		root = simplex_execution_begin(
			&p.execution, SIMPLEX_EXECUTION_SOLVE, 0, 0, 0);
		simplex_execution_set_parent(&p.execution, root);
		simplex_pan_service_context_init(&services, &p);
		simplex_pan_service_bind(&b.a, &services);
		if (mode >= 2 && mode <= 4)
			while (h.count <
				SIMPLEX_EXECUTION_CAPACITY - (mode - 2))
				(void)simplex_execution_begin(&p.execution,
					SIMPLEX_EXECUTION_SOLVE, 0, 0, 0);
		if (mode == 5)
			context.execution_serial = ULONG_MAX;
		x[0] = x[1] = .5;
		x[2] = 0.;
		work = 7;
		CHECK(!simplex_pan_service_face(&owner, &b, x, y, &o, &work,
			saved, r, scratch, SIMPLEX_PAN_FACE_INITIAL, 0, 0, NULL,
			&out));
		CHECK(work == 7 && out.valid && !out.grant &&
			out.stop == SIMPLEX_STOP_REASSESS &&
			out.output == SIMPLEX_PAN_OUTPUT_NONE &&
			out.face.result.action == SIMPLEX_PAN_FACE_NO_CHANGE &&
			out.face.result.remembered && owner.known);
		CHECK(report.services
				.calls[SIMPLEX_PAN_SERVICE_FACE_RECOVERY] == 1);
		CHECK(!report.services.work[SIMPLEX_PAN_SERVICE_FACE_RECOVERY]);
		if (mode == 1 || mode == 4) {
			const struct simplex_ExecutionEvent *e;
			CHECK(out.event > 0);
			e = h.events + out.event - 1;
			CHECK(e->parent == root && e->complete &&
				e->reason == SIMPLEX_REASON_FACE_INITIAL &&
				e->service_valid && !e->work && !e->grant);
			CHECK(e->cause > 0 &&
				h.events[e->cause - 1].kind ==
					SIMPLEX_EXECUTION_PAN_FACE_OPPORTUNITY &&
				h.events[e->cause - 1].complete &&
				h.events[e->cause - 1].parent == root &&
				!h.events[e->cause - 1].result_available);
		} else
			CHECK(!out.event);
		CHECK(simplex_pan_service_begin(&frame, &b.a,
			SIMPLEX_PAN_SERVICE_DUAL_REFINE, work,
			o.iteration_limit, 0));
		CHECK(simplex_pan_service_finish(&frame, &work, 0,
			lp_simplex_CondUnsatisfied, SIMPLEX_STOP_REASSESS,
			SIMPLEX_PAN_OUTPUT_NONE, &search));
		unrelated = simplex_execution_begin(&p.execution,
			SIMPLEX_EXECUTION_POINT_SELECTION, 0, 0, 0);
		simplex_execution_end(
			p.execution.history, unrelated, NULL, NULL);
		saved[0] = r[0] = scratch[0] = 123.;
		/* Poison output-only aggregates; the exact owner still decides. */
		report.services.calls[SIMPLEX_PAN_SERVICE_FACE_RECOVERY] =
			LONG_MAX;
		CHECK(!simplex_pan_service_face(&owner, &b, x, y, &o, &work,
			saved, r, scratch,
			SIMPLEX_PAN_FACE_AFTER_MULTIPLIER_SEARCH, 3, 2, &search,
			&out));
		CHECK(work == 7 && out.valid && !out.output &&
			out.stop == SIMPLEX_STOP_NOT_APPLICABLE &&
			out.face.result.action == SIMPLEX_PAN_FACE_EXACT_SKIP &&
			!out.face.result.changed && owner.skipped == 1);
		CHECK(saved[0] == 123. && r[0] == 123. && scratch[0] == 123.);
		if (mode == 1) {
			const struct simplex_ExecutionEvent *e, *cause;
			CHECK(out.event > 0);
			e = h.events + out.event - 1;
			CHECK(e->cause > 0);
			cause = h.events + e->cause - 1;
			CHECK(e->recovery_skipped && e->parent == root &&
				e->reason ==
					SIMPLEX_REASON_FACE_MULTIPLIER_SEARCH);
			CHECK(cause->kind ==
					SIMPLEX_EXECUTION_PAN_FACE_OPPORTUNITY &&
				cause->cause == search.event &&
				cause->cause != unrelated &&
				cause->face.steps == 3 &&
				cause->face.low_progress == 2 &&
				cause->observation_work == 7);
			/* Unrelated, unavailable or stale diagnostic sources cannot
			 * authorize a cause, nor suppress the numerical action. */
			for (mutation = 0; mutation < 6; mutation++) {
				stale = search;
				if (mutation == 0)
					stale.actor.serial++;
				else if (mutation == 1)
					stale.service =
						SIMPLEX_PAN_SERVICE_CONE;
				else if (mutation == 2)
					stale.input = NULL;
				else if (mutation == 3)
					stale.end = work + 1;
				else if (mutation == 4)
					stale.valid = 0;
				else
					stale.actor.context = NULL;
				CHECK(!simplex_pan_service_face(&owner, &b, x,
					y, &o, &work, saved, r, scratch,
					SIMPLEX_PAN_FACE_AFTER_MULTIPLIER_SEARCH,
					3, 2, &stale, &out));
				CHECK(out.valid && out.event > 0 && work == 7 &&
					out.face.result.action ==
						SIMPLEX_PAN_FACE_EXACT_SKIP);
				e = h.events + out.event - 1;
				CHECK(e->cause > 0 &&
					!h.events[e->cause - 1].cause);
			}
		}
		x[0] = 1.;
		x[1] = 1e-9;
		CHECK(simplex_pan_service_face(&owner, &b, x, y, &o, &work,
			      saved, r, scratch, SIMPLEX_PAN_FACE_LOW_PROGRESS,
			      4, 4, NULL, &out) == 1);
		CHECK(work == 7 && out.valid &&
			out.output ==
				SIMPLEX_PAN_OUTPUT_NATIVE_FACE_ADJUSTMENT &&
			out.face.result.action == SIMPLEX_PAN_FACE_ADJUSTED &&
			out.face.result.changed == 1 &&
			!out.face.result.remembered && x[0] == 1. &&
			x[1] == 0. && !y[0]);
		/* A face action cannot promote an adjustment to a KKT certificate. */
		CHECK(simplex_pan_service_begin(&frame, &b.a,
			SIMPLEX_PAN_SERVICE_FACE_RECOVERY, work,
			o.iteration_limit, 0));
		CHECK(!simplex_pan_service_finish(&frame, &work, 1,
			lp_simplex_Success, SIMPLEX_STOP_OPTIMAL,
			SIMPLEX_PAN_OUTPUT_NATIVE_KKT, &out));
		CHECK(!out.output && services.failed && work == 7);
		simplex_pan_face_recovery_destroy(&owner);
		simplex_solve_context_destroy(&context);
	}
	CHECK(!y[0]);
	simplex_pan_bounded_destroy(&b);
	simplex_problem_destroy(&p);
	lp_model_free(m);
	return 0;
}

static int native_producer(void)
{
	struct lp_Model *m = lp_model_create(10, 24);
	struct simplex_Problem p;
	struct simplex_PanConfig cfg;
	struct lp_simplex_Options o;
	struct lp_simplex_Result result;
	struct simplex_SolveContext context;
	struct simplex_Execution h;
	struct simplex_StageReport report;
	double start[24], out[24], reference[24];
	int i, j, mode, work = -1, faces, linked;
	CHECK(m);
	for (i = 0; i < 10; i++) {
		m->constraints[i].type = optm_CONS_T_EQ;
		m->constraints[i].rhs = 1.;
	}
	for (j = 0; j < 24; j++) {
		m->constraints[j / 3].coef[j] = 1.;
		m->objective[j] = 2 - j % 3;
		m->bounds[j].b_type = optm_BOUND_T_BS;
		m->bounds[j].ub = 1.;
		start[j] = j % 3 == 0 ? 1. : 0.;
		m->constraints[8].coef[j] = m->constraints[0].coef[j];
		m->constraints[9].coef[j] = m->constraints[1].coef[j];
	}
	CHECK(!simplex_problem_from_model(&p, m));
	p.report = &report;
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_PAN_BDA);
	o.iteration_limit = 20000;
	CHECK(simplex_pan_config_default(&cfg, SIMPLEX_PAN_NATIVE_BOUNDS,
		SIMPLEX_PAN_CONE_DIRECTION, SIMPLEX_PAN_EXPLICIT_ORACLE));
	for (mode = 0; mode < 3; mode++) {
		simplex_stage_reset(&report);
		simplex_solve_context_init(&context, &h, &p.execution);
		if (mode == 0)
			p.execution.history = NULL;
		if (mode == 2)
			while (h.count < SIMPLEX_EXECUTION_CAPACITY)
				(void)simplex_execution_begin(&p.execution,
					SIMPLEX_EXECUTION_SOLVE, 0, 0, 0);
		CHECK(!simplex_pan_bounded_solve_certified_configured(
			&p, &o, start, out, NULL, &result, &cfg));
		CHECK(result.status == lp_simplex_Success &&
			result.objective == 0. &&
			result.primal_infeasibility == 0.);
		for (j = 0; j < 24; j++)
			CHECK(out[j] == (j % 3 == 2 ? 1. : 0.));
		if (!mode) {
			work = result.iterations;
			memcpy(reference, out, sizeof(reference));
		} else
			CHECK(result.iterations == work &&
				!memcmp(reference, out, sizeof(reference)));
		CHECK(report.services.calls[SIMPLEX_PAN_SERVICE_FACE_RECOVERY] >
			1);
		faces = linked = 0;
		if (mode == 1) {
			for (i = 0; i < h.count; i++) {
				const struct simplex_ExecutionEvent *e =
					h.events + i;
				if (e->kind != SIMPLEX_EXECUTION_PAN_SERVICE ||
					e->detail !=
						SIMPLEX_PAN_SERVICE_FACE_RECOVERY)
					continue;
				faces++;
				CHECK(e->complete && e->service_valid &&
					!e->work && !e->grant);
				if (e->reason ==
					SIMPLEX_REASON_FACE_MULTIPLIER_SEARCH) {
					const struct simplex_ExecutionEvent
						*cause;
					CHECK(e->cause > 0);
					cause = h.events + e->cause - 1;
					CHECK(cause->cause > 0);
					CHECK(h.events[cause->cause - 1].kind ==
							SIMPLEX_EXECUTION_PAN_SERVICE &&
						h.events[cause->cause - 1]
								.detail ==
							SIMPLEX_PAN_SERVICE_DUAL_REFINE &&
						h.events[cause->cause - 1]
							.complete);
					linked++;
				}
			}
			CHECK(faces > 1 && linked > 0);
		}
		simplex_solve_context_destroy(&context);
	}
	simplex_problem_destroy(&p);
	lp_model_free(m);
	return 0;
}

int main(void)
{
	return lifecycle() || native_producer();
}
