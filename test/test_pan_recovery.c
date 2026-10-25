/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file test_pan_recovery.c
 * @brief Regression tests for Pan feasible-point recovery and acceptance
 * checks.
 */

#include "simplex_pan_recovery.h"
#include "simplex_pan_service.h"
#include "simplex_pan_config.h"
#include "simplex_pan_cone.h"
#include "simplex_pan_basis.h"
#include "simplex_solve_context.h"
#include "simplex_execution.h"
#include "simplex_stage.h"
#include "utils.h"

#include <lp_simplex/lp_simplex.h>

#include <math.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>

/* Only this numerical owner's initial workspace is refused. Basis and QR
 * allocations remain delegated and are not claimed as covered by this seam. */
static void *refuse_restore_storage(size_t bytes)
{
	(void)bytes;
	return NULL;
}

#define lp_simplex_malloc refuse_restore_storage
#define simplex_pan_progress_threshold refused_progress_threshold
#define simplex_pan_recovery_compress refused_recovery_compress
#define simplex_pan_restore_configured_reported refused_restore
#include "../src/bda/simplex_pan_recovery.c"
#undef simplex_pan_restore_configured_reported
#undef simplex_pan_recovery_compress
#undef simplex_pan_progress_threshold
#undef lp_simplex_malloc

static int restore_cases(void)
{
	struct simplex_PanStandard a;
	struct lp_simplex_Options o;
	struct simplex_PanConfig cfg;
	int starts[3] = {0, 1, 2}, rows[2] = {0, 1}, used, state;
	double values[2] = {1., 1e12}, norm[2] = {1., 1e12}, unit[2] = {1., 1.};
	double rhs[2] = {1., 1.}, cost[2] = {1., 0.}, point[2], saved[2];
	int stop;

	memset(&a, 0, sizeof(a));
	a.rows = 2;
	a.columns = 2;
	a.nonzeros = 2;
	a.column_start = starts;
	a.row_index = rows;
	a.value = values;
	a.normalized_value = unit;
	a.column_norm = norm;
	a.rhs = rhs;
	a.objective = cost;
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_PAN_BDA);
	simplex_pan_config_init(&cfg, &o);
	o.primal_tolerance = 1e-7;
	o.iteration_limit = 100;
	point[0] = 1. - 9e-8;
	point[1] = (1. + 8e-8) * 1e-12;
	memcpy(saved, point, sizeof(point));
	used = 7;
	if (refused_restore(&a, point, &o, 1.L, &used, &cfg, &stop) != -1 ||
		stop != SIMPLEX_STOP_RESOURCE || used != 7 ||
		memcmp(point, saved, sizeof(point)))
		return 0;
	used = 0;
	state = simplex_pan_service_restore_standard(
		&a, point, &o, 1.L, &used, &cfg, NULL, NULL);
	if (state != 1 || fabs(point[0] - 1.) > 1e-14 ||
		fabs(point[1] * 1e12 - 1.) > 1e-14 || used <= 0 || used > 100)
		return 0;
	/* Repair must use original costs even when working costs are shifted. */
	a.audit_objective = cost;
	{
		double shifted[2] = {-1e18, -1e18};

		a.objective = shifted;
		point[0] = 1. - 9e-8;
		point[1] = 1e-12;
		memcpy(saved, point, sizeof(point));
		used = 0;
		if (simplex_pan_service_restore_standard(&a, point, &o, .5L,
			    &used, &cfg, NULL, NULL) != -1 ||
			memcmp(point, saved, sizeof(point)))
			return 0;
	}
	a.objective = cost;
	a.audit_objective = NULL;
	/* No remaining budget and an invalid negative coordinate are failures
	 * with exact rollback; neither is silently clipped. */
	point[0] = 1. - 9e-8;
	point[1] = 1e-12;
	memcpy(saved, point, sizeof(point));
	used = 100;
	if (simplex_pan_service_restore_standard(
		    &a, point, &o, 1.L, &used, &cfg, NULL, NULL) != -1 ||
		used != 100 || memcmp(point, saved, sizeof(point)))
		return 0;
	point[1] = -1e-20;
	memcpy(saved, point, sizeof(point));
	used = 0;
	if (simplex_pan_service_restore_standard(
		    &a, point, &o, 1.L, &used, &cfg, NULL, NULL) != -1 ||
		memcmp(point, saved, sizeof(point)))
		return 0;
	/* An empty positive face cannot repair an inconsistent equality. */
	point[0] = point[1] = 0.;
	used = 0;
	if (simplex_pan_service_restore_standard(
		    &a, point, &o, 1.L, &used, &cfg, NULL, NULL) != -1 ||
		point[0] != 0. || point[1] != 0.)
		return 0;
	/* A correction crossing the nonnegative boundary is shortened. It must
	 * still reduce the full original residual enough, or roll back. */
	rhs[0] = -1e-6;
	point[0] = 1e-8;
	point[1] = 1e-12;
	memcpy(saved, point, sizeof(point));
	used = 0;
	return simplex_pan_service_restore_standard(
		       &a, point, &o, 1.L, &used, &cfg, NULL, NULL) == -1 &&
	       !memcmp(point, saved, sizeof(point));
}

static int expanding_face(void)
{
	struct simplex_PanStandard a;
	struct lp_simplex_Options o;
	struct simplex_PanConfig cfg;
	int starts[71], rows[70], j, used = 0;
	double values[70], norm[70], rhs[70], cost[70], point[70];

	memset(&a, 0, sizeof(a));
	for (j = 0; j < 70; j++) {
		starts[j] = rows[j] = j;
		values[j] = norm[j] = rhs[j] = point[j] = 1.;
		cost[j] = 0.;
	}
	starts[70] = 70;
	point[69] -= 9e-8;
	a.rows = a.columns = a.nonzeros = 70;
	a.column_start = starts;
	a.row_index = rows;
	a.value = a.normalized_value = values;
	a.column_norm = norm;
	a.rhs = rhs;
	a.objective = cost;
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_PAN_BDA);
	simplex_pan_config_init(&cfg, &o);
	o.primal_tolerance = 1e-7;
	o.iteration_limit = 200;
	return simplex_pan_service_restore_standard(
		       &a, point, &o, 0.L, &used, &cfg, NULL, NULL) == 1 &&
	       used > 65 && used <= 200 && fabs(point[69] - 1.) < 1e-14;
}

/* Exercise the real trial -> residual repair -> commit -> dual audit path.
 * The tangent direction preserves an initial residual near the tolerance;
 * repair must regain headroom without increasing the original objective. */
static int repair_in_solver(void)
{
	struct simplex_PanStandard a;
	struct simplex_PanBasis b;
	struct lp_simplex_Options o;
	struct simplex_PanConfig cfg;
	struct simplex_StageReport report;
	int starts[4] = {0, 1, 2, 4}, rows[4] = {0, 1, 0, 1};
	double values[4] = {1., 1., 1., 1.},
	       normalized[4] = {1., 1., .7071067811865475, .7071067811865475};
	double norms[3] = {1., 1., 1.4142135623730951}, rhs[2] = {1., 1.},
	       cost[3] = {1., 0., 0.};
	double point[3] = {1. - 9e-8, 1., 0.}, pi[2], dual, residual[2],
	       work[2];
	int used = 0, status = -1, result, good;

	memset(&a, 0, sizeof(a));
	a.rows = 2;
	a.columns = 3;
	a.nonzeros = 4;
	a.column_start = starts;
	a.row_index = rows;
	a.value = values;
	a.normalized_value = normalized;
	a.column_norm = norms;
	a.rhs = rhs;
	a.objective = cost;
	simplex_stage_reset(&report);
	a.report = &report;
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_PAN_BDA);
	simplex_pan_config_init(&cfg, &o);
	o.primal_tolerance = 1e-7;
	o.dual_tolerance = 1e-9;
	o.iteration_limit = 1000;
	if (simplex_pan_basis_create_configured(
		    &b, &a, o.pivot_tolerance, &cfg))
		return 0;
	result = simplex_pan_cone_phase2_configured(
		&a, &b, point, pi, &o, &used, &status, &dual, &cfg);
	simplex_pan_basis_destroy(&b);
	good = result == 0 && status == lp_simplex_Success && used <= 1000 &&
	       point[0] >= 0. && point[1] >= 0. && point[2] >= 0. &&
	       point[0] < 1e-9 &&
	       simplex_pan_residual(&a, point, residual, work) < 1e-12 &&
	       report.services.calls[SIMPLEX_PAN_SERVICE_RESTORE] > 0 &&
	       report.services.output[SIMPLEX_PAN_SERVICE_RESTORE] ==
		       SIMPLEX_PAN_OUTPUT_NONNEGATIVE_EQUALITY_POINT;
	/* Exhaustion during a trial/repair cannot commit an extra LP step. */
	point[0] = 1. - 9e-8;
	point[1] = 1.;
	point[2] = 0.;
	used = 0;
	o.iteration_limit = 2;
	if (simplex_pan_basis_create_configured(
		    &b, &a, o.pivot_tolerance, &cfg))
		return 0;
	result = simplex_pan_cone_phase2_configured(
		&a, &b, point, pi, &o, &used, &status, &dual, &cfg);
	simplex_pan_basis_destroy(&b);
	return good && result != 0 && status == lp_simplex_ExceedIterLimit &&
	       used <= 2 && point[0] == 1. - 9e-8 && point[1] == 1. &&
	       point[2] == 0.;
}

/* Real repair through the shared lease protocol, under independent diagnostic
 * capacities. Mathematical work/output must not depend on journal storage. */
static int restore_service_contracts(void)
{
	struct lp_Model *m = lp_model_create(2, 2);
	struct simplex_Problem p;
	struct simplex_PanStandard a;
	struct simplex_PanConfig cfg;
	struct simplex_PanServiceContext services;
	struct simplex_PanServiceOutcome out;
	struct simplex_PanServiceTrigger trigger;
	struct simplex_SolveContext context;
	struct simplex_Execution history;
	struct simplex_StageReport report;
	struct lp_simplex_Options o;
	double point[2], saved[2], residual[2], scratch[2], reference[2];
	int mode, scenario, work, code, i, cause, good = 1, reference_work = 0;

	if (!m)
		return 0;
	m->constraints[0].type = m->constraints[1].type = optm_CONS_T_EQ;
	m->constraints[0].rhs = m->constraints[1].rhs = 1.;
	m->constraints[0].coef[0] = 1.;
	m->constraints[1].coef[1] = 1e12;
	m->objective[0] = 1.;
	if (simplex_problem_from_model(&p, m)) {
		lp_model_free(m);
		return 0;
	}
	p.report = &report;
	if (simplex_pan_standard_create(&a, &p)) {
		simplex_problem_destroy(&p);
		lp_model_free(m);
		return 0;
	}
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_PAN_BDA);
	o.primal_tolerance = 1e-7;
	simplex_pan_config_init(&cfg, &o);
	for (mode = 0; mode < 4 && good; mode++)
		for (scenario = 0; scenario < 6 && good; scenario++) {
			simplex_stage_reset(&report);
			simplex_solve_context_init(&context,
				mode == 1 ? NULL : &history, &p.execution);
			simplex_pan_service_context_init(&services, &p);
			simplex_pan_service_bind(&a, &services);
			if (mode == 3)
				context.execution_serial = ULONG_MAX;
			if (mode == 2)
				for (i = history.count;
					i < SIMPLEX_EXECUTION_CAPACITY; i++)
					(void)simplex_execution_begin(
						&p.execution,
						SIMPLEX_EXECUTION_SOLVE, 0, 0,
						0);
			point[0] = scenario == 1 ? 1. : 1. - 9e-8;
			point[1] = scenario == 3 ? -1e-20 : 1e-12;
			memcpy(saved, point, sizeof(point));
			work = 7;
			o.iteration_limit = scenario == 2 ? 7 : 100;
			memset(&trigger, 0, sizeof(trigger));
			trigger.kind = SIMPLEX_PAN_TRIGGER_PRIMAL_MARGIN;
			trigger.work = work;
			trigger.epoch = 3;
			trigger.primal = simplex_pan_residual(
				&a, point, residual, scratch);
			trigger.previous_primal = HUGE_VAL;
			if (scenario == 5)
				simplex_problem_mark_changed(&p);
			code = simplex_pan_service_restore_standard(&a, point,
				&o, scenario == 4 ? .5L : 1.L, &work, &cfg,
				&trigger, &out);
			good = out.service == SIMPLEX_PAN_SERVICE_RESTORE &&
			       out.before == 7 && out.end == work &&
			       work >= 7 && work <= o.iteration_limit &&
			       out.raw_status == lp_simplex_CondUnsatisfied &&
			       report.services.calls
					       [SIMPLEX_PAN_SERVICE_RESTORE] ==
				       1 &&
			       report.services.work
					       [SIMPLEX_PAN_SERVICE_RESTORE] ==
				       work - 7;
			if (scenario == 0) {
				good = good && code == 1 && out.valid &&
				       out.stop == SIMPLEX_STOP_HANDOFF &&
				       out.output ==
					       SIMPLEX_PAN_OUTPUT_NONNEGATIVE_EQUALITY_POINT &&
				       work > 7 &&
				       simplex_pan_residual(&a, point, residual,
					       scratch) < 1e-12;
				if (!mode) {
					reference_work = work;
					memcpy(reference, point, sizeof(point));
				} else
					good = good && work == reference_work &&
					       !memcmp(reference, point,
						       sizeof(point));
			} else {
				good = good &&
				       !memcmp(saved, point, sizeof(point)) &&
				       out.output == SIMPLEX_PAN_OUTPUT_NONE;
				if (scenario == 1)
					good = good && !code && out.valid &&
					       work == 7 &&
					       out.stop ==
						       SIMPLEX_STOP_NOT_APPLICABLE;
				else if (scenario == 2)
					good = good && code == -1 &&
					       out.valid && work == 7 &&
					       out.stop == SIMPLEX_STOP_BUDGET;
				else
					good = good && code == -1 &&
					       out.stop ==
						       SIMPLEX_STOP_NUMERICAL &&
					       out.valid == (scenario != 5);
			}
			good = good && services.failed == (scenario == 5);
			if (!mode && scenario != 5 && out.event > 0) {
				cause = history.events[out.event - 1].cause;
				good = good && cause > 0 &&
				       cause != out.event &&
				       history.events[out.event - 1].reason ==
					       SIMPLEX_REASON_PRIMAL_MARGIN &&
				       history.events[cause - 1].kind ==
					       SIMPLEX_EXECUTION_PAN_TRIGGER &&
				       history.events[cause - 1].complete &&
				       !history.events[cause - 1]
						.result_available &&
				       history.events[cause - 1]
						       .observation_work == 7 &&
				       history.events[cause - 1].primal.value ==
					       trigger.primal;
			} else if (mode)
				good = good && !out.event;
			else if (scenario != 5)
				good = 0;
			simplex_solve_context_destroy(&context);
			if (!good)
				fprintf(stderr,
					"standard restore mode=%d case=%d code=%d stop=%d work=%d valid=%d output=%d\n",
					mode, scenario, code, out.stop, work,
					out.valid, out.output);
		}
	simplex_pan_standard_destroy(&a);
	simplex_problem_destroy(&p);
	lp_model_free(m);
	return good;
}

int main(void)
{
	struct simplex_PanRecovery s;
	double point[2] = {1., 2.};
	int good;

	memset(&s, 0, sizeof(s));
	s.low_progress = 2;
	good = simplex_pan_recovery_compress(&s, point, 2) == 1 &&
	       simplex_pan_recovery_compress(&s, point, 2) == 0 &&
	       s.low_progress == 2;
	point[1] = 3.;
	good = good && simplex_pan_recovery_compress(&s, point, 2) == 1 &&
	       s.compressions == 2 && s.compression_skips == 1;
	good = good && restore_cases() && expanding_face() &&
	       repair_in_solver() && restore_service_contracts();
	printf("Pan bounded recovery and original-model repair %s\n",
		good ? "passed" : "failed");
	return good ? 0 : 1;
}
