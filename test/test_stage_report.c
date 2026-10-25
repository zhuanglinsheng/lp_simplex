/**
 * @file test_stage_report.c
 * @brief Evidence, retention and checkpoint contracts, independent of
 * benchmarks.
 */

#include "simplex_stage.h"
#include "simplex_auto.h"
#include "simplex_problem.h"
#include "simplex_dual.h"
#include "simplex_progress.h"
#include <lp_simplex/lp_simplex.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#define CHECK(c) \
	do { \
		if (!(c)) { \
			fprintf(stderr, "stage line %d: %s\n", __LINE__, #c); \
			return 1; \
		} \
	} while (0)

static int calls, mode, initializer_mode;

static int reported_dual(const struct simplex_Problem *p,
	const struct lp_simplex_Options *o, double *x, double *y,
	struct lp_simplex_Result *r, int propagate)
{
	(void)o;
	(void)y;
	(void)propagate;
	calls++;
	x[0] = calls == 1 ? 0. : .75;
	x[1] = calls == 1 ? 0. : .25;
	r->iterations = 1;
	r->status = lp_simplex_ExceedIterLimit;
	if (mode == 1 && calls == 1) {
		r->iterations = -1;
		r->status = lp_simplex_Success;
		r->dual_infeasibility = 0.;
		x[0] = 1.;
		x[1] = 0.;
	}
	if (mode == 2 && calls == 1) {
		p->report->lower_bound.value = -1.;
		p->report->lower_bound.kind = SIMPLEX_OBSERVATION_LOWER_BOUND;
		p->report->lower_bound.space = SIMPLEX_SPACE_INPUT;
		p->report->lower_bound.quality = SIMPLEX_EVIDENCE_VERIFIED;
		p->report->lower_bound.work = 0; /* stale at work 1 */
	}
	/* Deliberately no stop report: one step does not imply a yield. */
	CHECK(p->report->stop == SIMPLEX_STOP_UNKNOWN);
	return -1;
}

static int local_resource(const struct simplex_Problem *p,
	const struct lp_simplex_Options *o, const double *seed, double *x,
	struct lp_simplex_Result *r)
{
	(void)o;
	(void)seed;
	x[0] = x[1] = .5;
	r->iterations = 1;
	r->status = lp_simplex_MemoryAllocError;
	p->report->stop = SIMPLEX_STOP_RESOURCE;
	if (mode == 2) {
		p->report->lower_bound.value = -2.;
		p->report->lower_bound.kind = SIMPLEX_OBSERVATION_LOWER_BOUND;
		p->report->lower_bound.space = SIMPLEX_SPACE_INPUT;
		p->report->lower_bound.quality = SIMPLEX_EVIDENCE_VERIFIED;
		p->report->lower_bound.work = 1;
	}
	return -1;
}

static void count_notifications(
	void *context, const struct lp_simplex_Progress *event)
{
	int *count = context;
	(void)event;
	(*count)++;
}

static int partial_initializer(void *context,
	const struct simplex_Problem *input, const struct lp_simplex_Options *o,
	double *x, int *work)
{
	(void)input;
	(void)context;
	(void)o;
	*work = 1;
	if (initializer_mode == 0)
		return 0; /* no output at all */
	x[0] = .5;
	x[1] = initializer_mode == 1 ? .5 : NAN;
	return initializer_mode == 2;
}

int main(void)
{
	struct simplex_StageReport report;
	struct simplex_AutoHistory history;
	struct simplex_AutoKernels kernels;
	struct simplex_AutoBindings bindings = {NULL, NULL, 0};
	struct simplex_ProgressTracker tracker;
	struct simplex_DualCheckpoint checkpoint;
	struct simplex_Problem p;
	struct lp_simplex_Options o;
	struct lp_simplex_Result result;
	struct lp_Model *m = lp_model_create(1, 2);
	int indices[1], status;
	unsigned char membership[3];
	double candidate[2] = {1., 0.}, worse[2] = {0., 1.};
	double invalid[2] = {NAN, 0.}, out[2];

	bindings.kernels = &kernels;
	CHECK(m);
	m->constraints[0].type = optm_CONS_T_EQ;
	m->constraints[0].rhs = 1.;
	m->constraints[0].coef[0] = m->constraints[0].coef[1] = 1.;
	m->objective[1] = 1.;
	CHECK(!simplex_problem_from_model(&p, m));
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_AUTO);
	{
		struct simplex_ObservationHistory progress;
		struct simplex_Observation sample;
		memset(&progress, 0, sizeof(progress));
		memset(&sample, 0, sizeof(sample));
		sample.value = 100.;
		sample.kind = SIMPLEX_OBSERVATION_PRIMAL_SUM;
		sample.space = SIMPLEX_SPACE_WORKING;
		sample.quality = SIMPLEX_EVIDENCE_RECOMPUTED;
		sample.epoch = 1;
		CHECK(simplex_stage_progress_observe(&progress, &sample,
			      sample.work, .01,
			      0.) == SIMPLEX_OBSERVATION_BASELINE);
		sample.value = 90.;
		sample.work = 128;
		CHECK(simplex_stage_progress_observe(&progress, &sample,
			      sample.work, .01,
			      0.) == SIMPLEX_OBSERVATION_IMPROVED);
		sample.value = 200.;
		sample.work = 256;
		CHECK(simplex_stage_progress_observe(&progress, &sample,
			      sample.work, .01,
			      0.) == SIMPLEX_OBSERVATION_UNCHANGED);
		sample.value = 150.;
		sample.work = 384;
		CHECK(simplex_stage_progress_observe(&progress, &sample,
			      sample.work, .01,
			      0.) == SIMPLEX_OBSERVATION_UNCHANGED);
		sample.value = 50.;
		sample.work = 300;
		CHECK(simplex_stage_progress_observe(&progress, &sample,
			      sample.work, .01,
			      0.) == SIMPLEX_OBSERVATION_STALE);
		sample.value = 0.;
		sample.work = 512;
		sample.quality = SIMPLEX_EVIDENCE_ESTIMATED;
		CHECK(simplex_stage_progress_observe(&progress, &sample,
			      sample.work, .01,
			      0.) == SIMPLEX_OBSERVATION_UNAVAILABLE);
		CHECK(progress.best.value == 90.);
		/* A recomputed sample can still be old relative to the live state. */
		sample.work = 500;
		sample.quality = SIMPLEX_EVIDENCE_RECOMPUTED;
		CHECK(simplex_stage_progress_observe(&progress, &sample, 512,
			      .01, 0.) == SIMPLEX_OBSERVATION_STALE);
		sample.work = 512;

		sample.quality = SIMPLEX_EVIDENCE_RECOMPUTED;
		sample.epoch = 2;
		CHECK(simplex_stage_progress_observe(&progress, &sample,
			      sample.work, .01,
			      0.) == SIMPLEX_OBSERVATION_BASELINE);
		sample.kind = SIMPLEX_OBSERVATION_LOWER_BOUND;
		sample.space = SIMPLEX_SPACE_INPUT;
		sample.quality = SIMPLEX_EVIDENCE_VERIFIED;
		sample.value = -100.;
		sample.work = 640;
		CHECK(simplex_stage_progress_observe(&progress, &sample,
			      sample.work, .01,
			      0.) == SIMPLEX_OBSERVATION_BASELINE);
		sample.value = -50.;
		sample.work = 768;
		CHECK(simplex_stage_progress_observe(&progress, &sample,
			      sample.work, .01,
			      0.) == SIMPLEX_OBSERVATION_IMPROVED);
		CHECK(simplex_stage_observation_current(&sample,
			SIMPLEX_SPACE_INPUT, 2, 768,
			SIMPLEX_EVIDENCE_VERIFIED));
		CHECK(!simplex_stage_observation_current(&sample,
			SIMPLEX_SPACE_INPUT, 2, 769,
			SIMPLEX_EVIDENCE_VERIFIED));
	}
	simplex_stage_reset(&report);
	CHECK(report.stop == SIMPLEX_STOP_UNKNOWN &&
		isnan(report.primal.value));
	CHECK(!simplex_stage_observation_valid(&report.primal,
		SIMPLEX_SPACE_INPUT, 0, SIMPLEX_EVIDENCE_ESTIMATED));
	report.primal.value = 0.;
	report.primal.quality = SIMPLEX_EVIDENCE_RECOMPUTED;
	report.primal.space = SIMPLEX_SPACE_WORKING;
	report.primal.epoch = 3;
	CHECK(!simplex_stage_observation_valid(&report.primal,
		SIMPLEX_SPACE_INPUT, 3, SIMPLEX_EVIDENCE_RECOMPUTED));
	CHECK(!simplex_stage_observation_valid(&report.primal,
		SIMPLEX_SPACE_WORKING, 4, SIMPLEX_EVIDENCE_RECOMPUTED));
	CHECK(simplex_stage_observation_valid(&report.primal,
		SIMPLEX_SPACE_WORKING, 3, SIMPLEX_EVIDENCE_RECOMPUTED));
	report.stop = SIMPLEX_STOP_NUMERICAL;
	CHECK(!simplex_stage_oracle_retry(&report));
	report.failure = SIMPLEX_FAILURE_FACTORIZATION;
	CHECK(!simplex_stage_oracle_retry(&report));
	report.failure = SIMPLEX_FAILURE_IMPLICIT_ORACLE;
	CHECK(simplex_stage_oracle_retry(&report));
	report.stop = SIMPLEX_STOP_BUDGET;
	CHECK(!simplex_stage_oracle_retry(&report));

	{
		struct simplex_Problem working = p;
		double rhs[1] = {0.}, zero[2] = {0., 0.};
		working.rhs = rhs;
		simplex_stage_observe(&p, &o, zero, NULL, 1, 10, &report);
		CHECK(report.primal.value == 1. &&
			report.feasibility_merit.value == 1.);
		CHECK(report.violation_count.value == 1.);
		CHECK(report.feasible_objective.quality ==
			SIMPLEX_EVIDENCE_UNKNOWN);
		CHECK(report.dual.quality == SIMPLEX_EVIDENCE_UNKNOWN);
		CHECK(report.lower_bound.quality == SIMPLEX_EVIDENCE_UNKNOWN);
		CHECK(simplex_problem_primal_infeasibility(&working, zero) ==
			0.);
		CHECK(simplex_problem_primal_infeasibility(&p, zero) == 1.);
		simplex_stage_observe(&p, &o, invalid, NULL, 1, 11, &report);
		CHECK(isnan(report.primal.value) &&
			isnan(report.violation_count.value));
		CHECK(report.primal.quality == SIMPLEX_EVIDENCE_UNKNOWN);
	}
	CHECK(simplex_stage_prefer_point(&p, &o, candidate, worse));
	CHECK(simplex_stage_prefer_point(&p, &o, candidate, invalid));
	CHECK(!simplex_stage_prefer_point(&p, &o, invalid, candidate));
	CHECK(!simplex_stage_prefer_point(&p, &o, worse, candidate));

	/* Archive exists with no callback, and keeps its independently checked
	 * point and conservative bound when a later observation is invalid. */
	CHECK(!o.progress_callback);
	CHECK(!simplex_progress_init(&tracker, m, &o));
	tracker.sink.report(
		tracker.sink.context, candidate, NULL, 1., "test", 1, 1, NULL);
	CHECK(tracker.have_primal && tracker.best_primal == 0.);
	tracker.sink.report(
		tracker.sink.context, invalid, NULL, 1., "failed", 2, 1, NULL);
	CHECK(tracker.have_primal && tracker.point[0] == 1.);
	simplex_progress_destroy(&tracker);
	{
		int notifications = 0;
		o.progress_callback = count_notifications;
		o.progress_context = &notifications;
		o.progress_interval = 1e9;
		CHECK(!simplex_progress_init(&tracker, m, &o));
		CHECK(notifications == 1);
		p.progress = &tracker.sink;
		simplex_progress_observe(
			&p, candidate, NULL, 1., "trusted-audit", 10);
		CHECK(tracker.have_primal && tracker.best_primal == 0.);
		CHECK(notifications == 1);
		p.progress = NULL;
		simplex_progress_destroy(&tracker);
		o.progress_callback = NULL;
		o.progress_context = NULL;
	}

	simplex_dual_checkpoint_init(&checkpoint, &p, indices, membership);
	o.iteration_limit = 0;
	p.report = &report;
	simplex_stage_reset(&report);
	status = simplex_dual_solve_problem_checkpoint(
		&p, &o, out, NULL, &result, 0, &checkpoint);
	CHECK(status != 0 && result.status == lp_simplex_ExceedIterLimit);
	CHECK(report.stop == SIMPLEX_STOP_BUDGET && checkpoint.basis.valid);
	CHECK(report.basis_available && !report.factor_available);
	CHECK(report.work.basis_counts_available &&
		report.work.factorizations > 0);
	o.iteration_limit = 100;
	simplex_stage_reset(&report);
	status = simplex_dual_solve_problem_checkpoint(
		&p, &o, out, NULL, &result, 0, &checkpoint);
	CHECK(status == 0 && result.status == lp_simplex_Success);
	CHECK(checkpoint.basis.accepted && report.stop == SIMPLEX_STOP_OPTIMAL);
	CHECK(simplex_stage_observation_valid(&report.primal,
		SIMPLEX_SPACE_INPUT, 0, SIMPLEX_EVIDENCE_RECOMPUTED));
	CHECK(report.primal.value <= o.primal_tolerance);
	/* A local resource failure cannot prevent later progress. Unknown short
	 * returns stay unknown, and failed routes can improve the saved point. */
	kernels.dual = reported_dual;
	kernels.hybrid = local_resource;
	calls = 0;
	CHECK(simplex_auto_solve(
		      &p, &o, out, &result, 0, &bindings, NULL, &history) != 0);
	CHECK(calls == 2 && history.count == 3 && result.iterations == 3);
	CHECK(history.stages[0].report.stop == SIMPLEX_STOP_UNKNOWN);
	CHECK(history.stages[1].report.stop == SIMPLEX_STOP_RESOURCE);
	CHECK(history.stages[2].report.stop == SIMPLEX_STOP_UNKNOWN);
	CHECK(history.have_point && history.best_objective == .25 &&
		out[1] == .25);
	CHECK(history.stages[1].action.cause == 0);
	CHECK(history.stages[2].action.cause == 1);
	CHECK(history.stages[2].action.failure_action ==
		SIMPLEX_AUTO_NO_RECOVERY);
	CHECK(history.stages[0].action.failure_action ==
		SIMPLEX_AUTO_NEXT_ROUTE);
	CHECK(history.stages[2].action.point_source ==
		-1); /* dual ignores primal seed */
	CHECK(history.stages[0].report.primal.quality ==
		SIMPLEX_EVIDENCE_UNKNOWN);
	CHECK(history.stages[0].primal_audit.value == 1.);
	for (mode = 1; mode <= 2; mode++) {
		calls = 0;
		CHECK(simplex_auto_solve(&p, &o, out, &result, 0, &bindings,
			      NULL, &history) != 0);
		if (mode == 1) {
			CHECK(history.stages[0].kernel_status ==
				lp_simplex_Success);
			CHECK(history.stages[0].reported_work == -1);
			CHECK(history.stages[0].work ==
				history.stages[0].action.grant);
			CHECK(history.stages[0].validation ==
				SIMPLEX_AUTO_TERMINAL_INVALID);
			CHECK(!history.stages[0].accepted_point);
			CHECK(history.stages[0].primal_audit.value == 0.);
			CHECK(history.stages[0].report.failure ==
				SIMPLEX_FAILURE_STAGE_CONTRACT);
			CHECK(history.best_objective == .25 && out[1] == .25);
		} else {
			CHECK(history.best_bound == -2. &&
				history.bound_source == 1);
		}
	}
	mode = 0;
	for (initializer_mode = 0; initializer_mode <= 2; initializer_mode++) {
		struct simplex_AutoInitializer initializer;
		initializer.context = NULL;
		initializer.run = partial_initializer;
		calls = 0;
		CHECK(simplex_auto_solve(&p, &o, out, &result, 0, &bindings,
			      &initializer, &history) != 0);
		CHECK(history.count == 4 && result.iterations == 4);
		CHECK(history.stages[1].action.stage ==
			SIMPLEX_AUTO_INITIALIZER);
		CHECK(history.stages[1].action.cause == 0);
		CHECK(history.stages[2].action.cause == 0);
		CHECK(history.stages[2].action.point_source ==
			(initializer_mode == 1 ? 1 : -1));
		CHECK(history.stages[2].action.purpose ==
			(initializer_mode == 1 ? SIMPLEX_AUTO_OPTIMIZE_POINT
					       : SIMPLEX_AUTO_SEEK_FEASIBLE));
		CHECK(history.stages[1].validation ==
			(initializer_mode == 1
					? SIMPLEX_AUTO_POINT_VALID
					: SIMPLEX_AUTO_VALIDATION_UNKNOWN));
		CHECK(out[1] == .25);
	}
	simplex_problem_destroy(&p);
	lp_model_free(m);
	/* Cover is a real production action, not an unrecorded prelude to the
	 * journal. Validate the full input and account for its work once. */
	{
		struct lp_Model *wide = lp_model_create(1, 40);
		struct simplex_Problem covering;
		double point[40];
		int j;
		CHECK(wide);
		wide->constraints[0].type = optm_CONS_T_GE;
		wide->constraints[0].rhs = 1.;
		for (j = 0; j < 40; j++) {
			wide->constraints[0].coef[j] = 1.;
			wide->objective[j] = j + 1.;
		}
		CHECK(!simplex_problem_from_model(&covering, wide));
		CHECK(!simplex_auto_solve(&covering, &o, point, &result, 0,
			NULL, NULL, &history));
		CHECK(history.count == 1);
		CHECK(history.stages[0].action.stage == SIMPLEX_AUTO_COVER);
		CHECK(history.stages[0].validation ==
			SIMPLEX_AUTO_TERMINAL_VALID);
		CHECK(history.stages[0].work == result.iterations);
		CHECK(history.stages[0].primal_audit.value <=
			o.primal_tolerance);
		simplex_problem_destroy(&covering);
		lp_model_free(wide);
	}
	return 0;
}
