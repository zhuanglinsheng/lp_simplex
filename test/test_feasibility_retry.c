/**
 * @file test_feasibility_retry.c
 * @brief Stable LP input identity across feasibility service stack/diagnostic
 * views.
 */

#include "simplex_stage.h"
#include "simplex_solve_context.h"
#include "simplex_execution.h"
#include "simplex_pan_service.h"
#include "simplex_dual.h"
#include "simplex_progress.h"
#include <lp_simplex/lp_simplex.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#define CHECK(c) \
	do { \
		if (!(c)) { \
			fprintf(stderr, "feasibility retry line %d: %s\n", \
				__LINE__, #c); \
			return 1; \
		} \
	} while (0)

static int nested(struct simplex_PanServiceContext *services,
	const struct simplex_Problem *p, const struct lp_simplex_Options *o,
	double *x, int *work, struct lp_simplex_Result *result, int depth)
{
	volatile unsigned char padding[1024];
	int code;
	padding[0] = 0;
	if (depth)
		code = nested(services, p, o, x, work, result, depth - 1);
	else
		code = simplex_pan_service_dual_feasible(
			services, p, o, x, work, result, NULL);
	/* Keep the caller frame alive; no tail-call or shared wrapper address
	 * can accidentally make an ephemeral child look like stable identity. */
	return code + padding[0];
}

static struct lp_Model *fixture(void)
{
	struct lp_Model *model = lp_model_create(1, 1);
	if (model) {
		model->constraints[0].type = optm_CONS_T_EQ;
		model->constraints[0].rhs = 1.;
		model->constraints[0].coef[0] = 1e-320;
		model->bounds[0].b_type = optm_BOUND_T_FR;
	}
	return model;
}

static int service_identity(int diagnostics)
{
	struct lp_Model *model = fixture();
	struct simplex_Problem p;
	struct simplex_SolveContext context;
	struct simplex_Execution history, alternate;
	struct simplex_PanServiceContext services;
	struct lp_simplex_Options options;
	struct lp_simplex_Result result;
	double x[1];
	int i, round, work = 0, failures = 0, before, backends;
	CHECK(model && !simplex_problem_from_model(&p, model));
	simplex_solve_context_init(
		&context, diagnostics ? &history : NULL, &p.execution);
	simplex_solve_context_configure(&context, 0);
	context.controls.disable_scaling = 1;
	simplex_pan_service_context_init(&services, &p);
	lp_simplex_default_options(&options, lp_simplex_ALGORITHM_PAN_BDA);
	options.iteration_limit = 100;
	if (diagnostics == 2)
		for (i = 0; i < SIMPLEX_EXECUTION_CAPACITY; i++)
			(void)simplex_execution_begin(
				&p.execution, SIMPLEX_EXECUTION_SOLVE, 0, 0, 0);
	for (round = 0; round < 3; round++) {
		if (round == 1 && diagnostics == 3) {
			memset(&alternate, 0, sizeof(alternate));
			p.execution.history = services.scope.history =
				&alternate;
		}
		before = p.execution.history ? p.execution.history->count : 0;
		CHECK(nested(&services, &p, &options, x, &work, &result,
			      round * 2) != 0);
		CHECK(result.status == lp_simplex_PrecisionError ||
			result.status == lp_simplex_Singularity);
		CHECK(work == 0 && result.iterations == 0 && !services.failed);
		if (round == 0)
			failures = context.retries.count;
		CHECK(failures > 0 && context.retries.count == failures);
		for (i = 0; i < failures; i++) {
			CHECK(context.retries.records[i].key.input.input == &p);
			CHECK(context.retries.records[i].key.feasibility_only);
		}
		backends = 0;
		if (p.execution.history)
			for (i = before; i < p.execution.history->count; i++)
				backends +=
					p.execution.history->events[i].kind ==
					SIMPLEX_EXECUTION_BACKEND;
		if (round > 0)
			CHECK(!backends && !isfinite(x[0]));
	}
	/* A new declared revision permits a genuinely new attempt. */
	simplex_problem_mark_changed(&p);
	simplex_pan_service_context_init(&services, &p);
	CHECK(nested(&services, &p, &options, x, &work, &result, 1) != 0);
	CHECK(context.retries.count == 2 * failures);
	simplex_solve_context_destroy(&context);
	simplex_problem_destroy(&p);
	lp_model_free(model);
	return 0;
}

static int diagnostic_binding(void)
{
	struct lp_Model *model = fixture();
	struct simplex_Problem p;
	struct simplex_SolveContext context;
	struct simplex_Execution owner, actual;
	struct simplex_ExecutionScope scope;
	struct simplex_StageReport report;
	struct lp_simplex_Options options;
	struct lp_simplex_Result result;
	double x[1];
	int failures, i, count;
	CHECK(model && !simplex_problem_from_model(&p, model));
	simplex_solve_context_init(&context, &owner, &p.execution);
	simplex_solve_context_configure(&context, 0);
	context.controls.disable_scaling = 1;
	memset(&actual, 0, sizeof(actual));
	scope = p.execution;
	scope.history = &actual;
	lp_simplex_default_options(&options, lp_simplex_ALGORITHM_DUAL_REVISED);
	CHECK(simplex_dual_find_feasible_reported(
		      &p, &options, x, &result, &scope, &report) != 0);
	failures = context.retries.count;
	CHECK(failures > 0 && !owner.count);
	for (i = 0; i < failures; i++) {
		const struct simplex_RetryRecord *record =
			context.retries.records + i;
		CHECK(record->key.input.input == &p &&
			record->history == &actual);
		CHECK(record->event_serial && record->event > 0);
	}
	CHECK(simplex_dual_find_feasible_reported(
		      &p, &options, x, &result, NULL, &report) != 0);
	CHECK(context.retries.count == failures &&
		report.recovery_skipped == failures);
	CHECK(!result.iterations && !isfinite(x[0]));
	/* A diagnostic scope is not permission to solve another mathematical
	 * model under the input's failure key. Reject before any backend. */
	count = actual.count;
	scope.model++;
	x[0] = 123.;
	report.point_available = report.basis_available =
		report.factor_available = 1;
	CHECK(simplex_dual_find_feasible_reported(
		      &p, &options, x, &result, &scope, &report) != 0);
	CHECK(result.status == lp_simplex_PrecisionError &&
		!result.iterations && x[0] == 123.);
	CHECK(report.failure == SIMPLEX_FAILURE_STAGE_CONTRACT &&
		!report.point_available && !report.basis_available &&
		!report.factor_available);
	CHECK(actual.count == count && context.retries.count == failures);
	scope = p.execution;
	scope.context = NULL;
	count = owner.count;
	CHECK(simplex_dual_find_feasible_reported(
		      &p, &options, x, &result, &scope, &report) != 0);
	CHECK(report.failure == SIMPLEX_FAILURE_STAGE_CONTRACT &&
		!result.iterations);
	CHECK(owner.count == count && context.retries.count == failures);
	simplex_solve_context_destroy(&context);
	simplex_problem_destroy(&p);
	lp_model_free(model);
	return 0;
}

struct mutation {
	struct simplex_Problem *input;
	int enabled, changed;
};

static int due(void *context)
{
	(void)context;
	return 1;
}

static void observe(void *context, const double *x, const double *y,
	double scale, const char *stage, int work, int force,
	const struct simplex_ProgressSource *source)
{
	struct mutation *mutation = context;
	(void)x;
	(void)y;
	(void)scale;
	(void)stage;
	(void)work;
	(void)force;
	(void)source;
	if (mutation->enabled && !mutation->changed) {
		simplex_problem_mark_changed(mutation->input);
		mutation->changed = 1;
	}
}

static int input_lifetime(void)
{
	struct lp_Model *model = lp_model_create(1, 2);
	struct simplex_Problem p;
	struct simplex_SolveContext context;
	struct simplex_StageReport report;
	struct simplex_ProgressSink sink;
	struct mutation mutation;
	struct lp_simplex_Options options;
	struct lp_simplex_Result result;
	double x[2];
	int change, code, work = -1;
	CHECK(model);
	model->constraints[0].type = optm_CONS_T_EQ;
	model->constraints[0].rhs = 1.;
	model->constraints[0].coef[0] = model->constraints[0].coef[1] = 1.;
	CHECK(!simplex_problem_from_model(&p, model));
	lp_simplex_default_options(&options, lp_simplex_ALGORITHM_DUAL_REVISED);
	sink.context = &mutation;
	sink.due = due;
	sink.report = observe;
	p.progress = &sink;
	for (change = 0; change < 2; change++) {
		simplex_solve_context_init(&context, NULL, &p.execution);
		simplex_solve_context_configure(&context, 0);
		context.controls.disable_scaling = 1;
		mutation.input = &p;
		mutation.enabled = change;
		mutation.changed = 0;
		code = simplex_dual_find_feasible_reported(
			&p, &options, x, &result, NULL, &report);
		if (!change) {
			CHECK(!code &&
				result.status == lp_simplex_CondUnsatisfied);
			CHECK(fabs(x[0] + x[1] - 1.) <=
				options.primal_tolerance);
			work = result.iterations;
		} else {
			CHECK(code != 0 && mutation.changed &&
				result.status == lp_simplex_PrecisionError);
			CHECK(result.iterations == work &&
				!report.point_available);
			CHECK(report.failure == SIMPLEX_FAILURE_STAGE_CONTRACT);
		}
		CHECK(!context.retries.count);
		simplex_solve_context_destroy(&context);
	}
	simplex_problem_destroy(&p);
	lp_model_free(model);
	return 0;
}

int main(void)
{
	int diagnostics;
	for (diagnostics = 0; diagnostics < 4; diagnostics++)
		if (service_identity(diagnostics))
			return 1;
	return diagnostic_binding() || input_lifetime();
}
