/**
 * @file test_auto_output.c
 * @brief Global exit, returned candidate and raw producer reports are distinct.
 */

#include "simplex_stage.h"
#include "simplex_solve_context.h"
#include "simplex_execution.h"
#include "simplex_auto.h"
#include <lp_simplex/lp_simplex.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>

#define CHECK(c) \
	do { \
		if (!(c)) { \
			fprintf(stderr, "AUTO output line %d mode=%d: %s\n", \
				__LINE__, mode, #c); \
			return 1; \
		} \
	} while (0)
static int mode, calls;
static const struct simplex_Problem *running;

static void diagnostics(const struct simplex_Problem *p,
	const struct lp_simplex_Options *o, const double *x,
	struct lp_simplex_Result *r)
{
	simplex_stage_observe(p, o, x, NULL, 1, r->iterations, p->report);
	p->report->stop = SIMPLEX_STOP_CYCLE;
	p->report->phase = 7;
	p->report->perturbed_bounds = p->report->perturbed_costs = 1;
	p->report->working_primal.quality = p->report->working_dual.quality =
		SIMPLEX_EVIDENCE_ESTIMATED;
	p->report->first_order_stop = SIMPLEX_STOP_STAGNATION;
	p->report->input_improvements = SIMPLEX_INPUT_OBJECTIVE_PROGRESS;
}

static int dual(const struct simplex_Problem *p,
	const struct lp_simplex_Options *o, double *x, double *y,
	struct lp_simplex_Result *r, int propagate)
{
	(void)y;
	(void)propagate;
	calls++;
	running = p;
	x[0] = calls == 1 ? .5 : .25;
	x[1] = calls == 1 ? .5 : .75;
	if (mode == 1 && calls > 1) {
		x[0] = 0.;
		x[1] = 2.;
	}
	if (mode == 4 || mode == 6 || mode == 8)
		x[0] = x[1] = 0.;
	r->iterations = mode == 2 && calls > 1 ? o->iteration_limit : 1;
	r->status = mode == 8 ? lp_simplex_Infeasibility
			      : lp_simplex_PrecisionError;
	r->objective = 99.;
	r->dual_infeasibility = .7;
	r->dual_residual_space = lp_simplex_RESIDUAL_ORIGINAL;
	diagnostics(p, o, x, r);
	if (mode == 5 && calls > 1) {
		r->iterations = INT_MIN;
		x[0] = 1.;
		x[1] = 0.;
	}
	return -1;
}

static int hybrid(const struct simplex_Problem *p,
	const struct lp_simplex_Options *o, const double *seed, double *x,
	struct lp_simplex_Result *r)
{
	(void)seed;
	x[0] = .75;
	x[1] = .25;
	if (mode == 3) {
		x[0] = 1.;
		x[1] = 0.;
	}
	if (mode == 6)
		x[0] = x[1] = NAN;
	r->iterations = 1;
	r->status = mode == 3	? lp_simplex_Success
		    : mode == 7 ? lp_simplex_Unboundedness
				: lp_simplex_PrecisionError;
	r->dual_infeasibility = mode == 3 ? 0. : .2;
	r->dual_residual_space = lp_simplex_RESIDUAL_ORIGINAL;
	diagnostics(p, o, x, r);
	return mode == 3 ? 0 : -1;
}

static int initialize(void *context, const struct simplex_Problem *input,
	const struct lp_simplex_Options *o, double *x, int *work)
{
	(void)input;
	(void)context;
	(void)o;
	x[0] = x[1] = .5;
	*work = 1;
	((struct simplex_Problem *)running)->columns = INT_MAX;
	return 1;
}

static int run(int tracing)
{
	struct lp_Model *m = lp_model_create(1, 2);
	struct simplex_Problem p;
	struct simplex_StageReport report;
	struct simplex_SolveContext context;
	struct simplex_Execution execution;
	struct simplex_AutoHistory history;
	struct simplex_AutoKernels kernels;
	struct simplex_AutoBindings bindings = {NULL, NULL, 0};
	struct simplex_AutoInitializer initializer;
	struct lp_simplex_Options o;
	struct lp_simplex_Result result;
	double x[2];
	int state, i, last;
	bindings.kernels = &kernels;
	CHECK(m);
	m->constraints[0].type = optm_CONS_T_EQ;
	m->constraints[0].rhs = 1.;
	m->constraints[0].coef[0] = m->constraints[0].coef[1] = 1.;
	m->objective[1] = 1.;
	CHECK(!simplex_problem_from_model(&p, m));
	simplex_solve_context_init(
		&context, tracing ? &execution : NULL, &p.execution);
	if (tracing == 2)
		for (i = 0; i < SIMPLEX_EXECUTION_CAPACITY; i++)
			(void)simplex_execution_begin(
				&p.execution, SIMPLEX_EXECUTION_SOLVE, 0, 0, 0);
	p.report = &report;
	kernels.dual = dual;
	kernels.hybrid = hybrid;
	initializer.context = NULL;
	initializer.run = initialize;
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_AUTO);
	o.iteration_limit = 16;
	calls = 0;
	state = simplex_auto_solve(&p, &o, x, &result, 0, &bindings,
		mode == 4 ? &initializer : NULL, &history);
	CHECK(history.outcome.status == result.status &&
		history.outcome.state == state);
	CHECK(history.outcome.work == result.iterations &&
		report.stop == history.outcome.stop);
	CHECK(report.stop == simplex_stage_stop_from_status(result.status));
	CHECK(report.phase == -1 && !report.perturbed_bounds &&
		!report.perturbed_costs);
	CHECK(!report.working_primal.quality && !report.working_dual.quality);
	CHECK(!report.dual.quality && !report.lower_bound.quality &&
		!report.input_improvements);
	CHECK(!report.basis_available && !report.factor_available &&
		!report.first_order_stop);
	CHECK(!report.work.basis_counts_available);
	last = history.count - 1;
	CHECK(history.outcome.cause.record == last);
	if (mode == 4) {
		CHECK(state && result.status == lp_simplex_PrecisionError);
		CHECK(history.outcome.exit == SIMPLEX_AUTO_EXIT_INPUT_CHANGED);
		CHECK(!history.outcome.point_available &&
			!report.primal.quality && !history.have_point);
		CHECK(isnan(result.objective) &&
			isinf(result.primal_infeasibility));
	} else if (mode == 6 || mode == 8) {
		CHECK(state && !history.outcome.point_available &&
			!report.point_available);
		CHECK(history.outcome.point.record == -1 &&
			!history.outcome.point.event);
		CHECK(history.outcome.exit ==
			(mode == 8 ? SIMPLEX_AUTO_EXIT_TERMINAL
				   : SIMPLEX_AUTO_EXIT_ROUTES));
	} else {
		double objective = mode == 3 ? 0. : .25;
		CHECK(result.objective == objective &&
			report.feasible_objective.value == objective);
		CHECK(report.primal.value == 0. &&
			report.feasibility_merit.value == 0. &&
			report.violation_count.value == 0.);
		CHECK(report.primal.work == result.iterations &&
			report.feasible_objective.work == result.iterations);
		CHECK(report.point_available &&
			history.outcome.point_available &&
			history.outcome.point.record == 1);
		CHECK(history.stages[last].report.stop == SIMPLEX_STOP_CYCLE ||
			mode == 5);
		if (mode == 3 || mode == 7) {
			CHECK((mode == 3 ? !state : state != 0) &&
				history.outcome.exit ==
					SIMPLEX_AUTO_EXIT_TERMINAL &&
				!history.outcome.restored);
		} else {
			CHECK(state && history.outcome.restored &&
				x[0] == .75 && x[1] == .25);
			CHECK(history.outcome.exit ==
				(mode == 2 || mode == 5
						? SIMPLEX_AUTO_EXIT_BUDGET
						: SIMPLEX_AUTO_EXIT_ROUTES));
		}
		if (mode == 0)
			CHECK(history.stages[last]
					.report.feasible_objective.value ==
				.75);
		if (mode == 1)
			CHECK(history.stages[last]
					.report.violation_count.value == 1.);
	}
	if (tracing == 1) {
		const struct simplex_ExecutionEvent *event =
			execution.events + execution.count - 1;
		CHECK(event->kind == SIMPLEX_EXECUTION_AUTO_FINISH &&
			event->detail == history.outcome.exit);
		CHECK(event->point_source == history.outcome.point.event &&
			event->cause == history.outcome.cause.event);
		CHECK(event->work == result.iterations &&
			event->stop == report.stop && event->complete);
		CHECK(history.outcome.execution_id == execution.count);
	} else
		CHECK(!history.outcome.execution_id);
	if (tracing == 2)
		CHECK(execution.dropped == history.count + 1);
	simplex_solve_context_destroy(&context);
	simplex_problem_destroy(&p);
	lp_model_free(m);
	return 0;
}

int main(void)
{
	int tracing;
	for (mode = 0; mode <= 8; mode++)
		for (tracing = 0; tracing < 3; tracing++)
			CHECK(!run(tracing));
	return 0;
}
