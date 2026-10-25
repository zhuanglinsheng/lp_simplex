/**
 * @file test_auto_contract.c
 * @brief Accounting validity gates every AUTO output, including feasibility
 * seeds.
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
			fprintf(stderr, \
				"AUTO contract line %d: %s (mode %d)\n", \
				__LINE__, #c, mode); \
			return 1; \
		} \
	} while (0)

static int mode, calls, initializations, seed_ok;

static int dual(const struct simplex_Problem *p,
	const struct lp_simplex_Options *o, double *x, double *y,
	struct lp_simplex_Result *r, int propagate)
{
	(void)p;
	(void)o;
	(void)y;
	(void)propagate;
	calls++;
	x[0] = x[1] = 0.;
	r->iterations = 1;
	r->status = lp_simplex_ExceedIterLimit;
	if (mode == 5 && calls == 1) {
		x[0] = 1.;
		r->status = lp_simplex_Success;
		r->iterations = INT_MIN;
		r->dual_infeasibility = 0.;
		return 0;
	}
	return -1;
}

static int initialize(void *context, const struct simplex_Problem *input,
	const struct lp_simplex_Options *o, double *x, int *work)
{
	(void)input;
	(void)context;
	initializations++;
	x[0] = x[1] = .5;
	*work = mode == 1   ? -1
		: mode == 2 ? o->iteration_limit + 1
		: mode == 6 ? 0
			    : 1;
	/* A return code alone neither grants nor revokes feasibility. */
	return 0;
}

static int hybrid(const struct simplex_Problem *p,
	const struct lp_simplex_Options *o, const double *seed, double *x,
	struct lp_simplex_Result *r)
{
	seed_ok = mode == 1 || mode == 2
			  ? !seed
			  : seed && seed[0] == .5 && seed[1] == .5;
	x[0] = .75;
	x[1] = .25;
	r->status = lp_simplex_PrecisionError;
	r->iterations = 1;
	if (mode == 3 || mode == 4) {
		r->status = lp_simplex_Success;
		r->dual_infeasibility = 0.;
		r->iterations = mode == 3 ? -1 : o->iteration_limit + 1;
		/* A malformed producer cannot enable an oracle retry or bound. */
		p->report->failure = SIMPLEX_FAILURE_IMPLICIT_ORACLE;
		p->report->lower_bound.value = 0.;
		p->report->lower_bound.kind = SIMPLEX_OBSERVATION_LOWER_BOUND;
		p->report->lower_bound.space = SIMPLEX_SPACE_INPUT;
		p->report->lower_bound.quality = SIMPLEX_EVIDENCE_VERIFIED;
		p->report->lower_bound.work = o->iteration_limit;
		return 0;
	}
	return -1;
}

int main(void)
{
	struct lp_Model *m = lp_model_create(1, 2);
	struct simplex_Problem p;
	struct simplex_AutoHistory history;
	struct simplex_AutoKernels kernels;
	struct simplex_AutoBindings bindings = {NULL, NULL, 0};
	struct simplex_AutoInitializer initializer;
	struct lp_simplex_Options o;
	struct lp_simplex_Result r;
	double x[2];
	int state, diagnostics;
	struct simplex_Execution journal;
	struct simplex_SolveContext context;
	bindings.kernels = &kernels;
	CHECK(m);
	m->constraints[0].type = optm_CONS_T_EQ;
	m->constraints[0].rhs = 1.;
	m->constraints[0].coef[0] = m->constraints[0].coef[1] = 1.;
	m->objective[1] = 1.;
	CHECK(!simplex_problem_from_model(&p, m));
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_AUTO);
	o.iteration_limit = 160;
	kernels.dual = dual;
	kernels.hybrid = hybrid;
	initializer.context = NULL;
	initializer.run = initialize;
	for (diagnostics = 0; diagnostics < 2; diagnostics++) {
		for (mode = 0; mode < 7; mode++) {
			simplex_solve_context_init(&context,
				diagnostics ? &journal : NULL, &p.execution);
			calls = initializations = seed_ok = 0;
			state = simplex_auto_solve(&p, &o, x, &r, 0, &bindings,
				&initializer, &history);
			CHECK(state && r.status != lp_simplex_Success &&
				r.iterations <= o.iteration_limit);
			CHECK(calls == 2 && initializations == 1 && seed_ok &&
				history.count == 4);
			CHECK(history.have_point &&
				simplex_problem_primal_infeasibility(&p, x) ==
					0.);
			if (mode == 1 || mode == 2) {
				CHECK(!history.stages[1].accepted_point);
				CHECK(history.stages[1].work ==
					history.stages[1].action.grant);
				CHECK(history.stages[1].status ==
					lp_simplex_PrecisionError);
				CHECK(history.stages[1].report.failure ==
					SIMPLEX_FAILURE_STAGE_CONTRACT);
				CHECK(history.stages[2].action.point_source ==
					-1);
			} else
				CHECK(history.stages[1].accepted_point);
			if (mode == 6)
				CHECK(history.stages[1].work == 0 &&
					r.iterations == 3);
			if (mode == 3 || mode == 4) {
				CHECK(history.stages[2].validation ==
					SIMPLEX_AUTO_TERMINAL_INVALID);
				CHECK(!history.stages[2].accepted_point &&
					history.point_source == 1);
				CHECK(history.stages[2].report.failure ==
					SIMPLEX_FAILURE_STAGE_CONTRACT);
				CHECK(x[0] == .5 && x[1] == .5 &&
					history.best_objective == .5);
				CHECK(isinf(history.best_bound) &&
					history.best_bound < 0.);
			} else
				CHECK(history.point_source == 2 && x[1] == .25);
			if (mode == 5) {
				CHECK(!history.stages[0].accepted_point &&
					history.stages[0].reported_work ==
						INT_MIN);
				CHECK(history.stages[2].action.reason ==
					SIMPLEX_REASON_STAGE_CONTRACT);
				CHECK(history.stages[2].action.cause == 0);
			}
			simplex_solve_context_destroy(&context);
		}
	}
	simplex_problem_destroy(&p);
	lp_model_free(m);
	puts("AUTO rejects points, bounds and terminal claims from invalid work reports");
	return 0;
}
