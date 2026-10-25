/**
 * @file test_auto_resource.c
 * @brief Production AUTO startup ownership under required/optional storage failure.
 */

#include "simplex_auto.h"
#include "simplex_auto_routes.h"
#include "simplex_dual.h"
#include "simplex_execution.h"
#include "simplex_solve_context.h"
#include <lp_simplex/lp_simplex.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

static void *live[3];
static int requests, fail_at, active, bad_free;

static void *auto_resource_malloc(size_t bytes)
{
	void *p;
	int i;
	requests++;
	if (requests == fail_at || bytes > 1024 * 1024)
		return NULL;
	p = malloc(bytes ? bytes : 1);
	if (!p)
		return NULL;
	memset(p, 0xa5, bytes);
	for (i = 0; i < 3; i++)
		if (!live[i]) {
			live[i] = p;
			active++;
			return p;
		}
	free(p);
	return NULL;
}

static void auto_resource_free(void *p)
{
	int i;
	if (!p)
		return;
	for (i = 0; i < 3; i++)
		if (live[i] == p) {
			live[i] = NULL;
			active--;
			free(p);
			return;
		}
	bad_free++;
}

#define lp_simplex_malloc auto_resource_malloc
#define lp_simplex_free auto_resource_free
#define simplex_auto_solve resource_auto_solve
#include "../src/core/simplex_auto.c"
#undef lp_simplex_malloc
#undef lp_simplex_free

#define CHECK(c) \
	do { \
		if (!(c)) { \
			fprintf(stderr, "AUTO resource line %d: %s\n", \
				__LINE__, #c); \
			return 1; \
		} \
	} while (0)

static int startup(void)
{
	struct lp_Model *m = lp_model_create(1, 2);
	struct simplex_Problem p, before;
	struct simplex_SolveContext context;
	struct simplex_Execution events;
	struct simplex_AutoHistory history;
	struct simplex_AutoKernels kernels;
	struct simplex_AutoRoutes routes;
	struct simplex_AutoBindings bindings = {NULL, NULL, 0};
	struct simplex_StageReport report;
	struct lp_simplex_Options o;
	struct lp_simplex_Result result;
	double x[2], reference[2];
	int mode, failure, i, state, reference_work = -1;
	bindings.kernels = &kernels;
	bindings.routes = &routes;
	CHECK(m);
	m->constraints[0].type = optm_CONS_T_EQ;
	m->constraints[0].rhs = 1.;
	m->constraints[0].coef[0] = m->constraints[0].coef[1] = 1.;
	m->objective[1] = 1.;
	CHECK(!simplex_problem_from_model(&p, m));
	p.report = &report;
	kernels.dual = simplex_dual_solve_problem;
	kernels.hybrid = NULL;
	routes.solve = NULL;
	CHECK(simplex_auto_routes_has_checkpoint(&kernels));
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_AUTO);
	o.iteration_limit = 100;
	for (mode = 0; mode < 3; mode++)
		for (failure = 0; failure <= 3; failure++) {
			CHECK(!active && !bad_free);
			requests = 0;
			fail_at = failure;
			simplex_stage_reset(&report);
			simplex_solve_context_init(
				&context, &events, &p.execution);
			if (mode == 0)
				p.execution.history = NULL;
			if (mode == 2)
				for (i = 0; i < SIMPLEX_EXECUTION_CAPACITY; i++)
					(void)simplex_execution_begin(
						&p.execution,
						SIMPLEX_EXECUTION_SOLVE, 0, 0,
						0);
			before = p;
			/* Even a preexisting feasible output buffer is not an inherited
			 * AUTO candidate when required startup storage is unavailable. */
			x[0] = failure == 1 && mode == 1 ? 1. : 123.;
			x[1] = failure == 1 && mode == 1 ? 0. : 123.;
			state = resource_auto_solve(&p, &o, x, &result, 0,
				&bindings, NULL, &history);
			CHECK(!active && !bad_free &&
				!memcmp(&p, &before, sizeof(p)));
			CHECK(p.rhs[0] == 1. && p.objective[0] == 0. &&
				p.objective[1] == 1.);
			CHECK(p.matrix.nonzeros == 2 &&
				p.matrix.value[0] == 1. &&
				p.matrix.value[1] == 1. &&
				p.bounds[0].lb == 0. && p.bounds[1].lb == 0.);
			if (failure == 1) {
				CHECK(requests == 1 && state &&
					result.status ==
						lp_simplex_MemoryAllocError &&
					!result.iterations && !history.count);
				CHECK(history.outcome.exit ==
						SIMPLEX_AUTO_EXIT_RESOURCE &&
					!history.outcome.point_available &&
					!report.point_available &&
					!report.basis_available &&
					report.stop == SIMPLEX_STOP_RESOURCE);
				CHECK(!isfinite(result.objective) &&
					x[0] == (mode == 1 ? 1. : 123.) &&
					x[1] == (mode == 1 ? 0. : 123.));
				if (mode == 1)
					for (i = 0; i < events.count; i++)
						CHECK(events.events[i].kind !=
								SIMPLEX_EXECUTION_BACKEND &&
							events.events[i].kind !=
								SIMPLEX_EXECUTION_PAN_SERVICE);
			} else {
				CHECK(requests == 3 && !state &&
					result.status == lp_simplex_Success &&
					result.objective == 0. &&
					result.primal_infeasibility == 0. &&
					history.count == 1 &&
					history.outcome.point_available);
				CHECK(x[0] == 1. && x[1] == 0. &&
					history.stages[0].action.grant ==
						o.iteration_limit);
				CHECK(!history.stages[0].basis_accepted);
				CHECK(!!history.stages[0].basis_exported ==
					!failure);
				if (!mode && !failure) {
					reference_work = result.iterations;
					memcpy(reference, x, sizeof(reference));
				} else
					CHECK(result.iterations ==
							reference_work &&
						!memcmp(reference, x,
							sizeof(reference)));
			}
			simplex_solve_context_destroy(&context);
		}
	simplex_problem_destroy(&p);
	lp_model_free(m);
	return 0;
}

int main(void)
{
	return startup();
}
