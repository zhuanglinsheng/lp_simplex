/**
 * @file test_tableau_report.c
 * @brief Raw Tableau phases and output permissions survive adapter boundaries.
 */

#include "simplex_stage.h"
#include "simplex_solve_context.h"
#include "simplex_tableau_solver.h"
#include "simplex_auto.h"
#include "simplex_auto_routes.h"
#include "simplex_execution.h"
#include <lp_simplex/lp_simplex.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#define CHECK(c) \
	do { \
		if (!(c)) { \
			fprintf(stderr, "tableau report %d: %s\n", __LINE__, \
				#c); \
			return 1; \
		} \
	} while (0)

static int reports(void)
{
	struct lp_Model *m = lp_model_create(1, 1);
	struct simplex_SolveContext context;
	struct simplex_Execution h;
	struct simplex_ExecutionScope scope;
	struct simplex_StageReport report;
	struct lp_simplex_Options o;
	double x[1], objective;
	int diagnostics, kind, i, status, work, state, event;
	CHECK(m);
	for (diagnostics = 0; diagnostics < 3; diagnostics++)
		for (kind = 0; kind < 5; kind++) {
			simplex_solve_context_init(
				&context, diagnostics ? &h : NULL, &scope);
			simplex_execution_set_parent(&scope,
				simplex_execution_begin(&scope,
					SIMPLEX_EXECUTION_SOLVE, 0, 0, 0));
			if (diagnostics == 2)
				for (i = h.count;
					i < SIMPLEX_EXECUTION_CAPACITY; i++)
					(void)simplex_execution_begin(&scope,
						SIMPLEX_EXECUTION_SOLVE, 0, 0,
						0);
			lp_simplex_default_options(
				&o, lp_simplex_ALGORITHM_TABLEAU);
			o.iteration_limit = kind < 3 ? kind + 1 : 100;
			m->constraints[0].type =
				kind == 4 ? optm_CONS_T_LE : optm_CONS_T_EQ;
			m->constraints[0].coef[0] = kind == 3 ? 0. : -1.;
			m->constraints[0].rhs = kind == 3 ? 1. : 0.;
			m->objective[0] = kind == 4 ? -1. : 1.;
			memset(&report, 0xff, sizeof(report));
			x[0] = NAN;
			work = 99;
			state = simplex_tableau_solve_model_reported(m, &o, x,
				&objective, &status, &work, &report, &scope);
			CHECK(report.tableau.available &&
				report.tableau.phase_one_work >= 1);
			CHECK(report.tableau.phase_one_work +
					report.tableau.cleanup_work +
					report.tableau.phase_two_work ==
				work);
			CHECK(!report.basis_available &&
				!report.factor_available &&
				!report.dual.quality &&
				!report.lower_bound.quality);
			CHECK(report.phase ==
				-1); /* Tableau phases have their own type. */
			if (kind < 3) {
				CHECK(work == kind + 1 &&
					report.tableau.phase_one_work == 1);
				CHECK(report.tableau.cleanup_work ==
					(kind ? 1 : 0));
				CHECK(report.tableau.phase ==
					(kind ? SIMPLEX_TABLEAU_PHASE_TWO
					      : SIMPLEX_TABLEAU_ARTIFICIAL_CLEANUP));
				CHECK(report.point_available == (kind == 2) &&
					report.tableau.point_available ==
						(kind == 2));
				CHECK(report.stop ==
					(kind == 2 ? SIMPLEX_STOP_OPTIMAL
						   : SIMPLEX_STOP_BUDGET));
				CHECK((state == 0) == (kind == 2));
				if (kind == 2)
					CHECK(x[0] == 0. && objective == 0.);
				else
					CHECK(isnan(x[0]));
			} else {
				CHECK(state != 0 && !report.point_available &&
					!report.tableau.point_available);
				CHECK(report.stop ==
					(kind == 3 ? SIMPLEX_STOP_INFEASIBLE
						   : SIMPLEX_STOP_UNBOUNDED));
				CHECK(report.tableau.phase ==
					(kind == 3 ? SIMPLEX_TABLEAU_PHASE_ONE
						   : SIMPLEX_TABLEAU_PHASE_TWO));
			}
			if (diagnostics == 1) {
				event = 0;
				for (i = scope.parent; i < h.count; i++)
					if (h.events[i].kind ==
						SIMPLEX_EXECUTION_TABLEAU)
						event = i + 1;
				CHECK(event > scope.parent &&
					h.events[event - 1].kind ==
						SIMPLEX_EXECUTION_TABLEAU);
				CHECK(h.events[event - 1].context == &context &&
					h.events[event - 1].model ==
						scope.model);
				CHECK(h.events[event - 1].parent ==
						scope.parent &&
					h.events[event - 1].grant ==
						o.iteration_limit);
				CHECK(h.events[event - 1].work == work &&
					h.events[event - 1].status == status &&
					h.events[event - 1].stop ==
						report.stop);
				CHECK(h.events[event - 1].tableau.phase ==
					report.tableau.phase);
			}
			simplex_solve_context_destroy(&context);
		}
	lp_model_free(m);
	return 0;
}

static int dual(const struct simplex_Problem *p,
	const struct lp_simplex_Options *o, double *x, double *y,
	struct lp_simplex_Result *r, int propagate)
{
	(void)o;
	(void)x;
	(void)y;
	(void)propagate;
	r->status = lp_simplex_PrecisionError;
	r->iterations = 1;
	p->report->stop = SIMPLEX_STOP_NUMERICAL;
	return -1;
}

static int hybrid(const struct simplex_Problem *p,
	const struct lp_simplex_Options *o, const double *seed, double *x,
	struct lp_simplex_Result *r)
{
	(void)seed;
	return dual(p, o, x, NULL, r, 0);
}

static int route(int kind, const struct simplex_Problem *p,
	const struct lp_simplex_Options *o, const double *seed, double *x,
	struct lp_simplex_Result *r)
{
	if (kind == SIMPLEX_AUTO_TABLEAU)
		return simplex_auto_route_solve(kind, p, o, seed, x, r);
	r->status = lp_simplex_CondUnsatisfied;
	r->iterations = 0;
	p->report->stop = SIMPLEX_STOP_NOT_APPLICABLE;
	return -1;
}

static int controller(void)
{
	struct lp_Model *m = lp_model_create(1, 1);
	struct simplex_Problem p;
	struct simplex_StageReport final;
	struct simplex_AutoHistory history;
	struct simplex_AutoKernels kernels;
	struct simplex_AutoRoutes routes;
	struct simplex_AutoBindings bindings = {NULL, NULL, 0};
	struct simplex_SolveContext context;
	struct simplex_Execution h;
	struct lp_simplex_Options o;
	struct lp_simplex_Result r;
	double x[1];
	int diagnostics, i, found, event;
	bindings.kernels = &kernels;
	bindings.routes = &routes;
	CHECK(m);
	m->constraints[0].type = optm_CONS_T_EQ;
	m->constraints[0].coef[0] = -1.;
	m->objective[0] = 1.;
	CHECK(!simplex_problem_from_model(&p, m));
	p.report = &final;
	kernels.dual = dual;
	kernels.hybrid = hybrid;
	routes.solve = route;
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_AUTO);
	o.iteration_limit = 100;
	for (diagnostics = 0; diagnostics < 3; diagnostics++) {
		simplex_solve_context_init(
			&context, diagnostics ? &h : NULL, &p.execution);
		if (diagnostics == 2)
			for (i = 0; i < SIMPLEX_EXECUTION_CAPACITY; i++)
				(void)simplex_execution_begin(&p.execution,
					SIMPLEX_EXECUTION_SOLVE, 0, 0, 0);
		CHECK(!simplex_auto_solve(
			&p, &o, x, &r, 0, &bindings, NULL, &history));
		CHECK(r.status == lp_simplex_Success && r.iterations == 5 &&
			x[0] == 0.);
		CHECK(isnan(r.dual_infeasibility) &&
			r.dual_residual_space ==
				lp_simplex_RESIDUAL_UNAVAILABLE);
		CHECK(final.stop == SIMPLEX_STOP_OPTIMAL &&
			final.point_available && !final.tableau.available);
		found = 0;
		for (i = 0; i < history.count; i++)
			if (history.stages[i].action.stage ==
				SIMPLEX_AUTO_TABLEAU) {
				const struct simplex_AutoStageRecord *s =
					history.stages + i;
				CHECK(s->report.stop == SIMPLEX_STOP_OPTIMAL &&
					s->report.point_available &&
					s->accepted_point);
				CHECK(s->work == 3 &&
					s->report.tableau.phase_one_work == 1 &&
					s->report.tableau.cleanup_work == 1 &&
					s->report.tableau.phase_two_work == 1);
				if (diagnostics == 1) {
					found = 0;
					for (event = 0; event < h.count;
						event++)
						if (h.events[event].kind ==
							SIMPLEX_EXECUTION_TABLEAU) {
							CHECK(h.events[event].parent ==
									s->action
										.execution_id &&
								h.events[event].work ==
									3);
							CHECK(h.events[event]
									.tableau
									.available &&
								h.events[event]
									.tableau
									.point_available);
							found++;
						}
					CHECK(found == 1);
				}
				found = 1;
			}
		CHECK(found);
		simplex_solve_context_destroy(&context);
	}
	simplex_problem_destroy(&p);
	lp_model_free(m);
	return 0;
}

static int public_entry(void)
{
	struct lp_Model *m = lp_model_create(1, 1);
	struct lp_simplex_Options o;
	struct lp_simplex_Result r;
	struct simplex_Execution h;
	double x[1];
	int diagnostics, i, found;
	CHECK(m);
	m->constraints[0].type = optm_CONS_T_EQ;
	m->constraints[0].coef[0] = -1.;
	m->constraints[0].rhs = -2.;
	m->objective[0] = 1.;
	m->bounds[0].b_type = optm_BOUND_T_BS;
	m->bounds[0].lb = 1.;
	m->bounds[0].ub = 3.;
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_TABLEAU);
	o.presolve = 0;
	for (diagnostics = 0; diagnostics < 2; diagnostics++) {
		CHECK(!simplex_solve_reported(
			m, &o, x, &r, diagnostics ? &h : NULL));
		CHECK(r.status == lp_simplex_Success &&
			fabs(x[0] - 2.) < 1e-10);
		CHECK(fabs(r.objective - 2.) < 1e-10 &&
			r.primal_infeasibility <= o.primal_tolerance);
		CHECK(isnan(r.dual_infeasibility) &&
			r.dual_residual_space ==
				lp_simplex_RESIDUAL_UNAVAILABLE);
		if (diagnostics) {
			found = 0;
			for (i = 0; i < h.count; i++)
				if (h.events[i].kind ==
					SIMPLEX_EXECUTION_TABLEAU) {
					const struct simplex_ExecutionEvent *e =
						h.events + i;
					CHECK(e->parent > 0 && e->model > 0 &&
						e->complete);
					CHECK(e->work == r.iterations &&
						e->stop ==
							SIMPLEX_STOP_OPTIMAL);
					CHECK(e->tableau.point_available &&
						e->tableau.phase ==
							SIMPLEX_TABLEAU_PHASE_TWO);
					CHECK(e->tableau.phase_one_work +
							e->tableau
								.cleanup_work +
							e->tableau
								.phase_two_work ==
						e->work);
					found++;
				}
			CHECK(found == 1);
		}
	}
	lp_model_free(m);
	return 0;
}

int main(void)
{
	return reports() || controller() || public_entry();
}
