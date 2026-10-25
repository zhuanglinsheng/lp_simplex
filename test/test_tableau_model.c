/**
 * @file test_tableau_model.c
 * @brief Tableau representation identity is runtime state, not journal state.
 */

#include "simplex_stage.h"
#include "simplex_solve_context.h"
#include "simplex_tableau_solver.h"
#include "simplex_execution.h"
#include <lp_simplex/lp_simplex.h>
#include <math.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>

#define CHECK(c) \
	do { \
		if (!(c)) { \
			fprintf(stderr, "tableau model %d: %s\n", __LINE__, \
				#c); \
			return 1; \
		} \
	} while (0)

static int transitions(void)
{
	struct lp_Model *m = lp_model_create(1, 1);
	struct optm_VariableBound *bounds;
	struct simplex_SolveContext context;
	struct simplex_Execution h;
	struct simplex_ExecutionScope scope;
	struct simplex_StageReport report;
	struct lp_simplex_Options o;
	double x[1], objective;
	int diagnostics, form, limited, i, before, count, status, work, state;
	int outer, standard, model, converted;
	CHECK(m);
	bounds = m->bounds;
	m->constraints[0].type = optm_CONS_T_EQ;
	m->constraints[0].coef[0] = 1.;
	m->constraints[0].rhs = 2.;
	m->objective[0] = 1.;
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_TABLEAU);
	for (diagnostics = 0; diagnostics < 5; diagnostics++)
		for (form = 0; form < 5; form++)
			for (limited = 0; limited < 2; limited++) {
				simplex_solve_context_init(&context,
					diagnostics ? &h : NULL, &scope);
				simplex_execution_set_parent(
					&scope, simplex_execution_begin(&scope,
							SIMPLEX_EXECUTION_SOLVE,
							0, 0, 0));
				if (diagnostics >= 2) {
					count = SIMPLEX_EXECUTION_CAPACITY -
						(diagnostics - 2);
					while (h.count < count)
						(void)simplex_execution_begin(
							&scope,
							SIMPLEX_EXECUTION_SOLVE,
							0, 0, 0);
				}
				before = diagnostics ? h.count : 0;
				converted = form != 0;
				m->bounds = converted ? bounds : NULL;
				bounds[0].b_type = form == 1   ? optm_BOUND_T_LO
						   : form == 2 ? optm_BOUND_T_UP
						   : form == 3
							   ? optm_BOUND_T_BS
							   : optm_BOUND_T_FR;
				bounds[0].lb = 1.;
				bounds[0].ub = 3.;
				o.iteration_limit = limited ? 1 : 100;
				x[0] = objective = NAN;
				memset(&report, 0xff, sizeof(report));
				state = simplex_tableau_solve_model_reported(m,
					&o, x, &objective, &status, &work,
					&report, &scope);
				CHECK(report.tableau.available &&
					report.tableau.converted == converted);
				CHECK(report.tableau.input_model == 1 &&
					report.tableau.standard_model ==
						1 + converted);
				CHECK(context.models == 1 + converted);
				CHECK(report.tableau.standard_rows ==
					1 + (form == 3));
				CHECK(report.tableau.standard_columns ==
					1 + (form == 4));
				CHECK(report.tableau.phase_one_work +
						report.tableau.cleanup_work +
						report.tableau.phase_two_work ==
					work);
				CHECK(report.point_available == !limited);
				CHECK(report.tableau.recovered_point ==
					(converted && !limited));
				CHECK(!report.basis_available &&
					!report.factor_available &&
					!report.dual.quality &&
					!report.lower_bound.quality);
				if (limited) {
					CHECK(state &&
						status ==
							lp_simplex_ExceedIterLimit &&
						report.stop ==
							SIMPLEX_STOP_BUDGET &&
						work == 1);
					CHECK(isnan(x[0]) && isnan(objective));
				} else {
					CHECK(!state &&
						status == lp_simplex_Success &&
						report.stop ==
							SIMPLEX_STOP_OPTIMAL);
					CHECK(fabs(x[0] - 2.) < 1e-12 &&
						fabs(objective - 2.) < 1e-12);
				}
				if (diagnostics) {
					outer = standard = model = 0;
					for (i = before; i < h.count; i++) {
						const struct
							simplex_ExecutionEvent
								*e = h.events +
								     i;
						CHECK(e->complete &&
							e->context == &context);
						if (e->kind ==
							SIMPLEX_EXECUTION_TABLEAU)
							outer = i + 1;
						if (e->kind ==
							SIMPLEX_EXECUTION_TABLEAU_MODEL)
							model = i + 1;
						if (e->kind ==
							SIMPLEX_EXECUTION_TABLEAU_STANDARD)
							standard = i + 1;
					}
					if (outer) {
						const struct
							simplex_ExecutionEvent
								*e = h.events +
								     outer - 1;
						CHECK(e->model == 1 &&
							e->parent ==
								scope.parent &&
							e->work == work &&
							e->grant ==
								o.iteration_limit);
						CHECK(e->point_source ==
							(limited ? 0
								 : standard));
						CHECK(e->tableau.standard_model ==
								1 + converted &&
							e->tableau.recovered_point ==
								(converted &&
									!limited));
						CHECK(e->source_model ==
							(limited ? 1
								 : 1 + converted));
					}
					if (model) {
						const struct
							simplex_ExecutionEvent
								*e = h.events +
								     model - 1;
						CHECK(converted &&
							e->model == 2 &&
							e->source_model == 1 &&
							e->parent == outer);
						CHECK(!e->result_available &&
							!e->work);
						CHECK(e->model_rows ==
								report.tableau
									.standard_rows &&
							e->model_columns ==
								report.tableau
									.standard_columns);
					}
					if (standard) {
						const struct
							simplex_ExecutionEvent
								*e = h.events +
								     standard -
								     1;
						CHECK(e->model ==
								1 + converted &&
							e->parent ==
								(converted ? model
									   : outer));
						CHECK(e->work == work &&
							e->status == status &&
							e->stop == report.stop);
						CHECK(e->tableau.point_available ==
								!limited &&
							!e->tableau
								.recovered_point);
						CHECK(e->tableau.standard_model ==
								1 + converted &&
							!e->point_source);
					}
					if (diagnostics == 1)
						CHECK(outer && standard &&
							(!!model == converted));
					if (diagnostics == 2)
						CHECK(!outer && !model &&
							!standard);
					if (diagnostics == 3)
						CHECK(outer && !model &&
							!standard);
					if (diagnostics == 4)
						CHECK(outer &&
							(converted ? (model &&
									     !standard)
								   : standard));
				}
				simplex_solve_context_destroy(&context);
			}
	m->bounds = bounds;
	lp_model_free(m);
	return 0;
}

static int lifetime(void)
{
	struct lp_Model *m = lp_model_create(1, 1);
	struct simplex_SolveContext context;
	struct simplex_Execution h;
	struct simplex_ExecutionScope scope;
	struct simplex_StageReport report;
	struct lp_simplex_Options o;
	double x[1], objective;
	int diagnostics, call, status, work;
	CHECK(m);
	m->constraints[0].type = optm_CONS_T_EQ;
	m->constraints[0].coef[0] = 1.;
	m->constraints[0].rhs = 2.;
	m->objective[0] = 1.;
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_TABLEAU);
	for (diagnostics = 0; diagnostics < 2; diagnostics++) {
		simplex_solve_context_init(
			&context, diagnostics ? &h : NULL, &scope);
		for (call = 0; call < 2; call++) {
			CHECK(!simplex_tableau_solve_model_reported(m, &o, x,
				&objective, &status, &work, &report, &scope));
			CHECK(report.tableau.input_model == 1 &&
				report.tableau.standard_model == call + 2);
			CHECK(context.models == call + 2 && x[0] == 2. &&
				objective == 2.);
		}
		context.models = INT_MAX;
		CHECK(!simplex_tableau_solve_model_reported(
			m, &o, x, &objective, &status, &work, &report, &scope));
		CHECK(!report.tableau.standard_model &&
			context.models == INT_MAX && report.point_available);
		simplex_solve_context_destroy(&context);
	}
	CHECK(!simplex_tableau_solve_model_reported(
		m, &o, x, &objective, &status, &work, &report, NULL));
	CHECK(!report.tableau.input_model && !report.tableau.standard_model &&
		report.tableau.recovered_point);
	lp_model_free(m);
	return 0;
}

/* Classification must survive each bounds representation without exporting
 * transformed coordinates, partial points or nonexistent dual/basis state. */
static int classifications(void)
{
	struct lp_Model *m = lp_model_create(2, 1);
	struct optm_VariableBound *bounds;
	struct simplex_SolveContext context;
	struct simplex_Execution h;
	struct simplex_ExecutionScope scope;
	struct simplex_StageReport report;
	struct lp_simplex_Options o;
	double x, objective;
	int diagnostics, form, kind, status, state, work, i, expected, phase;
	CHECK(m);
	bounds = m->bounds;
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_TABLEAU);
	for (diagnostics = 0; diagnostics < 3; diagnostics++)
		for (form = 0; form < 5; form++)
			for (kind = 0; kind < 2; kind++) {
				simplex_solve_context_init(&context,
					diagnostics ? &h : NULL, &scope);
				if (diagnostics == 2)
					while (h.count <
						SIMPLEX_EXECUTION_CAPACITY)
						(void)simplex_execution_begin(
							&scope,
							SIMPLEX_EXECUTION_SOLVE,
							0, 0, 0);
				m->bounds = form ? bounds : NULL;
				bounds[0].b_type = form == 1   ? optm_BOUND_T_LO
						   : form == 2 ? optm_BOUND_T_UP
						   : form == 3
							   ? optm_BOUND_T_BS
							   : optm_BOUND_T_FR;
				bounds[0].lb = -2.;
				bounds[0].ub = 3.;
				m->constraints[0].type =
					kind ? optm_CONS_T_LE : optm_CONS_T_EQ;
				m->constraints[1].type = optm_CONS_T_EQ;
				m->constraints[0].coef[0] =
					kind && form != 2 ? -1. : 1.;
				m->constraints[1].coef[0] = kind ? 0. : 1.;
				m->constraints[0].rhs = kind ? 0. : .5;
				m->constraints[1].rhs = kind ? 0. : 1.5;
				m->objective[0] = kind && form != 2 ? -1. : 1.;
				expected = !kind ? lp_simplex_Infeasibility
					   : form == 3
						   ? lp_simplex_Success
						   : lp_simplex_Unboundedness;
				phase = kind ? SIMPLEX_TABLEAU_PHASE_TWO
					     : SIMPLEX_TABLEAU_PHASE_ONE;
				x = objective = -123.;
				state = simplex_tableau_solve_model_reported(m,
					&o, &x, &objective, &status, &work,
					&report, &scope);
				CHECK(status == expected && work > 0 &&
					work <= o.iteration_limit);
				CHECK((state == 0) ==
					(expected == lp_simplex_Success));
				CHECK(report.tableau.phase == phase &&
					report.tableau.converted ==
						(form != 0));
				CHECK(report.tableau.input_model == 1 &&
					report.tableau.standard_model ==
						1 + (form != 0));
				CHECK(report.tableau.phase_one_work +
						report.tableau.cleanup_work +
						report.tableau.phase_two_work ==
					work);
				CHECK(!report.basis_available &&
					!report.factor_available &&
					!report.dual.quality &&
					!report.lower_bound.quality);
				if (expected == lp_simplex_Success) {
					CHECK(report.point_available &&
						report.tableau
							.recovered_point &&
						x == 3. && objective == -3.);
				} else {
					CHECK(!report.point_available &&
						!report.tableau
							.point_available &&
						!report.tableau
							.recovered_point &&
						x == -123. &&
						objective == -123.);
					CHECK(report.stop ==
						(kind ? SIMPLEX_STOP_UNBOUNDED
						      : SIMPLEX_STOP_INFEASIBLE));
				}
				if (diagnostics == 1)
					for (i = 0; i < h.count; i++) {
						const struct
							simplex_ExecutionEvent
								*e = h.events +
								     i;
						CHECK(e->complete);
						if (e->kind == SIMPLEX_EXECUTION_TABLEAU ||
							e->kind ==
								SIMPLEX_EXECUTION_TABLEAU_STANDARD) {
							CHECK(e->status ==
									expected &&
								e->work ==
									work &&
								e->stop ==
									report.stop);
							if (expected !=
								lp_simplex_Success)
								CHECK(!e->point_source);
						}
					}
				simplex_solve_context_destroy(&context);
			}
	m->bounds = bounds;
	lp_model_free(m);
	return 0;
}

int main(void)
{
	return transitions() || lifetime() || classifications();
}
