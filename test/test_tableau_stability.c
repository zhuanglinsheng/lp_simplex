/**
 * @file test_tableau_stability.c
 * @brief Stable ratio windows must preserve blockers and honest numerical
 * exits.
 */

#include "simplex_pivot.h"
#include "simplex_tableau_solver.h"
#include "simplex_stage.h"
#include "simplex_solve_context.h"
#include "simplex_execution.h"
#include <lp_simplex/lp_simplex.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#define CHECK(c) \
	do { \
		if (!(c)) { \
			fprintf(stderr, "tableau stability %d: %s\n", \
				__LINE__, #c); \
			return 1; \
		} \
	} while (0)

static int ratio_windows(void)
{
	struct lp_simplex_Options o;
	double t[12], saved[12];
	int basis[2], pricing, kind, work, code;
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_TABLEAU);
	o.iteration_limit = 1;
	for (pricing = 0; pricing < 2; pricing++)
		for (kind = 0; kind < 8; kind++) {
			o.pricing = pricing ? lp_simplex_PRICING_DANTZIG
					    : lp_simplex_PRICING_BLAND;
			memset(t, 0, sizeof(t));
			t[0] = 1.;
			t[4] = 1e-8;
			t[5] = 1.;
			t[8] = 1.;
			t[10] = 1.;
			basis[0] = 1;
			basis[1] = 2;
			if (kind == 1) {
				t[4] = 1e-6;
				t[11] = 1.;
			}
			if (kind == 2) {
				t[4] = 1e-12;
				t[11] = 1e8;
			}
			if (kind == 3) {
				t[4] = 1e-6;
				t[11] = .05;
			}
			if (kind == 4) {
				t[4] = 1.;
				t[7] = -5e-8;
				t[8] = -1e8;
			}
			if (kind == 5) {
				t[4] = 1e-8;
				t[7] = DBL_MAX;
				t[8] = 0.;
			}
			if (kind == 7) {
				t[4] = 1e-8;
				t[7] = -5e-8;
				t[8] = 0.;
			}
			if (kind == 6) {
				t[4] = 1.;
				t[8] = 1.;
			}
			memcpy(saved, t, sizeof(t));
			work = 0;
			code = simplex_run_pivots(
				&work, t, 4, basis, 2, 3, 3, &o);
			CHECK(work == 1);
			if (kind == 2 || kind == 4 || kind == 5 || kind == 7) {
				CHECK(code == 9 && basis[0] == 1 &&
					basis[1] == 2 &&
					!memcmp(t, saved, sizeof(t)));
			} else {
				CHECK(code == 0);
				if (kind == 1 || kind == 6)
					CHECK(basis[0] == 0 && basis[1] == 2);
				else
					CHECK(basis[0] == 1 && basis[1] == 0);
				CHECK(t[7] >= -o.primal_tolerance &&
					t[11] >= -o.primal_tolerance);
				if (kind == 3)
					CHECK(t[7] < 0. &&
						fabs(t[7] + 5e-8) < 1e-20);
			}
		}
	return 0;
}

static int canonical_columns(void)
{
	double t[12] = {1e16, 0., 0., 0., 49., 1., 0., 49., 1e16, 0., 1., 1e16};
	double partial[6] = {3., 0., 0., 2., 1., 2.};
	double saved[12];
	memcpy(saved, t, sizeof(t));
	simplex_apply_pivot(t, 4, 2, 3, 0, 0, 1, 1, 1);
	CHECK(t[0] == 0. && t[4] == 1. && t[8] == 0.);
	CHECK(fabs(t[9] / (-1e16 / 49.) - 1.) < 1e-15 && t[10] == 1.);
	memcpy(t, saved, sizeof(t));
	simplex_apply_pivot(t, 4, 2, 3, 0, 0, 1, 1, 0);
	CHECK(t[0] == saved[0] && t[4] == 1. && t[8] == 0.);
	memcpy(t, saved, sizeof(t));
	simplex_apply_pivot(t, 4, 2, 3, 0, 0, 0, 0, 0);
	CHECK(!memcmp(t, saved, sizeof(t)));
	/* Without normalization this is a raw row operation, not a unit column. */
	simplex_apply_pivot(partial, 3, 1, 2, 0, 0, 0, 0, 1);
	CHECK(partial[0] == -3. && partial[3] == 2.);
	return 0;
}

static int damaged_terminal(void)
{
	struct lp_simplex_Options o;
	double t[6];
	int basis[1] = {1}, kind, work;
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_TABLEAU);
	for (kind = 0; kind < 6; kind++) {
		memset(t, 0, sizeof(t));
		t[4] = 1.;
		if (kind == 0)
			t[0] = NAN;
		if (kind == 1)
			t[5] = NAN;
		if (kind == 2)
			t[5] = -2. * o.primal_tolerance;
		if (kind == 3)
			t[0] = HUGE_VAL;
		if (kind == 4)
			t[2] = NAN;
		if (kind == 5)
			t[2] = HUGE_VAL;
		work = 0;
		CHECK(simplex_run_pivots(&work, t, 3, basis, 1, 2, 2, &o) ==
				9 &&
			work == 1);
	}
	return 0;
}

/* Four consistent equations have a finite known feasible point. The sum of
 * their artificial objective RHS values overflows; it proves no contradiction.
 * A second case overflows the real objective only after entering Phase II. */
static int objective_overflow(void)
{
	struct lp_Model *m = lp_model_create(4, 2);
	struct optm_VariableBound *bounds;
	struct simplex_SolveContext context;
	struct simplex_Execution h;
	struct simplex_ExecutionScope scope;
	struct simplex_StageReport report;
	struct lp_simplex_Options o;
	struct lp_simplex_Result r;
	double x[2], value;
	int diagnostics, converted, pricing, i, status, work, state, kind;
	CHECK(m);
	bounds = m->bounds;
	for (i = 0; i < 4; i++) {
		m->constraints[i].type = optm_CONS_T_EQ;
		m->constraints[i].coef[i / 2] = 1.;
		m->constraints[i].rhs = 8e307;
	}
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_TABLEAU);
	o.presolve = 0;
	for (diagnostics = 0; diagnostics < 3; diagnostics++)
		for (converted = 0; converted < 2; converted++)
			for (pricing = 0; pricing < 2; pricing++)
				for (kind = 0; kind < 2; kind++) {
					for (i = 0; i < 4; i++)
						m->constraints[i].rhs =
							kind ? 2. : 8e307;
					m->objective[0] = kind ? DBL_MAX : 0.;
					m->bounds = converted ? bounds : NULL;
					o.pricing =
						pricing ? lp_simplex_PRICING_DANTZIG
							: lp_simplex_PRICING_BLAND;
					simplex_solve_context_init(&context,
						diagnostics ? &h : NULL,
						&scope);
					if (diagnostics == 2)
						while (h.count <
							SIMPLEX_EXECUTION_CAPACITY)
							(void)simplex_execution_begin(
								&scope,
								SIMPLEX_EXECUTION_SOLVE,
								0, 0, 0);
					x[0] = x[1] = value = -123.;
					state = simplex_tableau_solve_model_reported(
						m, &o, x, &value, &status,
						&work, &report, &scope);
					if (status != lp_simplex_PrecisionError)
						fprintf(stderr,
							"finite feasible overflow: status=%d work=%d phase=%d\n",
							status, work,
							report.tableau.phase);
					CHECK(state &&
						status ==
							lp_simplex_PrecisionError &&
						report.stop ==
							SIMPLEX_STOP_NUMERICAL);
					CHECK(report.tableau.phase ==
							(kind ? SIMPLEX_TABLEAU_PHASE_TWO
							      : SIMPLEX_TABLEAU_PHASE_ONE) &&
						work == (kind ? 4 : 1) &&
						report.tableau.phase_one_work ==
							(kind ? 3 : 1) &&
						report.tableau.phase_two_work ==
							kind);
					CHECK(!report.point_available &&
						!report.tableau
							.recovered_point &&
						x[0] == -123. &&
						x[1] == -123. &&
						value == -123.);
					if (diagnostics == 1)
						for (i = 0; i < h.count; i++)
							if (h.events[i].kind ==
									SIMPLEX_EXECUTION_TABLEAU ||
								h.events[i].kind ==
									SIMPLEX_EXECUTION_TABLEAU_STANDARD)
								CHECK(h.events[i]
										.complete &&
									h.events[i].stop ==
										SIMPLEX_STOP_NUMERICAL &&
									h.events[i].status ==
										lp_simplex_PrecisionError &&
									!h.events[i]
										.point_source);
					simplex_solve_context_destroy(&context);
					CHECK(lp_simplex_solve(m, &o, x, &r) &&
						r.status ==
							lp_simplex_PrecisionError &&
						r.iterations ==
							(kind ? 4 : 1) &&
						isnan(x[0]) && isnan(x[1]));
				}
	m->bounds = bounds;
	lp_model_free(m);
	return 0;
}

int main(void)
{
	return ratio_windows() || canonical_columns() || objective_overflow() ||
	       damaged_terminal();
}
