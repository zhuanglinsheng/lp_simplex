/**
 * @file test_numerical_policy.c
 * @brief Regression checks for numerical policy and scale-dependent recovery.
 */

#include "policy/pan.h"
#include "simplex_tableau.h"
#include "simplex_pivot.h"
#include "simplex_pan_pdhg.h"
#include "simplex_pan_bounded.h"
#include "simplex_pan_config.h"
#include "simplex_pan_recovery.h"
#include <lp_simplex/lp_simplex.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define CHECK(c) \
	do { \
		if (!(c)) { \
			fprintf(stderr, "line %d: %s\n", __LINE__, #c); \
			return 0; \
		} \
	} while (0)

static int solve_with_manual_pan_experiments(const struct lp_Model *m,
	const struct lp_simplex_Options *o, double *x,
	struct lp_simplex_Result *r)
{
	return o && o->algorithm == lp_simplex_ALGORITHM_PAN_BDA
		       ? lp_simplex_solve_experimental(m, o, x, r)
		       : lp_simplex_solve(m, o, x, r);
}

static int artificial_rows(void)
{
	double t[8], saved[8], scale;
	int basis[1], exponent, column, work = 0;
	/* An equality may be tiny in caller units without being redundant.
	 * Exercise coefficients both below and above the old absolute cutoff. */
	for (exponent = -48; exponent <= 48; exponent += 12) {
		for (column = 0; column < 2; column++) {
			memset(t, 0, sizeof(t));
			scale = ldexp(1., exponent);
			t[4 + column] = scale;
			t[6] = 1.;
			t[7] = scale * 2.;
			basis[0] = 2;
			CHECK(!simplex_tableau_remove_artificial_basis(
				t, 4, basis, 1, 2, 3, 1e-9, &work, 1000));
			CHECK(basis[0] == column && t[4 + column] == 1. &&
				t[7] == 2.);
		}
	}
	/* No real coefficient plus a nonzero RHS is unresolved, not redundant.
	 * Refusal must not silently zero either the row or RHS. */
	memset(t, 0, sizeof(t));
	t[6] = 1.;
	t[7] = 1e-12;
	basis[0] = 2;
	memcpy(saved, t, sizeof(t));
	CHECK(simplex_tableau_remove_artificial_basis(
		t, 4, basis, 1, 2, 3, 1e-9, &work, 1000));
	CHECK(!memcmp(saved, t, sizeof(t)) && basis[0] == 2);
	t[7] = 0.;
	CHECK(!simplex_tableau_remove_artificial_basis(
		t, 4, basis, 1, 2, 3, 1e-9, &work, 1000));
	CHECK(basis[0] == 2);
	/* Dividing by a subnormal risks overflow: refuse rather than erase. */
	t[4] = DBL_MIN / 2.;
	CHECK(simplex_tableau_remove_artificial_basis(
		t, 4, basis, 1, 2, 3, 1e-9, &work, 1000));
	return 1;
}

static int bland_index(void)
{
	struct lp_simplex_Options o;
	/* Smallest eligible reduced cost is below the old staged threshold. */
	double table[8] = {1e-5, 1., 0., 0., 1., 1., 1., 1.};
	int basis[1] = {2}, work = 0;
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_TABLEAU);
	o.iteration_limit = 1;
	(void)simplex_run_pivots(&work, table, 4, basis, 1, 3, 3, &o);
	CHECK(basis[0] == 0);
	return 1;
}

static int repair_progress(void)
{
	double previous = HUGE_VAL;
	CHECK(simplex_pan_pdhg_repair_progress(1., previous, 1e-7));
	previous =
		1.; /* Both failed and successful attempts bind this baseline. */
	CHECK(!simplex_pan_pdhg_repair_progress(1., previous, 1e-7));
	CHECK(!simplex_pan_pdhg_repair_progress(.75, previous, 1e-7));
	CHECK(!simplex_pan_pdhg_repair_progress(.5, previous, 1e-7));
	CHECK(simplex_pan_pdhg_repair_progress(
		nextafter(.5, 0.), previous, 1e-7));
	CHECK(!simplex_pan_pdhg_repair_progress(0., previous, 1e-7));
	CHECK(!simplex_pan_pdhg_repair_progress(NAN, previous, 1e-7));
	CHECK(!simplex_pan_pdhg_repair_progress(HUGE_VAL, previous, 1e-7));
	CHECK(!simplex_pan_pdhg_repair_progress(1., NAN, 1e-7));
	return 1;
}

static int progress_history(void)
{
	struct simplex_PdhgProgress h, primal;
	int k, stalled;
	simplex_pan_pdhg_progress_reset(&h);
	CHECK(!simplex_pan_pdhg_progress_stalled(&h, 0, 100.));
	/* Improvement followed by a rebound must keep the historical baseline. */
	CHECK(!simplex_pan_pdhg_progress_stalled(&h, 128, 50.));
	CHECK(!simplex_pan_pdhg_progress_stalled(&h, 256, 100.));
	CHECK(simplex_pan_pdhg_progress_stalled(
		&h, 128 + SIMPLEX_PDHG_STAGNANT_WORK, 60.));
	simplex_pan_pdhg_progress_reset(&h);
	for (k = 0; k < 100; k++)
		CHECK(!simplex_pan_pdhg_progress_stalled(
			&h, k * 128, ldexp(1., -k)));
	/* Productive work may exceed many restart windows without repair. */
	CHECK(!simplex_pan_pdhg_progress_stalled(&h, 12800, ldexp(1., -100)));
	CHECK(!simplex_pan_pdhg_progress_stalled(
		&h, 12800 + SIMPLEX_PDHG_STAGNANT_WORK - 1, ldexp(1., -100)));
	CHECK(simplex_pan_pdhg_progress_stalled(
		&h, 12800 + SIMPLEX_PDHG_STAGNANT_WORK, ldexp(1., -100)));
	CHECK(!simplex_pan_pdhg_progress_stalled(&h, 1024, NAN));
	CHECK(!simplex_pan_pdhg_progress_stalled(&h, 1024, -1.));
	/* A changed metric is explicitly reset; no false comparison with 100. */
	simplex_pan_pdhg_progress_reset(&h);
	CHECK(!simplex_pan_pdhg_progress_stalled(&h, 2000, 1e10));
	CHECK(!simplex_pan_pdhg_progress_stalled(&h, 2128, 1e10));
	CHECK(simplex_pan_pdhg_progress_stalled(
		&h, 2000 + SIMPLEX_PDHG_STAGNANT_WORK, 1e10));
	/* Recalibrating a joint metric cannot erase an unchanged physical
	 * feasibility component, even while the joint metric keeps improving. */
	simplex_pan_pdhg_progress_reset(&primal);
	for (k = 0; k <= SIMPLEX_PDHG_STAGNANT_WORK; k += 128) {
		simplex_pan_pdhg_progress_reset(&h);
		stalled = simplex_pan_pdhg_progress_stalled_pair(
			&h, &primal, k, 1. / (1. + k), 1.);
		CHECK(stalled == (k == SIMPLEX_PDHG_STAGNANT_WORK
						 ? SIMPLEX_PDHG_PRIMAL_STALLED
						 : SIMPLEX_PDHG_PROGRESSING));
	}
	/* Productive physical progress retains its matrix-free trajectory. */
	simplex_pan_pdhg_progress_reset(&primal);
	for (k = 0; k < 100; k++) {
		simplex_pan_pdhg_progress_reset(&h);
		CHECK(simplex_pan_pdhg_progress_stalled_pair(&h, &primal,
			      k * 128, ldexp(1., -k),
			      ldexp(1., -k)) == SIMPLEX_PDHG_PROGRESSING);
	}
	return 1;
}

static int objective_units(void)
{
	struct simplex_PanStandard a;
	double cost[2] = {1., -1.}, centered[2] = {100., 98.};
	double point[2] = {1., 1.};
	long double baseline, actual;
	int e;
	memset(&a, 0, sizeof(a));
	a.columns = 2;
	a.objective = centered;
	a.audit_objective = cost;
	baseline = simplex_pan_progress_threshold(&a, point, 0., 1e-7, 100.);
	CHECK(baseline > 0.); /* cancellation does not erase the error scale */
	for (e = -40; e <= 40; e += 10) {
		cost[0] = ldexp(1., e);
		cost[1] = -cost[0];
		a.objective_offset = 1e20;
		a.audit_objective_offset = -1e20;
		actual = simplex_pan_progress_threshold(
			&a, point, 0., 1e-7, 100.);
		CHECK(fabsl(actual / ldexp(1., e) - baseline) <=
			baseline * 1e-14);
	}
	return 1;
}

static int handoff_bounds(void)
{
	struct simplex_PanBounded b;
	int starts[2] = {0, 0};
	double cost[1] = {-1.}, lo[1] = {0.}, hi[1] = {HUGE_VAL};
	double x[1] = {0.}, y[1] = {0.};
	memset(&b, 0, sizeof(b));
	b.a.columns = 1;
	b.a.column_start = starts;
	b.a.objective = cost;
	b.lower = lo;
	b.upper = hi;
	/* The old gap omitted -infinity contributions and returned zero here. */
	CHECK(isinf(simplex_pan_pdhg_handoff_gap(&b, x, y, 1.)));
	lo[0] = -HUGE_VAL;
	cost[0] = 1.;
	CHECK(isinf(simplex_pan_pdhg_handoff_gap(&b, x, y, 1.)));
	cost[0] = 0.;
	/* Zero reduced cost contributes zero even with both bounds infinite. */
	CHECK(simplex_pan_pdhg_handoff_gap(&b, x, y, 1.) == 0.);
	return 1;
}

/* A pure activity-relative audit rejects this vanishing free coordinate until
 * underflow. Row scaling must permit ordinary finite-tolerance convergence. */
static int zero_activity_certificate(void)
{
	struct simplex_PanConfig controls;
	struct lp_Model *m = lp_model_create(1, 1);
	struct simplex_Problem p;
	struct simplex_PanBounded b;
	struct lp_simplex_Options o;
	double point[1] = {0.}, multiplier[1] = {1.};
	int work = 0, certificate;

	CHECK(m);
	m->constraints[0].type = optm_CONS_T_EQ;
	m->constraints[0].rhs = 0.;
	m->constraints[0].coef[0] = 1.;
	m->bounds[0].b_type = optm_BOUND_T_FR;
	CHECK(!simplex_problem_from_model(&p, m));
	CHECK(!simplex_pan_bounded_create(&b, &p));
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_AUTO);
	simplex_pan_config_init(&controls, &o);
	o.iteration_limit = SIMPLEX_PDHG_AUDIT_PERIOD + 1;
	certificate = simplex_pan_bounded_pdhg_configured_reported(&b, point,
		multiplier, &o, &work, o.iteration_limit, &controls, NULL,
		NULL);
	CHECK(certificate && work <= o.iteration_limit);
	CHECK(fabs(point[0]) <= o.primal_tolerance);
	simplex_pan_bounded_destroy(&b);
	simplex_problem_destroy(&p);
	lp_model_free(m);
	return 1;
}

static int transformed_case(int algorithm, int shift, int flip, int exponent)
{
	struct lp_Model *m = lp_model_create(2, 3);
	struct lp_simplex_Options o;
	struct lp_simplex_Result r;
	double coeff[2][3] = {{1., 0., 1.}, {0., 1., 1.}};
	double objective[3] = {2., 3., 10.}, optimum[3] = {1., 2., 0.};
	double x[3], variable_scale[3], scale, cost_scale = ldexp(1., exponent);
	int i, j, row, col, ok;

	CHECK(m);
	m->objective_offset = 17. * cost_scale;
	for (j = 0; j < 3; j++) {
		col = (j + shift) % 3;
		variable_scale[j] = ldexp(1., (col - 1) * (exponent / 4));
		m->objective[j] =
			objective[col] * cost_scale * variable_scale[j];
	}
	for (i = 0; i < 2; i++) {
		row = (i + flip) % 2;
		scale = ldexp(1., i ? exponent : -exponent);
		m->constraints[i].rhs = (row + 1.) * scale;
		for (j = 0; j < 3; j++) {
			col = (j + shift) % 3;
			m->constraints[i].coef[j] =
				coeff[row][col] * scale * variable_scale[j];
		}
	}
	CHECK(!lp_model_build_sparse(m));
	lp_simplex_default_options(&o, algorithm);
	o.presolve = 0;
	ok = !solve_with_manual_pan_experiments(m, &o, x, &r);
	if (!ok)
		fprintf(stderr,
			"algo=%d shift=%d flip=%d exponent=%d status=%d\n",
			algorithm, shift, flip, exponent, r.status);
	for (j = 0; j < 3; j++)
		ok = ok && fabs(x[j] * variable_scale[j] -
				   optimum[(j + shift) % 3]) < 1e-6;
	ok = ok && fabs(r.objective / cost_scale - 25.) < 1e-6 &&
	     r.primal_infeasibility <= o.primal_tolerance;
	lp_model_free(m);
	CHECK(ok);
	return 1;
}

static int transformed_models(void)
{
	int algorithm, shift, flip, exponent;
	for (algorithm = 0; algorithm <= 3; algorithm++)
		for (shift = 0; shift < 3; shift++)
			for (flip = 0; flip < 2; flip++)
				for (exponent = -12; exponent <= 12;
					exponent += 12)
					CHECK(transformed_case(algorithm, shift,
						flip, exponent));
	return 1;
}

static int parser(void)
{
	CHECK(simplex_pan_integer_setting_value("128") == 128);
	CHECK(simplex_pan_integer_setting_value("999999999999999999999999") ==
		INT_MAX);
	CHECK(simplex_pan_integer_setting_value("12garbage") == 0);
	CHECK(simplex_pan_integer_setting_value("-10") == 0);
	CHECK(simplex_pan_integer_setting_value("nan") == 0);
	CHECK(simplex_pan_integer_setting_value(NULL) == 0);
	return 1;
}

int main(void)
{
	return artificial_rows() && bland_index() && repair_progress() &&
			       progress_history() && objective_units() &&
			       handoff_bounds() &&
			       zero_activity_certificate() &&
			       transformed_models() && parser()
		       ? 0
		       : 1;
}
