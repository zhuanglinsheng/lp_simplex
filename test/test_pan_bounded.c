/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file test_pan_bounded.c
 * @brief Regression tests for native-bound Pan models and original-model
 * solution checks.
 */

#include "simplex_pan_config.h"
#include "simplex_pan_bounded.h"
#include "simplex_pan.h"

#include <lp_simplex/lp_simplex.h>

#include <math.h>
#include <stdio.h>
#include <string.h>

static int representation(void)
{
	struct lp_Model *m = lp_model_create(3, 5);
	struct simplex_Problem p;
	struct simplex_PanBounded b;
	double point[7], r[3], work[3], x[5] = {-2., 3., -1., .25, 1.};
	int j, good;

	if (!m)
		return 0;
	for (j = 0; j < 5; j++)
		m->constraints[0].coef[j] = 1.;
	m->constraints[0].type = optm_CONS_T_EQ;
	m->constraints[0].rhs = 1.25;
	m->constraints[1].coef[0] = 1.;
	m->constraints[1].type = optm_CONS_T_LE;
	m->constraints[1].rhs = 0.;
	m->constraints[2].coef[1] = 1.;
	m->constraints[2].type = optm_CONS_T_GE;
	m->constraints[2].rhs = 2.;
	for (j = 0; j < 5; j++) {
		m->bounds[j].b_type = optm_BOUND_T_BS;
		m->bounds[j].lb = -2.;
		m->bounds[j].ub = 3.;
	}
	m->bounds[2].b_type = optm_BOUND_T_FR;
	m->bounds[3].lb = m->bounds[3].ub = .25;
	if (simplex_problem_from_model(&p, m))
		return 0;
	if (simplex_pan_bounded_create(&b, &p))
		return 0;
	good = b.a.rows == 3 && b.a.columns == 7 && b.a.original_columns == 5 &&
	       simplex_pan_bounded_encode(&b, x, point, r, work, 1e-7) == 0 &&
	       point[0] == -2. && point[2] == -1. && point[3] == .25 &&
	       point[5] == 2. && point[6] == 1.;
	x[0] = -3.;
	good = good &&
	       simplex_pan_bounded_encode(&b, x, point, r, work, 1e-7) != 0;
	simplex_pan_bounded_destroy(&b);
	simplex_problem_destroy(&p);
	lp_model_free(m);
	return good;
}

static int kkt(void)
{
	struct simplex_PanBounded b;
	int starts[5] = {0, 0, 0, 0, 0};
	double cost[4] = {1., -2., 0., 123.},
	       lower[4] = {-2., -2., -HUGE_VAL, .25};
	double upper[4] = {3., 3., HUGE_VAL, .25}, x[4] = {-2., 3., -7., .25},
	       y = 0., v;

	memset(&b, 0, sizeof(b));
	b.a.columns = 4;
	b.a.column_start = starts;
	b.a.objective = cost;
	b.lower = lower;
	b.upper = upper;
	if (!simplex_pan_bounded_check_dual(&b, x, &y, 1e-8, &v))
		return 0;
	cost[0] = -1.;
	if (simplex_pan_bounded_check_dual(&b, x, &y, 1e-8, &v) || v < .9)
		return 0;
	cost[0] = 1.;
	cost[1] = 2.;
	if (simplex_pan_bounded_check_dual(&b, x, &y, 1e-8, &v))
		return 0;
	cost[1] = -2.;
	cost[2] = 1.;
	return !simplex_pan_bounded_check_dual(&b, x, &y, 1e-8, &v);
}

static unsigned int random_number(unsigned int *state)
{
	*state = *state * 1664525U + 1013904223U;
	return *state;
}

/* Native box KKT witnesses: g positive at lower, negative at upper, zero
 * inside/free. Fixed columns allow arbitrary g. Includes duplicate rows,
 * negative bounds and power-of-two column scaling. No reference solver. */
static int known(unsigned int seed, int warm)
{
	struct simplex_PanConfig controls;
	struct lp_Model *m = lp_model_create(10, 24);
	struct simplex_Problem p;
	struct lp_simplex_Options o;
	struct lp_simplex_Result result;
	double expected[24], start[24], out[24], y[10], objective = 0.,
							error = 0.;
	unsigned int state = seed;
	int i, j, code, good;

	if (!m)
		return 0;
	for (i = 0; i < 10; i++) {
		m->constraints[i].type = optm_CONS_T_EQ;
		y[i] = (int)(random_number(&state) % 5) - 2.;
	}
	for (j = 0; j < 24; j++) {
		double scale = ldexp(1., (int)(random_number(&state) % 17) - 8),
		       g = 0.;

		m->bounds[j].b_type = optm_BOUND_T_BS;
		m->bounds[j].lb = -2. / scale;
		m->bounds[j].ub = 3. / scale;
		switch (j % 6) {
		case 0:
			expected[j] = -2. / scale;
			g = scale;
			break;
		case 1:
			expected[j] = 3. / scale;
			g = -scale;
			break;
		case 2:
			expected[j] = .5 / scale;
			break;
		case 3:
			expected[j] = -1. / scale;
			m->bounds[j].b_type = optm_BOUND_T_FR;
			break;
		case 4:
			expected[j] = .25 / scale;
			m->bounds[j].lb = m->bounds[j].ub = expected[j];
			g = 2. * scale;
			break;
		default:
			expected[j] = -2. / scale;
			m->bounds[j].b_type = optm_BOUND_T_LO;
			g = scale;
			break;
		}
		for (i = 0; i < 8; i++)
			m->constraints[i].coef[j] =
				j < 8 ? (i == j ? scale : 0.)
				      : (random_number(&state) % 3 == 0
							? scale * ((int)(random_number(
										 &state) %
									   5) -
									  2.)
							: 0.);
		m->constraints[8].coef[j] = m->constraints[0].coef[j];
		m->constraints[9].coef[j] = m->constraints[1].coef[j];
		m->objective[j] = g;
		for (i = 0; i < 10; i++)
			m->objective[j] -= y[i] * m->constraints[i].coef[j];
		objective += m->objective[j] * expected[j];
	}
	for (i = 0; i < 10; i++)
		for (j = 0; j < 24; j++)
			m->constraints[i].rhs +=
				m->constraints[i].coef[j] * expected[j];
	for (j = 0; j < 24; j++) {
		start[j] = expected[j];
		/* Feasible to roundoff, but exact comparisons see spurious free
		 * coordinates at both lower and upper bounds. */
		if (warm == 2 && (j % 6 == 0 || j % 6 == 5))
			start[j] += 2e-15 * (1. + fabs(start[j]));
		if (warm == 2 && j % 6 == 1)
			start[j] -= 2e-15 * (1. + fabs(start[j]));
	}
	if (simplex_problem_from_model(&p, m))
		return 0;
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_PAN_BDA);
	simplex_pan_config_resolve(&controls, &p, &o);
	controls.model_form = SIMPLEX_PAN_NATIVE_BOUNDS;
	controls.direction_method = SIMPLEX_PAN_CONE_DIRECTION;
	if (controls.oracle == SIMPLEX_PAN_NO_ORACLE)
		controls.oracle = SIMPLEX_PAN_EXPLICIT_ORACLE;
	o.iteration_limit = 20000;
	o.presolve = 0;
	code = simplex_pan_bounded_solve_certified_configured(
		&p, &o, warm ? start : NULL, out, NULL, &result, &controls);
	for (i = 0; i < 10; i++) {
		long double activity = 0.;

		for (j = 0; j < 24; j++)
			activity +=
				(long double)m->constraints[i].coef[j] * out[j];
		error = fmax(error,
			fabs((double)(activity - m->constraints[i].rhs)));
	}
	for (j = 0; j < 24; j++) {
		if (!isfinite(out[j]))
			error = HUGE_VAL;
		if (m->bounds[j].b_type != optm_BOUND_T_FR)
			error = fmax(error, m->bounds[j].lb - out[j]);
		if (m->bounds[j].b_type == optm_BOUND_T_BS)
			error = fmax(error, out[j] - m->bounds[j].ub);
	}
	good = code == 0 && result.status == lp_simplex_Success &&
	       error < 1e-7 &&
	       fabs(result.objective - objective) <
		       1e-7 * (1. + fabs(objective));
	if (!good)
		printf("native seed=%u warm=%d status=%d iterations=%d obj=%.12g expected=%.12g primal=%.6g\n",
			seed, warm, result.status, result.iterations,
			result.objective, objective, error);
	if (!good)
		for (j = 0; j < 24; j++)
			printf("coordinate %d out=%.17g expected=%.17g lower=%.17g upper=%.17g\n",
				j, out[j], expected[j], m->bounds[j].lb,
				m->bounds[j].ub);
	simplex_problem_destroy(&p);
	lp_model_free(m);
	return good;
}

static int free_ray(void)
{
	struct simplex_PanConfig controls;
	struct lp_Model *m = lp_model_create(1, 2);
	struct simplex_Problem p;
	struct lp_simplex_Options o;
	struct lp_simplex_Result result;
	double start[2] = {0., 0.}, out[2];
	int code, good;

	if (!m)
		return 0;
	m->constraints[0].type = optm_CONS_T_EQ;
	m->constraints[0].rhs = 0.;
	m->constraints[0].coef[0] = 1.;
	m->constraints[0].coef[1] = -1.;
	m->objective[0] = m->objective[1] = -1.;
	m->bounds[0].b_type = m->bounds[1].b_type = optm_BOUND_T_FR;
	if (simplex_problem_from_model(&p, m))
		return 0;
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_PAN_BDA);
	simplex_pan_config_resolve(&controls, &p, &o);
	controls.model_form = SIMPLEX_PAN_NATIVE_BOUNDS;
	controls.direction_method = SIMPLEX_PAN_CONE_DIRECTION;
	if (controls.oracle == SIMPLEX_PAN_NO_ORACLE)
		controls.oracle = SIMPLEX_PAN_EXPLICIT_ORACLE;
	code = simplex_pan_bounded_solve_certified_configured(
		&p, &o, start, out, NULL, &result, &controls);
	good = code != 0 && result.status == lp_simplex_Unboundedness;
	simplex_problem_destroy(&p);
	lp_model_free(m);
	return good;
}

static int bound_exchange(void)
{
	struct simplex_PanConfig controls;
	struct lp_Model *m = lp_model_create(1, 2);
	struct simplex_Problem p;
	struct lp_simplex_Options o;
	struct lp_simplex_Result result;
	double start[2] = {3., -2.}, out[2];
	int j, good;

	if (!m)
		return 0;
	m->constraints[0].type = optm_CONS_T_EQ;
	m->constraints[0].rhs = 1.;
	m->constraints[0].coef[0] = m->constraints[0].coef[1] = 1.;
	m->objective[0] = 1.;
	for (j = 0; j < 2; j++) {
		m->bounds[j].b_type = optm_BOUND_T_BS;
		m->bounds[j].lb = -2.;
		m->bounds[j].ub = 3.;
	}
	if (simplex_problem_from_model(&p, m))
		return 0;
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_PAN_BDA);
	simplex_pan_config_resolve(&controls, &p, &o);
	controls.model_form = SIMPLEX_PAN_NATIVE_BOUNDS;
	controls.direction_method = SIMPLEX_PAN_CONE_DIRECTION;
	if (controls.oracle == SIMPLEX_PAN_NO_ORACLE)
		controls.oracle = SIMPLEX_PAN_EXPLICIT_ORACLE;
	good = simplex_pan_bounded_solve_certified_configured(
		       &p, &o, start, out, NULL, &result, &controls) == 0 &&
	       fabs(out[0] + 2.) < 1e-10 && fabs(out[1] - 3.) < 1e-10 &&
	       fabs(result.objective + 2.) < 1e-10;
	o.iteration_limit = 0;
	good = good &&
	       simplex_pan_bounded_solve_certified_configured(
		       &p, &o, start, out, NULL, &result, &controls) != 0 &&
	       result.status == lp_simplex_ExceedIterLimit && out[0] == 3. &&
	       out[1] == -2.;
	simplex_problem_destroy(&p);
	lp_model_free(m);
	return good;
}

static int scaled_coordinate(void)
{
	struct simplex_PanConfig controls;
	struct lp_Model *m = lp_model_create(1, 1);
	struct simplex_Problem p;
	struct lp_simplex_Options o;
	struct lp_simplex_Result result;
	double start[1] = {1e-12}, out[1];
	int good;

	if (!m)
		return 0;
	m->constraints[0].type = optm_CONS_T_EQ;
	m->constraints[0].rhs = 1.;
	m->constraints[0].coef[0] = 1e12;
	m->bounds[0].b_type = optm_BOUND_T_BS;
	m->bounds[0].lb = 0.;
	m->bounds[0].ub = 1e-12;
	m->objective[0] = -1.;
	if (simplex_problem_from_model(&p, m))
		return 0;
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_PAN_BDA);
	simplex_pan_config_resolve(&controls, &p, &o);
	controls.model_form = SIMPLEX_PAN_NATIVE_BOUNDS;
	controls.direction_method = SIMPLEX_PAN_CONE_DIRECTION;
	if (controls.oracle == SIMPLEX_PAN_NO_ORACLE)
		controls.oracle = SIMPLEX_PAN_EXPLICIT_ORACLE;
	good = simplex_pan_bounded_solve_certified_configured(
		       &p, &o, start, out, NULL, &result, &controls) == 0 &&
	       out[0] == 1e-12;
	simplex_problem_destroy(&p);
	/* A small nonzero equality column cannot certify a recession ray. */
	m->constraints[0].rhs = 0.;
	m->constraints[0].coef[0] = 1e-12;
	m->bounds[0].b_type = optm_BOUND_T_LO;
	start[0] = 0.;
	if (simplex_problem_from_model(&p, m))
		return 0;
	good = good &&
	       simplex_pan_bounded_solve_certified_configured(
		       &p, &o, start, out, NULL, &result, &controls) == 0 &&
	       fabs(out[0]) < 1e-12;
	simplex_problem_destroy(&p);
	lp_model_free(m);
	return good;
}

static int local_restore(void)
{
	struct simplex_PanConfig controls;
	struct lp_Model *m = lp_model_create(1, 2);
	struct simplex_Problem p;
	struct simplex_PanBounded b;
	struct lp_simplex_Options o;
	double x[2] = {.4, .60001}, r[1], work[1];
	int good = 0, iterations = 0, j;

	if (!m)
		return 0;
	m->constraints[0].type = optm_CONS_T_EQ;
	m->constraints[0].rhs = 1.;
	for (j = 0; j < 2; j++) {
		m->constraints[0].coef[j] = 1.;
		m->bounds[j].b_type = optm_BOUND_T_BS;
		m->bounds[j].lb = 0.;
		m->bounds[j].ub = 1.;
		m->objective[j] = 1.;
	}
	if (simplex_problem_from_model(&p, m))
		return 0;
	if (simplex_pan_bounded_create(&b, &p))
		return 0;
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_PAN_BDA);
	simplex_pan_config_resolve(&controls, b.problem, &o);
	good = simplex_pan_bounded_restore_configured_reported(
		       &b, x, &o, &iterations, 2., &controls, NULL) == 1 &&
	       iterations > 0 &&
	       simplex_pan_bounded_primal(&b, x, r, work) < 1e-12 &&
	       x[0] >= 0. && x[1] <= 1.;
	/* No free columns: the correction must release an inward bound. */
	x[0] = x[1] = 0.;
	iterations = 0;
	good = good &&
	       simplex_pan_bounded_restore_configured_reported(
		       &b, x, &o, &iterations, 2., &controls, NULL) == 1 &&
	       simplex_pan_bounded_primal(&b, x, r, work) < 1e-8 &&
	       x[0] >= 0. && x[1] <= 1.;
	/* A ceiling that excludes every feasible point must roll back exactly. */
	x[0] = .4;
	x[1] = .60001;
	iterations = 0;
	good = good &&
	       simplex_pan_bounded_restore_configured_reported(
		       &b, x, &o, &iterations, .5, &controls, NULL) == -1 &&
	       x[0] == .4 && x[1] == .60001;
	o.iteration_limit = 0;
	iterations = 0;
	good = good &&
	       simplex_pan_bounded_restore_configured_reported(
		       &b, x, &o, &iterations, 2., &controls, NULL) == -1 &&
	       iterations == 0 && x[0] == .4 && x[1] == .60001;
	simplex_pan_bounded_destroy(&b);
	simplex_problem_destroy(&p);
	lp_model_free(m);
	return good;
}

static int multiplier_repair(void)
{
	struct simplex_PanBounded b;
	int starts[3] = {0, 1, 2}, rows[2] = {0, 0}, used = 0;
	double values[2] = {1., 1.}, rhs[1] = {0.};
	double cost[2] = {-2e-7, 0.}, lower[2] = {0., 0.};
	double upper[2] = {HUGE_VAL, HUGE_VAL}, x[2] = {0., 0.};
	double y[1] = {0.}, violation;
	int good;

	memset(&b, 0, sizeof(b));
	b.a.rows = 1;
	b.a.columns = 2;
	b.a.column_start = starts;
	b.a.row_index = rows;
	b.a.value = values;
	b.a.objective = cost;
	b.a.rhs = rhs;
	b.lower = lower;
	b.upper = upper;
	/* Auditing cannot alter a multiplier or consume solver work. */
	good = !simplex_pan_bounded_check_dual(&b, x, y, 1e-7, &violation) &&
	       y[0] == 0.;
	good = good &&
	       !simplex_pan_bounded_repair_dual(
		       &b, x, y, 1e-7, 0, &used, &violation) &&
	       used == 0 && y[0] == 0.;
	good = good &&
	       simplex_pan_bounded_repair_dual(
		       &b, x, y, 1e-7, 2, &used, &violation) &&
	       used == 1 &&
	       simplex_pan_bounded_check_dual(&b, x, y, 1e-7, &violation);
	/* Inconsistent stationarity equations: failed correction must roll
	 * back exactly, while still charging the exhausted search budget. */
	lower[0] = lower[1] = -HUGE_VAL;
	cost[0] = -1.;
	cost[1] = 1.;
	y[0] = .125;
	used = 0;
	good = good &&
	       !simplex_pan_bounded_repair_dual(
		       &b, x, y, 1e-7, 2, &used, &violation) &&
	       used == 2 && y[0] == .125;
	/* A material objective gap is an optimization problem, not polishing. */
	x[0] = 1.;
	rhs[0] = 1.;
	y[0] = 0.;
	used = 0;
	good = good &&
	       !simplex_pan_bounded_repair_dual(
		       &b, x, y, 1e-7, 2, &used, &violation) &&
	       used == 0 && y[0] == 0.;
	return good;
}

static int pdhg_short_budgets(void)
{
	struct simplex_PanConfig controls;
	struct lp_Model *m = lp_model_create(1, 2);
	struct simplex_Problem p;
	struct simplex_PanBounded b;
	struct lp_simplex_Options o;
	double point[2], dual[1], residual[1], work[1], violation;
	int j, budget, used, certificate, good = 1;
	if (!m)
		return 0;
	m->constraints[0].type = optm_CONS_T_EQ;
	m->constraints[0].rhs = 1.;
	for (j = 0; j < 2; j++) {
		m->constraints[0].coef[j] = 1.;
		m->objective[j] = j + 1.;
		m->bounds[j].b_type = optm_BOUND_T_BS;
		m->bounds[j].lb = 0.;
		m->bounds[j].ub = 1.;
	}
	if (simplex_problem_from_model(&p, m)) {
		lp_model_free(m);
		return 0;
	}
	if (simplex_pan_bounded_create(&b, &p)) {
		simplex_problem_destroy(&p);
		lp_model_free(m);
		return 0;
	}
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_PAN_BDA);
	simplex_pan_config_init(&controls, &o);
	/* Both buffer parities and the zero-work exit retain an audited feasible
	 * incumbent. An unfinished iterate must never become an optimum claim. */
	for (budget = 0; budget <= 3; budget++) {
		point[0] = point[1] = .5;
		dual[0] = 0.;
		used = 0;
		certificate = simplex_pan_bounded_pdhg_configured_reported(&b,
			point, dual, &o, &used, budget, &controls, NULL, NULL);
		good = good && used <= budget &&
		       simplex_pan_bounded_primal(&b, point, residual, work) <=
			       o.primal_tolerance &&
		       point[0] + 2. * point[1] <= 1.5 + o.dual_tolerance;
		if (certificate)
			good = good &&
			       simplex_pan_bounded_check_dual(&b, point, dual,
				       o.dual_tolerance, &violation);
	}
	simplex_pan_bounded_destroy(&b);
	simplex_problem_destroy(&p);
	lp_model_free(m);
	return good;
}

int main(void)
{
	unsigned int seed;

	if (!representation() || !kkt() || !free_ray() || !bound_exchange() ||
		!scaled_coordinate() || !local_restore() ||
		!multiplier_repair() || !pdhg_short_budgets())
		return 1;
	for (seed = 1; seed <= 100; seed++)
		if (!known(seed, 0) || !known(seed, 1) || !known(seed, 2))
			return 1;
	printf("Native bounded Pan: model, box KKT, ray, 100 cold/warm/bound-dust witnesses passed\n");
	return 0;
}
