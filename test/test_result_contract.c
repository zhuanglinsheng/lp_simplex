/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file test_result_contract.c
 * @brief Public result units, model space and success gate across every entry.
 */

#include <lp_simplex/lp_simplex.h>

#include "simplex_problem.h"
#include "simplex_dual.h"
#include "simplex_execution.h"
#include <string.h>

#include <math.h>
#include <stdio.h>

static int solve_with_manual_pan_experiments(const struct lp_Model *m,
	const struct lp_simplex_Options *o, double *x,
	struct lp_simplex_Result *r)
{
	return o && o->algorithm == lp_simplex_ALGORITHM_PAN_BDA
		       ? lp_simplex_solve_experimental(m, o, x, r)
		       : lp_simplex_solve(m, o, x, r);
}

static int check_case(int algorithm, int pre, int maximize, int budget)
{
	struct lp_Model *m = lp_model_create(1, 2);
	struct lp_simplex_Options o;
	struct lp_simplex_Result r;
	double x[2], expected, primal;
	int state, ok, point;

	if (!m)
		return 0;
	m->constraints[0].type = optm_CONS_T_EQ;
	m->constraints[0].rhs = 1.;
	m->constraints[0].coef[0] = m->constraints[0].coef[1] = 1.;
	m->objective_sense = maximize ? -1 : 1;
	m->objective[0] = maximize ? -1. : 1.;
	m->objective[1] = maximize ? -2. : 2.;
	m->objective_offset = maximize ? -3. : 3.;
	lp_simplex_default_options(&o, algorithm);
	o.presolve = pre;
	o.iteration_limit = budget;
	state = solve_with_manual_pan_experiments(m, &o, x, &r);
	point = isfinite(x[0]) && isfinite(x[1]);
	ok = r.iterations >= 0 && r.iterations <= budget;
	if (point) {
		expected = x[0] + 2. * x[1] + 3.;
		primal = fmax(
			fabs(x[0] + x[1] - 1.), fmax(0., fmax(-x[0], -x[1])));
		ok = ok && isfinite(r.objective) &&
		     fabs(r.objective - expected) < 1e-10 &&
		     fabs(r.primal_infeasibility - primal) < 1e-10;
	} else
		ok = ok && isnan(r.objective) && isinf(r.primal_infeasibility);
	if (r.status == lp_simplex_Success)
		ok = ok && state == 0 && point &&
		     r.primal_infeasibility <= o.primal_tolerance;
	else
		ok = ok && state != 0;
	if (budget > 1)
		ok = ok && state == 0 &&
		     fabs(r.objective - (maximize ? 5. : 4.)) < 1e-7;
	if (r.dual_residual_space == lp_simplex_RESIDUAL_UNAVAILABLE)
		ok = ok && isnan(r.dual_infeasibility);
	else
		ok = ok && !isnan(r.dual_infeasibility);
	if (algorithm == lp_simplex_ALGORITHM_TABLEAU)
		ok = ok &&
		     r.dual_residual_space == lp_simplex_RESIDUAL_UNAVAILABLE;
	if (!ok)
		printf("result contract algorithm=%d pre=%d max=%d budget=%d status=%d state=%d x=%g,%g obj=%g primal=%g dual=%g space=%d work=%d\n",
			algorithm, pre, maximize, budget, r.status, state, x[0],
			x[1], r.objective, r.primal_infeasibility,
			r.dual_infeasibility, r.dual_residual_space,
			r.iterations);
	lp_model_free(m);
	return ok;
}

/* No point can satisfy both rows within 1e-7: a 10x success gate is invalid. */
static int strict_primal_gate(int algorithm, int pre)
{
	struct lp_Model *m = lp_model_create(2, 1);
	struct lp_simplex_Options o;
	struct lp_simplex_Result r;
	double x[1];
	int state, ok;

	if (!m)
		return 0;
	m->constraints[0].type = m->constraints[1].type = optm_CONS_T_EQ;
	m->constraints[0].coef[0] = m->constraints[1].coef[0] = 1.;
	m->constraints[0].rhs = 1.;
	m->constraints[1].rhs = 1. + 5e-7;
	lp_simplex_default_options(&o, algorithm);
	o.presolve = pre;
	o.iteration_limit = 2000;
	state = solve_with_manual_pan_experiments(m, &o, x, &r);
	ok = state != 0 && r.status != lp_simplex_Success;
	lp_model_free(m);
	return ok;
}

static int dual_units(void)
{
	struct lp_Model *m = lp_model_create(1, 1);
	struct simplex_Problem p;
	struct lp_simplex_Options o;
	struct lp_simplex_Result r;
	double x[1], y[1], bad[1];
	int ok, state;

	if (!m)
		return 0;
	m->objective[0] = -3.;
	m->constraints[0].type = optm_CONS_T_LE;
	m->constraints[0].coef[0] = 1048576.;
	m->constraints[0].rhs = 1048576.;
	if (simplex_problem_from_model(&p, m))
		return 0;
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_DUAL_REVISED);
	o.presolve = 0;
	state = simplex_dual_solve_problem(&p, &o, x, y, &r, 0);
	ok = state == 0 && fabs(x[0] - 1.) < 1e-7 &&
	     fabs(y[0] + 3. / 1048576.) < 1e-12 &&
	     r.dual_residual_space == lp_simplex_RESIDUAL_ORIGINAL &&
	     r.dual_infeasibility < 1e-7;
	bad[0] = -1. / 1048576.;
	ok = ok &&
	     fabs(simplex_problem_dual_infeasibility(&p, x, bad, 1e-7, -1) -
		     2.) < 1e-12;
	x[0] = .5;
	ok = ok && simplex_problem_dual_infeasibility(&p, x, y, 1e-7, -1) > 0.;
	state = simplex_dual_find_feasible(&p, &o, x, &r);
	ok = ok && state == 0 && r.status == lp_simplex_CondUnsatisfied &&
	     r.dual_residual_space == lp_simplex_RESIDUAL_AUXILIARY;
	simplex_problem_destroy(&p);
	lp_model_free(m);
	return ok;
}

static int reduced_space(void)
{
	struct lp_Model *m = lp_model_create(2, 4);
	struct lp_simplex_Options o;
	struct lp_simplex_Result r;
	double x[4];
	int j, state, ok;

	if (!m)
		return 0;
	m->constraints[0].type = m->constraints[1].type = optm_CONS_T_LE;
	m->constraints[0].rhs = 4.;
	m->constraints[1].rhs = 12.;
	for (j = 0; j < 4; j++) {
		m->objective[j] = -1. - j;
		m->constraints[0].coef[j] = 1.;
		m->constraints[1].coef[j] = j == 0 ? 2. : 2. * j + 1.;
	}
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_DUAL_REVISED);
	state = solve_with_manual_pan_experiments(m, &o, x, &r);
	ok = state == 0 &&
	     r.dual_residual_space == lp_simplex_RESIDUAL_PRESOLVED;
	if (!ok)
		printf("reduced scope status=%d scope=%d\n", r.status,
			r.dual_residual_space);
	lp_model_free(m);
	return ok;
}

static int invalid_options(void)
{
	struct lp_Model *m = lp_model_create(1, 1);
	struct lp_simplex_Options o;
	struct lp_simplex_Result r;
	double x[1] = {7.};
	int state, ok;

	if (!m)
		return 0;
	m->constraints[0].type = optm_CONS_T_LE;
	m->constraints[0].rhs = 1.;
	m->constraints[0].coef[0] = 1.;
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_AUTO);
	o.primal_tolerance = NAN;
	state = solve_with_manual_pan_experiments(m, &o, x, &r);
	ok = state != 0 && isnan(r.objective) &&
	     isinf(r.primal_infeasibility) && isnan(r.dual_infeasibility) &&
	     r.dual_residual_space == lp_simplex_RESIDUAL_UNAVAILABLE &&
	     x[0] == 7.;
	lp_model_free(m);
	return ok;
}

static int invalid_callbacks;

static void invalid_progress(void *context, const struct lp_simplex_Progress *p)
{
	(void)context;
	(void)p;
	invalid_callbacks++;
}

static int invalid_model_data(void)
{
	struct lp_Model *m;
	struct optm_VariableBound *saved_bounds;
	struct lp_simplex_Options o;
	struct lp_simplex_Result r;
	struct simplex_Execution h;
	double x;
	int algorithm, pre, form, kind, code, valid;
	for (algorithm = lp_simplex_ALGORITHM_TABLEAU;
		algorithm <= lp_simplex_ALGORITHM_AUTO; algorithm++)
		for (pre = 0; pre < 2; pre++)
			for (form = 0; form < 3; form++)
				for (kind = 0; kind < 19; kind++) {
					m = form == 2 ? lp_model_create_sparse(
								1, 1, 1)
						      : lp_model_create(1, 1);
					if (!m)
						return 0;
					m->constraints[0].type = optm_CONS_T_EQ;
					m->constraints[0].rhs =
						m->objective[0] = 1.;
					if (form == 2) {
						m->nnz = 1;
						m->column_start[0] =
							m->row_start[0] = 0;
						m->column_start[1] =
							m->row_start[1] = 1;
						m->row_index[0] =
							m->column_index[0] = 0;
						m->value[0] = m->row_value[0] =
							1.;
					} else {
						m->constraints[0].coef[0] = 1.;
						if (form &&
							lp_model_build_sparse(
								m)) {
							lp_model_free(m);
							return 0;
						}
					}
					saved_bounds = m->bounds;
					switch (kind) {
					case 0:
						m->objective[0] = NAN;
						break;
					case 1:
						m->objective[0] = HUGE_VAL;
						break;
					case 2:
						m->constraints[0].rhs = NAN;
						break;
					case 3:
						m->constraints[0].rhs =
							-HUGE_VAL;
						break;
					case 4:
					case 5:
						if (form)
							m->value[0] =
								kind == 4
									? NAN
									: HUGE_VAL;
						else
							m->constraints[0]
								.coef[0] =
								kind == 4
									? NAN
									: HUGE_VAL;
						break;
					case 6:
					case 7:
						if (form)
							m->row_value[0] =
								kind == 6
									? NAN
									: -HUGE_VAL;
						else
							m->constraints[0]
								.coef[0] =
								kind == 6
									? NAN
									: -HUGE_VAL;
						break;
					case 8:
						m->objective_offset = NAN;
						break;
					case 9:
						m->objective_offset = HUGE_VAL;
						break;
					case 10:
						m->objective_sense = 0.;
						break;
					case 11:
						m->objective_sense = NAN;
						break;
					case 12:
						m->bounds[0].lb = NAN;
						break;
					case 13:
						m->bounds[0].b_type =
							optm_BOUND_T_BS;
						m->bounds[0].ub = NAN;
						break;
					case 14:
						m->bounds[0].lb = HUGE_VAL;
						break;
					case 15:
						m->bounds[0].b_type =
							optm_BOUND_T_UP;
						m->bounds[0].ub = -HUGE_VAL;
						break;
					case 16:
						if (form == 2)
							m->objective_sense = 2.;
						else
							m->constraints[0]
								.coef[0] = NAN;
						break;
					case 17:
						m->bounds = NULL;
						m->objective[0] = NAN;
						break;
					case 18:
						m->bounds = NULL;
						m->constraints[0].rhs =
							HUGE_VAL;
						break;
					}
					lp_simplex_default_options(
						&o, algorithm);
					o.presolve = pre;
					o.progress_callback = invalid_progress;
					invalid_callbacks = 0;
					x = 17.;
					memset(&r, 0xff, sizeof(r));
					code = simplex_solve_reported(
						m, &o, &x, &r, &h);
					valid = code &&
						r.status ==
							lp_simplex_CondUnsatisfied &&
						!r.iterations && x == 17. &&
						isnan(r.objective) &&
						isnan(r.dual_infeasibility) &&
						isinf(r.primal_infeasibility) &&
						r.dual_residual_space ==
							lp_simplex_RESIDUAL_UNAVAILABLE &&
						!invalid_callbacks &&
						h.count == 1;
					m->bounds = saved_bounds;
					lp_model_free(m);
					if (!valid) {
						fprintf(stderr,
							"invalid model accepted: alg=%d pre=%d form=%d kind=%d status=%d callbacks=%d\n",
							algorithm, pre, form,
							kind, r.status,
							invalid_callbacks);
						return 0;
					}
				}
	return 1;
}

static int inactive_bound_fields(void)
{
	struct lp_Model *m = lp_model_create(1, 1);
	struct lp_simplex_Options o;
	struct lp_simplex_Result r;
	double x;
	int algorithm, form, state, valid;
	if (!m)
		return 0;
	m->constraints[0].type = optm_CONS_T_EQ;
	m->constraints[0].coef[0] = m->constraints[0].rhs = m->objective[0] =
		1.;
	for (algorithm = lp_simplex_ALGORITHM_TABLEAU;
		algorithm <= lp_simplex_ALGORITHM_AUTO; algorithm++)
		for (form = 0; form < 3; form++) {
			m->bounds[0].b_type = form == 0	  ? optm_BOUND_T_LO
					      : form == 1 ? optm_BOUND_T_UP
							  : optm_BOUND_T_FR;
			m->bounds[0].lb = form == 0 ? 0. : NAN;
			m->bounds[0].ub = form == 1 ? 2. : NAN;
			lp_simplex_default_options(&o, algorithm);
			state = lp_simplex_solve(m, &o, &x, &r);
			valid = !state && r.status == lp_simplex_Success &&
				x == 1. && r.objective == 1.;
			if (!valid) {
				lp_model_free(m);
				return 0;
			}
		}
	lp_model_free(m);
	return 1;
}

int main(void)
{
	int algorithm, pre, maximize;

	for (algorithm = lp_simplex_ALGORITHM_TABLEAU;
		algorithm <= lp_simplex_ALGORITHM_AUTO; algorithm++)
		for (pre = 0; pre <= 1; pre++) {
			for (maximize = 0; maximize <= 1; maximize++)
				if (!check_case(
					    algorithm, pre, maximize, 10000) ||
					!check_case(
						algorithm, pre, maximize, 1))
					return 1;
			if (!strict_primal_gate(algorithm, pre))
				return 1;
		}
	if (!inactive_bound_fields() || !invalid_model_data() ||
		!dual_units() || !reduced_space() || !invalid_options())
		return 1;
	puts("Original objective/primal, strict gate, partial points and dual model spaces passed");
	return 0;
}
