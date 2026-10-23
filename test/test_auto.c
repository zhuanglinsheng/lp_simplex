/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file test_auto.c
 * @brief AUTO contracts: terminal classification, shared budget, and mixed
 * bounds.
 */

#include <lp_simplex/lp_simplex.h>

#include "simplex_auto.h"
#include "simplex_degeneracy.h"
#include "simplex_pan_config.h"
#include "simplex_pan_bounded.h"
#include "simplex_problem.h"

#include <math.h>
#include <stdio.h>
#include <string.h>


static int solve_status(int kind, int presolve, int budget)
{
	struct lp_Model *m = lp_model_create(2, 2);
	struct lp_simplex_Options o;
	struct lp_simplex_Result r;
	double x[2];
	int state, ok;

	if (!m)
		return 0;
	m->constraints[0].type = optm_CONS_T_EQ;
	m->constraints[0].coef[0] = 1.;
	m->constraints[0].coef[1] = 1.;
	m->constraints[0].rhs = 1.;
	m->constraints[1].type = optm_CONS_T_EQ;
	m->constraints[1].coef[0] = 1.;
	m->constraints[1].coef[1] = 1.;
	m->constraints[1].rhs = kind == 1 ? 2. : 1.;
	m->objective[0] = kind == 2 ? -1. : 1.;
	m->objective[1] = 2.;
	if (kind == 2) {
		m->constraints[0].coef[0] = 0.;
		m->constraints[1].coef[0] = 0.;
	}
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_AUTO);
	o.presolve = presolve;
	o.iteration_limit = budget;
	state = lp_simplex_solve(m, &o, x, &r);
	ok = r.iterations >= 0 && r.iterations <= budget;
	if (budget == 1)
		ok = ok && (r.status != lp_simplex_Success ||
				   (state == 0 && isfinite(r.objective) &&
					   r.primal_infeasibility <=
						   o.primal_tolerance));
	else if (kind == 0)
		ok = ok && state == 0 && r.status == lp_simplex_Success &&
		     fabs(r.objective - 1.) < 1e-7;
	else
		ok = ok && state != 0 &&
		     r.status == (kind == 1 ? lp_simplex_Infeasibility
					    : lp_simplex_Unboundedness);
	if (!ok)
		printf("AUTO kind=%d presolve=%d budget=%d state=%d status=%d work=%d obj=%g\n",
			kind, presolve, budget, state, r.status, r.iterations,
			r.objective);
	lp_model_free(m);
	return ok;
}


static int mixed_hybrid(void)
{
	struct lp_Model *m = lp_model_create(2, 3);
	struct simplex_Problem p;
	struct lp_simplex_Options o;
	struct lp_simplex_Result r;
	double x[3], start[3] = {2., -1., 1.};
	int state;

	if (!m)
		return 0;
	m->constraints[0].type = optm_CONS_T_EQ;
	m->constraints[0].rhs = 1.;
	m->constraints[0].coef[0] = 1.;
	m->constraints[0].coef[1] = 1.;
	m->constraints[1].type = optm_CONS_T_EQ;
	m->constraints[1].rhs = 1.;
	m->constraints[1].coef[2] = 1.;
	m->bounds[0].b_type = optm_BOUND_T_BS;
	m->bounds[0].lb = 0.;
	m->bounds[0].ub = 2.;
	m->bounds[1].b_type = optm_BOUND_T_FR;
	m->bounds[2].b_type = optm_BOUND_T_BS;
	m->bounds[2].lb = 1.;
	m->bounds[2].ub = 1.;
	m->objective[0] = 1.;
	m->objective[1] = 2.;
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_AUTO);
	o.iteration_limit = 10000;
	if (simplex_problem_from_model(&p, m)) {
		lp_model_free(m);
		return 0;
	}
	state = simplex_pan_bounded_solve(&p, &o, start, x, &r);
	simplex_problem_destroy(&p);
	lp_model_free(m);
	return state == 0 && r.status == lp_simplex_Success &&
	       fabs(r.objective) < 1e-7 && r.iterations <= o.iteration_limit;
}


static int primal_check(void)
{
	struct lp_Model *m = lp_model_create(1, 3);
	struct simplex_Problem p;
	double x[3] = {1., 1., 1.};
	int ok;

	if (!m)
		return 0;
	m->constraints[0].type = optm_CONS_T_EQ;
	m->constraints[0].rhs = 1.;
	m->constraints[0].coef[0] = 1e16;
	m->constraints[0].coef[1] = 1.;
	m->constraints[0].coef[2] = -1e16;
	if (simplex_problem_from_model(&p, m)) {
		lp_model_free(m);
		return 0;
	}
	ok = simplex_problem_primal_infeasibility(&p, x) == 0.;
	x[0] = HUGE_VAL;
	ok = ok && !isfinite(simplex_problem_primal_infeasibility(&p, x));
	x[0] = 1.;
	x[1] = -1.;
	ok = ok && simplex_problem_primal_infeasibility(&p, x) >= 1.;
	simplex_problem_destroy(&p);
	lp_model_free(m);
	return ok;
}


int main(void)
{
	int kind, pre;
	struct lp_simplex_Options o;
	struct simplex_PanConfig config;

	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_AUTO);
	if (o.pricing != lp_simplex_PRICING_DUAL_STEEPEST_EDGE || !o.presolve)
		return 1;
	simplex_pan_config_init(&config, &o);
	if (!config.implicit_cone || !config.implicit_updates ||
		!config.restore || config.dual_polish_iterations != 8192 ||
		config.pdhg_iterations != 500000 || config.predict ||
		config.implicit_bulk || config.full_refinement ||
		config.dual_center || config.newton_qr ||
		config.nnls_block_limit != 1 || config.block_crash ||
		!config.singleton_crash)
		return 1;
	for (pre = 0; pre <= 1; pre++)
		for (kind = 0; kind < 3; kind++)
			if (!solve_status(kind, pre, 10000) ||
				!solve_status(kind, pre, 1))
				return 1;
	if (!mixed_hybrid() || !primal_check())
		return 1;
	puts("AUTO terminal, budget and mixed-bound hybrid contracts passed");
	return 0;
}
