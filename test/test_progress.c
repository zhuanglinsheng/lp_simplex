/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file test_progress.c
 * @brief Regression tests for progress.
 */

#include <lp_simplex/lp_simplex.h>
#include "simplex_progress.h"
#include <math.h>
#include <float.h>
#include <stdio.h>
#include <string.h>

#define CHECK(x) \
	do { \
		if (!(x)) { \
			fprintf(stderr, "progress test line %d: %s\n", \
				__LINE__, #x); \
			return 0; \
		} \
	} while (0)

struct Audit {
	const struct lp_Model *model;
	double optimum, lower, upper;
	int calls, points, invalid, pdhg;
	struct lp_simplex_Progress last;
};

static int solve_with_manual_pan_experiments(const struct lp_Model *m,
	const struct lp_simplex_Options *o, double *x,
	struct lp_simplex_Result *r)
{
	return o && o->algorithm == lp_simplex_ALGORITHM_PAN_BDA
		       ? lp_simplex_solve_experimental(m, o, x, r)
		       : lp_simplex_solve(m, o, x, r);
}

static void observe(void *context, const struct lp_simplex_Progress *p)
{
	struct Audit *a = context;
	int i, j;
	a->calls++;
	if (strstr(p->stage, "pdhg"))
		a->pdhg++;
	if (p->lower_bound < a->lower || p->upper_bound > a->upper ||
		p->lower_bound > a->optimum + 1e-7 ||
		p->upper_bound < a->optimum - 1e-7 || isnan(p->relative_gap))
		a->invalid++;
	if ((!isfinite(p->lower_bound) || !isfinite(p->upper_bound)) &&
		isfinite(p->relative_gap))
		a->invalid++;
	if (isfinite(p->relative_gap)) {
		double gap = p->upper_bound - p->lower_bound;
		if (gap < 0. || fabs(gap - p->absolute_gap) > 1e-12 ||
			fabs(p->relative_gap -
				gap / fmax(1., fmax(fabs(p->lower_bound),
						       fabs(p->upper_bound)))) >
				1e-12)
			a->invalid++;
	}
	if (p->point) {
		long double cost = a->model->objective_offset;
		double error = 0.;
		a->points++;
		for (j = 0; j < a->model->n; j++) {
			if (!isfinite(p->point[j]))
				a->invalid++;
			cost += (long double)a->model->objective[j] *
				p->point[j];
			error = fmax(
				error, fmax(-p->point[j], p->point[j] - 1.));
		}
		for (i = 0; i < a->model->m; i++) {
			long double activity = 0.;
			for (j = 0; j < a->model->n; j++)
				activity +=
					(long double)a->model->constraints[i]
						.coef[j] *
					p->point[j];
			if (a->model->constraints[i].type == optm_CONS_T_EQ)
				error = fmax(error,
					(double)fabsl(
						activity -
						a->model->constraints[i].rhs));
			else
				error = fmax(error,
					(double)(a->model->constraints[i].rhs -
						 activity));
		}
		if (error > 1e-7 || p->primal_infeasibility > 1e-7 ||
			fabs((double)cost * a->model->objective_sense -
				(a->model->objective_sense == 1.
						? p->upper_bound
						: p->lower_bound)) > 1e-7)
			a->invalid++;
	}
	a->lower = p->lower_bound;
	a->upper = p->upper_bound;
	a->last = *p;
}

static int integration(void)
{
	int mode, presolve, algorithm;
	for (mode = 0; mode < 3; mode++)
		for (presolve = 0; presolve < 2; presolve++)
			for (algorithm = 0; algorithm < 4; algorithm++) {
				struct lp_Model *m = lp_model_create(3, 3);
				struct lp_simplex_Options o;
				struct lp_simplex_Result r;
				struct Audit a;
				double x[3], cost = mode == 1 ? 4096. : 1.;
				int i, state;
				CHECK(m);
				for (i = 0; i < 3; i++) {
					double scale =
						mode == 1
							? (i == 0 ? 1048576.
								  : i == 1
									  ? .0009765625
									  : 32.)
							: 1.;
					m->constraints[i].type = optm_CONS_T_GE;
					m->constraints[i].rhs = scale;
					m->constraints[i].coef[i] = scale;
					m->constraints[i].coef[(i + 1) % 3] =
						scale;
					m->bounds[i].b_type = optm_BOUND_T_BS;
					m->bounds[i].ub = 1.;
					m->objective[i] =
						mode == 2 ? -cost : cost;
				}
				m->objective_offset = mode == 2 ? -7.25 : 7.25;
				m->objective_sense = mode == 2 ? -1. : 1.;
				memset(&a, 0, sizeof(a));
				a.model = m;
				a.optimum =
					mode == 2 ? 10.25 : 1.5 * cost + 7.25;
				a.lower = -HUGE_VAL;
				a.upper = HUGE_VAL;
				lp_simplex_default_options(&o, algorithm);
				o.presolve = presolve;
				o.progress_callback = observe;
				o.progress_context = &a;
				o.progress_interval = 1e-12;
				state = solve_with_manual_pan_experiments(
					m, &o, x, &r);
				if (state)
					fprintf(stderr,
						"progress integration mode=%d presolve=%d algorithm=%d status=%d primal=%g dual=%g work=%d\n",
						mode, presolve, algorithm,
						r.status,
						r.primal_infeasibility,
						r.dual_infeasibility,
						r.iterations);
				CHECK(state == 0);
				CHECK(fabs(r.objective - a.optimum) < 1e-7);
				CHECK(a.calls >= 2 && a.points >= 1 &&
					!a.invalid);
				CHECK(strcmp(a.last.stage, "finish") == 0);
				/* Tableau does not yet export its row multipliers. */
				if (algorithm != lp_simplex_ALGORITHM_TABLEAU) {
					if (!(a.last.relative_gap < 1e-7))
						fprintf(stderr,
							"mode=%d presolve=%d algorithm=%d lower=%.17g upper=%.17g\n",
							mode, presolve,
							algorithm,
							a.last.lower_bound,
							a.last.upper_bound);
					CHECK(a.last.relative_gap < 1e-7);
				}
				lp_model_free(m);
			}
	return 1;
}

static int bounds(void)
{
	struct lp_Model *m = lp_model_create(1, 2);
	struct simplex_Problem p;
	double y[1], bound;
	int j;
	CHECK(m);
	m->constraints[0].type = optm_CONS_T_GE;
	m->constraints[0].rhs = 1.;
	for (j = 0; j < 2; j++) {
		m->constraints[0].coef[j] = 1.;
		m->objective[j] = j + 1.;
		m->bounds[j].b_type = optm_BOUND_T_BS;
		m->bounds[j].lb = -2.;
		m->bounds[j].ub = 3.;
	}
	CHECK(simplex_problem_from_model(&p, m) == 0);
	y[0] = -1.;
	bound = simplex_progress_dual_bound(&p, y, 1.);
	CHECK(bound <= -1. && bound > -1. - 1e-12);
	y[0] = 100.; /* Wrong GE sign is projected, giving the box bound. */
	bound = simplex_progress_dual_bound(&p, y, 1.);
	CHECK(bound <= -6. && bound > -6. - 1e-12);
	p.row_type[0] = optm_CONS_T_LE;
	y[0] = -100.;
	CHECK(simplex_progress_dual_bound(&p, y, 1.) <= -6.);
	p.row_type[0] = optm_CONS_T_EQ;
	p.bounds[0].b_type = p.bounds[1].b_type = optm_BOUND_T_LO;
	p.bounds[0].lb = p.bounds[1].lb = 0.;
	y[0] = -1.;
	bound = simplex_progress_dual_bound(&p, y, 1.);
	CHECK(isfinite(bound) && bound <= 1. && bound > 1. - 1e-12);
	y[0] = -3.;
	CHECK(simplex_progress_dual_bound(&p, y, 1.) == -HUGE_VAL);
	p.bounds[0].b_type = optm_BOUND_T_FR;
	y[0] = -1.;
	CHECK(isfinite(simplex_progress_dual_bound(&p, y, 1.)));
	y[0] = 0.;
	CHECK(simplex_progress_dual_bound(&p, y, 1.) == -HUGE_VAL);
	y[0] = NAN;
	CHECK(simplex_progress_dual_bound(&p, y, 1.) == -HUGE_VAL);
	/* The product lies just below DBL_MIN; a double-width FMA
	 * residual underflows to negative zero. Never round the bound up. */
	p.bounds[0].b_type = optm_BOUND_T_BS;
	p.bounds[0].lb = p.bounds[0].ub = 1.;
	p.bounds[1].b_type = optm_BOUND_T_BS;
	p.bounds[1].lb = p.bounds[1].ub = 0.;
	p.rhs[0] = p.objective[0] = 0.;
	p.matrix.value[p.matrix.column_start[0]] = DBL_MIN * (1. + DBL_EPSILON);
	y[0] = 1. - DBL_EPSILON;
	CHECK(simplex_progress_dual_bound(&p, y, 1.) < DBL_MIN);

	simplex_problem_destroy(&p);
	lp_model_free(m);
	return 1;
}

/* Presolve substitution changes costs and removes a row. Extending the
 * remaining multipliers by zero may give a weaker bound, but never an
 * optimistic reduced-model bound. The incumbent must be fully postsolved. */
static int postsolve(void)
{
	struct lp_Model *m = lp_model_create(3, 3);
	struct lp_simplex_Options o;
	struct lp_simplex_Result r;
	struct Audit a;
	double x[3];
	int i;
	CHECK(m);
	for (i = 0; i < 3; i++) {
		m->objective[i] = 1.;
		m->bounds[i].b_type = optm_BOUND_T_BS;
		m->bounds[i].ub = 1.;
		m->constraints[i].type = optm_CONS_T_GE;
	}
	m->constraints[0].type = optm_CONS_T_EQ;
	m->constraints[0].coef[0] = m->constraints[0].coef[1] = 1.;
	m->constraints[0].rhs = 1.;
	m->constraints[1].coef[1] = m->constraints[1].coef[2] = 1.;
	m->constraints[1].rhs = 1.;
	m->constraints[2].coef[0] = m->constraints[2].coef[2] = 1.;
	m->constraints[2].rhs = 1.;
	m->objective_offset = 7.25;
	memset(&a, 0, sizeof(a));
	a.model = m;
	a.optimum = 8.75;
	a.lower = -HUGE_VAL;
	a.upper = HUGE_VAL;
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_AUTO);
	o.progress_callback = observe;
	o.progress_context = &a;
	o.progress_interval = 1e-12;
	CHECK(solve_with_manual_pan_experiments(m, &o, x, &r) == 0);
	CHECK(a.points && !a.invalid);
	lp_model_free(m);
	return 1;
}

/* Dyadic data permit an independent exact analytical box minimization.
 * Exercise all row signs and negative bounds/costs, not just cover matrices. */
static int random_bounds(void)
{
	unsigned long seed = 731;
	int sample, i, j;
	for (sample = 0; sample < 200; sample++) {
		struct lp_Model *m = lp_model_create(3, 5);
		struct simplex_Problem p;
		double y[3], projected[3], bound;
		long double expected = 0.;
		CHECK(m);
		for (i = 0; i < 3; i++) {
			seed = seed * 1664525UL + 1013904223UL;
			y[i] = ((int)(seed % 33) - 16) / 4.;
			m->constraints[i].type = i == 0	  ? optm_CONS_T_GE
						 : i == 1 ? optm_CONS_T_LE
							  : optm_CONS_T_EQ;
			projected[i] = i == 0	? fmin(0., y[i])
				       : i == 1 ? fmax(0., y[i])
						: y[i];
			m->constraints[i].rhs = i - 1.;
			expected -= (long double)m->constraints[i].rhs *
				    projected[i];
		}
		for (j = 0; j < 5; j++) {
			long double reduced;
			seed = seed * 1664525UL + 1013904223UL;
			m->objective[j] = ((int)(seed % 17) - 8) / 2.;
			reduced = m->objective[j];
			m->bounds[j].b_type = optm_BOUND_T_BS;
			m->bounds[j].lb = j - 3.;
			m->bounds[j].ub = j - 1.;
			for (i = 0; i < 3; i++) {
				seed = seed * 1664525UL + 1013904223UL;
				m->constraints[i].coef[j] =
					((int)(seed % 17) - 8) / 2.;
				reduced +=
					(long double)m->constraints[i].coef[j] *
					projected[i];
			}
			expected += reduced * (reduced >= 0. ? m->bounds[j].lb
							     : m->bounds[j].ub);
		}
		CHECK(simplex_problem_from_model(&p, m) == 0);
		bound = simplex_progress_dual_bound(&p, y, 1.);
		CHECK(isfinite(bound) && bound <= expected &&
			expected - bound < 1e-10);
		simplex_problem_destroy(&p);
		lp_model_free(m);
	}
	return 1;
}

static int incomplete(void)
{
	struct lp_Model *m = lp_model_create(3, 3);
	struct lp_simplex_Options o;
	struct lp_simplex_Result r;
	struct Audit a;
	double x[3];
	int i;
	CHECK(m);
	for (i = 0; i < 3; i++) {
		m->objective[i] = 1.;
		m->bounds[i].b_type = optm_BOUND_T_BS;
		m->bounds[i].ub = 1.;
		m->constraints[i].type = optm_CONS_T_GE;
		m->constraints[i].rhs = 1.;
		m->constraints[i].coef[i] =
			m->constraints[i].coef[(i + 1) % 3] = 1.;
	}
	memset(&a, 0, sizeof(a));
	a.model = m;
	a.optimum = 1.5;
	a.lower = -HUGE_VAL;
	a.upper = HUGE_VAL;
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_DUAL_REVISED);
	o.presolve = 0;
	o.iteration_limit = 1;
	o.progress_callback = observe;
	o.progress_context = &a;
	CHECK(solve_with_manual_pan_experiments(m, &o, x, &r) != 0);
	CHECK(r.status == lp_simplex_ExceedIterLimit);
	CHECK(a.calls >= 2 && !a.invalid && !a.points);
	CHECK(!isfinite(a.last.relative_gap) && !isfinite(a.last.upper_bound));
	a.calls = 0;
	o.progress_interval = 0.;
	x[0] = 42.;
	CHECK(solve_with_manual_pan_experiments(m, &o, x, &r) != 0);
	CHECK(r.status == lp_simplex_CondUnsatisfied && a.calls == 0 &&
		x[0] == 42.);
	lp_model_free(m);
	return 1;
}

int main(void)
{
	return bounds() && random_bounds() && integration() && postsolve() &&
			       incomplete()
		       ? 0
		       : 1;
}
