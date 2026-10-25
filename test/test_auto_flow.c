/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file test_auto_flow.c
 * @brief Failure injection tests the scheduler's promises, not numerical
 * internals.
 */

#include "simplex_auto.h"

#include <lp_simplex/lp_simplex.h>

#include <math.h>
#include <stdio.h>

static int scenario, dual_calls, hybrid_calls, seed_seen, allocated, spent;

static int dual_mock(const struct simplex_Problem *p,
	const struct lp_simplex_Options *o, double *x, double *y,
	struct lp_simplex_Result *r, int propagate)
{
	(void)p;
	(void)y;
	(void)propagate;
	dual_calls++;
	allocated += o->iteration_limit;
	r->dual_residual_space = lp_simplex_RESIDUAL_ORIGINAL;
	r->iterations = o->iteration_limit;
	r->dual_infeasibility = 1.;
	r->status = lp_simplex_ExceedIterLimit;
	x[0] = x[1] = 0.;
	if (scenario == 8)
		r->iterations = 2;
	if (scenario == 0) {
		r->iterations = 1;
		r->status = lp_simplex_Infeasibility;
	}
	if (scenario == 1 || scenario == 6 || scenario == 7 || scenario == 12) {
		r->iterations = 2;
		r->status = lp_simplex_PrecisionError;
		x[0] = x[1] = .5;
	}
	if (scenario == 2 && dual_calls == 2) {
		x[0] = .25;
		x[1] = .75;
	}
	if (scenario == 3) {
		r->iterations = dual_calls == 1 ? 3 : 2;
		r->status = dual_calls == 1 ? lp_simplex_Singularity
					    : lp_simplex_Success;
		if (dual_calls == 2) {
			r->dual_infeasibility = 0.;
			x[0] = 1.;
			x[1] = 0.;
		}
	}
	if (scenario == 5) {
		x[0] = HUGE_VAL;
		r->status = lp_simplex_Success;
		r->dual_infeasibility = 0.;
		r->iterations = 1;
	}
	if (scenario == 9 || scenario == 10 || scenario == 11) {
		r->iterations = 2;
		x[0] = x[1] = .5;
		r->status = scenario == 11 ? lp_simplex_PrecisionError
					   : lp_simplex_Success;
		r->dual_infeasibility = scenario == 10 ? -1. : 1.;
	}
	spent += r->iterations;
	return r->status == lp_simplex_Success ? 0 : -1;
}

static int hybrid_mock(const struct simplex_Problem *p,
	const struct lp_simplex_Options *o, const double *start, double *x,
	struct lp_simplex_Result *r)
{
	(void)p;

	hybrid_calls++;
	allocated += o->iteration_limit;
	r->dual_residual_space = lp_simplex_RESIDUAL_STANDARD;
	seed_seen = start != NULL && fabs(start[0] - .5) < 1e-12 &&
		    fabs(start[1] - .5) < 1e-12;
	r->iterations = o->iteration_limit;
	r->dual_infeasibility = 1.;
	r->status = lp_simplex_ExceedIterLimit;
	x[0] = .75;
	x[1] = .25;
	if (scenario == 1 || scenario == 8) {
		x[0] = 1.;
		x[1] = 0.;
		r->iterations = 3;
		r->status = lp_simplex_Success;
		r->dual_infeasibility = 0.;
	}
	if (scenario == 3) {
		r->iterations = 5;
		r->status = lp_simplex_PrecisionError;
	}
	if (scenario == 5) {
		x[0] = HUGE_VAL;
		r->iterations = 1;
		r->status = lp_simplex_PrecisionError;
	}
	if (scenario == 6 || scenario == 7) {
		x[0] = HUGE_VAL;
		r->iterations = 1;
		r->status = scenario == 6 ? lp_simplex_Infeasibility
					  : lp_simplex_MemoryAllocError;
	}
	if (scenario == 11) {
		r->status = lp_simplex_Success;
		r->dual_infeasibility = NAN;
	}
	if (scenario == 12) {
		x[0] = .25;
		x[1] = .75;
		r->status = lp_simplex_Unboundedness;
	}
	spent += r->iterations;
	return r->status == lp_simplex_Success ? 0 : -1;
}

static int initializer_calls, initializer_invalid;

static int initialize_mock(void *context, const struct simplex_Problem *input,
	const struct lp_simplex_Options *o, double *x, int *work)
{
	int *called_after_dual = (int *)context;
	(void)input;

	initializer_calls++;
	*called_after_dual = dual_calls == 1 && hybrid_calls == 0;
	*work = o->iteration_limit > 4 ? 4 : o->iteration_limit;
	if (initializer_invalid == 2)
		*work = o->iteration_limit;
	spent += *work;
	x[0] = .5;
	x[1] = initializer_invalid == 1 ? HUGE_VAL : .5;
	return 1;
}

int main(void)
{
	struct lp_Model *m = lp_model_create(1, 2);
	struct simplex_Problem p;
	struct lp_simplex_Options o;
	struct lp_simplex_Result r;
	struct simplex_AutoKernels k;
	struct simplex_AutoBindings bindings = {NULL, NULL, 0};
	double x[2];
	int state, ok;

	bindings.kernels = &k;
	if (!m)
		return 1;
	m->constraints[0].type = optm_CONS_T_EQ;
	m->constraints[0].rhs = 1.;
	m->constraints[0].coef[0] = m->constraints[0].coef[1] = 1.;
	m->objective[1] = 1.;
	if (simplex_problem_from_model(&p, m))
		return 1;
	k.dual = dual_mock;
	k.hybrid = hybrid_mock;
	for (scenario = 0; scenario <= 12; scenario++) {
		if (scenario == 8)
			continue;
		dual_calls = hybrid_calls = seed_seen = allocated = spent = 0;
		lp_simplex_default_options(&o, lp_simplex_ALGORITHM_AUTO);
		o.iteration_limit = scenario == 4 ? 1 : 160;
		state = simplex_auto_solve(
			&p, &o, x, &r, 0, &bindings, NULL, NULL);
		ok = r.iterations == spent && spent <= o.iteration_limit;
		if (scenario == 0)
			ok = ok && state != 0 &&
			     r.status == lp_simplex_Infeasibility &&
			     hybrid_calls == 0 && dual_calls == 1;
		if (scenario == 1)
			ok = ok && state == 0 &&
			     r.status == lp_simplex_Success && seed_seen &&
			     r.iterations == 5;
		if (scenario == 2)
			ok = ok && state != 0 &&
			     r.status == lp_simplex_ExceedIterLimit &&
			     fabs(x[1] - .25) < 1e-12 && dual_calls == 2 &&
			     r.dual_residual_space ==
				     lp_simplex_RESIDUAL_STANDARD;
		if (scenario == 3)
			ok = ok && state == 0 &&
			     r.status == lp_simplex_Success &&
			     r.iterations == 10 && dual_calls == 2;
		if (scenario == 4)
			ok = ok && state != 0 &&
			     r.status == lp_simplex_ExceedIterLimit &&
			     hybrid_calls == 0;
		if (scenario == 5)
			ok = ok && state != 0 && r.status != lp_simplex_Success;
		if (scenario == 6)
			ok = ok && state != 0 &&
			     r.status != lp_simplex_Infeasibility &&
			     dual_calls == 2 && fabs(x[1] - .5) < 1e-12;
		if (scenario == 7)
			ok = ok && state != 0 &&
			     r.status == lp_simplex_PrecisionError &&
			     dual_calls == 2 && fabs(x[1] - .5) < 1e-12 &&
			     r.dual_residual_space ==
				     lp_simplex_RESIDUAL_ORIGINAL;
		/* Public residuals may use unscaled input units. A kernel's valid
		 * success cannot be rejected by comparing that diagnostic to a
		 * scaled internal search tolerance. Negative/NaN remain invalid. */
		if (scenario == 9)
			ok = ok && state == 0 &&
			     r.status == lp_simplex_Success &&
			     r.dual_infeasibility == 1.;
		if (scenario == 10 || scenario == 11)
			ok = ok && state != 0 && r.status != lp_simplex_Success;
		if (scenario == 12)
			ok = ok && state != 0 &&
			     r.status == lp_simplex_Unboundedness &&
			     fabs(x[1] - .5) < 1e-12;
		if (!ok) {
			printf("AUTO flow scenario=%d state=%d status=%d spent=%d reported=%d calls=%d/%d seed=%d x1=%g\n",
				scenario, state, r.status, spent, r.iterations,
				dual_calls, hybrid_calls, seed_seen, x[1]);
			return 1;
		}
	}
	/* An alternate initializer runs after dual yield, consumes shared work,
	 * and may hand off only an independently validated point. */
	for (initializer_invalid = 0; initializer_invalid <= 1;
		initializer_invalid++) {
		struct simplex_AutoInitializer initializer;
		int ordered = 0;

		scenario = 8;
		dual_calls = hybrid_calls = seed_seen = spent = 0;
		initializer_calls = 0;
		o.iteration_limit = 160;
		initializer.context = &ordered;
		initializer.run = initialize_mock;
		state = simplex_auto_solve(
			&p, &o, x, &r, 0, &bindings, &initializer, NULL);
		if (state || r.status != lp_simplex_Success ||
			initializer_calls != 1 || !ordered ||
			seed_seen != !initializer_invalid ||
			r.iterations != spent || spent != 9) {
			puts("AUTO alternate initializer contract failed");
			return 1;
		}
	}
	/* A successful feasibility initializer can use its full lease without
	 * starving the optimizer of the work needed to use the resulting seed. */
	{
		struct simplex_AutoInitializer initializer;
		int ordered = 0;
		scenario = 8;
		initializer_invalid = 2;
		dual_calls = hybrid_calls = seed_seen = spent =
			initializer_calls = 0;
		o.iteration_limit = 160;
		initializer.context = &ordered;
		initializer.run = initialize_mock;
		state = simplex_auto_solve(
			&p, &o, x, &r, 0, &bindings, &initializer, NULL);
		if (state || r.status != lp_simplex_Success || !seed_seen ||
			initializer_calls != 1 || r.iterations != spent ||
			spent > 160) {
			puts("AUTO initializer starved optimization");
			return 1;
		}
	}
	simplex_problem_destroy(&p);
	lp_model_free(m);
	puts("AUTO stage transitions, seed, incumbent and cumulative budget passed");
	return 0;
}
