/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file test_pan_known_optima.c
 * @brief Independent LP optimality oracle: c=A^T*y+s with s>=0 and s*x*=0.
 *
 * Duplicate rows and zero slacks create degeneracy; binary column scaling
 * varies numerical magnitude without changing the known optimum.
 */

#include <lp_simplex/lp_simplex.h>

#include <math.h>
#include <stdio.h>


static unsigned int next_number(unsigned int *seed)
{
	*seed = *seed * 1664525U + 1013904223U;
	return *seed;
}


static int check_case(unsigned int seed, int pricing)
{
	const int rows = 12, columns = 28;
	struct lp_Model *model = lp_model_create(rows, columns);
	struct lp_simplex_Options options;
	struct lp_simplex_Result result;
	double known[28], x[28], dual[12], expected = 0., maximum = 0.;
	unsigned int generator = seed;
	int i, j, state, valid;

	if (!model)
		return 0;
	for (i = 0; i < rows; i++) {
		model->constraints[i].type = optm_CONS_T_EQ;
		dual[i] = (int)(next_number(&generator) % 5U) - 2.;
	}
	for (j = 0; j < columns; j++) {
		double scale =
			ldexp(1., (int)(next_number(&generator) % 17U) - 8);

		known[j] = j < rows - 2 && (next_number(&generator) % 3U == 0)
				   ? 1. / scale
				   : 0.;
		for (i = 0; i < rows - 2; i++)
			model->constraints[i].coef[j] =
				j < rows - 2
					? (i == j ? scale : 0.)
					: (next_number(&generator) % 4U == 0
							  ? scale * (1. + next_number(
										  &generator) %
										    3U)
							  : 0.);
		if (j >= rows - 2)
			model->constraints[j % (rows - 2)].coef[j] += scale;
		model->constraints[10].coef[j] = model->constraints[0].coef[j];
		model->constraints[11].coef[j] = model->constraints[1].coef[j];
		model->objective[j] =
			known[j] > 0. ? 0.
				      : scale * (next_number(&generator) % 3U);
		for (i = 0; i < rows; i++)
			model->objective[j] +=
				model->constraints[i].coef[j] * dual[i];
	}
	for (i = 0; i < rows; i++) {
		for (j = 0; j < columns; j++)
			model->constraints[i].rhs +=
				model->constraints[i].coef[j] * known[j];
		expected += model->constraints[i].rhs * dual[i];
	}
#ifdef TEST_AUTO
	lp_simplex_default_options(&options, lp_simplex_ALGORITHM_AUTO);
#else
	lp_simplex_default_options(&options, lp_simplex_ALGORITHM_PAN_BDA);
#endif
	options.presolve = 0;
	options.iteration_limit = 5000;
#ifndef TEST_AUTO
	options.pricing = pricing;
#endif
	state = lp_simplex_solve(model, &options, x, &result);
	for (i = 0; i < rows; i++) {
		long double activity = 0.;

		for (j = 0; j < columns; j++)
			activity += (long double)model->constraints[i].coef[j] *
				    x[j];
		maximum = fmax(maximum,
			fabs((double)(activity - model->constraints[i].rhs)));
	}
	for (j = 0; j < columns; j++)
		maximum = fmax(maximum, -x[j]);
	valid = state == lp_simplex_EXIT_SUCCESS &&
		result.status == lp_simplex_Success &&
		isfinite(result.objective) &&
		fabs(result.objective - expected) <
			1e-7 * (1. + fabs(expected)) &&
		maximum < 1e-7;
	if (!valid)
		printf("known LP seed=%u pricing=%d status=%d iterations=%d obj=%.12g expected=%.12g primal=%.6g\n",
			seed, pricing, result.status, result.iterations,
			result.objective, expected, maximum);
	lp_model_free(model);
	return valid;
}


int main(void)
{
	unsigned int seed;

	for (seed = 1; seed <= 100; seed++)
#ifdef TEST_AUTO
		if (!check_case(seed, lp_simplex_PRICING_DUAL_STEEPEST_EDGE))
			return 1;
	printf("AUTO 100 LPs with independent primal-dual optimality witnesses passed\n");
#else
		if (!check_case(seed, lp_simplex_PRICING_DANTZIG) ||
			!check_case(seed, lp_simplex_PRICING_PAN_NORMALIZED))
			return 1;
	printf("Pan 200 degenerate LPs with independent primal-dual optimality witnesses passed\n");
#endif
	return 0;
}
