/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file test_pan_pdhg_kernels.c
 * @brief Compare fused PDHG residuals against the compensated reference.
 */

#include "simplex_pan_pdhg.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

static int same_number(double left, double right)
{
	if (isnan(left) || isnan(right))
		return isnan(left) && isnan(right);
	return memcmp(&left, &right, sizeof(double)) == 0;
}

static int compare_residual(const struct simplex_PanStandard *a,
	const struct simplex_PanPdhgColumns *plan, const double *next,
	const double *previous)
{
	double extrapolated[5], expected[3], correction[3];
	double actual[3], actual_correction[3];
	int i;

	for (i = 0; i < a->columns; i++)
		extrapolated[i] = 2. * next[i] - previous[i];
	(void)simplex_pan_residual(a, extrapolated, expected, correction);
	simplex_pan_pdhg_residual(
		a, plan, next, previous, actual, actual_correction);
	for (i = 0; i < a->rows; i++)
		if (!same_number(expected[i], actual[i]) ||
			!same_number(correction[i], actual_correction[i]))
			return 0;
	return 1;
}

static int residual_cases(void)
{
	struct simplex_PanStandard a;
	struct simplex_PanPdhgColumns plan, fallback;
	int starts[6] = {0, 1, 3, 5, 8, 8};
	int rows[8] = {0, 1, 2, 0, 0, 0, 1, 2};
	double values[8] = {1., 1., 1., 1., 1., .25, -3., 1e-12};
	double rhs[3] = {1., -2., 3.};
	double next[5] = {1e16, .5, -1e16, .25, 0.};
	double previous[5] = {0., 1., -1e16, -.25, 0.};
	int signs, i, good = 1;

	memset(&a, 0, sizeof(a));
	memset(&fallback, 0, sizeof(fallback));
	a.rows = 3;
	a.columns = 5;
	a.nonzeros = 8;
	a.column_start = starts;
	a.row_index = rows;
	a.value = values;
	a.rhs = rhs;
	/* Include all short-column signs, repeated row indices, a long generic
	 * column and an empty column. Both plans must retain summation order. */
	for (signs = 0; signs < 32; signs++) {
		for (i = 0; i < 5; i++)
			values[i] = signs & (1 << i) ? -1. : 1.;
		simplex_pan_pdhg_columns_create(&a, &plan);
		good = good && plan.compiled == 3 &&
		       compare_residual(&a, &plan, next, previous) &&
		       compare_residual(&a, &fallback, next, previous) &&
		       compare_residual(&a, &plan, next, next);
		simplex_pan_pdhg_columns_destroy(&plan);
	}
	/* Subnormal products, signed zero and nonfinite extrapolation use the
	 * same reference semantics, including the generic allocation fallback. */
	next[0] = nextafter(0., 1.);
	next[1] = -0.;
	previous[0] = previous[1] = 0.;
	simplex_pan_pdhg_columns_create(&a, &plan);
	good = good && compare_residual(&a, &plan, next, previous) &&
	       compare_residual(&a, &fallback, next, previous);
	next[0] = HUGE_VAL;
	good = good && compare_residual(&a, &plan, next, previous) &&
	       compare_residual(&a, &fallback, next, previous);
	next[0] = NAN;
	good = good && compare_residual(&a, &plan, next, previous);
	simplex_pan_pdhg_columns_destroy(&plan);
	return good;
}

int main(void)
{
	if (!residual_cases()) {
		puts("PDHG residual kernel regression failed");
		return 1;
	}
	puts("PDHG residual kernels preserve compensated reference results");
	return 0;
}
