/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file linalg.c
 * @brief BLAS/LAPACK wrappers and compensated sparse arithmetic.
 */

#include "linalg.h"

#include <assert.h>
#include <math.h>

void lp_simplex_linalg_daxpy(
	int n, double a, double *x, int incx, double *y, int incy)
{
#if USE_BLAS
	extern void daxpy_(int *, double *, double *, int *, double *, int *);

	daxpy_(&n, &a, x, &incx, y, &incy);
#else
	int i = 0, j = 0;

	assert(x != NULL);
	assert(y != NULL);
	while (i < n && j < n) {
		y[i] += a * x[j];
		i += incy;
		j += incx;
	}
#endif
}

double lp_simplex_linalg_ddot(
	int n, const double *x, int incx, const double *y, int incy)
{
#if USE_BLAS
	extern double ddot_(int *, double *, int *, double *, int *);

	return ddot_(&n, (double *)x, &incx, (double *)y, &incy);
#else
	int i = 0, j = 0;
	double result = 0.;

	assert(x != NULL);
	assert(y != NULL);
	while (i < n && j < n) {
		result += x[i] * y[j];
		i += incx;
		j += incy;
	}
	return result;
#endif
}

void lp_simplex_linalg_dscal(int n, double scale, double *array, int increment)
{
#if USE_BLAS
	extern void dscal_(int *, double *, double *, int *);
	dscal_(&n, &scale, array, &increment);
#else
	int i;

	assert(array != NULL);
	for (i = 0; i < n; i += increment)
		array[i] *= scale;
#endif
}

void lp_simplex_linalg_dlarfg(
	int n, double *alpha, double *x, int incx, double *tau)
{
#if USE_LAPACK
	extern void dlarfg_(int *, double *, double *, int *, double *);

	dlarfg_(&n, alpha, x, &incx, tau);
#else
	(void)n;
	(void)alpha;
	(void)x;
	(void)incx;
	(void)tau;
#endif
}

double lp_simplex_linalg_sparse_dot(int first, int end, const int *index,
	const double *value, const double *vector)
{
	double sum = 0., correction = 0.;
	int k;
	int count = end - first;

	/* Signed unit products are exact, including subnormals. For one or two
	 * such terms the compensated expression rounds to the same finite sum.
	 * Keep arbitrary coefficients and longer supports on the full path: a
	 * rounded FMA remainder can change even a one-term result near underflow. */
	if (count == 1 && fabs(value[first]) == 1.) {
		double product = value[first] * vector[index[first]];
		if (isfinite(product))
			return product == 0. ? 0. : product;
	}
	if (count == 2 && fabs(value[first]) == 1. &&
		fabs(value[first + 1]) == 1.) {
		double pair = value[first] * vector[index[first]] +
			      value[first + 1] * vector[index[first + 1]];
		if (isfinite(pair))
			return pair == 0. ? 0. : pair;
	}

	for (k = first; k < end; k++) {
		double a = value[k], x = vector[index[k]];
		double product = a * x;
		double next = sum + product;

		correction +=
			(fabs(sum) >= fabs(product) ? (sum - next) + product
						    : (product - next) + sum) +
			fma(a, x, -product);
		sum = next;
	}
	return sum + correction;
}
