/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#include "linalg.h"

#include <assert.h>


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
