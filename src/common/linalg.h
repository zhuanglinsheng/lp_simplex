/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file linalg.h
 * @brief BLAS/LAPACK wrappers and compensated sparse arithmetic.
 */

#ifndef LP_SIMPLEX_LINALG_INTERNAL_H
#define LP_SIMPLEX_LINALG_INTERNAL_H

void lp_simplex_linalg_daxpy(
	int n, double a, double *x, int incx, double *y, int incy);

double lp_simplex_linalg_ddot(
	int n, const double *x, int incx, const double *y, int incy);

void lp_simplex_linalg_dscal(int n, double x, double *arr, int inc);

void lp_simplex_linalg_dlarfg(
	int n, double *alpha, double *x, int incx, double *tau);

/* Compensated sparse dot over [first, end), including FMA product residuals.
 * Shared arithmetic for model-independent column pricing and certification. */
double lp_simplex_linalg_sparse_dot(int first, int end, const int *index,
	const double *value, const double *vector);

#endif
