/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_pan_basis_qr.c
 * @brief Explicit numerical recovery on immutable active columns; no LP
 * staging.
 */

#include <lp_simplex/status.h>
#include "simplex_pan_basis.h"
#include "policy/resources.h"
#include "policy/controls.h"
#include "utils.h"
#include <math.h>
#include <float.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>

/* Fill (B D^-1)^T and D^-1 right, without forming B^T B. Row scaling
 * preserves the exact equations and therefore their minimum-norm solution. */
static int pan_minimum_norm_qr_input(const struct simplex_PanBasis *basis,
	const double *right, double *matrix, double *rhs, int *pivot)
{
	int i, k, rows = basis->standard->rows, count = basis->count;

	lp_simplex_memset(matrix, 0, (size_t)rows * count * sizeof(double));
	lp_simplex_memset(rhs, 0, (size_t)rows * sizeof(double));
	lp_simplex_memset(pivot, 0, (size_t)rows * sizeof(int));
	for (i = 0; i < count; i++) {
		int column = basis->column[i];
		double norm = simplex_pan_column_norm(basis->standard, column);

		if (!isfinite(norm) || norm <= 0. || !isfinite(right[i]))
			return lp_simplex_EXIT_FAILURE;
		rhs[i] = right[i] / norm;
		if (!isfinite(rhs[i]))
			return lp_simplex_EXIT_FAILURE;
		for (k = basis->standard->column_start[column];
			k < basis->standard->column_start[column + 1]; k++) {
			int position =
				i + basis->standard->row_index[k] * count;
			double value = basis->standard->value[k] / norm;

			if (!isfinite(value))
				return lp_simplex_EXIT_FAILURE;
			matrix[position] += value;
			if (!isfinite(matrix[position]))
				return lp_simplex_EXIT_FAILURE;
		}
	}
	return lp_simplex_EXIT_SUCCESS;
}

int simplex_pan_basis_minimum_norm_qr(struct simplex_PanBasis *basis,
	const double *right, double *solution, double *residual_norm)
{
	int count = basis->count, rows = basis->standard->rows;
	int one = 1, rank = 0, info = 0, lwork = -1;
	int *pivot = NULL;
	double *matrix = NULL, *rhs = NULL, *work = NULL;
	/* The working basis has already passed its admission test. Resolve its
	 * equations to arithmetic precision; do not apply the user pivot policy
	 * again as DGELSY's global condition-number cutoff. Original residuals
	 * and the caller's unchanged dual tolerance still decide acceptance. */
	double work_size = 0., rcond = DBL_EPSILON, error;
	long double bytes;
	int i;

	extern void dgelsy_(int *, int *, int *, double *, int *, double *,
		int *, int *, double *, int *, double *, int *, int *);
	basis->minimum_norm_qr_calls++;
	*residual_norm = HUGE_VAL;
	if (count < 0 || count > rows || rows < 0)
		return lp_simplex_EXIT_FAILURE;
	if (count == 0) {
		lp_simplex_memset(solution, 0, (size_t)rows * sizeof(double));
		*residual_norm = 0.;
		return lp_simplex_EXIT_SUCCESS;
	}
	/* LAPACK's minimum LWORK includes 3*N+M+1. Reject dimensions before
	 * its integer arithmetic, and reject scratch beyond this local lease. */
	if (rows > (INT_MAX - 1) / 4 ||
		(size_t)rows > (size_t)-1 / sizeof(double) / count)
		return lp_simplex_EXIT_FAILURE;
	bytes = ((long double)rows * count + rows) * sizeof(double) +
		(long double)rows * sizeof(int);
	if (bytes > SIMPLEX_PAN_MINIMUM_NORM_QR_BYTES)
		return lp_simplex_EXIT_FAILURE;
	matrix = (double *)lp_simplex_malloc(
		(size_t)rows * count * sizeof(double));
	rhs = (double *)lp_simplex_malloc((size_t)rows * sizeof(double));
	pivot = (int *)lp_simplex_malloc((size_t)rows * sizeof(int));
	if (!matrix || !rhs || !pivot ||
		pan_minimum_norm_qr_input(basis, right, matrix, rhs, pivot))
		goto failure;
	dgelsy_(&count, &rows, &one, matrix, &count, rhs, &rows, pivot, &rcond,
		&rank, &work_size, &lwork, &info);
	if (info != 0 || !isfinite(work_size) || work_size < 1. ||
		work_size > INT_MAX)
		goto failure;
	lwork = (int)ceil(work_size);
	if (bytes + (long double)lwork * sizeof(double) >
		SIMPLEX_PAN_MINIMUM_NORM_QR_BYTES)
		goto failure;
	work = (double *)lp_simplex_malloc((size_t)lwork * sizeof(double));
	if (!work ||
		pan_minimum_norm_qr_input(basis, right, matrix, rhs, pivot))
		goto failure;
	dgelsy_(&count, &rows, &one, matrix, &count, rhs, &rows, pivot, &rcond,
		&rank, work, &lwork, &info);
	if (getenv(SIMPLEX_ENV_PROFILE))
		fprintf(stderr,
			"pan-bda multiplier QR: info=%d rank=%d/%d rcond=%.6g\n",
			info, rank, count, rcond);
	if (info != 0 || rank != count)
		goto failure;
	for (i = 0; i < rows; i++)
		if (!isfinite(rhs[i]))
			goto failure;
	error = 0.;
	for (i = 0; i < count; i++) {
		double residual =
			right[i] - simplex_pan_column_dot(basis->standard,
					   basis->column[i], rhs);
		if (!isfinite(residual))
			goto failure;
		error = hypot(error, residual);
	}
	if (!isfinite(error))
		goto failure;
	if (getenv(SIMPLEX_ENV_PROFILE))
		fprintf(stderr,
			"pan-bda multiplier QR: original-residual=%.17g\n",
			error);
	lp_simplex_memcpy(solution, rhs, (size_t)rows * sizeof(double));
	*residual_norm = error;
	lp_simplex_free(matrix);
	lp_simplex_free(rhs);
	lp_simplex_free(pivot);
	lp_simplex_free(work);
	return lp_simplex_EXIT_SUCCESS;
failure:
	lp_simplex_free(matrix);
	lp_simplex_free(rhs);
	lp_simplex_free(pivot);
	lp_simplex_free(work);
	return lp_simplex_EXIT_FAILURE;
}
