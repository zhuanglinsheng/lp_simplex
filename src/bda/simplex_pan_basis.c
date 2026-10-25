/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_pan_basis.c
 * @brief Dynamic deficient-basis storage and numerical operations for Pan/BDA.
 */

#include "simplex_pan_basis.h"
#include "simplex_pan_config.h"
#include "utils.h"
#include "policy/resources.h"

#include <lp_simplex/status.h>

#include <math.h>
#include <limits.h>
#include <stdlib.h>

#ifndef LP_SIMPLEX_HAVE_SPQR

int simplex_pan_basis_select_columns(struct simplex_PanBasis *basis,
	const int *column, int count, const double *required)
{
	int i;

	if (basis->count != 0 || count < 0 || count > basis->standard->rows)
		return lp_simplex_EXIT_FAILURE;
	for (i = 0; i < count; i++)
		if (simplex_pan_basis_add(basis, column[i]))
			return lp_simplex_EXIT_FAILURE;
	return simplex_pan_basis_repair_rank(basis, required);
}

int simplex_pan_dual_newton_qr(const struct simplex_PanStandard *a,
	const double *point, const unsigned char *free_mask, const double *map,
	const double *multiplier, double scale, double *correction)
{
	(void)a;
	(void)point;
	(void)free_mask;
	(void)map;
	(void)multiplier;
	(void)scale;
	(void)correction;
	return lp_simplex_EXIT_FAILURE;
}

int simplex_pan_basis_pricing_probes(struct simplex_PanBasis *basis,
	const double *value, double *span_probe, double *edge_probe)
{
	(void)basis;
	(void)value;
	(void)span_probe;
	(void)edge_probe;
	return lp_simplex_EXIT_FAILURE;
}

int simplex_pan_basis_stabilize(struct simplex_PanBasis *basis,
	const double *value, const int protected_column)
{
	(void)basis;
	(void)value;
	(void)protected_column;
	return lp_simplex_EXIT_FAILURE;
}

int simplex_pan_basis_complete(
	struct simplex_PanBasis *basis, const double *value)
{
	(void)basis;
	(void)value;
	return lp_simplex_EXIT_FAILURE;
}

int simplex_pan_basis_repair_rank(
	struct simplex_PanBasis *basis, const double *value)
{
	(void)value;

	return simplex_pan_basis_rebuild(basis);
}

int simplex_pan_basis_rebuild(struct simplex_PanBasis *basis)
{
	return simplex_pan_basis_factorize(basis);
}

static double pan_sparse_column_dot(const struct simplex_PanStandard *standard,
	const int first, const int second)
{
	int a = standard->column_start[first];
	int b = standard->column_start[second];
	int a_end = standard->column_start[first + 1];
	int b_end = standard->column_start[second + 1];
	double dot = 0.;

	while (a < a_end && b < b_end) {
		int row_a = standard->row_index[a];
		int row_b = standard->row_index[b];

		if (row_a < row_b)
			a++;
		else if (row_b < row_a)
			b++;
		else {
			dot += standard->value[a] * standard->value[b];
			a++;
			b++;
		}
	}
	return dot;
}

double simplex_pan_column_norm(
	const struct simplex_PanStandard *standard, const int column)
{
	return standard->column_norm[column];
}

static void pan_column_axpy(const struct simplex_PanStandard *standard,
	const int column, const double scale, double *dense)
{
	int k;

	for (k = standard->column_start[column];
		k < standard->column_start[column + 1]; k++)
		dense[standard->row_index[k]] += scale * standard->value[k];
}

static double pan_column_dense_dot(const struct simplex_PanStandard *standard,
	const int column, const double *dense)
{
	return simplex_pan_column_dot(standard, column, dense);
}

static double pan_dense_norm(const double *vector, const int count)
{
	double scale = 0.;
	double sum = 1.;
	int i;

	for (i = 0; i < count; i++) {
		double value = __lp_simplex_ABS__(vector[i]);

		if (!isfinite(value))
			return HUGE_VAL;
		if (value == 0.)
			continue;
		if (scale < value) {
			double ratio = scale / value;

			sum = 1. + sum * ratio * ratio;
			scale = value;
		} else {
			double ratio = value / scale;

			sum += ratio * ratio;
		}
	}
	return scale == 0. ? 0. : scale * sqrt(sum);
}

static int pan_basis_reserve(struct simplex_PanBasis *basis, const int capacity)
{
	double *factor;
	double *norm;
	double *work;
	int *columns;

	if (capacity <= basis->capacity)
		return lp_simplex_EXIT_SUCCESS;
	columns = (int *)lp_simplex_malloc((size_t)capacity * sizeof(int));
	factor = (double *)lp_simplex_malloc(
		(size_t)capacity * capacity * sizeof(double));
	norm = (double *)lp_simplex_malloc((size_t)capacity * sizeof(double));
	work = (double *)lp_simplex_malloc((size_t)capacity * sizeof(double));
	if (columns == NULL || factor == NULL || norm == NULL || work == NULL) {
		lp_simplex_free(columns);
		lp_simplex_free(factor);
		lp_simplex_free(norm);
		lp_simplex_free(work);
		return lp_simplex_EXIT_FAILURE;
	}
	if (basis->count > 0)
		lp_simplex_memcpy(columns, basis->column,
			(size_t)basis->count * sizeof(int));
	lp_simplex_free(basis->column);
	lp_simplex_free(basis->factor);
	lp_simplex_free(basis->column_norm);
	lp_simplex_free(basis->basis_work);
	basis->column = columns;
	basis->factor = factor;
	basis->column_norm = norm;
	basis->basis_work = work;
	basis->capacity = capacity;
	return lp_simplex_EXIT_SUCCESS;
}

int simplex_pan_basis_create_configured(struct simplex_PanBasis *basis,
	const struct simplex_PanStandard *standard, const double rank_tolerance,
	const struct simplex_PanConfig *config)
{

	lp_simplex_memset(basis, 0, sizeof(*basis));
	basis->standard = standard;
	simplex_pan_basis_configure(basis, config);
	basis->rank_tolerance = rank_tolerance;
	basis->position = (int *)lp_simplex_malloc(
		(size_t)(standard->columns > 0 ? standard->columns : 1) *
		sizeof(int));
	basis->row_work = (double *)lp_simplex_malloc(
		(size_t)standard->rows * sizeof(double));
	basis->right_work = (double *)lp_simplex_malloc(
		(size_t)standard->rows * sizeof(double));
	if (basis->position == NULL || basis->row_work == NULL ||
		basis->right_work == NULL ||
		pan_basis_reserve(
			basis, standard->rows > 8 ? 8 : standard->rows) ==
			lp_simplex_EXIT_FAILURE) {
		simplex_pan_basis_destroy(basis);
		return lp_simplex_EXIT_FAILURE;
	}
	if (standard->columns > 0)
		lp_simplex_memset(basis->position, 0xff,
			(size_t)standard->columns * sizeof(int));
	return lp_simplex_EXIT_SUCCESS;
}

void simplex_pan_basis_destroy(struct simplex_PanBasis *basis)
{
	if (basis == NULL)
		return;
	lp_simplex_free(basis->column);
	lp_simplex_free(basis->position);
	lp_simplex_free(basis->factor);
	lp_simplex_free(basis->column_norm);
	lp_simplex_free(basis->basis_work);
	lp_simplex_free(basis->row_work);
	lp_simplex_free(basis->right_work);
	lp_simplex_memset(basis, 0, sizeof(*basis));
}

int simplex_pan_basis_add(struct simplex_PanBasis *basis, const int column)
{
	int capacity;

	if (column < 0 || column >= basis->standard->columns ||
		basis->position[column] >= 0 ||
		basis->count >= basis->standard->rows)
		return lp_simplex_EXIT_FAILURE;
	if (basis->count == basis->capacity) {
		capacity = basis->capacity > 0 ? 2 * basis->capacity : 8;
		if (capacity > basis->standard->rows)
			capacity = basis->standard->rows;
		if (pan_basis_reserve(basis, capacity) ==
			lp_simplex_EXIT_FAILURE)
			return lp_simplex_EXIT_FAILURE;
	}
	basis->position[column] = basis->count;
	basis->column[basis->count++] = column;
	return lp_simplex_EXIT_SUCCESS;
}

void simplex_pan_basis_remove_marked(
	struct simplex_PanBasis *basis, const unsigned char *remove)
{
	int destination = 0;
	int i;

	for (i = 0; i < basis->count; i++) {
		int column = basis->column[i];

		if (remove[i]) {
			basis->position[column] = -1;
			continue;
		}
		basis->column[destination] = column;
		basis->position[column] = destination++;
	}
	basis->count = destination;
}

/* Normal equations square the condition number. A failed Gram factor is
 * not evidence that the admitted original columns are dependent. Validate
 * them with bounded QR scratch and retain that backend for this working set. */
static int pan_basis_recover_qr(struct simplex_PanBasis *basis)
{
	double error;

	lp_simplex_memset(basis->right_work, 0,
		(size_t)basis->standard->rows * sizeof(double));
	basis->factorizations++;
	if (simplex_pan_basis_least_squares_qr(
		    basis, basis->right_work, basis->basis_work, &error))
		return lp_simplex_EXIT_FAILURE;
	basis->dense_qr = 1;
	return lp_simplex_EXIT_SUCCESS;
}

int simplex_pan_basis_factorize(struct simplex_PanBasis *basis)
{
	char lower = 'L';
	int info = 0;
	int i, j, count = basis->count;

	extern void dpotrf_(char *, int *, double *, int *, int *);
	basis->dense_qr = 0;
	if (count == 0)
		return lp_simplex_EXIT_SUCCESS;
	for (i = 0; i < count; i++) {
		basis->column_norm[i] = simplex_pan_column_norm(
			basis->standard, basis->column[i]);
		if (basis->column_norm[i] == 0.)
			return lp_simplex_EXIT_FAILURE;
	}
	for (j = 0; j < count; j++)
		for (i = 0; i <= j; i++) {
			double dot =
				pan_sparse_column_dot(basis->standard,
					basis->column[i], basis->column[j]) /
				(basis->column_norm[i] * basis->column_norm[j]);

			basis->factor[i + j * count] = dot;
			basis->factor[j + i * count] = dot;
		}
	dpotrf_(&lower, &count, basis->factor, &count, &info);
	basis->factorizations++;
	if (info < 0)
		return lp_simplex_EXIT_FAILURE;
	if (info > 0)
		return pan_basis_recover_qr(basis);
	for (i = 0; i < count; i++) {
		double diagonal =
			__lp_simplex_ABS__(basis->factor[i + i * count]);

		if (!isfinite(diagonal) || diagonal <= basis->rank_tolerance)
			return pan_basis_recover_qr(basis);
	}
	return lp_simplex_EXIT_SUCCESS;
}

int simplex_pan_basis_replace_factorized(
	struct simplex_PanBasis *basis, const int position, const int column)
{
	int leaving;

	if (position < 0 || position >= basis->count || column < 0 ||
		column >= basis->standard->columns ||
		basis->position[column] >= 0)
		return lp_simplex_EXIT_FAILURE;
	leaving = basis->column[position];
	basis->position[leaving] = -1;
	basis->column[position] = column;
	basis->position[column] = position;
	return simplex_pan_basis_factorize(basis);
}

static int pan_basis_solve_factor(
	const struct simplex_PanBasis *basis, double *right)
{
	char lower = 'L';
	int count = basis->count;
	int one = 1;
	int info = 0;

	extern void dpotrs_(
		char *, int *, int *, double *, int *, double *, int *, int *);
	if (count == 0)
		return lp_simplex_EXIT_SUCCESS;
	dpotrs_(&lower, &count, &one, basis->factor, &count, right, &count,
		&info);
	return info == 0 ? lp_simplex_EXIT_SUCCESS : lp_simplex_EXIT_FAILURE;
}

static void pan_basis_multiply(const struct simplex_PanBasis *basis,
	const double *coefficient, double *result)
{
	int i;
	lp_simplex_memset(
		result, 0, (size_t)basis->standard->rows * sizeof(double));
	for (i = 0; i < basis->count; i++)

		pan_column_axpy(basis->standard, basis->column[i],
			coefficient[i], result);
}

int simplex_pan_basis_least_squares(struct simplex_PanBasis *basis,
	const double *right, double *solution, double *residual_norm)
{
	int refinement, i;

	if (basis->dense_qr && basis->count)
		return simplex_pan_basis_least_squares_qr(
			basis, right, solution, residual_norm);
	if (basis->count == 0) {
		lp_simplex_memcpy(basis->row_work, right,
			(size_t)basis->standard->rows * sizeof(double));
		*residual_norm =
			pan_dense_norm(basis->row_work, basis->standard->rows);
		return lp_simplex_EXIT_SUCCESS;
	}
	for (i = 0; i < basis->count; i++)
		solution[i] = pan_column_dense_dot(basis->standard,
				      basis->column[i], right) /
			      basis->column_norm[i];
	if (pan_basis_solve_factor(basis, solution) == lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	for (i = 0; i < basis->count; i++)
		solution[i] /= basis->column_norm[i];
	for (refinement = 0; refinement < 2; refinement++) {
		pan_basis_multiply(basis, solution, basis->row_work);
		for (i = 0; i < basis->standard->rows; i++)
			basis->row_work[i] = right[i] - basis->row_work[i];
		for (i = 0; i < basis->count; i++)
			basis->basis_work[i] =
				pan_column_dense_dot(basis->standard,
					basis->column[i], basis->row_work) /
				basis->column_norm[i];
		if (pan_basis_solve_factor(basis, basis->basis_work) ==
			lp_simplex_EXIT_FAILURE)
			return lp_simplex_EXIT_FAILURE;
		for (i = 0; i < basis->count; i++)
			solution[i] +=
				basis->basis_work[i] / basis->column_norm[i];
		basis->refinements++;
	}
	pan_basis_multiply(basis, solution, basis->row_work);
	for (i = 0; i < basis->standard->rows; i++)
		basis->row_work[i] = right[i] - basis->row_work[i];
	*residual_norm = pan_dense_norm(basis->row_work, basis->standard->rows);
	return lp_simplex_EXIT_SUCCESS;
}

int simplex_pan_basis_least_squares_qr(struct simplex_PanBasis *basis,
	const double *right, double *solution, double *residual_norm)
{
	int rows = basis->standard->rows;
	int columns = basis->count;
	int one = 1;
	int leading_a = rows;
	int leading_b = rows > columns ? rows : columns;
	int rank = 0;
	int info = 0;
	int lwork = -1;
	int *pivot = NULL;
	double *matrix = NULL;
	double *right_work = NULL;
	double *work = NULL;
	double work_size = 0.;
	double rcond = basis->rank_tolerance;
	long double bytes;
	int i, j, k;

	extern void dgelsy_(int *, int *, int *, double *, int *, double *,
		int *, int *, double *, int *, double *, int *, int *);
	if (columns == 0)
		return simplex_pan_basis_least_squares(
			basis, right, solution, residual_norm);
	*residual_norm = HUGE_VAL;
	if (rows < 1 || columns < 0 || columns > rows ||
		rows > (INT_MAX - 1) / 4 ||
		(size_t)rows > (size_t)-1 / sizeof(double) / columns)
		return lp_simplex_EXIT_FAILURE;
	bytes = ((long double)rows * columns + leading_b) * sizeof(double) +
		(long double)columns * sizeof(int);
	if (bytes > SIMPLEX_PAN_MINIMUM_NORM_QR_BYTES)
		return lp_simplex_EXIT_FAILURE;
	matrix = (double *)lp_simplex_malloc(
		(size_t)rows * columns * sizeof(double));
	right_work =
		(double *)lp_simplex_malloc((size_t)leading_b * sizeof(double));
	pivot = (int *)lp_simplex_malloc((size_t)columns * sizeof(int));
	if (matrix == NULL || right_work == NULL || pivot == NULL)
		goto failure;
	lp_simplex_memset(matrix, 0, (size_t)rows * columns * sizeof(double));
	lp_simplex_memset(pivot, 0, (size_t)columns * sizeof(int));
	lp_simplex_memset(right_work, 0, (size_t)leading_b * sizeof(double));
	lp_simplex_memcpy(right_work, right, (size_t)rows * sizeof(double));
	for (j = 0; j < columns; j++) {
		int column = basis->column[j];
		double norm = simplex_pan_column_norm(basis->standard, column);

		if (norm == 0.)
			goto failure;
		for (k = basis->standard->column_start[column];
			k < basis->standard->column_start[column + 1]; k++)
			matrix[basis->standard->row_index[k] + j * rows] =
				basis->standard->value[k] / norm;
	}
	dgelsy_(&rows, &columns, &one, matrix, &leading_a, right_work,
		&leading_b, pivot, &rcond, &rank, &work_size, &lwork, &info);
	if (info != 0 || !isfinite(work_size) || work_size < 1. ||
		work_size > INT_MAX)
		goto failure;
	lwork = (int)ceil(work_size);
	if (bytes + (long double)lwork * sizeof(double) >
		SIMPLEX_PAN_MINIMUM_NORM_QR_BYTES)
		goto failure;
	work = (double *)lp_simplex_malloc((size_t)lwork * sizeof(double));
	if (work == NULL)
		goto failure;
	/* The workspace query may overwrite A and B on some LAPACK builds. */
	lp_simplex_memset(matrix, 0, (size_t)rows * columns * sizeof(double));
	lp_simplex_memset(pivot, 0, (size_t)columns * sizeof(int));
	lp_simplex_memset(right_work, 0, (size_t)leading_b * sizeof(double));
	lp_simplex_memcpy(right_work, right, (size_t)rows * sizeof(double));
	for (j = 0; j < columns; j++) {
		int column = basis->column[j];
		double norm = simplex_pan_column_norm(basis->standard, column);

		for (k = basis->standard->column_start[column];
			k < basis->standard->column_start[column + 1]; k++)
			matrix[basis->standard->row_index[k] + j * rows] =
				basis->standard->value[k] / norm;
	}
	dgelsy_(&rows, &columns, &one, matrix, &leading_a, right_work,
		&leading_b, pivot, &rcond, &rank, work, &lwork, &info);
	if (info != 0 || rank != columns)
		goto failure;
	for (i = 0; i < columns; i++)
		solution[i] =
			right_work[i] / simplex_pan_column_norm(basis->standard,
						basis->column[i]);
	pan_basis_multiply(basis, solution, basis->row_work);
	for (i = 0; i < rows; i++)
		basis->row_work[i] = right[i] - basis->row_work[i];
	*residual_norm = pan_dense_norm(basis->row_work, rows);
	lp_simplex_free(matrix);
	lp_simplex_free(right_work);
	lp_simplex_free(pivot);
	lp_simplex_free(work);
	return lp_simplex_EXIT_SUCCESS;
failure:
	lp_simplex_free(matrix);
	lp_simplex_free(right_work);
	lp_simplex_free(pivot);
	lp_simplex_free(work);
	return lp_simplex_EXIT_FAILURE;
}

int simplex_pan_basis_minimum_norm(struct simplex_PanBasis *basis,
	const double *right, double *solution, double *residual_norm)
{
	int refinement, i;

	if (basis->dense_qr && basis->count)
		return simplex_pan_basis_minimum_norm_qr(
			basis, right, solution, residual_norm);
	if (basis->count == 0) {
		lp_simplex_memset(solution, 0,
			(size_t)basis->standard->rows * sizeof(double));
		*residual_norm = 0.;
		return lp_simplex_EXIT_SUCCESS;
	}
	for (i = 0; i < basis->count; i++)
		basis->basis_work[i] = right[i] / basis->column_norm[i];
	if (pan_basis_solve_factor(basis, basis->basis_work) ==
		lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	for (i = 0; i < basis->count; i++)
		basis->basis_work[i] /= basis->column_norm[i];
	pan_basis_multiply(basis, basis->basis_work, solution);
	for (refinement = 0; refinement < 2; refinement++) {
		for (i = 0; i < basis->count; i++)
			basis->basis_work[i] =
				right[i] - pan_column_dense_dot(basis->standard,
						   basis->column[i], solution);
		for (i = 0; i < basis->count; i++)
			basis->basis_work[i] /= basis->column_norm[i];
		if (pan_basis_solve_factor(basis, basis->basis_work) ==
			lp_simplex_EXIT_FAILURE)
			return lp_simplex_EXIT_FAILURE;
		for (i = 0; i < basis->count; i++)
			pan_column_axpy(basis->standard, basis->column[i],
				basis->basis_work[i] / basis->column_norm[i],
				solution);
		basis->refinements++;
	}
	for (i = 0; i < basis->count; i++)
		basis->basis_work[i] =
			right[i] - pan_column_dense_dot(basis->standard,
					   basis->column[i], solution);
	*residual_norm = pan_dense_norm(basis->basis_work, basis->count);
	return lp_simplex_EXIT_SUCCESS;
}

int simplex_pan_basis_column_projection(struct simplex_PanBasis *basis,
	const int column, double *coefficient, double *residual_norm)
{
	int refinement, i, k;

	if (basis->dense_qr && basis->count) {
		int state;

		lp_simplex_memset(basis->right_work, 0,
			(size_t)basis->standard->rows * sizeof(double));
		pan_column_axpy(basis->standard, column, 1., basis->right_work);
		*residual_norm = HUGE_VAL;
		state = simplex_pan_basis_least_squares_qr(
			basis, basis->right_work, coefficient, residual_norm);
		basis->projection_solve_error = *residual_norm;
		return state;
	}
	if (basis->count == 0) {
		*residual_norm =
			simplex_pan_column_norm(basis->standard, column);
		return lp_simplex_EXIT_SUCCESS;
	}
	for (i = 0; i < basis->count; i++)
		coefficient[i] = pan_sparse_column_dot(basis->standard,
					 basis->column[i], column) /
				 basis->column_norm[i];
	if (pan_basis_solve_factor(basis, coefficient) ==
		lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	for (i = 0; i < basis->count; i++)
		coefficient[i] /= basis->column_norm[i];
	for (refinement = 0; refinement < 2; refinement++) {
		pan_basis_multiply(basis, coefficient, basis->row_work);
		for (i = 0; i < basis->standard->rows; i++)
			basis->row_work[i] = -basis->row_work[i];
		for (k = basis->standard->column_start[column];
			k < basis->standard->column_start[column + 1]; k++)
			basis->row_work[basis->standard->row_index[k]] +=
				basis->standard->value[k];
		for (i = 0; i < basis->count; i++)
			basis->basis_work[i] =
				pan_column_dense_dot(basis->standard,
					basis->column[i], basis->row_work) /
				basis->column_norm[i];
		if (pan_basis_solve_factor(basis, basis->basis_work) ==
			lp_simplex_EXIT_FAILURE)
			return lp_simplex_EXIT_FAILURE;
		for (i = 0; i < basis->count; i++)
			coefficient[i] +=
				basis->basis_work[i] / basis->column_norm[i];
		basis->refinements++;
	}
	pan_basis_multiply(basis, coefficient, basis->row_work);
	for (i = 0; i < basis->standard->rows; i++)
		basis->row_work[i] = -basis->row_work[i];
	for (k = basis->standard->column_start[column];
		k < basis->standard->column_start[column + 1]; k++)
		basis->row_work[basis->standard->row_index[k]] +=
			basis->standard->value[k];
	*residual_norm = pan_dense_norm(basis->row_work, basis->standard->rows);
	basis->projection_solve_error = *residual_norm;
	return lp_simplex_EXIT_SUCCESS;
}

#endif /* !LP_SIMPLEX_HAVE_SPQR */
