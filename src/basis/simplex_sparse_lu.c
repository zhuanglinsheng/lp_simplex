/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_sparse_lu.c
 * @brief Sparse row-oriented LU with partial pivoting for simplex basis
 * matrices.
 */

#include "policy/factorization.h"
#include "simplex_sparse_lu.h"
#include "utils.h"
#include "simplex_storage.h"
#include <limits.h>

#include <lp_simplex/status.h>

/* Order columns by increasing initial degree, breaking ties by the original
 * position.  The old insertion sort made every reinversion quadratic when a
 * changing basis happened to present columns in reverse degree order. */
static int sparse_lu_degree_greater(
	const int left, const int right, const int *degree)
{
	return degree[left] > degree[right] ||
	       (degree[left] == degree[right] && left > right);
}

static void sparse_lu_sift_degree_heap(
	int *order, const int count, int root, const int *degree)
{
	for (;;) {
		int child;
		int largest = root;
		int value;
		if (root >= count / 2)
			return;
		child = 2 * root + 1;
		if (child < count && sparse_lu_degree_greater(order[child],
					     order[largest], degree))
			largest = child;
		if (child + 1 < count &&
			sparse_lu_degree_greater(
				order[child + 1], order[largest], degree))
			largest = child + 1;
		if (largest == root)
			return;
		value = order[root];
		order[root] = order[largest];
		order[largest] = value;
		root = largest;
	}
}

static void sparse_lu_sort_columns_by_degree(
	int *order, const int count, const int *degree)
{
	int adjacent_inversions = 0;
	int i;
	/* Basis updates often preserve most of the previous structural order.  In
	 * that common case insertion sort is linear in practice and moves less
	 * memory than a heap.  Use the bounded-complexity path once disorder is
	 * visible across a material number of adjacent pairs. */
	for (i = 1; i < count; i++)
		if (sparse_lu_degree_greater(order[i - 1], order[i], degree))
			adjacent_inversions++;
	if (adjacent_inversions <= __lp_simplex_MAX__(8, count / 64)) {
		for (i = 1; i < count; i++) {
			int original = order[i];
			int position = i;

			while (position > 0 &&
				sparse_lu_degree_greater(order[position - 1],
					original, degree)) {
				order[position] = order[position - 1];
				position--;
			}
			order[position] = original;
		}
		return;
	}
	for (i = count / 2; i > 0; i--)
		sparse_lu_sift_degree_heap(order, count, i - 1, degree);
	for (i = count - 1; i > 0; i--) {
		int value = order[0];

		order[0] = order[i];
		order[i] = value;
		sparse_lu_sift_degree_heap(order, i, 0, degree);
	}
}

static void sparse_row_destroy(struct simplex_SparseRow *row)
{
	lp_simplex_free(row->column);
	lp_simplex_free(row->value);
	row->column = NULL;
	row->value = NULL;
	row->count = 0;
	row->capacity = 0;
}

static void sparse_column_rows_destroy(struct simplex_SparseColumnRows *column)
{
	lp_simplex_free(column->row);
	column->row = NULL;
	column->count = 0;
	column->capacity = 0;
}

static int sparse_column_rows_append(struct simplex_SparseLu *factor,
	struct simplex_SparseColumnRows *column, const int row)
{
	int *grown_rows;
	int capacity;
	size_t bytes;

	if (column->count == INT_MAX)
		goto resource_failure;
	if (column->count < column->capacity) {
		column->row[column->count++] = row;
		return lp_simplex_EXIT_SUCCESS;
	}
	capacity = column->capacity > INT_MAX / 2 ? column->count + 1
		   : column->capacity > 0	  ? 2 * column->capacity
						  : 8;
	if (!simplex_storage_bytes(
		    (size_t)capacity, sizeof(int), (size_t)-1, &bytes))
		goto resource_failure;
	grown_rows = (int *)lp_simplex_malloc(bytes);
	if (grown_rows == NULL)
		goto resource_failure;
	if (column->count > 0)
		lp_simplex_memcpy(grown_rows, column->row,
			(size_t)column->count * sizeof(int));
	lp_simplex_free(column->row);
	column->row = grown_rows;
	column->capacity = capacity;
	column->row[column->count++] = row;
	return lp_simplex_EXIT_SUCCESS;
resource_failure:
	factor->storage_failed = 1;
	return lp_simplex_EXIT_FAILURE;
}

static int sparse_row_reserve(struct simplex_SparseLu *factor,
	struct simplex_SparseRow *row, const int capacity)
{
	int *column;
	double *value;
	int grown;
	size_t column_bytes, value_bytes;

	if (capacity < 0)
		goto resource_failure;
	if (capacity <= row->capacity)
		return lp_simplex_EXIT_SUCCESS;
	grown = row->capacity > 0 ? row->capacity : 4;
	while (grown < capacity) {
		int increment = grown < 1024 ? grown : grown / 2;
		if (grown > INT_MAX - increment) {
			grown = capacity;
			break;
		}
		grown += increment;
	}
	if (!simplex_storage_bytes(
		    (size_t)grown, sizeof(int), (size_t)-1, &column_bytes) ||
		!simplex_storage_bytes((size_t)grown, sizeof(double),
			(size_t)-1, &value_bytes))
		goto resource_failure;
	column = (int *)lp_simplex_malloc(column_bytes);
	value = (double *)lp_simplex_malloc(value_bytes);
	if (column == NULL || value == NULL) {
		lp_simplex_free(column);
		lp_simplex_free(value);
		goto resource_failure;
	}
	if (row->count > 0) {
		lp_simplex_memcpy(
			column, row->column, (size_t)row->count * sizeof(int));
		lp_simplex_memcpy(
			value, row->value, (size_t)row->count * sizeof(double));
	}
	lp_simplex_free(row->column);
	lp_simplex_free(row->value);
	row->column = column;
	row->value = value;
	row->capacity = grown;
	return lp_simplex_EXIT_SUCCESS;
resource_failure:
	factor->storage_failed = 1;
	return lp_simplex_EXIT_FAILURE;
}

static int sparse_row_find_counted(struct simplex_SparseLu *factor,
	const struct simplex_SparseRow *row, const int column)
{
	int left = 0, right = row->count - 1;

	while (left <= right) {
		int middle = left + (right - left) / 2;

		factor->factor_work++;
		if (row->column[middle] < column)
			left = middle + 1;
		else if (row->column[middle] > column)
			right = middle - 1;
		else
			return middle;
	}
	return -1;
}

struct sparse_LuStorage {
	size_t integer, numeric, rows, columns, etas, order, pivots, cache;
};

static int sparse_lu_storage(
	int dimension, size_t limit, struct sparse_LuStorage *layout)
{
	struct sparse_LuStorage next = {0};
	*layout = next;
	if (dimension < 0 || limit < sizeof(int) ||
		!simplex_storage_bytes((size_t)dimension, 6 * sizeof(int),
			limit - sizeof(int), &next.integer) ||
		!simplex_storage_bytes((size_t)dimension, 5 * sizeof(double),
			limit, &next.numeric) ||
		!simplex_storage_bytes((size_t)dimension,
			sizeof(struct simplex_SparseRow), limit, &next.rows) ||
		!simplex_storage_bytes((size_t)dimension,
			sizeof(struct simplex_SparseColumnRows), limit,
			&next.columns) ||
		!simplex_storage_bytes(SPARSE_LU_FT_UPDATE_LIMIT,
			sizeof(struct simplex_SparseRow), limit, &next.etas) ||
		!simplex_storage_bytes((size_t)dimension, 2 * sizeof(int),
			limit, &next.order) ||
		!simplex_storage_bytes(SPARSE_LU_FT_UPDATE_LIMIT, sizeof(int),
			limit, &next.pivots) ||
		!simplex_storage_bytes((size_t)dimension, 3 * sizeof(double),
			limit, &next.cache))
		return 0;
	next.integer += sizeof(int);
	*layout = next;
	return 1;
}

static void *sparse_lu_zero_storage(size_t bytes)
{
	void *storage = lp_simplex_malloc(bytes);
	if (storage != NULL)
		lp_simplex_memset(storage, 0, bytes);
	return storage;
}

int simplex_sparse_lu_create(
	struct simplex_SparseLu *factor, const int dimension)
{
	struct sparse_LuStorage layout;
	int *integer;
	double *numeric;
	lp_simplex_memset(factor, 0, sizeof(*factor));
	if (!sparse_lu_storage(dimension, (size_t)-1, &layout)) {
		factor->storage_failed = 1;
		return lp_simplex_EXIT_FAILURE;
	}
	factor->dimension = dimension;
	factor->integer_storage = (int *)lp_simplex_malloc(layout.integer);
	factor->numeric_storage = (double *)lp_simplex_malloc(layout.numeric);
	integer = factor->integer_storage;
	numeric = factor->numeric_storage;
	if (integer != NULL) {
		factor->permutation = integer;
		factor->column_permutation = integer + dimension;
		factor->work_column = integer + (size_t)2 * dimension;
		factor->pivot_row = integer + (size_t)3 * dimension;
		factor->diagonal_position = integer + (size_t)4 * dimension;
		factor->packed_start = integer + (size_t)5 * dimension;
	}
	if (numeric != NULL) {
		factor->column_scale = numeric;
		factor->row_scale = numeric + dimension;
		factor->diagonal_value = numeric + (size_t)2 * dimension;
		factor->work_value = numeric + (size_t)3 * dimension;
		factor->solve_work = numeric + (size_t)4 * dimension;
	}
	/* Every allocated descriptor is destructible before another allocation
	 * can fail. Numerical row contents remain lazily allocated. */
	factor->row = sparse_lu_zero_storage(layout.rows);
	factor->column_rows = sparse_lu_zero_storage(layout.columns);
	factor->ft_column = sparse_lu_zero_storage(layout.rows);
	factor->ft_row_eta = sparse_lu_zero_storage(layout.etas);
	factor->ft_order = (int *)lp_simplex_malloc(layout.order);
	factor->ft_pivot_position =
		factor->ft_order != NULL ? factor->ft_order + dimension : NULL;
	factor->ft_row_pivot = (int *)lp_simplex_malloc(layout.pivots);
	factor->ft_update_capacity = SPARSE_LU_FT_UPDATE_LIMIT;
	if (factor->integer_storage == NULL ||
		factor->numeric_storage == NULL || factor->row == NULL ||
		factor->column_rows == NULL || factor->ft_column == NULL ||
		factor->ft_row_eta == NULL || factor->ft_order == NULL ||
		factor->ft_row_pivot == NULL) {
		simplex_sparse_lu_destroy(factor);
		factor->storage_failed = 1;
		return lp_simplex_EXIT_FAILURE;
	}
	return lp_simplex_EXIT_SUCCESS;
}

void simplex_sparse_lu_destroy(struct simplex_SparseLu *factor)
{
	int i;

	if (factor == NULL)
		return;
	if (factor->row != NULL) {
		for (i = 0; i < factor->dimension; i++)
			sparse_row_destroy(factor->row + i);
	}
	if (factor->column_rows != NULL) {
		for (i = 0; i < factor->dimension; i++)
			sparse_column_rows_destroy(factor->column_rows + i);
	}
	if (factor->ft_column != NULL)
		for (i = 0; i < factor->dimension; i++)
			sparse_row_destroy(factor->ft_column + i);
	if (factor->ft_row_eta != NULL)
		for (i = 0; i < factor->ft_update_capacity; i++)
			sparse_row_destroy(factor->ft_row_eta + i);
	lp_simplex_free(factor->row);
	lp_simplex_free(factor->column_rows);
	lp_simplex_free(factor->ft_column);
	lp_simplex_free(factor->ft_row_eta);
	lp_simplex_free(factor->ft_order);
	lp_simplex_free(factor->ft_row_pivot);
	lp_simplex_free(factor->ft_spike_cache);
	lp_simplex_free(factor->integer_storage);
	lp_simplex_free(factor->numeric_storage);
	lp_simplex_free(factor->packed_column);
	lp_simplex_free(factor->packed_value);
	lp_simplex_memset(factor, 0, sizeof(*factor));
}

static int sparse_lu_assemble(struct simplex_SparseLu *factor,
	const struct simplex_CscMatrix *matrix, const int structural_columns,
	const int *basis)
{
	int i, j, k, *count;
	int n = factor->dimension;

	count = factor->work_column;
	lp_simplex_memset(count, 0, (size_t)n * sizeof(int));
	for (i = 0; i < n; i++) {
		factor->permutation[i] = i;
		factor->row[i].count = 0;
	}
	for (j = 0; j < n; j++) {
		int variable = basis[j];
		double largest = 0.;

		factor->column_permutation[j] = j;
		count[j] = variable >= structural_columns
				   ? 1
				   : matrix->column_start[variable + 1] -
					     matrix->column_start[variable];
		if (variable >= structural_columns)
			largest = 1.;
		else
			for (k = matrix->column_start[variable];
				k < matrix->column_start[variable + 1]; k++)
				largest = __lp_simplex_MAX__(largest,
					__lp_simplex_ABS__(matrix->value[k]));
		factor->column_scale[j] = largest > 0. ? 1. / largest : 1.;
	}
	lp_simplex_memset(factor->row_scale, 0, (size_t)n * sizeof(double));
	for (j = 0; j < n; j++) {
		int variable = basis[j];
		double column_scale = factor->column_scale[j];

		if (variable >= structural_columns) {
			int row_index = variable - structural_columns;

			factor->row_scale[row_index] = __lp_simplex_MAX__(
				factor->row_scale[row_index], column_scale);
		} else {
			for (k = matrix->column_start[variable];
				k < matrix->column_start[variable + 1]; k++) {
				int row_index = matrix->row_index[k];

				factor->row_scale[row_index] =
					__lp_simplex_MAX__(
						factor->row_scale[row_index],
						__lp_simplex_ABS__(
							matrix->value[k] *
							column_scale));
			}
		}
	}
	for (i = 0; i < n; i++)
		factor->row_scale[i] = factor->row_scale[i] > 0.
					       ? 1. / factor->row_scale[i]
					       : 1.;
	sparse_lu_sort_columns_by_degree(factor->column_permutation, n, count);
	lp_simplex_memset(count, 0, (size_t)n * sizeof(int));
	for (j = 0; j < n; j++) {
		int variable = basis[factor->column_permutation[j]];

		if (variable >= structural_columns)
			count[variable - structural_columns]++;
		else
			for (k = matrix->column_start[variable];
				k < matrix->column_start[variable + 1]; k++)
				count[matrix->row_index[k]]++;
	}
	for (i = 0; i < n; i++) {
		if (sparse_row_reserve(factor, factor->row + i, count[i]) ==
			lp_simplex_EXIT_FAILURE)
			return lp_simplex_EXIT_FAILURE;
	}
	for (j = 0; j < n; j++) {
		int variable = basis[factor->column_permutation[j]];

		if (variable >= structural_columns) {
			int row_index = variable - structural_columns;
			struct simplex_SparseRow *row = factor->row + row_index;

			row->column[row->count] = j;
			row->value[row->count++] =
				-factor->column_scale[factor
						->column_permutation[j]] *
				factor->row_scale[row_index];
		} else {
			for (k = matrix->column_start[variable];
				k < matrix->column_start[variable + 1]; k++) {
				int row_index = matrix->row_index[k];
				struct simplex_SparseRow *row =
					factor->row + row_index;

				row->column[row->count] = j;
				row->value[row->count++] =
					matrix->value[k] *
					factor->column_scale[factor
							->column_permutation
								[j]] *
					factor->row_scale[row_index];
			}
		}
	}
	return lp_simplex_EXIT_SUCCESS;
}

static int sparse_lu_assemble_submatrix(struct simplex_SparseLu *factor,
	const struct simplex_CscMatrix *matrix, const int *columns,
	const int *row_to_core)
{
	int i, j, k, *count;
	int n = factor->dimension;

	count = factor->work_column;
	lp_simplex_memset(count, 0, (size_t)n * sizeof(int));
	for (i = 0; i < n; i++) {
		factor->permutation[i] = i;
		factor->row[i].count = 0;
	}
	for (j = 0; j < n; j++) {
		double largest = 0.;

		factor->column_permutation[j] = j;
		for (k = matrix->column_start[columns[j]];
			k < matrix->column_start[columns[j] + 1]; k++) {
			int row = row_to_core[matrix->row_index[k]];

			if (row < 0)
				continue;
			count[j]++;
			largest = __lp_simplex_MAX__(
				largest, __lp_simplex_ABS__(matrix->value[k]));
		}
		factor->column_scale[j] = largest > 0. ? 1. / largest : 1.;
	}
	lp_simplex_memset(factor->row_scale, 0, (size_t)n * sizeof(double));
	for (j = 0; j < n; j++) {
		double column_scale = factor->column_scale[j];

		for (k = matrix->column_start[columns[j]];
			k < matrix->column_start[columns[j] + 1]; k++) {
			int row = row_to_core[matrix->row_index[k]];

			if (row >= 0)
				factor->row_scale[row] = __lp_simplex_MAX__(
					factor->row_scale[row],
					__lp_simplex_ABS__(matrix->value[k] *
							   column_scale));
		}
	}
	for (i = 0; i < n; i++)
		factor->row_scale[i] = factor->row_scale[i] > 0.
					       ? 1. / factor->row_scale[i]
					       : 1.;
	sparse_lu_sort_columns_by_degree(factor->column_permutation, n, count);
	lp_simplex_memset(count, 0, (size_t)n * sizeof(int));
	for (j = 0; j < n; j++) {
		int column = columns[factor->column_permutation[j]];

		for (k = matrix->column_start[column];
			k < matrix->column_start[column + 1]; k++) {
			int row = row_to_core[matrix->row_index[k]];

			if (row >= 0)
				count[row]++;
		}
	}
	for (i = 0; i < n; i++)
		if (sparse_row_reserve(factor, factor->row + i, count[i]) ==
			lp_simplex_EXIT_FAILURE)
			return lp_simplex_EXIT_FAILURE;
	for (j = 0; j < n; j++) {
		int original = factor->column_permutation[j];
		int column = columns[original];

		for (k = matrix->column_start[column];
			k < matrix->column_start[column + 1]; k++) {
			int row_index = row_to_core[matrix->row_index[k]];

			if (row_index >= 0) {
				struct simplex_SparseRow *row =
					factor->row + row_index;

				row->column[row->count] = j;
				row->value[row->count++] =
					matrix->value[k] *
					factor->column_scale[original] *
					factor->row_scale[row_index];
			}
		}
	}
	return lp_simplex_EXIT_SUCCESS;
}

static int sparse_lu_eliminate_row(struct simplex_SparseLu *factor,
	const int target_index, const int pivot_index, const int column,
	const double multiplier)
{
	int a = 0, b = 0, count = 0;
	struct simplex_SparseRow *target = factor->row + target_index;
	const struct simplex_SparseRow *pivot = factor->row + pivot_index;

	while (a < target->count || b < pivot->count) {
		int target_column = a < target->count ? target->column[a]
						      : factor->dimension;
		int pivot_column =
			b < pivot->count ? pivot->column[b] : factor->dimension;
		int result_column;
		double result_value;

		factor->factor_work++;
		if (pivot_column <= column) {
			b++;
			continue;
		}
		if (target_column < pivot_column) {
			result_column = target_column;
			result_value = target->value[a++];
		} else if (pivot_column < target_column) {
			result_column = pivot_column;
			result_value = -multiplier * pivot->value[b++];
		} else {
			result_column = target_column;
			result_value = target->value[a++] -
				       multiplier * pivot->value[b++];
		}
		if (result_column == column)
			result_value = multiplier;
		if (result_value != 0.) {
			if (target_column > pivot_column &&
				factor->column_rows != NULL &&
				sparse_column_rows_append(factor,
					factor->column_rows + result_column,
					target_index) ==
					lp_simplex_EXIT_FAILURE)
				return lp_simplex_EXIT_FAILURE;
			factor->work_column[count] = result_column;
			factor->work_value[count] = result_value;
			count++;
		}
	}
	if (sparse_row_reserve(factor, target, count) ==
		lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	lp_simplex_memcpy(target->column, factor->work_column,
		(size_t)count * sizeof(int));
	lp_simplex_memcpy(target->value, factor->work_value,
		(size_t)count * sizeof(double));
	target->count = count;
	return lp_simplex_EXIT_SUCCESS;
}

static int sparse_lu_build_column_rows(struct simplex_SparseLu *factor)
{
	int i, k;

	for (i = 0; i < factor->dimension; i++)
		factor->column_rows[i].count = 0;
	for (i = 0; i < factor->dimension; i++) {
		const struct simplex_SparseRow *row = factor->row + i;

		for (k = 0; k < row->count; k++) {
			factor->factor_work++;
			if (sparse_column_rows_append(factor,
				    factor->column_rows + row->column[k],
				    i) == lp_simplex_EXIT_FAILURE)
				return lp_simplex_EXIT_FAILURE;
		}
	}
	return lp_simplex_EXIT_SUCCESS;
}

static int sparse_lu_pack_rows(struct simplex_SparseLu *factor)
{
	int i, next = 0, nonzeros = 0;
	int *column;
	double *value;
	size_t column_bytes, value_bytes;

	for (i = 0; i < factor->dimension; i++) {
		if (factor->row[i].count > INT_MAX - nonzeros)
			goto resource_failure;
		nonzeros += factor->row[i].count;
	}
	if (!simplex_storage_bytes(
		    (size_t)nonzeros, sizeof(int), (size_t)-1, &column_bytes) ||
		!simplex_storage_bytes((size_t)nonzeros, sizeof(double),
			(size_t)-1, &value_bytes))
		goto resource_failure;
	column = nonzeros > 0 ? (int *)lp_simplex_malloc(column_bytes) : NULL;
	value = nonzeros > 0 ? (double *)lp_simplex_malloc(value_bytes) : NULL;
	if (nonzeros > 0 && (column == NULL || value == NULL)) {
		lp_simplex_free(column);
		lp_simplex_free(value);
		goto resource_failure;
	}
	for (i = 0; i < factor->dimension; i++) {
		const struct simplex_SparseRow *row = factor->row + i;

		factor->packed_start[i] = next;
		if (row->count > 0) {
			lp_simplex_memcpy(column + next, row->column,
				(size_t)row->count * sizeof(int));
			lp_simplex_memcpy(value + next, row->value,
				(size_t)row->count * sizeof(double));
			next += row->count;
		}
	}
	factor->packed_start[factor->dimension] = next;
	lp_simplex_free(factor->packed_column);
	lp_simplex_free(factor->packed_value);
	factor->packed_column = column;
	factor->packed_value = value;
	factor->packed_nonzeros = nonzeros;
	return lp_simplex_EXIT_SUCCESS;
resource_failure:
	factor->storage_failed = 1;
	return lp_simplex_EXIT_FAILURE;
}

/* Build the initial U eta file.  For an upper triangular U, its raw column
 * etas have the factor product
 *
 *     U = E_{n-1} ... E_1 E_0.
 *
 * Forrest--Tomlin updates preserve this product form while changing the eta
 * order; no triangular pattern is assumed after the first update. */
static int sparse_lu_ft_initialize(struct simplex_SparseLu *factor)
{
	int i, k;
	int n = factor->dimension;
	long nonzeros = 0;
	/* KLU-backed bases never use these caches.  Internal LU allocates them
	 * once, retaining a verified spike fast path without per-pivot allocation. */
	if (factor->ft_spike_cache == NULL) {
		size_t bytes;
		if (!simplex_storage_bytes((size_t)n, 3 * sizeof(double),
			    (size_t)-1, &bytes)) {
			factor->storage_failed = 1;
			return lp_simplex_EXIT_FAILURE;
		}
		factor->ft_spike_cache = (double *)lp_simplex_malloc(bytes);
		if (factor->ft_spike_cache == NULL) {
			factor->storage_failed = 1;
			return lp_simplex_EXIT_FAILURE;
		}
		factor->ft_btran_cache = factor->ft_spike_cache + n;
		factor->ft_direction_cache =
			factor->ft_spike_cache + (size_t)2 * n;
	}
	for (i = 0; i < n; i++) {
		factor->ft_column[i].count = 0;
		/* U = E_{n-1} ... E_1 E_0.  The eta file stores that factor
		 * order; applying the inverse to a vector follows the same list. */
		factor->ft_order[i] = n - 1 - i;
		factor->ft_pivot_position[n - 1 - i] = i;
	}
	for (i = 0; i < factor->ft_update_capacity; i++)
		factor->ft_row_eta[i].count = 0;
	for (i = 0; i < n; i++) {
		int diagonal =
			factor->packed_start[i] + factor->diagonal_position[i];
		int end = factor->packed_start[i + 1];

		for (k = diagonal; k < end; k++) {
			int column = factor->packed_column[k];
			struct simplex_SparseRow *eta =
				factor->ft_column + column;

			if (sparse_row_reserve(factor, eta, eta->count + 1) ==
				lp_simplex_EXIT_FAILURE)
				return lp_simplex_EXIT_FAILURE;
			eta->column[eta->count] = i;
			eta->value[eta->count++] = factor->packed_value[k];
			nonzeros++;
		}
	}
	factor->ft_update_count = 0;
	factor->ft_active = 1;
	factor->ft_spike_valid = 0;
	factor->ft_btran_pivot = -1;
	factor->ft_initial_nonzeros = nonzeros;
	factor->ft_nonzeros = nonzeros;
	factor->ft_row_nonzeros = 0;
	return lp_simplex_EXIT_SUCCESS;
}

static double sparse_lu_ft_pivot(
	const struct simplex_SparseLu *factor, const int pivot)
{
	const struct simplex_SparseRow *eta = factor->ft_column + pivot;
	int k;

	for (k = 0; k < eta->count; k++)
		if (eta->column[k] == pivot)
			return eta->value[k];
	return 0.;
}

static int sparse_lu_ft_u_inverse(const struct simplex_SparseLu *factor,
	double *work, const int transpose)
{
	int q;
	int n = factor->dimension;

	if (!transpose) {
		for (q = 0; q < n; q++) {
			int k;
			int p = factor->ft_order[q];
			const struct simplex_SparseRow *eta =
				factor->ft_column + p;
			double pivot = sparse_lu_ft_pivot(factor, p);
			double multiplier;

			if (__lp_simplex_ABS__(pivot) <=
				SPARSE_LU_PIVOT_TOLERANCE)
				return lp_simplex_EXIT_FAILURE;
			if (work[p] == 0.)
				continue;
			multiplier = work[p] / pivot;
			for (k = 0; k < eta->count; k++)
				if (eta->column[k] != p)
					work[eta->column[k]] -=
						multiplier * eta->value[k];
			work[p] = multiplier;
		}
	} else {
		for (q = n - 1; q >= 0; q--) {
			int k;
			int p = factor->ft_order[q];
			const struct simplex_SparseRow *eta =
				factor->ft_column + p;
			double pivot = sparse_lu_ft_pivot(factor, p);
			double value = work[p];

			if (__lp_simplex_ABS__(pivot) <=
				SPARSE_LU_PIVOT_TOLERANCE)
				return lp_simplex_EXIT_FAILURE;
			for (k = 0; k < eta->count; k++)
				if (eta->column[k] != p)
					value -= eta->value[k] *
						 work[eta->column[k]];
			work[p] = value / pivot;
		}
	}
	return lp_simplex_EXIT_SUCCESS;
}

static void sparse_lu_ft_u_product(
	const struct simplex_SparseLu *factor, double *work)
{
	int q;

	for (q = factor->dimension - 1; q >= 0; q--) {
		int k;
		int p = factor->ft_order[q];
		const struct simplex_SparseRow *eta = factor->ft_column + p;
		double value = work[p];

		for (k = 0; k < eta->count; k++)
			if (eta->column[k] != p)
				work[eta->column[k]] += eta->value[k] * value;
		work[p] = sparse_lu_ft_pivot(factor, p) * value;
	}
}

static void sparse_lu_ft_apply_rows(const struct simplex_SparseLu *factor,
	double *work, const int transpose)
{
	int q;

	if (!transpose)
		for (q = 0; q < factor->ft_update_count; q++) {
			int k;
			const struct simplex_SparseRow *row =
				factor->ft_row_eta + q;
			int p = factor->ft_row_pivot[q];
			double dot = 0.;

			for (k = 0; k < row->count; k++)
				dot += row->value[k] * work[row->column[k]];
			work[p] -= dot;
		}
	else
		for (q = factor->ft_update_count - 1; q >= 0; q--) {
			int k;
			const struct simplex_SparseRow *row =
				factor->ft_row_eta + q;
			int p = factor->ft_row_pivot[q];
			double value = work[p];

			for (k = 0; k < row->count; k++)
				work[row->column[k]] -= row->value[k] * value;
		}
}

static int sparse_lu_append_row_columns(
	struct simplex_SparseLu *factor, const int row_index)
{
	int k;
	const struct simplex_SparseRow *row = factor->row + row_index;

	for (k = 0; k < row->count; k++)
		if (sparse_column_rows_append(factor,
			    factor->column_rows + row->column[k],
			    row_index) == lp_simplex_EXIT_FAILURE)
			return lp_simplex_EXIT_FAILURE;
	return lp_simplex_EXIT_SUCCESS;
}

static void sparse_row_swap_columns(struct simplex_SparseLu *factor,
	struct simplex_SparseRow *row, const int left, const int right)
{
	int left_position = sparse_row_find_counted(factor, row, left);
	int right_position = sparse_row_find_counted(factor, row, right);

	if (left_position >= 0 && right_position >= 0) {
		double value = row->value[left_position];

		row->value[left_position] = row->value[right_position];
		row->value[right_position] = value;
	} else if (left_position >= 0) {
		double value = row->value[left_position];

		while (left_position + 1 < row->count &&
			row->column[left_position + 1] < right) {
			factor->factor_work++;
			row->column[left_position] =
				row->column[left_position + 1];
			row->value[left_position] =
				row->value[left_position + 1];
			left_position++;
		}
		row->column[left_position] = right;
		row->value[left_position] = value;
	} else if (right_position >= 0) {
		double value = row->value[right_position];

		while (right_position > 0 &&
			row->column[right_position - 1] > left) {
			factor->factor_work++;
			row->column[right_position] =
				row->column[right_position - 1];
			row->value[right_position] =
				row->value[right_position - 1];
			right_position--;
		}
		row->column[right_position] = left;
		row->value[right_position] = value;
	}
}

static void sparse_lu_swap_columns(
	struct simplex_SparseLu *factor, const int left, const int right)
{
	int i, original;

	if (left == right)
		return;
	for (i = 0; i < factor->dimension; i++)
		sparse_row_swap_columns(factor, factor->row + i, left, right);
	original = factor->column_permutation[left];
	factor->column_permutation[left] = factor->column_permutation[right];
	factor->column_permutation[right] = original;
}

static int sparse_lu_choose_column(
	struct simplex_SparseLu *factor, const int first)
{
	int i, k, column, best = -1;
	int n = factor->dimension;
	double global_maximum = 0.;
	double best_score = 0.;

	for (column = first; column < n; column++) {
		factor->work_column[column] = 0;
		factor->work_value[column] = 0.;
		factor->pivot_row[column] = -1;
	}
	for (i = first; i < n; i++) {
		const struct simplex_SparseRow *row = factor->row + i;

		for (k = 0; k < row->count; k++) {
			factor->factor_work++;
			column = row->column[k];
			if (column >= first) {
				double magnitude =
					__lp_simplex_ABS__(row->value[k]);

				factor->work_column[column]++;
				if (magnitude > factor->work_value[column]) {
					factor->work_value[column] = magnitude;
					factor->pivot_row[column] = i;
				}
				global_maximum = __lp_simplex_MAX__(
					global_maximum, magnitude);
			}
		}
	}
	for (column = first; column < n; column++) {
		int row_nonzeros = 0;
		double score;
		const struct simplex_SparseRow *row;

		if (factor->pivot_row[column] < 0 ||
			factor->work_value[column] <
				SPARSE_LU_MARKOWITZ_THRESHOLD * global_maximum)
			continue;
		row = factor->row + factor->pivot_row[column];
		for (k = 0; k < row->count; k++) {
			factor->factor_work++;
			if (row->column[k] >= first)
				row_nonzeros++;
		}
		score = (double)(factor->work_column[column] - 1) *
			(double)(row_nonzeros - 1);
		if (best < 0 || score < best_score ||
			(score == best_score &&
				factor->work_value[column] >
					factor->work_value[best])) {
			best = column;
			best_score = score;
		}
	}
	return best;
}

static int sparse_lu_factorize_assembled(struct simplex_SparseLu *factor)
{
	int i, k;
	int n = factor->dimension;

	if (n >= SPARSE_LU_GLOBAL_MARKOWITZ_LIMIT) {
		if (sparse_lu_build_column_rows(factor) ==
			lp_simplex_EXIT_FAILURE)
			return lp_simplex_EXIT_FAILURE;
		for (i = 0; i < n; i++)
			factor->pivot_row[i] = -1;
	}
	for (k = 0; k < n; k++) {
		int pivot_row = -1;
		double largest = 0.;
		/* A full active-matrix scan per pivot becomes more expensive than
		 * the extra fill of the preordered columns on large bases. */
		int pivot_column = n >= SPARSE_LU_GLOBAL_MARKOWITZ_LIMIT
					   ? k
					   : sparse_lu_choose_column(factor, k);

		if (pivot_column < 0)
			return lp_simplex_EXIT_FAILURE;
		sparse_lu_swap_columns(factor, k, pivot_column);
		if (n >= SPARSE_LU_GLOBAL_MARKOWITZ_LIMIT) {
			const struct simplex_SparseColumnRows *rows =
				factor->column_rows + k;

			for (i = 0; i < rows->count; i++) {
				int row_index = rows->row[i];
				int position;
				double magnitude;

				factor->factor_work++;
				if (row_index < k)
					continue;
				position = sparse_row_find_counted(
					factor, factor->row + row_index, k);
				if (position < 0)
					continue;
				magnitude = __lp_simplex_ABS__(
					factor->row[row_index].value[position]);
				if (magnitude > largest) {
					largest = magnitude;
					pivot_row = row_index;
				}
			}
		} else
			for (i = k; i < n; i++) {
				int position = sparse_row_find_counted(
					factor, factor->row + i, k);

				factor->factor_work++;
				if (position >= 0) {
					double magnitude = __lp_simplex_ABS__(
						factor->row[i].value[position]);

					if (magnitude > largest) {
						largest = magnitude;
						pivot_row = i;
					}
				}
			}
		if (pivot_row < 0 || largest <= SPARSE_LU_PIVOT_TOLERANCE)
			return lp_simplex_EXIT_FAILURE;
		if (pivot_row != k) {
			struct simplex_SparseRow swapped = factor->row[k];
			int permutation = factor->permutation[k];

			factor->row[k] = factor->row[pivot_row];
			factor->row[pivot_row] = swapped;
			factor->permutation[k] = factor->permutation[pivot_row];
			factor->permutation[pivot_row] = permutation;
			if (n >= SPARSE_LU_GLOBAL_MARKOWITZ_LIMIT &&
				(sparse_lu_append_row_columns(factor, k) ==
						lp_simplex_EXIT_FAILURE ||
					sparse_lu_append_row_columns(
						factor, pivot_row) ==
						lp_simplex_EXIT_FAILURE))
				return lp_simplex_EXIT_FAILURE;
		}
		{
			int diagonal = sparse_row_find_counted(
				factor, factor->row + k, k);
			double pivot = factor->row[k].value[diagonal];

			if (n >= SPARSE_LU_GLOBAL_MARKOWITZ_LIMIT) {
				const struct simplex_SparseColumnRows *rows =
					factor->column_rows + k;

				for (i = 0; i < rows->count; i++) {
					int row_index = rows->row[i];
					int position;
					double multiplier;

					if (row_index <= k ||
						factor->pivot_row[row_index] ==
							k)
						continue;
					factor->pivot_row[row_index] = k;
					position = sparse_row_find_counted(
						factor, factor->row + row_index,
						k);
					if (position < 0)
						continue;
					multiplier = factor->row[row_index]
							     .value[position] /
						     pivot;
					if (sparse_lu_eliminate_row(factor,
						    row_index, k, k,
						    multiplier) ==
						lp_simplex_EXIT_FAILURE)
						return lp_simplex_EXIT_FAILURE;
				}
			} else
				for (i = k + 1; i < n; i++) {
					int position = sparse_row_find_counted(
						factor, factor->row + i, k);

					if (position >= 0) {
						double multiplier =
							factor->row[i].value
								[position] /
							pivot;

						if (sparse_lu_eliminate_row(
							    factor, i, k, k,
							    multiplier) ==
							lp_simplex_EXIT_FAILURE)
							return lp_simplex_EXIT_FAILURE;
					}
				}
		}
	}
	for (i = 0; i < n; i++) {
		int diagonal =
			sparse_row_find_counted(factor, factor->row + i, i);

		if (diagonal < 0 ||
			__lp_simplex_ABS__(factor->row[i].value[diagonal]) <=
				SPARSE_LU_PIVOT_TOLERANCE)
			return lp_simplex_EXIT_FAILURE;
		factor->diagonal_position[i] = diagonal;
		factor->diagonal_value[i] = factor->row[i].value[diagonal];
	}
	/* The lists used during elimination contain stale duplicate entries.
	 * Rebuild exact column adjacency for reach-based triangular solves. */
	if (sparse_lu_build_column_rows(factor) == lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	if (sparse_lu_pack_rows(factor) == lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	return sparse_lu_ft_initialize(factor);
}

int simplex_sparse_lu_factorize(struct simplex_SparseLu *factor,
	const struct simplex_CscMatrix *matrix, const int structural_columns,
	const int *basis)
{
	int result;
	factor->factor_work = 0;
	factor->factor_valid = factor->ft_active = factor->storage_failed = 0;

	if (sparse_lu_assemble(factor, matrix, structural_columns, basis) ==
		lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	result = sparse_lu_factorize_assembled(factor);
	factor->factor_valid = result == lp_simplex_EXIT_SUCCESS;
	return result;
}

int simplex_sparse_lu_factorize_submatrix(struct simplex_SparseLu *factor,
	const struct simplex_CscMatrix *matrix, const int *columns,
	const int *row_to_core)
{
	int result;
	factor->factor_work = 0;
	factor->factor_valid = factor->ft_active = factor->storage_failed = 0;

	if (sparse_lu_assemble_submatrix(factor, matrix, columns,
		    row_to_core) == lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	result = sparse_lu_factorize_assembled(factor);
	factor->factor_valid = result == lp_simplex_EXIT_SUCCESS;
	return result;
}

static int sparse_lu_ft_solve(
	struct simplex_SparseLu *factor, double *vector, const int transpose)
{
	int i, k;
	int n = factor->dimension;
	double *work = factor->solve_work;

	if (!transpose) {
		for (i = 0; i < n; i++) {
			work[i] = vector[factor->permutation[i]] *
				  factor->row_scale[factor->permutation[i]];
		}
		/* The initial L is immutable. */
		for (i = 0; i < n; i++) {
			int start = factor->packed_start[i];
			int diagonal = start + factor->diagonal_position[i];
			double value = work[i];

			for (k = start; k < diagonal; k++)
				value -= factor->packed_value[k] *
					 work[factor->packed_column[k]];
			work[i] = value;
		}
		sparse_lu_ft_apply_rows(factor, work, 0);
		lp_simplex_memcpy(factor->ft_spike_cache, work,
			(size_t)n * sizeof(double));
		factor->ft_spike_valid = 1;
		if (sparse_lu_ft_u_inverse(factor, work, 0) ==
			lp_simplex_EXIT_FAILURE)
			return lp_simplex_EXIT_FAILURE;
		for (i = 0; i < n; i++) {
			int column = factor->column_permutation[i];

			vector[column] = work[i] * factor->column_scale[column];
			factor->ft_direction_cache[column] = vector[column];
		}
	} else {
		int singleton = -1;
		double singleton_value = 0.;

		for (i = 0; i < n; i++)
			if (vector[i] != 0.) {
				if (singleton >= 0) {
					singleton = -2;
					break;
				}
				singleton = i;
				singleton_value = vector[i];
			}
		for (i = 0; i < n; i++)
			work[i] = vector[factor->column_permutation[i]] *
				  factor->column_scale[factor
						  ->column_permutation[i]];
		if (sparse_lu_ft_u_inverse(factor, work, 1) ==
			lp_simplex_EXIT_FAILURE)
			return lp_simplex_EXIT_FAILURE;
		factor->ft_btran_pivot = -1;
		if (singleton >= 0 && singleton_value != 0.)
			for (i = 0; i < n; i++)
				if (factor->column_permutation[i] ==
					singleton) {
					double normalization =
						singleton_value *
						factor->column_scale[singleton];

					for (k = 0; k < n; k++)
						factor->ft_btran_cache[k] =
							work[k] / normalization;
					factor->ft_btran_pivot = i;
					break;
				}
		sparse_lu_ft_apply_rows(factor, work, 1);
		for (i = n - 1; i >= 0; i--) {
			int start = factor->packed_start[i];
			int diagonal = start + factor->diagonal_position[i];
			double value = work[i];

			if (value == 0.)
				continue;
			for (k = start; k < diagonal; k++)
				work[factor->packed_column[k]] -=
					factor->packed_value[k] * value;
		}
		for (i = 0; i < n; i++)
			vector[factor->permutation[i]] =
				work[i] *
				factor->row_scale[factor->permutation[i]];
	}
	return lp_simplex_EXIT_SUCCESS;
}

int simplex_sparse_lu_solve(
	struct simplex_SparseLu *factor, double *vector, const int transpose)
{
	int i, k;
	int n = factor->dimension;
	int *active = factor->work_column;
	int active_count = 0;
	int use_reach;
	double *work = factor->solve_work;

	if (!factor->factor_valid)
		return lp_simplex_EXIT_FAILURE;
	if (factor->ft_active && factor->ft_update_count > 0)
		return sparse_lu_ft_solve(factor, vector, transpose);
	if (!transpose) {
		for (i = 0; i < n; i++) {
			work[i] = vector[factor->permutation[i]] *
				  factor->row_scale[factor->permutation[i]];
			active[i] = work[i] != 0.;
			active_count += active[i];
		}
		use_reach = n > 0 && active_count <= (n - 1) / 8;
		for (i = 0; i < n; i++) {
			int start = factor->packed_start[i];
			int diagonal = start + factor->diagonal_position[i];
			double value = work[i];

			if (use_reach && !active[i])
				continue;
			for (k = start; k < diagonal; k++)
				value -= factor->packed_value[k] *
					 work[factor->packed_column[k]];
			work[i] = value;
			if (use_reach && value != 0.) {
				const struct simplex_SparseColumnRows *rows =
					factor->column_rows + i;

				for (k = 0; k < rows->count; k++)
					if (rows->row[k] > i)
						active[rows->row[k]] = 1;
			}
		}
		if (use_reach)
			for (i = 0; i < n; i++)
				active[i] = work[i] != 0.;
		for (i = n - 1; i >= 0; i--) {
			int start = factor->packed_start[i];
			int diagonal = start + factor->diagonal_position[i];
			int end = factor->packed_start[i + 1];
			double value = work[i];

			if (use_reach && !active[i])
				continue;
			for (k = diagonal + 1; k < end; k++)
				value -= factor->packed_value[k] *
					 work[factor->packed_column[k]];
			value /= factor->diagonal_value[i];
			work[i] = value;
			if (use_reach && value != 0.) {
				const struct simplex_SparseColumnRows *rows =
					factor->column_rows + i;

				for (k = 0; k < rows->count; k++)
					if (rows->row[k] < i)
						active[rows->row[k]] = 1;
			}
		}
		for (i = 0; i < n; i++)
			vector[factor->column_permutation[i]] =
				work[i] *
				factor->column_scale[factor
						->column_permutation[i]];
	} else {
		for (i = 0; i < n; i++)
			work[i] = vector[factor->column_permutation[i]] *
				  factor->column_scale[factor
						  ->column_permutation[i]];
		for (i = 0; i < n; i++) {
			int diagonal = factor->packed_start[i] +
				       factor->diagonal_position[i];
			int end = factor->packed_start[i + 1];

			if (work[i] == 0.)
				continue;
			work[i] /= factor->diagonal_value[i];
			for (k = diagonal + 1; k < end; k++)
				work[factor->packed_column[k]] -=
					factor->packed_value[k] * work[i];
		}
		for (i = n - 1; i >= 0; i--) {
			int start = factor->packed_start[i];
			int diagonal = start + factor->diagonal_position[i];

			if (work[i] == 0.)
				continue;
			for (k = start; k < diagonal; k++)
				work[factor->packed_column[k]] -=
					factor->packed_value[k] * work[i];
		}
		for (i = 0; i < n; i++)
			vector[factor->permutation[i]] =
				work[i] *
				factor->row_scale[factor->permutation[i]];
	}
	return lp_simplex_EXIT_SUCCESS;
}

int simplex_sparse_lu_solve_pair(
	struct simplex_SparseLu *factor, double *first, double *second)
{
	int i, k;
	int n = factor->dimension;
	int *active = factor->work_column;
	int active_count = 0;
	int use_reach;
	double *a = factor->solve_work;
	double *b = factor->work_value;

	if (!factor->factor_valid)
		return lp_simplex_EXIT_FAILURE;
	if (factor->ft_active && factor->ft_update_count > 0) {
		if (sparse_lu_ft_solve(factor, first, 0) ==
			lp_simplex_EXIT_FAILURE)
			return lp_simplex_EXIT_FAILURE;
		lp_simplex_memcpy(factor->work_value, factor->ft_spike_cache,
			(size_t)n * sizeof(double));
		if (sparse_lu_ft_solve(factor, second, 0) ==
			lp_simplex_EXIT_FAILURE)
			return lp_simplex_EXIT_FAILURE;
		lp_simplex_memcpy(factor->ft_spike_cache, factor->work_value,
			(size_t)n * sizeof(double));
		lp_simplex_memcpy(factor->ft_direction_cache, first,
			(size_t)n * sizeof(double));
		factor->ft_spike_valid = 1;
		return lp_simplex_EXIT_SUCCESS;
	}
	for (i = 0; i < n; i++) {
		int row = factor->permutation[i];
		double scale = factor->row_scale[row];

		a[i] = first[row] * scale;
		b[i] = second[row] * scale;
		active[i] = a[i] != 0. || b[i] != 0.;
		active_count += active[i];
	}
	use_reach = n > 0 && active_count <= (n - 1) / 8;
	for (i = 0; i < n; i++) {
		int start = factor->packed_start[i];
		int diagonal = start + factor->diagonal_position[i];
		double av = a[i], bv = b[i];

		if (use_reach && !active[i])
			continue;
		for (k = start; k < diagonal; k++) {
			int column = factor->packed_column[k];
			double value = factor->packed_value[k];

			av -= value * a[column];
			bv -= value * b[column];
		}
		a[i] = av;
		b[i] = bv;
		if (use_reach && (av != 0. || bv != 0.)) {
			const struct simplex_SparseColumnRows *rows =
				factor->column_rows + i;

			for (k = 0; k < rows->count; k++)
				if (rows->row[k] > i)
					active[rows->row[k]] = 1;
		}
	}
	if (use_reach)
		for (i = 0; i < n; i++)
			active[i] = a[i] != 0. || b[i] != 0.;
	for (i = n - 1; i >= 0; i--) {
		int diagonal =
			factor->packed_start[i] + factor->diagonal_position[i];
		int end = factor->packed_start[i + 1];
		double av = a[i], bv = b[i];

		if (use_reach && !active[i])
			continue;
		for (k = diagonal + 1; k < end; k++) {
			int column = factor->packed_column[k];
			double value = factor->packed_value[k];

			av -= value * a[column];
			bv -= value * b[column];
		}
		a[i] = av / factor->diagonal_value[i];
		b[i] = bv / factor->diagonal_value[i];
		if (use_reach && (a[i] != 0. || b[i] != 0.)) {
			const struct simplex_SparseColumnRows *rows =
				factor->column_rows + i;

			for (k = 0; k < rows->count; k++)
				if (rows->row[k] < i)
					active[rows->row[k]] = 1;
		}
	}
	for (i = 0; i < n; i++) {
		int column = factor->column_permutation[i];
		double scale = factor->column_scale[column];

		first[column] = a[i] * scale;
		second[column] = b[i] * scale;
	}
	return lp_simplex_EXIT_SUCCESS;
}

static void sparse_lu_ft_remove_row(
	struct simplex_SparseRow *column, const int row)
{
	int k, out = 0;

	for (k = 0; k < column->count; k++)
		if (column->column[k] != row) {
			column->column[out] = column->column[k];
			column->value[out++] = column->value[k];
		}
	column->count = out;
}

int simplex_sparse_lu_ft_update(struct simplex_SparseLu *factor,
	const int leaving_position, const double *direction)
{
	struct simplex_SparseRow *row_eta;
	struct simplex_SparseRow *new_column;
	double *spike = factor->work_value;
	double *partial_btran = factor->solve_work;
	double pivot, dot = 0.;
	double column_scale;
	int i, k, p = -1, position;
	int new_nonzeros = 0;
	long removed_row = 0;
	long projected_nonzeros;
	int n = factor->dimension;

	factor->storage_failed = 0;
	if (!factor->factor_valid)
		return lp_simplex_EXIT_FAILURE;
	if (!factor->ft_active ||
		factor->ft_update_count >= factor->ft_update_capacity)
		return 1;
	for (i = 0; i < n; i++)
		if (factor->column_permutation[i] == leaving_position) {
			p = i;
			break;
		}
	if (p < 0)
		return lp_simplex_EXIT_FAILURE;
	column_scale = factor->column_scale[leaving_position];
	if (__lp_simplex_ABS__(column_scale) <= SPARSE_LU_PIVOT_TOLERANCE)
		return lp_simplex_EXIT_FAILURE;
	/* Reuse a spike only for the exact direction it produced.  Refinement,
	 * paired solves, or an auxiliary RHS can replace the last FTRAN cache.
	 * A mismatched direction is reconstructed from the current U instead. */
	if (factor->ft_spike_valid &&
		lp_simplex_memcmp(direction, factor->ft_direction_cache,
			(size_t)n * sizeof(double)) == 0)
		lp_simplex_memcpy(spike, factor->ft_spike_cache,
			(size_t)n * sizeof(double));
	else {
		for (i = 0; i < n; i++)
			spike[i] = direction[factor->column_permutation[i]] /
				   factor->column_scale[factor
						   ->column_permutation[i]];
		sparse_lu_ft_u_product(factor, spike);
	}
	for (i = 0; i < n; i++)
		spike[i] *= column_scale;

	/* Tomlin's shortcut: r^T = e_p^T-u_pp e_p^T U^{-1}. */
	if (factor->ft_btran_pivot == p)
		lp_simplex_memcpy(partial_btran, factor->ft_btran_cache,
			(size_t)n * sizeof(double));
	else {
		for (i = 0; i < n; i++)
			partial_btran[i] = 0.;
		partial_btran[p] = 1.;
		if (sparse_lu_ft_u_inverse(factor, partial_btran, 1) ==
			lp_simplex_EXIT_FAILURE)
			return 1;
	}
	pivot = sparse_lu_ft_pivot(factor, p);
	row_eta = factor->ft_row_eta + factor->ft_update_count;
	row_eta->count = 0;
	for (i = 0; i < n; i++) {
		double value = i == p ? 0. : -pivot * partial_btran[i];

		if (value != 0.) {
			if (sparse_row_reserve(
				    factor, row_eta, row_eta->count + 1) ==
				lp_simplex_EXIT_FAILURE)
				return lp_simplex_EXIT_FAILURE;
			row_eta->column[row_eta->count] = i;
			row_eta->value[row_eta->count++] = value;
			dot += value * spike[i];
		}
	}
	spike[p] -= dot;
	if (__lp_simplex_ABS__(spike[p]) <= SPARSE_LU_PIVOT_TOLERANCE)
		return 1;
	for (i = 0; i < n; i++) {
		if (spike[i] != 0.)
			new_nonzeros++;
		if (i != p) {
			const struct simplex_SparseRow *column =
				factor->ft_column + i;

			for (k = 0; k < column->count; k++)
				if (column->column[k] == p) {
					removed_row++;
					break;
				}
		}
	}
	projected_nonzeros = factor->ft_nonzeros - factor->ft_column[p].count -
			     removed_row + new_nonzeros;
	/* FT fill can change abruptly on hyper-sparse models.  Reinvert before
	 * committing an update whose U file or accumulated row file would cost
	 * more to traverse than a fresh sparse factorization. */
	if (factor->ft_update_count > 0 &&
		(projected_nonzeros >
				__lp_simplex_MAX__(
					2 * factor->ft_initial_nonzeros,
					factor->ft_initial_nonzeros + 8L * n) ||
			factor->ft_row_nonzeros + row_eta->count > 8L * n))
		return 1;

	/* Delete the old pivot column and zero the leaving row in every other
	 * eta.  In factor-product order the new eta is prepended (equivalently it
	 * is appended to the inverse-application file used in the literature). */
	new_column = factor->ft_column + p;
	/* All potentially failing allocation precedes the mathematical commit.
	 * Uncommitted row_eta scratch and enlarged capacities may be retained. */
	if (sparse_row_reserve(factor, new_column, new_nonzeros) ==
		lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	for (i = 0; i < n; i++)
		if (i != p)
			sparse_lu_ft_remove_row(factor->ft_column + i, p);
	new_column->count = 0;
	for (i = 0; i < n; i++)
		if (spike[i] != 0.) {
			new_column->column[new_column->count] = i;
			new_column->value[new_column->count++] = spike[i];
		}
	position = factor->ft_pivot_position[p];
	for (k = position; k > 0; k--) {
		factor->ft_order[k] = factor->ft_order[k - 1];
		factor->ft_pivot_position[factor->ft_order[k]] = k;
	}
	factor->ft_order[0] = p;
	factor->ft_pivot_position[p] = 0;
	factor->ft_row_pivot[factor->ft_update_count] = p;
	factor->ft_nonzeros = projected_nonzeros;
	factor->ft_row_nonzeros += row_eta->count;
	factor->ft_update_count++;
	factor->ft_spike_valid = 0;
	factor->ft_btran_pivot = -1;
	return lp_simplex_EXIT_SUCCESS;
}

int simplex_sparse_lu_ft_active(const struct simplex_SparseLu *factor)
{
	return factor != NULL && factor->ft_active;
}
