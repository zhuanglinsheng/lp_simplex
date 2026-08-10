/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
/* Sparse row-oriented LU with partial pivoting for simplex basis matrices. */
#include "simplex_sparse_lu.h"
#include "utils.h"

#include <lp_simplex/status.h>


#define SPARSE_LU_PIVOT_TOLERANCE 1e-13
#define SPARSE_LU_MARKOWITZ_THRESHOLD 1e-1
#define SPARSE_LU_GLOBAL_MARKOWITZ_LIMIT 1024


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


static int sparse_column_rows_append(
		struct simplex_SparseColumnRows *column, const int row)
{
	int *grown_rows;
	int capacity;
	if (column->count < column->capacity) {
		column->row[column->count++] = row;
		return lp_simplex_EXIT_SUCCESS;
	}
	capacity = column->capacity > 0 ? 2 * column->capacity : 8;
	grown_rows = (int *)lp_simplex_malloc((size_t)capacity * sizeof(int));
	if (grown_rows == NULL)
		return lp_simplex_EXIT_FAILURE;
	if (column->count > 0)
		lp_simplex_memcpy(grown_rows, column->row,
			(size_t)column->count * sizeof(int));
	lp_simplex_free(column->row);
	column->row = grown_rows;
	column->capacity = capacity;
	column->row[column->count++] = row;
	return lp_simplex_EXIT_SUCCESS;
}


static int sparse_row_reserve(struct simplex_SparseRow *row, const int capacity)
{
	int *column;
	double *value;
	int grown;
	if (capacity <= row->capacity)
		return lp_simplex_EXIT_SUCCESS;
	grown = row->capacity > 0 ? row->capacity : 4;
	while (grown < capacity)
		grown = grown < 1024 ? 2 * grown : grown + grown / 2;
	column = (int *)lp_simplex_malloc((size_t)grown * sizeof(int));
	value = (double *)lp_simplex_malloc((size_t)grown * sizeof(double));
	if (column == NULL || value == NULL) {
		lp_simplex_free(column);
		lp_simplex_free(value);
		return lp_simplex_EXIT_FAILURE;
	}
	if (row->count > 0) {
		lp_simplex_memcpy(column, row->column,
			(size_t)row->count * sizeof(int));
		lp_simplex_memcpy(value, row->value,
			(size_t)row->count * sizeof(double));
	}
	lp_simplex_free(row->column);
	lp_simplex_free(row->value);
	row->column = column;
	row->value = value;
	row->capacity = grown;
	return lp_simplex_EXIT_SUCCESS;
}


static int sparse_row_find_counted(
		struct simplex_SparseLu *factor,
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


int simplex_sparse_lu_create(struct simplex_SparseLu *factor, const int dimension)
{
	int *integer;
	double *numeric;
	lp_simplex_memset(factor, 0, sizeof(*factor));
	factor->dimension = dimension;
	factor->integer_storage = (int *)lp_simplex_malloc(
		((size_t)6 * dimension + 1) * sizeof(int));
	factor->numeric_storage = (double *)lp_simplex_malloc(
		(size_t)5 * dimension * sizeof(double));
	integer = factor->integer_storage;
	numeric = factor->numeric_storage;
	if (integer != NULL) {
		factor->permutation = integer;
		factor->column_permutation = factor->permutation + dimension;
		factor->work_column = factor->column_permutation + dimension;
		factor->pivot_row = factor->work_column + dimension;
		factor->diagonal_position = factor->pivot_row + dimension;
		factor->packed_start = factor->diagonal_position + dimension;
	}
	if (numeric != NULL) {
		factor->column_scale = numeric;
		factor->row_scale = factor->column_scale + dimension;
		factor->diagonal_value = factor->row_scale + dimension;
		factor->work_value = factor->diagonal_value + dimension;
		factor->solve_work = factor->work_value + dimension;
	}
	factor->row = (struct simplex_SparseRow *)lp_simplex_malloc(
		(size_t)dimension * sizeof(*factor->row));
	factor->column_rows = (struct simplex_SparseColumnRows *)lp_simplex_malloc(
		(size_t)dimension * sizeof(*factor->column_rows));
	if (factor->integer_storage == NULL || factor->numeric_storage == NULL ||
	    factor->row == NULL || factor->column_rows == NULL) {
		simplex_sparse_lu_destroy(factor);
		return lp_simplex_EXIT_FAILURE;
	}
	lp_simplex_memset(factor->row, 0,
		(size_t)dimension * sizeof(*factor->row));
	lp_simplex_memset(factor->column_rows, 0,
		(size_t)dimension * sizeof(*factor->column_rows));
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
	lp_simplex_free(factor->row);
	lp_simplex_free(factor->column_rows);
	lp_simplex_free(factor->integer_storage);
	lp_simplex_free(factor->numeric_storage);
	lp_simplex_free(factor->packed_column);
	lp_simplex_free(factor->packed_value);
	lp_simplex_memset(factor, 0, sizeof(*factor));
}


static int sparse_lu_assemble(
		struct simplex_SparseLu *factor,
		const struct simplex_CscMatrix *matrix,
		const int structural_columns, const int *basis)
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
		count[j] = variable >= structural_columns ? 1 :
			matrix->column_start[variable + 1] -
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
				factor->row_scale[row_index] = __lp_simplex_MAX__(
					factor->row_scale[row_index],
					__lp_simplex_ABS__(matrix->value[k] * column_scale));
			}
		}
	}
	for (i = 0; i < n; i++)
		factor->row_scale[i] = factor->row_scale[i] > 0.
			? 1. / factor->row_scale[i] : 1.;
	for (j = 1; j < n; j++) {
		int original = factor->column_permutation[j];
		int position = j;
		while (position > 0 &&
		       count[factor->column_permutation[position - 1]] > count[original]) {
			factor->column_permutation[position] =
				factor->column_permutation[position - 1];
			position--;
		}
		factor->column_permutation[position] = original;
	}
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
		if (sparse_row_reserve(factor->row + i, count[i]) ==
		    lp_simplex_EXIT_FAILURE)
			return lp_simplex_EXIT_FAILURE;
	}
	for (j = 0; j < n; j++) {
		int variable = basis[factor->column_permutation[j]];
		if (variable >= structural_columns) {
			int row_index = variable - structural_columns;
			struct simplex_SparseRow *row = factor->row + row_index;
			row->column[row->count] = j;
			row->value[row->count++] = -factor->column_scale[
				factor->column_permutation[j]] *
				factor->row_scale[row_index];
		} else {
			for (k = matrix->column_start[variable];
			     k < matrix->column_start[variable + 1]; k++) {
				int row_index = matrix->row_index[k];
				struct simplex_SparseRow *row = factor->row + row_index;
				row->column[row->count] = j;
				row->value[row->count++] = matrix->value[k] *
					factor->column_scale[factor->column_permutation[j]] *
					factor->row_scale[row_index];
			}
		}
	}
	return lp_simplex_EXIT_SUCCESS;
}


static int sparse_lu_assemble_submatrix(
		struct simplex_SparseLu *factor,
		const struct simplex_CscMatrix *matrix,
		const int *columns, const int *row_to_core)
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
			largest = __lp_simplex_MAX__(largest,
				__lp_simplex_ABS__(matrix->value[k]));
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
					__lp_simplex_ABS__(matrix->value[k] * column_scale));
		}
	}
	for (i = 0; i < n; i++)
		factor->row_scale[i] = factor->row_scale[i] > 0.
			? 1. / factor->row_scale[i] : 1.;
	for (j = 1; j < n; j++) {
		int original = factor->column_permutation[j];
		int position = j;
		while (position > 0 &&
		       count[factor->column_permutation[position - 1]] > count[original]) {
			factor->column_permutation[position] =
				factor->column_permutation[position - 1];
			position--;
		}
		factor->column_permutation[position] = original;
	}
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
		if (sparse_row_reserve(factor->row + i, count[i]) ==
		    lp_simplex_EXIT_FAILURE)
			return lp_simplex_EXIT_FAILURE;
	for (j = 0; j < n; j++) {
		int original = factor->column_permutation[j];
		int column = columns[original];
		for (k = matrix->column_start[column];
		     k < matrix->column_start[column + 1]; k++) {
			int row_index = row_to_core[matrix->row_index[k]];
			if (row_index >= 0) {
				struct simplex_SparseRow *row = factor->row + row_index;
				row->column[row->count] = j;
				row->value[row->count++] = matrix->value[k] *
					factor->column_scale[original] *
					factor->row_scale[row_index];
			}
		}
	}
	return lp_simplex_EXIT_SUCCESS;
}


static int sparse_lu_eliminate_row(
		struct simplex_SparseLu *factor, const int target_index,
		const int pivot_index, const int column, const double multiplier)
{
	int a = 0, b = 0, count = 0;
	struct simplex_SparseRow *target = factor->row + target_index;
	const struct simplex_SparseRow *pivot = factor->row + pivot_index;
	while (a < target->count || b < pivot->count) {
		int target_column = a < target->count ? target->column[a] : factor->dimension;
		int pivot_column = b < pivot->count ? pivot->column[b] : factor->dimension;
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
			result_value = target->value[a++] - multiplier * pivot->value[b++];
		}
		if (result_column == column)
			result_value = multiplier;
		if (result_value != 0.) {
			if (target_column > pivot_column &&
			    factor->column_rows != NULL &&
			    sparse_column_rows_append(
				factor->column_rows + result_column,
				target_index) == lp_simplex_EXIT_FAILURE)
				return lp_simplex_EXIT_FAILURE;
			factor->work_column[count] = result_column;
			factor->work_value[count] = result_value;
			count++;
		}
	}
	if (sparse_row_reserve(target, count) == lp_simplex_EXIT_FAILURE)
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
			if (sparse_column_rows_append(
				factor->column_rows + row->column[k], i) ==
			    lp_simplex_EXIT_FAILURE)
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
	for (i = 0; i < factor->dimension; i++)
		nonzeros += factor->row[i].count;
	column = nonzeros > 0 ? (int *)lp_simplex_malloc(
		(size_t)nonzeros * sizeof(int)) : NULL;
	value = nonzeros > 0 ? (double *)lp_simplex_malloc(
		(size_t)nonzeros * sizeof(double)) : NULL;
	if (nonzeros > 0 && (column == NULL || value == NULL)) {
		lp_simplex_free(column);
		lp_simplex_free(value);
		return lp_simplex_EXIT_FAILURE;
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
}


static int sparse_lu_append_row_columns(
		struct simplex_SparseLu *factor, const int row_index)
{
	int k;
	const struct simplex_SparseRow *row = factor->row + row_index;
	for (k = 0; k < row->count; k++)
		if (sparse_column_rows_append(
			factor->column_rows + row->column[k], row_index) ==
		    lp_simplex_EXIT_FAILURE)
			return lp_simplex_EXIT_FAILURE;
	return lp_simplex_EXIT_SUCCESS;
}


static void sparse_row_swap_columns(
		struct simplex_SparseLu *factor, struct simplex_SparseRow *row,
		const int left, const int right)
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
			row->column[left_position] = row->column[left_position + 1];
			row->value[left_position] = row->value[left_position + 1];
			left_position++;
		}
		row->column[left_position] = right;
		row->value[left_position] = value;
	} else if (right_position >= 0) {
		double value = row->value[right_position];
		while (right_position > 0 && row->column[right_position - 1] > left) {
			factor->factor_work++;
			row->column[right_position] = row->column[right_position - 1];
			row->value[right_position] = row->value[right_position - 1];
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
				double magnitude = __lp_simplex_ABS__(row->value[k]);
				factor->work_column[column]++;
				if (magnitude > factor->work_value[column]) {
					factor->work_value[column] = magnitude;
					factor->pivot_row[column] = i;
				}
				global_maximum = __lp_simplex_MAX__(global_maximum, magnitude);
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
		    (score == best_score && factor->work_value[column] >
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
		if (sparse_lu_build_column_rows(factor) == lp_simplex_EXIT_FAILURE)
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
			? k : sparse_lu_choose_column(factor, k);
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
		} else for (i = k; i < n; i++) {
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
			     sparse_lu_append_row_columns(factor, pivot_row) ==
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
					if (row_index <= k || factor->pivot_row[row_index] == k)
						continue;
					factor->pivot_row[row_index] = k;
					position = sparse_row_find_counted(
						factor, factor->row + row_index, k);
					if (position < 0)
						continue;
					multiplier = factor->row[row_index].value[position] / pivot;
					if (sparse_lu_eliminate_row(factor, row_index, k, k,
							multiplier) == lp_simplex_EXIT_FAILURE)
						return lp_simplex_EXIT_FAILURE;
				}
			} else for (i = k + 1; i < n; i++) {
				int position = sparse_row_find_counted(
					factor, factor->row + i, k);
				if (position >= 0) {
					double multiplier = factor->row[i].value[position] / pivot;
					if (sparse_lu_eliminate_row(factor, i, k, k,
							multiplier) == lp_simplex_EXIT_FAILURE)
						return lp_simplex_EXIT_FAILURE;
				}
			}
		}
	}
	for (i = 0; i < n; i++) {
		int diagonal = sparse_row_find_counted(
			factor, factor->row + i, i);
		if (diagonal < 0 || __lp_simplex_ABS__(
				factor->row[i].value[diagonal]) <=
		    SPARSE_LU_PIVOT_TOLERANCE)
			return lp_simplex_EXIT_FAILURE;
		factor->diagonal_position[i] = diagonal;
		factor->diagonal_value[i] = factor->row[i].value[diagonal];
	}
	/* The lists used during elimination contain stale duplicate entries.
	 * Rebuild exact column adjacency for reach-based triangular solves. */
	if (sparse_lu_build_column_rows(factor) == lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	return sparse_lu_pack_rows(factor);
}


int simplex_sparse_lu_factorize(
		struct simplex_SparseLu *factor,
		const struct simplex_CscMatrix *matrix,
		const int structural_columns, const int *basis)
{
	factor->factor_work = 0;
	if (sparse_lu_assemble(factor, matrix, structural_columns, basis) ==
	    lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	return sparse_lu_factorize_assembled(factor);
}


int simplex_sparse_lu_factorize_submatrix(
		struct simplex_SparseLu *factor,
		const struct simplex_CscMatrix *matrix,
		const int *columns, const int *row_to_core)
{
	factor->factor_work = 0;
	if (sparse_lu_assemble_submatrix(factor, matrix, columns, row_to_core) ==
	    lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	return sparse_lu_factorize_assembled(factor);
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
	if (!transpose) {
		for (i = 0; i < n; i++) {
			work[i] = vector[factor->permutation[i]] *
				factor->row_scale[factor->permutation[i]];
			active[i] = work[i] != 0.;
			active_count += active[i];
		}
		use_reach = active_count * 8 < n;
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
			vector[factor->column_permutation[i]] = work[i] *
				factor->column_scale[factor->column_permutation[i]];
	} else {
		for (i = 0; i < n; i++)
			work[i] = vector[factor->column_permutation[i]] *
				factor->column_scale[factor->column_permutation[i]];
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
			vector[factor->permutation[i]] = work[i] *
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
	for (i = 0; i < n; i++) {
		int row = factor->permutation[i];
		double scale = factor->row_scale[row];
		a[i] = first[row] * scale;
		b[i] = second[row] * scale;
		active[i] = a[i] != 0. || b[i] != 0.;
		active_count += active[i];
	}
	use_reach = active_count * 8 < n;
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
		int diagonal = factor->packed_start[i] +
			factor->diagonal_position[i];
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
