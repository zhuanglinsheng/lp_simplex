/* Sparse row-oriented LU with partial pivoting for simplex basis matrices. */
#include "simplex_sparse_lu.h"
#include "utils.h"
#include <lp_simplex/status.h>

#define SPARSE_LU_PIVOT_TOLERANCE 1e-13
#define SPARSE_LU_MARKOWITZ_THRESHOLD 1e-1
#define SPARSE_LU_GLOBAL_MARKOWITZ_LIMIT 2048


static void sparse_row_destroy(struct simplex_SparseRow *row)
{
	lp_simplex_free(row->column);
	lp_simplex_free(row->value);
	row->column = NULL;
	row->value = NULL;
	row->count = 0;
	row->capacity = 0;
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


static int sparse_row_find(const struct simplex_SparseRow *row, const int column)
{
	int left = 0, right = row->count - 1;
	while (left <= right) {
		int middle = left + (right - left) / 2;
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
	lp_simplex_memset(factor, 0, sizeof(*factor));
	factor->dimension = dimension;
	factor->permutation = (int *)lp_simplex_malloc(
		(size_t)dimension * sizeof(int));
	factor->column_permutation = (int *)lp_simplex_malloc(
		(size_t)dimension * sizeof(int));
	factor->column_scale = (double *)lp_simplex_malloc(
		(size_t)dimension * sizeof(double));
	factor->row_scale = (double *)lp_simplex_malloc(
		(size_t)dimension * sizeof(double));
	factor->work_column = (int *)lp_simplex_malloc(
		(size_t)dimension * sizeof(int));
	factor->pivot_row = (int *)lp_simplex_malloc(
		(size_t)dimension * sizeof(int));
	factor->work_value = (double *)lp_simplex_malloc(
		(size_t)dimension * sizeof(double));
	factor->solve_work = (double *)lp_simplex_malloc(
		(size_t)dimension * sizeof(double));
	factor->row = (struct simplex_SparseRow *)lp_simplex_malloc(
		(size_t)dimension * sizeof(*factor->row));
	if (factor->permutation == NULL || factor->column_permutation == NULL ||
	    factor->column_scale == NULL || factor->row_scale == NULL ||
	    factor->work_column == NULL || factor->pivot_row == NULL ||
	    factor->work_value == NULL || factor->solve_work == NULL ||
	    factor->row == NULL) {
		simplex_sparse_lu_destroy(factor);
		return lp_simplex_EXIT_FAILURE;
	}
	lp_simplex_memset(factor->row, 0,
		(size_t)dimension * sizeof(*factor->row));
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
	lp_simplex_free(factor->row);
	lp_simplex_free(factor->permutation);
	lp_simplex_free(factor->column_permutation);
	lp_simplex_free(factor->column_scale);
	lp_simplex_free(factor->row_scale);
	lp_simplex_free(factor->work_column);
	lp_simplex_free(factor->pivot_row);
	lp_simplex_free(factor->work_value);
	lp_simplex_free(factor->solve_work);
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


static void sparse_row_swap_columns(
		struct simplex_SparseRow *row, const int left, const int right)
{
	int left_position = sparse_row_find(row, left);
	int right_position = sparse_row_find(row, right);
	if (left_position >= 0 && right_position >= 0) {
		double value = row->value[left_position];
		row->value[left_position] = row->value[right_position];
		row->value[right_position] = value;
	} else if (left_position >= 0) {
		double value = row->value[left_position];
		while (left_position + 1 < row->count &&
		       row->column[left_position + 1] < right) {
			row->column[left_position] = row->column[left_position + 1];
			row->value[left_position] = row->value[left_position + 1];
			left_position++;
		}
		row->column[left_position] = right;
		row->value[left_position] = value;
	} else if (right_position >= 0) {
		double value = row->value[right_position];
		while (right_position > 0 && row->column[right_position - 1] > left) {
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
		sparse_row_swap_columns(factor->row + i, left, right);
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
		for (k = 0; k < row->count; k++)
			if (row->column[k] >= first)
				row_nonzeros++;
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


int simplex_sparse_lu_factorize(
		struct simplex_SparseLu *factor,
		const struct simplex_CscMatrix *matrix,
		const int structural_columns, const int *basis)
{
	int i, k;
	int n = factor->dimension;
	if (sparse_lu_assemble(factor, matrix, structural_columns, basis) ==
	    lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
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
		for (i = k; i < n; i++) {
			int position = sparse_row_find(factor->row + i, k);
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
		}
		{
			int diagonal = sparse_row_find(factor->row + k, k);
			double pivot = factor->row[k].value[diagonal];
			for (i = k + 1; i < n; i++) {
				int position = sparse_row_find(factor->row + i, k);
				if (position >= 0) {
					double multiplier = factor->row[i].value[position] / pivot;
					if (sparse_lu_eliminate_row(factor, i, k, k,
							multiplier) == lp_simplex_EXIT_FAILURE)
						return lp_simplex_EXIT_FAILURE;
				}
			}
		}
	}
	return lp_simplex_EXIT_SUCCESS;
}


int simplex_sparse_lu_solve(
		struct simplex_SparseLu *factor, double *vector, const int transpose)
{
	int i, k;
	int n = factor->dimension;
	double *work = factor->solve_work;
	if (!transpose) {
		for (i = 0; i < n; i++)
			work[i] = vector[factor->permutation[i]] *
				factor->row_scale[factor->permutation[i]];
		for (i = 0; i < n; i++) {
			const struct simplex_SparseRow *row = factor->row + i;
			double value = work[i];
			for (k = 0; k < row->count && row->column[k] < i; k++)
				value -= row->value[k] * work[row->column[k]];
			work[i] = value;
		}
		for (i = n - 1; i >= 0; i--) {
			const struct simplex_SparseRow *row = factor->row + i;
			double diagonal = 0.;
			double value = work[i];
			for (k = 0; k < row->count; k++) {
				if (row->column[k] == i)
					diagonal = row->value[k];
				else if (row->column[k] > i)
					value -= row->value[k] * work[row->column[k]];
			}
			if (__lp_simplex_ABS__(diagonal) <= SPARSE_LU_PIVOT_TOLERANCE)
				return lp_simplex_EXIT_FAILURE;
			work[i] = value / diagonal;
		}
		for (i = 0; i < n; i++)
			vector[factor->column_permutation[i]] = work[i] *
				factor->column_scale[factor->column_permutation[i]];
	} else {
		for (i = 0; i < n; i++)
			work[i] = vector[factor->column_permutation[i]] *
				factor->column_scale[factor->column_permutation[i]];
		for (i = 0; i < n; i++) {
			const struct simplex_SparseRow *row = factor->row + i;
			double diagonal = 0.;
			for (k = 0; k < row->count; k++)
				if (row->column[k] == i)
					diagonal = row->value[k];
			if (__lp_simplex_ABS__(diagonal) <= SPARSE_LU_PIVOT_TOLERANCE)
				return lp_simplex_EXIT_FAILURE;
			work[i] /= diagonal;
			for (k = 0; k < row->count; k++)
				if (row->column[k] > i)
					work[row->column[k]] -= row->value[k] * work[i];
		}
		for (i = n - 1; i >= 0; i--) {
			const struct simplex_SparseRow *row = factor->row + i;
			for (k = 0; k < row->count && row->column[k] < i; k++)
				work[row->column[k]] -= row->value[k] * work[i];
		}
		for (i = 0; i < n; i++)
			vector[factor->permutation[i]] = work[i] *
				factor->row_scale[factor->permutation[i]];
	}
	return lp_simplex_EXIT_SUCCESS;
}
