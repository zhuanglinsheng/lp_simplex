/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_sparse_vector.c
 * @brief Sparse-vector storage and indexed numerical operations.
 */

#include "simplex_sparse_vector.h"
#include "utils.h"
#include "simplex_storage.h"

#include <lp_simplex/status.h>
#include <math.h>

int simplex_sparse_vector_create(
	struct simplex_SparseVector *vector, const int dimension)
{
	size_t indices, values, marks;
	if (vector == NULL)
		return lp_simplex_EXIT_FAILURE;
	lp_simplex_memset(vector, 0, sizeof(*vector));
	if (dimension < 0 ||
		!simplex_storage_bytes((size_t)dimension, 2 * sizeof(int),
			(size_t)-1, &indices) ||
		!simplex_storage_bytes((size_t)dimension, 2 * sizeof(double),
			(size_t)-1, &values) ||
		!simplex_storage_bytes((size_t)dimension, sizeof(unsigned int),
			(size_t)-1, &marks))
		return lp_simplex_EXIT_FAILURE;
	vector->dimension = dimension;
	vector->count = 0;
	vector->capacity = dimension;
	/* Index/slot and value/dense have identical lifetimes.  Store each pair in
	 * one allocation so a packed workspace has one clear owner per type. */
	vector->index =
		dimension > 0 ? (int *)lp_simplex_malloc(indices) : NULL;
	vector->slot = dimension > 0 && vector->index != NULL
			       ? vector->index + dimension
			       : NULL;
	vector->value =
		dimension > 0 ? (double *)lp_simplex_malloc(values) : NULL;
	vector->dense = dimension > 0 && vector->value != NULL
				? vector->value + dimension
				: NULL;
	vector->mark =
		dimension > 0 ? (unsigned int *)lp_simplex_malloc(marks) : NULL;
	vector->generation = 1;
	if (dimension > 0 &&
		(vector->index == NULL || vector->value == NULL ||
			vector->dense == NULL || vector->mark == NULL ||
			vector->slot == NULL)) {
		simplex_sparse_vector_destroy(vector);
		return lp_simplex_EXIT_FAILURE;
	}
	if (dimension > 0) {
		lp_simplex_memset(
			vector->dense, 0, (size_t)dimension * sizeof(double));
		lp_simplex_memset(vector->mark, 0,
			(size_t)dimension * sizeof(unsigned int));
	}
	return lp_simplex_EXIT_SUCCESS;
}

void simplex_sparse_vector_destroy(struct simplex_SparseVector *vector)
{
	if (vector == NULL)
		return;
	lp_simplex_free(vector->index);
	lp_simplex_free(vector->value);
	lp_simplex_free(vector->mark);
	vector->index = NULL;
	vector->value = NULL;
	vector->dense = NULL;
	vector->mark = NULL;
	vector->slot = NULL;
	vector->dimension = 0;
	vector->count = 0;
	vector->capacity = 0;
	vector->generation = 0;
}

void simplex_sparse_vector_clear(struct simplex_SparseVector *vector)
{
	int k;

	if (vector == NULL)
		return;
	for (k = 0; k < vector->count; k++)
		vector->dense[vector->index[k]] = 0.;
	vector->count = 0;
	vector->generation++;
	if (vector->generation == 0) {
		lp_simplex_memset(vector->mark, 0,
			(size_t)vector->dimension * sizeof(unsigned int));
		vector->generation = 1;
	}
}

void simplex_sparse_vector_compact(struct simplex_SparseVector *vector)
{
	int i, count = 0;
	for (i = 0; i < vector->count; i++) {
		int row = vector->index[i];
		double value = vector->dense[row];
		if (value == 0.) {
			vector->mark[row] = 0;
			continue;
		}
		if (count != i) {
			vector->index[count] = row;
			vector->slot[row] = count;
		}
		vector->value[count++] = value;
	}
	vector->count = count;
}

/* Slot mappings are unread during sorting and are rebuilt below. Borrow that
 * existing integer array as radix scratch: no allocation and 256 counters of
 * bounded stack storage, independent of the ambient vector dimension. */
static void sparse_sort_indices(struct simplex_SparseVector *vector)
{
	int i, ordered = 1;
	unsigned int maximum = 0, shift;
	int *source = vector->index, *target = vector->slot;
	for (i = 0; i < vector->count; i++) {
		unsigned int value = (unsigned int)source[i];
		if (value > maximum)
			maximum = value;
		if (i > 0 && source[i - 1] > source[i])
			ordered = 0;
	}
	if (ordered)
		return;
	if (vector->count <= 32) {
		for (i = 1; i < vector->count; i++) {
			int j = i, value = source[i];
			while (j > 0 && source[j - 1] > value) {
				source[j] = source[j - 1];
				j--;
			}
			source[j] = value;
		}
		return;
	}
	for (shift = 0;
		shift < sizeof(unsigned int) * 8 && (maximum >> shift) != 0;
		shift += 8) {
		int bucket[256] = {0}, sum = 0, b;
		int *swap;
		for (i = 0; i < vector->count; i++)
			bucket[((unsigned int)source[i] >> shift) & 255U]++;
		for (b = 0; b < 256; b++) {
			int count = bucket[b];
			bucket[b] = sum;
			sum += count;
		}
		for (i = 0; i < vector->count; i++) {
			unsigned int b =
				((unsigned int)source[i] >> shift) & 255U;
			target[bucket[b]++] = source[i];
		}
		swap = source;
		source = target;
		target = swap;
	}
	if (source != vector->index)
		lp_simplex_memcpy(vector->index, source,
			(size_t)vector->count * sizeof(int));
}

void simplex_sparse_vector_sort(struct simplex_SparseVector *vector)
{
	int i;
	if (vector->count < 2)
		return;
	sparse_sort_indices(vector);
	for (i = 0; i < vector->count; i++) {
		int row = vector->index[i];
		vector->value[i] = vector->dense[row];
		vector->slot[row] = i;
	}
}

int simplex_sparse_vector_set(struct simplex_SparseVector *vector,
	const int index, const double value)
{
	if (vector == NULL || index < 0 || index >= vector->dimension)
		return lp_simplex_EXIT_FAILURE;

	if (vector->mark[index] != vector->generation) {
		if (value == 0.)
			return lp_simplex_EXIT_SUCCESS;
		if (vector->count >= vector->capacity)
			return lp_simplex_EXIT_FAILURE;
		vector->mark[index] = vector->generation;
		vector->slot[index] = vector->count;
		vector->index[vector->count] = index;
		vector->value[vector->count] = value;
		vector->dense[index] = value;
		vector->count++;
	} else {
		int k = vector->slot[index];

		vector->dense[index] = value;
		vector->value[k] = value;
	}
	return lp_simplex_EXIT_SUCCESS;
}

int simplex_sparse_vector_add(struct simplex_SparseVector *vector,
	const int index, const double value)
{
	if (vector == NULL || index < 0 || index >= vector->dimension)
		return lp_simplex_EXIT_FAILURE;

	if (vector->mark[index] != vector->generation)
		return simplex_sparse_vector_set(vector, index, value);
	return simplex_sparse_vector_set(
		vector, index, vector->dense[index] + value);
}

static int sparse_vector_pack(struct simplex_SparseVector *vector,
	const double *dense, const double tolerance, double *maximum)
{
	int i;

	if (vector == NULL || dense == NULL || tolerance < 0.)
		return lp_simplex_EXIT_FAILURE;
	if (maximum != NULL)
		*maximum = 0.;
	simplex_sparse_vector_clear(vector);
	for (i = 0; i < vector->dimension; i++) {
		double value = dense[i];

		if (__lp_simplex_ABS__(value) <= tolerance)
			continue;
		vector->mark[i] = vector->generation;
		vector->slot[i] = vector->count;
		vector->index[vector->count] = i;
		vector->value[vector->count] = value;
		vector->dense[i] = value;
		vector->count++;
		if (maximum != NULL)
			*maximum = isfinite(value)
					   ? __lp_simplex_MAX__(*maximum,
						     __lp_simplex_ABS__(value))
					   : HUGE_VAL;
	}
	return lp_simplex_EXIT_SUCCESS;
}

int simplex_sparse_vector_pack(struct simplex_SparseVector *vector,
	const double *dense, const double tolerance)
{
	return sparse_vector_pack(vector, dense, tolerance, NULL);
}

int simplex_sparse_vector_pack_maximum(struct simplex_SparseVector *vector,
	const double *dense, const double tolerance, double *maximum)
{
	return sparse_vector_pack(vector, dense, tolerance, maximum);
}

double simplex_sparse_vector_get(
	const struct simplex_SparseVector *vector, const int index)
{
	if (vector == NULL || index < 0 || index >= vector->dimension ||
		vector->mark[index] != vector->generation)
		return 0.;

	return vector->dense[index];
}
