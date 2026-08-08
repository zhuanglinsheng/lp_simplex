#include "simplex_sparse_vector.h"
#include "utils.h"
#include <lp_simplex/status.h>


int simplex_sparse_vector_create(
		struct simplex_SparseVector *vector, const int dimension)
{
	if (vector == NULL || dimension < 0)
		return lp_simplex_EXIT_FAILURE;
	vector->dimension = dimension;
	vector->count = 0;
	vector->capacity = dimension;
	vector->index = dimension > 0 ? (int *)lp_simplex_malloc(
		(size_t)dimension * sizeof(int)) : NULL;
	vector->value = dimension > 0 ? (double *)lp_simplex_malloc(
		(size_t)dimension * sizeof(double)) : NULL;
	vector->dense = dimension > 0 ? (double *)lp_simplex_malloc(
		(size_t)dimension * sizeof(double)) : NULL;
	vector->mark = dimension > 0 ? (unsigned int *)lp_simplex_malloc(
		(size_t)dimension * sizeof(unsigned int)) : NULL;
	vector->slot = dimension > 0 ? (int *)lp_simplex_malloc(
		(size_t)dimension * sizeof(int)) : NULL;
	vector->generation = 1;
	if (dimension > 0 && (vector->index == NULL || vector->value == NULL ||
	    vector->dense == NULL || vector->mark == NULL || vector->slot == NULL)) {
		simplex_sparse_vector_destroy(vector);
		return lp_simplex_EXIT_FAILURE;
	}
	if (dimension > 0) {
		lp_simplex_memset(vector->dense, 0,
			(size_t)dimension * sizeof(double));
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
	lp_simplex_free(vector->dense);
	lp_simplex_free(vector->mark);
	lp_simplex_free(vector->slot);
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


int simplex_sparse_vector_set(
		struct simplex_SparseVector *vector, const int index,
		const double value)
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


int simplex_sparse_vector_add(
		struct simplex_SparseVector *vector, const int index,
		const double value)
{
	if (vector == NULL || index < 0 || index >= vector->dimension)
		return lp_simplex_EXIT_FAILURE;
	if (vector->mark[index] != vector->generation)
		return simplex_sparse_vector_set(vector, index, value);
	return simplex_sparse_vector_set(vector, index,
		vector->dense[index] + value);
}


int simplex_sparse_vector_pack(
		struct simplex_SparseVector *vector, const double *dense,
		const double tolerance)
{
	int i;
	if (vector == NULL || dense == NULL || tolerance < 0.)
		return lp_simplex_EXIT_FAILURE;
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
	}
	return lp_simplex_EXIT_SUCCESS;
}


double simplex_sparse_vector_get(
		const struct simplex_SparseVector *vector, const int index)
{
	if (vector == NULL || index < 0 || index >= vector->dimension ||
	    vector->mark[index] != vector->generation)
		return 0.;
	return vector->dense[index];
}
