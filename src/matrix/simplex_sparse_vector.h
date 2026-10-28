/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_sparse_vector.h
 * @brief Sparse-vector storage and indexed numerical operations.
 */

#ifndef LP_SIMPLEX_SPARSE_VECTOR_INTERNAL_H
#define LP_SIMPLEX_SPARSE_VECTOR_INTERNAL_H

/* Packed vector with an optional dense scatter and O(1) membership marks. */
struct simplex_SparseVector {
	int dimension;
	int count;
	int capacity;
	int *index;
	double *value;
	double *dense;
	unsigned int *mark;
	int *slot;
	unsigned int generation;
};

int simplex_sparse_vector_create(
	struct simplex_SparseVector *vector, int dimension);

void simplex_sparse_vector_destroy(struct simplex_SparseVector *vector);

void simplex_sparse_vector_clear(struct simplex_SparseVector *vector);
/* Restore increasing coordinate order without changing values or support. */
void simplex_sparse_vector_sort(struct simplex_SparseVector *vector);
/* Remove only exact zeros, retaining every nonzero and membership invariant. */
void simplex_sparse_vector_compact(struct simplex_SparseVector *vector);

int simplex_sparse_vector_set(
	struct simplex_SparseVector *vector, int index, double value);

int simplex_sparse_vector_add(
	struct simplex_SparseVector *vector, int index, double value);

int simplex_sparse_vector_pack(struct simplex_SparseVector *vector,
	const double *dense, double tolerance);

/* Pack and measure in one traversal. A nonfinite retained entry yields an
 * infinite maximum; it is still packed for the caller's error handling. */
int simplex_sparse_vector_pack_maximum(struct simplex_SparseVector *vector,
	const double *dense, double tolerance, double *maximum);

double simplex_sparse_vector_get(
	const struct simplex_SparseVector *vector, int index);

#endif
