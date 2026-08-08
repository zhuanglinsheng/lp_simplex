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
int simplex_sparse_vector_set(
		struct simplex_SparseVector *vector, int index, double value);
int simplex_sparse_vector_add(
		struct simplex_SparseVector *vector, int index, double value);
int simplex_sparse_vector_pack(
		struct simplex_SparseVector *vector, const double *dense,
		double tolerance);
double simplex_sparse_vector_get(
		const struct simplex_SparseVector *vector, int index);

#endif
