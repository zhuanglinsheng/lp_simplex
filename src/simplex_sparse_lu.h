#ifndef LP_SIMPLEX_SPARSE_LU_INTERNAL_H
#define LP_SIMPLEX_SPARSE_LU_INTERNAL_H

#include "simplex_csc.h"

struct simplex_SparseRow {
	int count;
	int capacity;
	int *column;
	double *value;
};

struct simplex_SparseLu {
	int dimension;
	int *permutation;
	int *column_permutation;
	double *column_scale;
	double *row_scale;
	int *work_column;
	int *pivot_row;
	double *work_value;
	double *solve_work;
	struct simplex_SparseRow *row;
};

int simplex_sparse_lu_create(struct simplex_SparseLu *factor, int dimension);
void simplex_sparse_lu_destroy(struct simplex_SparseLu *factor);
int simplex_sparse_lu_factorize(
		struct simplex_SparseLu *factor,
		const struct simplex_CscMatrix *matrix,
		int structural_columns, const int *basis);
int simplex_sparse_lu_solve(
		struct simplex_SparseLu *factor, double *vector, int transpose);

#endif
