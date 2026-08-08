/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#ifndef LP_SIMPLEX_SPARSE_LU_INTERNAL_H
#define LP_SIMPLEX_SPARSE_LU_INTERNAL_H

#include "simplex_csc.h"


struct simplex_SparseRow {
	int count;
	int capacity;
	int *column;
	double *value;
};

struct simplex_SparseColumnRows {
	int count;
	int capacity;
	int *row;
};

struct simplex_SparseLu {
	int dimension;
	int *permutation;
	int *column_permutation;
	double *column_scale;
	double *row_scale;
	int *work_column;
	int *pivot_row;
	int *diagonal_position;
	double *diagonal_value;
	double *work_value;
	double *solve_work;
	int *packed_start;
	int *packed_column;
	double *packed_value;
	int packed_nonzeros;
	struct simplex_SparseRow *row;
	struct simplex_SparseColumnRows *column_rows;
};

int simplex_sparse_lu_create(struct simplex_SparseLu *factor, int dimension);

void simplex_sparse_lu_destroy(struct simplex_SparseLu *factor);

int simplex_sparse_lu_factorize(
		struct simplex_SparseLu *factor,
		const struct simplex_CscMatrix *matrix,
		int structural_columns, const int *basis);

int simplex_sparse_lu_factorize_submatrix(
		struct simplex_SparseLu *factor,
		const struct simplex_CscMatrix *matrix,
		const int *columns, const int *row_to_core);

int simplex_sparse_lu_solve(
		struct simplex_SparseLu *factor, double *vector, int transpose);

int simplex_sparse_lu_solve_pair(
		struct simplex_SparseLu *factor, double *first, double *second);

#endif
