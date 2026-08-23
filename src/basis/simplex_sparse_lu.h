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
	int *integer_storage;
	double *numeric_storage;
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
	long factor_work;
	struct simplex_SparseRow *row;
	struct simplex_SparseColumnRows *column_rows;
	/* Forrest--Tomlin representation of the mutable U factor.  ft_column is
	 * an eta file indexed by pivot row; ft_order is the order in which the
	 * column etas are applied for U^{-1}.  Each update moves one pivot to the
	 * end and appends a sparse row transformation. */
	struct simplex_SparseRow *ft_column;
	struct simplex_SparseRow *ft_row_eta;
	int *ft_order;
	int *ft_pivot_position;
	int *ft_row_pivot;
	double *ft_spike_cache;
	double *ft_btran_cache;
	int ft_spike_valid;
	int ft_btran_pivot;
	int ft_update_count;
	int ft_update_capacity;
	int ft_active;
	long ft_initial_nonzeros;
	long ft_nonzeros;
	long ft_row_nonzeros;
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

int simplex_sparse_lu_ft_update(
		struct simplex_SparseLu *factor,
		int leaving_position, const double *direction);

int simplex_sparse_lu_ft_active(const struct simplex_SparseLu *factor);

#endif
