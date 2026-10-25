/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_csc.h
 * @brief Sparse column and row views shared by presolve and simplex kernels.
 */

#ifndef LP_SIMPLEX_CSC_INTERNAL_H
#define LP_SIMPLEX_CSC_INTERNAL_H

#include <lp_simplex/model.h>

/* Ownership modes; value-only views borrow the immutable sparsity pattern
 * and own one contiguous value/row_value allocation rooted at value. */
#define SIMPLEX_CSC_OWNS_ALL 1
#define SIMPLEX_CSC_OWNS_VALUES 2

struct simplex_CscMatrix {
	int rows;
	int columns;
	int nonzeros;
	int *column_start;
	int *row_index;
	double *value;
	int *row_start;
	int *column_index;
	double *row_value;
	int owns_storage; /* 0 borrowed, or SIMPLEX_CSC_OWNS_* */
};

/* Initialize a fresh view. Complete model caches are borrowed; dense inputs
 * produce owned CSC/CSR arrays. Failure leaves no owned buffers. */
int simplex_csc_from_model(
	const struct lp_Model *model, struct simplex_CscMatrix *matrix);

void simplex_csc_destroy(struct simplex_CscMatrix *matrix);

double simplex_csc_column_dot(const struct simplex_CscMatrix *matrix,
	int column, const double *vector);

void simplex_csc_column_to_dense(const struct simplex_CscMatrix *matrix,
	int structural_columns, int variable, double *dense);

void simplex_csc_column_axpy(const struct simplex_CscMatrix *matrix,
	int structural_columns, int variable, double scale, double *dense);

#endif
