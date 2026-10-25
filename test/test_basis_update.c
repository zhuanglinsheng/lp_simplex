/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file test_basis_update.c
 * @brief Regression tests for basis factorization and update consistency.
 */

#include "simplex_basis.h"
#include "simplex_csc.h"
#include "simplex_sparse_lu.h"

#include <lp_simplex/status.h>

#include <stdio.h>
#include <math.h>

static int test_forrest_tomlin_sequence(const int interleaved)
{
	/* Dense columns make the expected matrix-vector residual independent of
	 * sparsity assumptions and exercise row/column pivot permutations. */
	int column_start[7] = {0, 4, 8, 12, 16, 20, 24};
	int row_index[24] = {0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0,
		1, 2, 3, 0, 1, 2, 3};
	double value[24] = {4., 2., 1., 0.5, 1., 5., 2., 1., 2., 1., 6., 2., 1.,
		2., 1., 7., 3., -1., 2., 4., -2., 3., 1., 5.};
	int basis_index[4] = {0, 1, 2, 3};
	double right[4] = {2., -1., 3., 4.};
	double direction[4];
	struct simplex_CscMatrix matrix;
	struct simplex_SparseLu factor;
	int update, i, k;

	matrix.rows = 4;
	matrix.columns = 6;
	matrix.nonzeros = 24;
	matrix.column_start = column_start;
	matrix.row_index = row_index;
	matrix.value = value;
	matrix.row_start = NULL;
	matrix.column_index = NULL;
	matrix.row_value = NULL;
	matrix.owns_storage = 0;
	if (simplex_sparse_lu_create(&factor, 4) == lp_simplex_EXIT_FAILURE)
		return 0;
	if (simplex_sparse_lu_factorize(&factor, &matrix, 6, basis_index) ==
		lp_simplex_EXIT_FAILURE) {
		simplex_sparse_lu_destroy(&factor);
		return 0;
	}
	for (i = 0; i < 4; i++)
		direction[i] = right[i];
	if (simplex_sparse_lu_solve(&factor, direction, 0) !=
		lp_simplex_EXIT_SUCCESS) {
		simplex_sparse_lu_destroy(&factor);
		return 0;
	}
	for (i = 0; i < 4; i++) {
		double residual = -right[i];

		for (k = 0; k < 4; k++)
			residual += value[column_start[basis_index[k]] + i] *
				    direction[k];
		if (fabs(residual) > 1e-9) {
			simplex_sparse_lu_destroy(&factor);
			return 0;
		}
	}
	for (update = 0; update < 2; update++) {
		int entering = 4 + update;
		int leaving = update == 0 ? 1 : 3;

		for (i = 0; i < 4; i++)
			direction[i] = value[column_start[entering] + i];
		if (simplex_sparse_lu_solve(&factor, direction, 0) !=
			lp_simplex_EXIT_SUCCESS) {
			simplex_sparse_lu_destroy(&factor);
			return 0;
		}
		/* Refinement and batched work may solve another RHS before update.
		 * The supplied entering direction must still define the new basis. */
		if (interleaved) {
			double auxiliary[4];

			for (i = 0; i < 4; i++)
				auxiliary[i] = right[i];
			if (simplex_sparse_lu_solve(&factor, auxiliary, 0) !=
				lp_simplex_EXIT_SUCCESS) {
				simplex_sparse_lu_destroy(&factor);
				return 0;
			}
		}
		if (simplex_sparse_lu_ft_update(&factor, leaving, direction) !=
			lp_simplex_EXIT_SUCCESS) {
			simplex_sparse_lu_destroy(&factor);
			return 0;
		}
		basis_index[leaving] = entering;
		for (i = 0; i < 4; i++)
			direction[i] = right[i];
		if (simplex_sparse_lu_solve(&factor, direction, 0) !=
			lp_simplex_EXIT_SUCCESS) {
			simplex_sparse_lu_destroy(&factor);
			return 0;
		}
		for (i = 0; i < 4; i++) {
			double residual = -right[i];

			for (k = 0; k < 4; k++)
				residual += value[column_start[basis_index[k]] +
						    i] *
					    direction[k];
			if (fabs(residual) > 1e-9) {
				simplex_sparse_lu_destroy(&factor);
				return 0;
			}
		}
		for (i = 0; i < 4; i++)
			direction[i] = right[i];
		if (simplex_sparse_lu_solve(&factor, direction, 1) !=
			lp_simplex_EXIT_SUCCESS) {
			simplex_sparse_lu_destroy(&factor);
			return 0;
		}
		for (k = 0; k < 4; k++) {
			double residual = -right[k];

			for (i = 0; i < 4; i++)
				residual += value[column_start[basis_index[k]] +
						    i] *
					    direction[i];
			if (fabs(residual) > 1e-9) {
				simplex_sparse_lu_destroy(&factor);
				return 0;
			}
		}
	}
	simplex_sparse_lu_destroy(&factor);
	return 1;
}

static int test_weak_update_requests_reinversion(const int backend)
{
	int column_start[2] = {0, 2};
	int row_index[2] = {0, 1};
	double value[2] = {1., 1.};
	int row_start[3] = {0, 1, 2};
	int column_index[2] = {0, 0};
	double row_value[2] = {1., 1.};
	int basis_index[2] = {1, 2};
	int packed_index[2] = {0, 1};
	double direction[2] = {1e-8, 1.};
	struct simplex_CscMatrix matrix;
	struct simplex_Basis basis;
	int result;

	matrix.rows = 2;
	matrix.columns = 1;
	matrix.nonzeros = 2;
	matrix.column_start = column_start;
	matrix.row_index = row_index;
	matrix.value = value;
	matrix.row_start = row_start;
	matrix.column_index = column_index;
	matrix.row_value = row_value;
	matrix.owns_storage = 0;
	if (simplex_basis_create(&basis, &matrix, 1, basis_index) ==
		lp_simplex_EXIT_FAILURE)
		return 0;
	simplex_basis_set_backend(&basis, backend);
	if (simplex_basis_factorize(&basis) == lp_simplex_EXIT_FAILURE) {
		simplex_basis_destroy(&basis);
		return 0;
	}
	/* update() is called after the basis mapping has atomically exchanged. */
	basis_index[0] = 0;
	result = simplex_basis_update_packed(
		&basis, 0, direction, packed_index, 2);
	if (result == 1)
		result = simplex_basis_factorize(&basis);
	result = result == lp_simplex_EXIT_SUCCESS &&
		 simplex_basis_update_count(&basis) == 0;
	simplex_basis_destroy(&basis);
	return result;
}

static int test_current_basis_refinement(const int backend, const int product)
{
	int column_start[2] = {0, 2}, row_index[2] = {0, 1};
	int row_start[3] = {0, 1, 2}, column_index[2] = {0, 0};
	int basis_index[2] = {1, 2};
	double value[2] = {1., 1.}, direction[2] = {-1., -1.};
	double transpose[2] = {5., -3.}, forward[2] = {2., -1.};
	struct simplex_CscMatrix matrix;
	struct simplex_Basis basis;
	int update, valid;

	matrix.rows = 2;
	matrix.columns = 1;
	matrix.nonzeros = 2;
	matrix.column_start = column_start;
	matrix.row_index = row_index;
	matrix.value = value;
	matrix.row_start = row_start;
	matrix.column_index = column_index;
	matrix.row_value = value;
	matrix.owns_storage = 0;
	if (simplex_basis_create(&basis, &matrix, 1, basis_index) !=
		lp_simplex_EXIT_SUCCESS)
		return 0;
	simplex_basis_set_backend(&basis, backend);
	valid = simplex_basis_factorize(&basis) == lp_simplex_EXIT_SUCCESS;
	if (product)
		simplex_basis_use_product_updates(&basis);
	basis_index[0] = 0;
	update = valid ? simplex_basis_update(&basis, 0, direction) : -1;
	if (update == 1)
		update = simplex_basis_factorize(&basis);
	valid = valid && update == lp_simplex_EXIT_SUCCESS &&
		simplex_basis_btran_refined(&basis, transpose) ==
			lp_simplex_EXIT_SUCCESS &&
		simplex_basis_ftran_refined(&basis, forward) ==
			lp_simplex_EXIT_SUCCESS &&
		fabs(transpose[0] - 2.) < 1e-12 &&
		fabs(transpose[1] - 3.) < 1e-12 &&
		fabs(forward[0] - 2.) < 1e-12 && fabs(forward[1] - 3.) < 1e-12;
	forward[0] = transpose[0] = HUGE_VAL;
	forward[1] = transpose[1] = 0.;
	valid = valid &&
		simplex_basis_ftran_refined(&basis, forward) ==
			lp_simplex_EXIT_FAILURE &&
		simplex_basis_btran_refined(&basis, transpose) ==
			lp_simplex_EXIT_FAILURE;
	simplex_basis_destroy(&basis);
	return valid;
}

int main(void)
{
	if (!test_forrest_tomlin_sequence(0) ||
		!test_forrest_tomlin_sequence(1) ||
		!test_current_basis_refinement(
			SIMPLEX_BASIS_BACKEND_DEFAULT, 0) ||
		!test_current_basis_refinement(
			SIMPLEX_BASIS_BACKEND_INTERNAL, 0) ||
		!test_current_basis_refinement(
			SIMPLEX_BASIS_BACKEND_INTERNAL, 1)) {
		printf("Forrest-Tomlin sequence regression failed\n");
		return 1;
	}
	if (!test_weak_update_requests_reinversion(
		    SIMPLEX_BASIS_BACKEND_DEFAULT) ||
		!test_weak_update_requests_reinversion(
			SIMPLEX_BASIS_BACKEND_INTERNAL)) {
		printf("basis update regression failed\n");
		return 1;
	}
	printf("basis update regression passed\n");
	return 0;
}
