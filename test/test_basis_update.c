#include "simplex_basis.h"
#include "simplex_csc.h"
#include "simplex_sparse_lu.h"

#include <lp_simplex/status.h>

#include <stdio.h>
#include <math.h>


static int test_forrest_tomlin_sequence(void)
{
	/* Dense columns make the expected matrix-vector residual independent of
	 * sparsity assumptions and exercise row/column pivot permutations. */
	int column_start[7] = {0, 4, 8, 12, 16, 20, 24};
	int row_index[24] = {
		0,1,2,3, 0,1,2,3, 0,1,2,3,
		0,1,2,3, 0,1,2,3, 0,1,2,3};
	double value[24] = {
		4.,2.,1.,0.5, 1.,5.,2.,1., 2.,1.,6.,2., 1.,2.,1.,7.,
		3.,-1.,2.,4., -2.,3.,1.,5.};
	int basis_index[4] = {0,1,2,3};
	double right[4] = {2.,-1.,3.,4.};
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
	if (simplex_sparse_lu_solve(&factor, direction, 0) != lp_simplex_EXIT_SUCCESS) {
		simplex_sparse_lu_destroy(&factor);
		return 0;
	}
	for (i = 0; i < 4; i++) {
		double residual = -right[i];
		for (k = 0; k < 4; k++)
			residual += value[column_start[basis_index[k]] + i] * direction[k];
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
		    lp_simplex_EXIT_SUCCESS ||
		    simplex_sparse_lu_ft_update(&factor, leaving, direction) !=
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
				residual += value[column_start[basis_index[k]] + i] *
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
				residual += value[column_start[basis_index[k]] + i] *
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


static int test_weak_update_requests_reinversion(void)
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


int main(void)
{
	if (!test_forrest_tomlin_sequence()) {
		printf("Forrest-Tomlin sequence regression failed\n");
		return 1;
	}
	if (!test_weak_update_requests_reinversion()) {
		printf("basis update regression failed\n");
		return 1;
	}
	printf("basis update regression passed\n");
	return 0;
}
