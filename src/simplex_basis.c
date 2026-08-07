/* Basis factorization backend.  Revised-simplex code depends only on this API. */
#include "simplex_basis.h"
#include "linalg.h"
#include "utils.h"
#include <lp_simplex/status.h>
#include <stdlib.h>
#include <time.h>

#define SIMPLEX_BASIS_UPDATE_LIMIT 128


int simplex_basis_create(
		struct simplex_Basis *basis,
		const struct simplex_CscMatrix *matrix,
		const int structural_columns, int *index)
{
	basis->rows = matrix->rows;
	basis->structural_columns = structural_columns;
	basis->matrix = matrix;
	basis->index = index;
#ifdef LP_SIMPLEX_HAVE_KLU
	basis->base_column_start = (int *)lp_simplex_malloc(
		(size_t)(matrix->rows + 1) * sizeof(int));
	basis->base_row_index = NULL;
	basis->base_value = NULL;
	basis->symbolic = NULL;
	basis->numeric = NULL;
	klu_defaults(&basis->common);
#else
	lp_simplex_memset(&basis->sparse, 0, sizeof(basis->sparse));
#endif
	basis->eta = (double *)lp_simplex_malloc(
		(size_t)SIMPLEX_BASIS_UPDATE_LIMIT * matrix->rows * sizeof(double));
	basis->eta_pivot = (int *)lp_simplex_malloc(
		(size_t)SIMPLEX_BASIS_UPDATE_LIMIT * sizeof(int));
	basis->update_count = 0;
	basis->update_limit = SIMPLEX_BASIS_UPDATE_LIMIT;
	basis->profile_enabled = getenv("LP_SIMPLEX_PROFILE") != NULL;
	basis->profile_factor_seconds = 0.;
	basis->profile_ftran_seconds = 0.;
	basis->profile_btran_seconds = 0.;
	basis->profile_factor_calls = 0;
	basis->profile_ftran_calls = 0;
	basis->profile_btran_calls = 0;
	if (
#ifdef LP_SIMPLEX_HAVE_KLU
	    basis->base_column_start == NULL ||
#endif
	    basis->eta == NULL || basis->eta_pivot == NULL) {
		simplex_basis_destroy(basis);
		return lp_simplex_EXIT_FAILURE;
	}
#ifndef LP_SIMPLEX_HAVE_KLU
	if (simplex_sparse_lu_create(&basis->sparse, matrix->rows) ==
	    lp_simplex_EXIT_FAILURE) {
		simplex_basis_destroy(basis);
		return lp_simplex_EXIT_FAILURE;
	}
#endif
	return lp_simplex_EXIT_SUCCESS;
}


void simplex_basis_destroy(struct simplex_Basis *basis)
{
	if (basis == NULL)
		return;
#ifdef LP_SIMPLEX_HAVE_KLU
	if (basis->numeric != NULL)
		klu_free_numeric(&basis->numeric, &basis->common);
	if (basis->symbolic != NULL)
		klu_free_symbolic(&basis->symbolic, &basis->common);
	lp_simplex_free(basis->base_column_start);
	lp_simplex_free(basis->base_row_index);
	lp_simplex_free(basis->base_value);
	basis->base_column_start = NULL;
	basis->base_row_index = NULL;
	basis->base_value = NULL;
#else
	simplex_sparse_lu_destroy(&basis->sparse);
#endif
	lp_simplex_free(basis->eta);
	lp_simplex_free(basis->eta_pivot);
	basis->eta = NULL;
	basis->eta_pivot = NULL;
	basis->update_count = 0;
}


static int simplex_basis_factorize_impl(struct simplex_Basis *basis)
{
#ifdef LP_SIMPLEX_HAVE_KLU
	int column, k;
	int m = basis->rows;
	int nonzeros = 0;
	if (basis->numeric != NULL)
		klu_free_numeric(&basis->numeric, &basis->common);
	if (basis->symbolic != NULL)
		klu_free_symbolic(&basis->symbolic, &basis->common);
	lp_simplex_free(basis->base_row_index);
	lp_simplex_free(basis->base_value);
	basis->base_row_index = NULL;
	basis->base_value = NULL;
	for (column = 0; column < m; column++) {
		int variable = basis->index[column];
		if (variable >= basis->structural_columns)
			nonzeros++;
		else
			nonzeros += basis->matrix->column_start[variable + 1] -
				basis->matrix->column_start[variable];
	}
	basis->base_row_index = (int *)lp_simplex_malloc(
		(size_t)nonzeros * sizeof(int));
	basis->base_value = (double *)lp_simplex_malloc(
		(size_t)nonzeros * sizeof(double));
	if (basis->base_row_index == NULL || basis->base_value == NULL)
		return lp_simplex_EXIT_FAILURE;
	nonzeros = 0;
	for (column = 0; column < m; column++) {
		int variable = basis->index[column];
		basis->base_column_start[column] = nonzeros;
		if (variable >= basis->structural_columns) {
			basis->base_row_index[nonzeros] =
				variable - basis->structural_columns;
			basis->base_value[nonzeros] = -1.;
			nonzeros++;
		} else {
			for (k = basis->matrix->column_start[variable];
			     k < basis->matrix->column_start[variable + 1]; k++) {
				basis->base_row_index[nonzeros] =
					basis->matrix->row_index[k];
				basis->base_value[nonzeros] = basis->matrix->value[k];
				nonzeros++;
			}
		}
	}
	basis->base_column_start[m] = nonzeros;
	basis->symbolic = klu_analyze(m, basis->base_column_start,
		basis->base_row_index, &basis->common);
	if (basis->symbolic == NULL)
		return lp_simplex_EXIT_FAILURE;
	basis->numeric = klu_factor(basis->base_column_start,
		basis->base_row_index, basis->base_value,
		basis->symbolic, &basis->common);
	if (basis->numeric == NULL)
		return lp_simplex_EXIT_FAILURE;
	basis->update_count = 0;
	return lp_simplex_EXIT_SUCCESS;
#else
	if (simplex_sparse_lu_factorize(&basis->sparse, basis->matrix,
		basis->structural_columns, basis->index) == lp_simplex_EXIT_SUCCESS) {
		basis->update_count = 0;
		return lp_simplex_EXIT_SUCCESS;
	}
	return lp_simplex_EXIT_FAILURE;
#endif
}


int simplex_basis_factorize(struct simplex_Basis *basis)
{
	clock_t started = 0;
	int result;
	if (basis->profile_enabled)
		started = clock();
	result = simplex_basis_factorize_impl(basis);
	if (basis->profile_enabled) {
		basis->profile_factor_seconds +=
			(double)(clock() - started) / (double)CLOCKS_PER_SEC;
		basis->profile_factor_calls++;
	}
	return result;
}


static int simplex_basis_solve(
		const struct simplex_Basis *basis, double *vector, const char trans)
{
#ifdef LP_SIMPLEX_HAVE_KLU
	int success;
	if (trans == 'N')
		success = klu_solve(basis->symbolic, basis->numeric,
			basis->rows, 1, vector, (klu_common *)&basis->common);
	else
		success = klu_tsolve(basis->symbolic, basis->numeric,
			basis->rows, 1, vector, (klu_common *)&basis->common);
	return success ? lp_simplex_EXIT_SUCCESS : lp_simplex_EXIT_FAILURE;
#else
	return simplex_sparse_lu_solve(
		&((struct simplex_Basis *)basis)->sparse, vector, trans == 'T');
#endif
}


static int simplex_basis_apply_eta(
		const struct simplex_Basis *basis, const int update, double *vector)
{
	int p = basis->eta_pivot[update];
	const double *direction = basis->eta + (size_t)update * basis->rows;
	double pivot = direction[p];
	double value;
	if (__lp_simplex_ABS__(pivot) <= 1e-14)
		return lp_simplex_EXIT_FAILURE;
	value = vector[p] / pivot;
	lp_simplex_linalg_daxpy(basis->rows, -value, (double *)direction, 1,
		vector, 1);
	vector[p] = value;
	return lp_simplex_EXIT_SUCCESS;
}


static int simplex_basis_apply_eta_transpose(
		const struct simplex_Basis *basis, const int update, double *vector)
{
	int p = basis->eta_pivot[update];
	const double *direction = basis->eta + (size_t)update * basis->rows;
	double pivot = direction[p];
	double pivot_value = vector[p];
	double value;
	if (__lp_simplex_ABS__(pivot) <= 1e-14)
		return lp_simplex_EXIT_FAILURE;
	value = pivot_value - lp_simplex_linalg_ddot(
		basis->rows, direction, 1, vector, 1) + pivot * pivot_value;
	vector[p] = value / pivot;
	return lp_simplex_EXIT_SUCCESS;
}


static int simplex_basis_ftran_raw(
		const struct simplex_Basis *basis, double *vector)
{
	int update;
	if (simplex_basis_solve(basis, vector, 'N') == lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	for (update = 0; update < basis->update_count; update++) {
		if (simplex_basis_apply_eta(basis, update, vector) ==
		    lp_simplex_EXIT_FAILURE)
			return lp_simplex_EXIT_FAILURE;
	}
	return lp_simplex_EXIT_SUCCESS;
}


static int simplex_basis_btran_raw(
		const struct simplex_Basis *basis, double *vector)
{
	int update;
	for (update = basis->update_count - 1; update >= 0; update--) {
		if (simplex_basis_apply_eta_transpose(basis, update, vector) ==
		    lp_simplex_EXIT_FAILURE)
			return lp_simplex_EXIT_FAILURE;
	}
	return simplex_basis_solve(basis, vector, 'T');
}


int simplex_basis_ftran(const struct simplex_Basis *basis, double *vector)
{
	struct simplex_Basis *mutable = (struct simplex_Basis *)basis;
	clock_t started = 0;
	int result;
	if (basis->profile_enabled)
		started = clock();
	result = simplex_basis_ftran_raw(basis, vector);
	if (basis->profile_enabled) {
		mutable->profile_ftran_seconds +=
			(double)(clock() - started) / (double)CLOCKS_PER_SEC;
		mutable->profile_ftran_calls++;
	}
	return result;
}


int simplex_basis_btran(const struct simplex_Basis *basis, double *vector)
{
	struct simplex_Basis *mutable = (struct simplex_Basis *)basis;
	clock_t started = 0;
	int result;
	if (basis->profile_enabled)
		started = clock();
	result = simplex_basis_btran_raw(basis, vector);
	if (basis->profile_enabled) {
		mutable->profile_btran_seconds +=
			(double)(clock() - started) / (double)CLOCKS_PER_SEC;
		mutable->profile_btran_calls++;
	}
	return result;
}


int simplex_basis_edge_weights(
		const struct simplex_Basis *basis, double *weights)
{
	int i, j;
	double *column;
	if (basis->update_count != 0)
		return lp_simplex_EXIT_FAILURE;
	column = (double *)lp_simplex_malloc(
		(size_t)basis->rows * sizeof(double));
	if (column == NULL)
		return lp_simplex_EXIT_FAILURE;
	for (j = 0; j < basis->rows; j++) {
		double weight = 0.;
		lp_simplex_memset(column, 0, (size_t)basis->rows * sizeof(double));
		column[j] = 1.;
		if (simplex_basis_btran(basis, column) == lp_simplex_EXIT_FAILURE) {
			lp_simplex_free(column);
			return lp_simplex_EXIT_FAILURE;
		}
		for (i = 0; i < basis->rows; i++)
			weight += column[i] * column[i];
		weights[j] = __lp_simplex_MAX__(weight, 1e-12);
	}
	lp_simplex_free(column);
	return lp_simplex_EXIT_SUCCESS;
}


int simplex_basis_update(
		struct simplex_Basis *basis, const int leaving_position,
		const double *direction)
{
	double *destination;
	if (basis->update_count >= basis->update_limit)
		return 1;
	if (__lp_simplex_ABS__(direction[leaving_position]) <= 1e-14)
		return lp_simplex_EXIT_FAILURE;
	destination = basis->eta + (size_t)basis->update_count * basis->rows;
	lp_simplex_memcpy(destination, direction,
		(size_t)basis->rows * sizeof(double));
	basis->eta_pivot[basis->update_count] = leaving_position;
	basis->update_count++;
	return lp_simplex_EXIT_SUCCESS;
}
