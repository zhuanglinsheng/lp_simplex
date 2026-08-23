/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 *
 * Basis factorization backend.
 * Revised-simplex code depends only on this API.
 */
#include "simplex_basis.h"
#include "linalg.h"
#include "simplex_sparse_lu.h"
#include "utils.h"

#ifdef LP_SIMPLEX_HAVE_KLU
#include <klu.h>
#endif

#include <lp_simplex/status.h>

#include <float.h>
#include <stdlib.h>
#include <time.h>


#define SIMPLEX_BASIS_UPDATE_LIMIT 512
#define SIMPLEX_BASIS_STABILITY_PIVOT 1e-7


struct simplex_BasisImpl {
	int rows;
	int structural_columns;
	const struct simplex_CscMatrix *matrix;
	int *index;
	int core_size;
	int factor_size;
	int compact_active;
	int compact_ever_active;
	int compact_requested;
	int allow_sparse_eta;
	int *core_basis;
	int *core_position;
	int *core_row;
	int *row_to_core;
	int *logical_position;
	int *base_index;
	double *base_work;
	double *core_work;
	double *refine_work;
	double *correction_work;
	double *certify_right;
	double *certify_residual;
	double *certify_correction;
#ifdef LP_SIMPLEX_HAVE_KLU
	int *base_column_start;
	int *base_row_index;
	double *base_value;
	klu_symbolic *symbolic;
	klu_numeric *numeric;
	klu_common common;
#else
	struct simplex_SparseLu sparse;
#endif
	double *eta_value;
	double *dense_eta;
	double *eta_pivot_value;
	int *eta_start;
	int *eta_index;
	int *eta_pivot;
	int *eta_dense_slot;
	int eta_capacity;
	int dense_eta_capacity;
	int dense_eta_count;
	int update_count;
	int update_limit;
	long update_work;
	long factor_work;
	long eta_apply_work;
	double minimum_relative_pivot;
	int compact_ftran_validation_countdown;
	int compact_btran_validation_countdown;
	int profile_enabled;
	double profile_factor_seconds;
	double profile_ftran_seconds;
	double profile_btran_seconds;
	long profile_factor_calls;
	long profile_ftran_calls;
	long profile_btran_calls;
	long profile_compact_calls;
	int profile_compact_min;
	int profile_compact_max;
	long profile_eta_nonzeros;
	long profile_eta_slots;
	long profile_fill_reinversions;
	long profile_stability_reinversions;
	long profile_compact_ftran_validations;
	long profile_compact_btran_validations;
	long profile_compact_ftran_refinements;
	long profile_compact_btran_refinements;
};


int simplex_basis_create(
		struct simplex_Basis *basis,
		const struct simplex_CscMatrix *matrix,
		const int structural_columns, int *index)
{
	if (basis == NULL || matrix == NULL)
		return lp_simplex_EXIT_FAILURE;
	lp_simplex_memset(basis, 0, sizeof(*basis));
	basis->impl = (struct simplex_BasisImpl *)lp_simplex_malloc(
		sizeof(*basis->impl));
	if (basis->impl == NULL)
		return lp_simplex_EXIT_FAILURE;
	lp_simplex_memset(basis->impl, 0, sizeof(*basis->impl));
	basis->impl->rows = matrix->rows;
	basis->impl->structural_columns = structural_columns;
	basis->impl->matrix = matrix;
	basis->impl->index = index;
	basis->impl->core_size = 0;
	basis->impl->factor_size = 0;
	basis->impl->compact_active = 0;
	basis->impl->compact_ever_active = 0;
	basis->impl->compact_requested = 0;
	basis->impl->allow_sparse_eta = 1;
	basis->impl->core_basis = (int *)lp_simplex_malloc(
		(size_t)6 * matrix->rows * sizeof(int));
	basis->impl->core_position = basis->impl->core_basis != NULL
		? basis->impl->core_basis + matrix->rows : NULL;
	basis->impl->core_row = basis->impl->core_basis != NULL
		? basis->impl->core_basis + 2 * matrix->rows : NULL;
	basis->impl->row_to_core = basis->impl->core_basis != NULL
		? basis->impl->core_basis + 3 * matrix->rows : NULL;
	basis->impl->logical_position = basis->impl->core_basis != NULL
		? basis->impl->core_basis + 4 * matrix->rows : NULL;
	basis->impl->base_index = basis->impl->core_basis != NULL
		? basis->impl->core_basis + 5 * matrix->rows : NULL;
	basis->impl->base_work = (double *)lp_simplex_malloc(
		(size_t)7 * matrix->rows * sizeof(double));
	basis->impl->core_work = basis->impl->base_work != NULL
		? basis->impl->base_work + matrix->rows : NULL;
	basis->impl->refine_work = basis->impl->base_work != NULL
		? basis->impl->base_work + 2 * matrix->rows : NULL;
	basis->impl->correction_work = basis->impl->base_work != NULL
		? basis->impl->base_work + 3 * matrix->rows : NULL;
	basis->impl->certify_right = basis->impl->base_work != NULL
		? basis->impl->base_work + 4 * matrix->rows : NULL;
	basis->impl->certify_residual = basis->impl->base_work != NULL
		? basis->impl->base_work + 5 * matrix->rows : NULL;
	basis->impl->certify_correction = basis->impl->base_work != NULL
		? basis->impl->base_work + 6 * matrix->rows : NULL;
#ifdef LP_SIMPLEX_HAVE_KLU
	basis->impl->base_column_start = (int *)lp_simplex_malloc(
		(size_t)(matrix->rows + 1) * sizeof(int));
	basis->impl->base_row_index = NULL;
	basis->impl->base_value = NULL;
	basis->impl->symbolic = NULL;
	basis->impl->numeric = NULL;
	klu_defaults(&basis->impl->common);
#else
	lp_simplex_memset(&basis->impl->sparse, 0, sizeof(basis->impl->sparse));
#endif
	basis->impl->eta_capacity = __lp_simplex_MAX__(matrix->rows * 4, 16);
	basis->impl->eta_value = (double *)lp_simplex_malloc(
		(size_t)basis->impl->eta_capacity * sizeof(double));
	/* Dense eta columns are grown by actual use.  Reserving update_limit * rows
	 * here makes the first dense update consume hundreds of megabytes on a
	 * large model even if reinversion follows immediately. */
	basis->impl->dense_eta = NULL;
	basis->impl->dense_eta_capacity = 0;
	basis->impl->dense_eta_count = 0;
	basis->impl->eta_index = (int *)lp_simplex_malloc(
		(size_t)basis->impl->eta_capacity * sizeof(int));
	basis->impl->eta_start = (int *)lp_simplex_malloc(
		(size_t)(3 * SIMPLEX_BASIS_UPDATE_LIMIT + 1) * sizeof(int));
	basis->impl->eta_pivot = basis->impl->eta_start != NULL
		? basis->impl->eta_start + SIMPLEX_BASIS_UPDATE_LIMIT + 1 : NULL;
	basis->impl->eta_pivot_value = (double *)lp_simplex_malloc(
		(size_t)SIMPLEX_BASIS_UPDATE_LIMIT * sizeof(double));
	basis->impl->eta_dense_slot = basis->impl->eta_pivot != NULL
		? basis->impl->eta_pivot + SIMPLEX_BASIS_UPDATE_LIMIT : NULL;
	basis->impl->update_count = 0;
	basis->impl->update_work = 0;
	basis->impl->factor_work = 0;
	basis->impl->eta_apply_work = 0;
	basis->impl->minimum_relative_pivot = 1.;
	basis->impl->compact_ftran_validation_countdown = 0;
	basis->impl->compact_btran_validation_countdown = 0;
	if (basis->impl->eta_start != NULL)
		basis->impl->eta_start[0] = 0;
	basis->impl->update_limit = SIMPLEX_BASIS_UPDATE_LIMIT;
	basis->impl->profile_enabled = getenv("LP_SIMPLEX_PROFILE") != NULL;
	basis->impl->profile_factor_seconds = 0.;
	basis->impl->profile_ftran_seconds = 0.;
	basis->impl->profile_btran_seconds = 0.;
	basis->impl->profile_factor_calls = 0;
	basis->impl->profile_ftran_calls = 0;
	basis->impl->profile_btran_calls = 0;
	basis->impl->profile_compact_calls = 0;
	basis->impl->profile_compact_min = matrix->rows;
	basis->impl->profile_compact_max = 0;
	basis->impl->profile_eta_nonzeros = 0;
	basis->impl->profile_eta_slots = 0;
	basis->impl->profile_fill_reinversions = 0;
	basis->impl->profile_stability_reinversions = 0;
	basis->impl->profile_compact_ftran_validations = 0;
	basis->impl->profile_compact_btran_validations = 0;
	basis->impl->profile_compact_ftran_refinements = 0;
	basis->impl->profile_compact_btran_refinements = 0;
	if (
#ifdef LP_SIMPLEX_HAVE_KLU
	    basis->impl->base_column_start == NULL ||
#endif
	    basis->impl->core_basis == NULL || basis->impl->core_position == NULL ||
	    basis->impl->core_row == NULL || basis->impl->row_to_core == NULL ||
	    basis->impl->logical_position == NULL || basis->impl->base_work == NULL ||
	    basis->impl->base_index == NULL || basis->impl->core_work == NULL ||
	    basis->impl->refine_work == NULL ||
	    basis->impl->correction_work == NULL ||
	    basis->impl->certify_right == NULL ||
	    basis->impl->certify_residual == NULL ||
	    basis->impl->certify_correction == NULL ||
	    basis->impl->eta_value == NULL ||
	    basis->impl->eta_index == NULL || basis->impl->eta_start == NULL ||
	    basis->impl->eta_pivot == NULL || basis->impl->eta_pivot_value == NULL ||
	    basis->impl->eta_dense_slot == NULL) {
		simplex_basis_destroy(basis);
		return lp_simplex_EXIT_FAILURE;
	}
#ifndef LP_SIMPLEX_HAVE_KLU
	/* The compact factor is allocated after the first basis is known. */
#endif
	return lp_simplex_EXIT_SUCCESS;
}


void simplex_basis_destroy(struct simplex_Basis *basis)
{
	if (basis == NULL || basis->impl == NULL)
		return;
#ifdef LP_SIMPLEX_HAVE_KLU
	if (basis->impl->numeric != NULL)
		klu_free_numeric(&basis->impl->numeric, &basis->impl->common);
	if (basis->impl->symbolic != NULL)
		klu_free_symbolic(&basis->impl->symbolic, &basis->impl->common);
	lp_simplex_free(basis->impl->base_column_start);
	lp_simplex_free(basis->impl->base_row_index);
	lp_simplex_free(basis->impl->base_value);
	basis->impl->base_column_start = NULL;
	basis->impl->base_row_index = NULL;
	basis->impl->base_value = NULL;
#else
	simplex_sparse_lu_destroy(&basis->impl->sparse);
#endif
	lp_simplex_free(basis->impl->eta_value);
	lp_simplex_free(basis->impl->dense_eta);
	lp_simplex_free(basis->impl->eta_index);
	lp_simplex_free(basis->impl->eta_start);
	lp_simplex_free(basis->impl->eta_pivot_value);
	lp_simplex_free(basis->impl->core_basis);
	lp_simplex_free(basis->impl->base_work);
	basis->impl->eta_value = NULL;
	basis->impl->dense_eta = NULL;
	basis->impl->eta_index = NULL;
	basis->impl->eta_start = NULL;
	basis->impl->eta_pivot = NULL;
	basis->impl->eta_pivot_value = NULL;
	basis->impl->eta_dense_slot = NULL;
	basis->impl->eta_capacity = 0;
	basis->impl->dense_eta_capacity = 0;
	basis->impl->dense_eta_count = 0;
	basis->impl->core_basis = NULL;
	basis->impl->core_position = NULL;
	basis->impl->core_row = NULL;
	basis->impl->row_to_core = NULL;
	basis->impl->logical_position = NULL;
	basis->impl->base_index = NULL;
	basis->impl->base_work = NULL;
	basis->impl->core_work = NULL;
	basis->impl->refine_work = NULL;
	basis->impl->correction_work = NULL;
	basis->impl->certify_right = NULL;
	basis->impl->certify_residual = NULL;
	basis->impl->certify_correction = NULL;
	basis->impl->core_size = 0;
	basis->impl->factor_size = 0;
	basis->impl->compact_active = 0;
	basis->impl->compact_ever_active = 0;
	basis->impl->compact_requested = 0;
	basis->impl->update_count = 0;
	basis->impl->update_work = 0;
	basis->impl->factor_work = 0;
	basis->impl->eta_apply_work = 0;
	basis->impl->minimum_relative_pivot = 1.;
	basis->impl->compact_ftran_validation_countdown = 0;
	basis->impl->compact_btran_validation_countdown = 0;
	lp_simplex_free(basis->impl);
	basis->impl = NULL;
}


static int simplex_basis_build_core(struct simplex_Basis *basis)
{
	int position, row;
	int structural_count = 0;
	int core_row_count = 0;
	int m = basis->impl->rows;
	for (row = 0; row < m; row++) {
		basis->impl->row_to_core[row] = -1;
		basis->impl->logical_position[row] = -1;
	}
	for (position = 0; position < m; position++) {
		int variable = basis->impl->index[position];
		if (variable < basis->impl->structural_columns) {
			basis->impl->core_basis[structural_count] = variable;
			basis->impl->core_position[structural_count] = position;
			structural_count++;
		} else {
			row = variable - basis->impl->structural_columns;
			if (row < 0 || row >= m || basis->impl->logical_position[row] >= 0)
				return lp_simplex_EXIT_FAILURE;
			basis->impl->logical_position[row] = position;
		}
	}
	for (row = 0; row < m; row++) {
		if (basis->impl->logical_position[row] < 0) {
			basis->impl->core_row[core_row_count] = row;
			basis->impl->row_to_core[row] = core_row_count;
			core_row_count++;
		}
	}
	if (core_row_count != structural_count)
		return lp_simplex_EXIT_FAILURE;
	lp_simplex_memcpy(basis->impl->base_index, basis->impl->index,
		(size_t)m * sizeof(int));
	basis->impl->core_size = structural_count;
	return lp_simplex_EXIT_SUCCESS;
}


static void simplex_basis_reset_updates(struct simplex_Basis *basis)
{
	basis->impl->update_count = 0;
	basis->impl->update_work = 0;
	basis->impl->minimum_relative_pivot = 1.;
	basis->impl->eta_start[0] = 0;
	basis->impl->dense_eta_count = 0;
}


static int simplex_basis_factorize_impl(struct simplex_Basis *basis)
{
	int core_size;
#ifdef LP_SIMPLEX_HAVE_KLU
	int column, k;
	int nonzeros = 0;
	int factor_size;
#endif
	if (simplex_basis_build_core(basis) == lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	core_size = basis->impl->core_size;
	basis->impl->compact_active = basis->impl->compact_requested &&
		basis->impl->rows >= 4096 &&
		core_size < basis->impl->rows &&
		core_size * 3 <= basis->impl->rows;
	basis->impl->factor_size = basis->impl->compact_active ? core_size : basis->impl->rows;
	basis->impl->compact_ftran_validation_countdown = 0;
	basis->impl->compact_btran_validation_countdown = 0;
	basis->impl->update_limit = SIMPLEX_BASIS_UPDATE_LIMIT;
	if (basis->impl->compact_active) {
		basis->impl->compact_ever_active = 1;
		basis->impl->profile_compact_calls++;
		basis->impl->profile_compact_min = __lp_simplex_MIN__(
			basis->impl->profile_compact_min, core_size);
		basis->impl->profile_compact_max = __lp_simplex_MAX__(
			basis->impl->profile_compact_max, core_size);
	}
#ifdef LP_SIMPLEX_HAVE_KLU
	factor_size = basis->impl->factor_size;
	if (basis->impl->numeric != NULL)
		klu_free_numeric(&basis->impl->numeric, &basis->impl->common);
	if (basis->impl->symbolic != NULL)
		klu_free_symbolic(&basis->impl->symbolic, &basis->impl->common);
	lp_simplex_free(basis->impl->base_row_index);
	lp_simplex_free(basis->impl->base_value);
	basis->impl->base_row_index = NULL;
	basis->impl->base_value = NULL;
	if (factor_size == 0) {
		simplex_basis_reset_updates(basis);
		return lp_simplex_EXIT_SUCCESS;
	}
	for (column = 0; column < factor_size; column++) {
		int variable = basis->impl->compact_active ? basis->impl->core_basis[column]
			: basis->impl->base_index[column];
		if (!basis->impl->compact_active && variable >= basis->impl->structural_columns)
			nonzeros++;
		else
			for (k = basis->impl->matrix->column_start[variable];
			     k < basis->impl->matrix->column_start[variable + 1]; k++)
				if (!basis->impl->compact_active ||
				    basis->impl->row_to_core[basis->impl->matrix->row_index[k]] >= 0)
					nonzeros++;
	}
	basis->impl->base_row_index = (int *)lp_simplex_malloc(
		(size_t)nonzeros * sizeof(int));
	basis->impl->base_value = (double *)lp_simplex_malloc(
		(size_t)nonzeros * sizeof(double));
	if (basis->impl->base_row_index == NULL || basis->impl->base_value == NULL)
		return lp_simplex_EXIT_FAILURE;
	nonzeros = 0;
	for (column = 0; column < factor_size; column++) {
		int variable = basis->impl->compact_active ? basis->impl->core_basis[column]
			: basis->impl->base_index[column];
		basis->impl->base_column_start[column] = nonzeros;
		if (!basis->impl->compact_active && variable >= basis->impl->structural_columns) {
			basis->impl->base_row_index[nonzeros] =
				variable - basis->impl->structural_columns;
			basis->impl->base_value[nonzeros++] = -1.;
		} else for (k = basis->impl->matrix->column_start[variable];
		     k < basis->impl->matrix->column_start[variable + 1]; k++) {
			int row = basis->impl->compact_active
				? basis->impl->row_to_core[basis->impl->matrix->row_index[k]]
				: basis->impl->matrix->row_index[k];
			if (row >= 0) {
				basis->impl->base_row_index[nonzeros] =
					row;
				basis->impl->base_value[nonzeros] = basis->impl->matrix->value[k];
				nonzeros++;
			}
		}
	}
	basis->impl->base_column_start[factor_size] = nonzeros;
	basis->impl->symbolic = klu_analyze(factor_size, basis->impl->base_column_start,
		basis->impl->base_row_index, &basis->impl->common);
	if (basis->impl->symbolic == NULL)
		return lp_simplex_EXIT_FAILURE;
	basis->impl->numeric = klu_factor(basis->impl->base_column_start,
		basis->impl->base_row_index, basis->impl->base_value,
		basis->impl->symbolic, &basis->impl->common);
	if (basis->impl->numeric == NULL)
		return lp_simplex_EXIT_FAILURE;
	klu_flops(basis->impl->symbolic, basis->impl->numeric, &basis->impl->common);
	simplex_basis_reset_updates(basis);
	return lp_simplex_EXIT_SUCCESS;
#else
	if (basis->impl->sparse.dimension != basis->impl->factor_size) {
		simplex_sparse_lu_destroy(&basis->impl->sparse);
		if (basis->impl->factor_size > 0 && simplex_sparse_lu_create(
				&basis->impl->sparse, basis->impl->factor_size) == lp_simplex_EXIT_FAILURE)
			return lp_simplex_EXIT_FAILURE;
	}
	if (basis->impl->factor_size == 0 ||
	    (basis->impl->compact_active && simplex_sparse_lu_factorize_submatrix(
		&basis->impl->sparse, basis->impl->matrix, basis->impl->core_basis,
		basis->impl->row_to_core) == lp_simplex_EXIT_SUCCESS) ||
	    (!basis->impl->compact_active && simplex_sparse_lu_factorize(
		&basis->impl->sparse, basis->impl->matrix, basis->impl->structural_columns,
		basis->impl->base_index) == lp_simplex_EXIT_SUCCESS)) {
		simplex_basis_reset_updates(basis);
		return lp_simplex_EXIT_SUCCESS;
	}
	return lp_simplex_EXIT_FAILURE;
#endif
}


int simplex_basis_factorize(struct simplex_Basis *basis)
{
	clock_t started = 0;
	double elapsed = 0.;
	int result;
	if (basis->impl->profile_enabled)
		started = clock();
	result = simplex_basis_factorize_impl(basis);
	if (result == lp_simplex_EXIT_SUCCESS) {
		if (basis->impl->profile_enabled)
			elapsed = (double)(clock() - started) /
				(double)CLOCKS_PER_SEC;
#ifdef LP_SIMPLEX_HAVE_KLU
		basis->impl->factor_work = (long)(basis->impl->common.flops +
			basis->impl->common.work + basis->impl->factor_size +
			(basis->impl->base_column_start != NULL
			 ? basis->impl->base_column_start[basis->impl->factor_size] : 0) +
			(basis->impl->numeric != NULL ? basis->impl->numeric->lnz +
			 basis->impl->numeric->unz + basis->impl->numeric->nzoff : 0));
#else
		basis->impl->factor_work = basis->impl->sparse.factor_work;
#endif
		basis->impl->eta_apply_work = 0;
	}
	if (basis->impl->profile_enabled) {
		basis->impl->profile_factor_seconds += elapsed;
		basis->impl->profile_factor_calls++;
	}
	return result;
}


static int simplex_basis_core_solve(
		const struct simplex_Basis *basis, double *vector, const char trans)
{
#ifdef LP_SIMPLEX_HAVE_KLU
	int success;
#endif
	if (basis->impl->factor_size == 0)
		return lp_simplex_EXIT_SUCCESS;
#ifdef LP_SIMPLEX_HAVE_KLU
	if (trans == 'N')
		success = klu_solve(basis->impl->symbolic, basis->impl->numeric,
			basis->impl->factor_size, 1, vector, (klu_common *)&basis->impl->common);
	else
		success = klu_tsolve(basis->impl->symbolic, basis->impl->numeric,
			basis->impl->factor_size, 1, vector, (klu_common *)&basis->impl->common);
	return success ? lp_simplex_EXIT_SUCCESS : lp_simplex_EXIT_FAILURE;
#else
	return simplex_sparse_lu_solve(
		&((struct simplex_Basis *)basis)->impl->sparse,
		vector, trans == 'T');
#endif
}


static int simplex_basis_core_solve_pair(
		const struct simplex_Basis *basis, double *first, double *second)
{
	if (basis->impl->factor_size == 0)
		return lp_simplex_EXIT_SUCCESS;
#ifdef LP_SIMPLEX_HAVE_KLU
	if (simplex_basis_core_solve(basis, first, 'N') ==
	    lp_simplex_EXIT_FAILURE ||
	    simplex_basis_core_solve(basis, second, 'N') ==
	    lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	return lp_simplex_EXIT_SUCCESS;
#else
	return simplex_sparse_lu_solve_pair(
		&((struct simplex_Basis *)basis)->impl->sparse, first, second);
#endif
}


static int simplex_basis_base_ftran_once(
		const struct simplex_Basis *basis, const double *right, double *vector)
{
	struct simplex_Basis *mutable = (struct simplex_Basis *)basis;
	const struct simplex_CscMatrix *matrix = basis->impl->matrix;
	int i, k, row, position;
	lp_simplex_memset(vector, 0, (size_t)basis->impl->rows * sizeof(double));
	for (i = 0; i < basis->impl->core_size; i++)
		mutable->impl->core_work[i] = right[basis->impl->core_row[i]];
	if (simplex_basis_core_solve(basis, mutable->impl->core_work, 'N') ==
	    lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	for (i = 0; i < basis->impl->core_size; i++)
		vector[basis->impl->core_position[i]] = mutable->impl->core_work[i];
	for (row = 0; row < basis->impl->rows; row++) {
		position = basis->impl->logical_position[row];
		if (position >= 0)
			vector[position] = -right[row];
	}
	for (i = 0; i < basis->impl->core_size; i++) {
		int column = basis->impl->core_basis[i];
		double value = mutable->impl->core_work[i];
		for (k = matrix->column_start[column];
		     k < matrix->column_start[column + 1]; k++) {
			row = matrix->row_index[k];
			position = basis->impl->logical_position[row];
			if (position >= 0)
				vector[position] += matrix->value[k] * value;
		}
	}
	return lp_simplex_EXIT_SUCCESS;
}


static int simplex_basis_base_ftran_pair_once(
		const struct simplex_Basis *basis,
		const double *first_right, const double *second_right,
		double *first, double *second)
{
	struct simplex_Basis *mutable = (struct simplex_Basis *)basis;
	const struct simplex_CscMatrix *matrix = basis->impl->matrix;
	int i, k, row, position;
	lp_simplex_memset(first, 0, (size_t)basis->impl->rows * sizeof(double));
	lp_simplex_memset(second, 0, (size_t)basis->impl->rows * sizeof(double));
	for (i = 0; i < basis->impl->core_size; i++) {
		row = basis->impl->core_row[i];
		mutable->impl->core_work[i] = first_right[row];
		mutable->impl->correction_work[i] = second_right[row];
	}
	if (simplex_basis_core_solve_pair(basis, mutable->impl->core_work,
			mutable->impl->correction_work) == lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	for (i = 0; i < basis->impl->core_size; i++) {
		position = basis->impl->core_position[i];
		first[position] = mutable->impl->core_work[i];
		second[position] = mutable->impl->correction_work[i];
	}
	for (row = 0; row < basis->impl->rows; row++) {
		position = basis->impl->logical_position[row];
		if (position >= 0) {
			first[position] = -first_right[row];
			second[position] = -second_right[row];
		}
	}
	for (i = 0; i < basis->impl->core_size; i++) {
		int column = basis->impl->core_basis[i];
		double first_value = mutable->impl->core_work[i];
		double second_value = mutable->impl->correction_work[i];
		for (k = matrix->column_start[column];
		     k < matrix->column_start[column + 1]; k++) {
			row = matrix->row_index[k];
			position = basis->impl->logical_position[row];
			if (position >= 0) {
				first[position] += matrix->value[k] * first_value;
				second[position] += matrix->value[k] * second_value;
			}
		}
	}
	return lp_simplex_EXIT_SUCCESS;
}


static void simplex_basis_ftran_residual(
		const struct simplex_Basis *basis, const double *right,
		const double *solution, double *residual, double *maximum)
{
	const struct simplex_CscMatrix *matrix = basis->impl->matrix;
	int i, k;
	double scale = 1.;
	for (i = 0; i < basis->impl->rows; i++) {
		residual[i] = right[i];
		scale = __lp_simplex_MAX__(scale, __lp_simplex_ABS__(right[i]));
	}
	for (i = 0; i < basis->impl->rows; i++) {
		int variable = basis->impl->base_index[i];
		double value = solution[i];
		if (variable >= basis->impl->structural_columns)
			residual[variable - basis->impl->structural_columns] += value;
		else
			for (k = matrix->column_start[variable];
			     k < matrix->column_start[variable + 1]; k++)
				residual[matrix->row_index[k]] -= matrix->value[k] * value;
	}
	*maximum = 0.;
	for (i = 0; i < basis->impl->rows; i++)
		*maximum = __lp_simplex_MAX__(*maximum,
			__lp_simplex_ABS__(residual[i]) / scale);
}


static int simplex_basis_base_ftran(
		const struct simplex_Basis *basis, double *vector)
{
	struct simplex_Basis *mutable = (struct simplex_Basis *)basis;
	int i, refinement, refined = 0;
	double residual;
	if (!basis->impl->compact_active)
		return simplex_basis_core_solve(basis, vector, 'N');
	lp_simplex_memcpy(mutable->impl->base_work, vector,
		(size_t)basis->impl->rows * sizeof(double));
	if (simplex_basis_base_ftran_once(basis, mutable->impl->base_work, vector) ==
	    lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	if (mutable->impl->compact_ftran_validation_countdown > 0) {
		mutable->impl->compact_ftran_validation_countdown--;
		return lp_simplex_EXIT_SUCCESS;
	}
	mutable->impl->profile_compact_ftran_validations++;
	for (refinement = 0; refinement < 2; refinement++) {
		simplex_basis_ftran_residual(basis, mutable->impl->base_work, vector,
			mutable->impl->refine_work, &residual);
		if (residual <= 1e-11)
			break;
		refined = 1;
		mutable->impl->profile_compact_ftran_refinements++;
		if (simplex_basis_base_ftran_once(basis, mutable->impl->refine_work,
				mutable->impl->correction_work) == lp_simplex_EXIT_FAILURE)
			return lp_simplex_EXIT_FAILURE;
		for (i = 0; i < basis->impl->rows; i++)
			vector[i] += mutable->impl->correction_work[i];
	}
	mutable->impl->compact_ftran_validation_countdown = refined ? 0 : 31;
	return lp_simplex_EXIT_SUCCESS;
}


static int simplex_basis_base_ftran_pair(
		const struct simplex_Basis *basis, double *first, double *second)
{
	struct simplex_Basis *mutable = (struct simplex_Basis *)basis;
	if (!basis->impl->compact_active)
		return simplex_basis_core_solve_pair(basis, first, second);
	/* Keep the certified path unchanged; the other 31 calls share mapping,
	 * sparse-LU traversal, and recovery for both right-hand sides. */
	if (mutable->impl->compact_ftran_validation_countdown <= 0) {
		if (simplex_basis_base_ftran(basis, first) == lp_simplex_EXIT_FAILURE ||
		    simplex_basis_base_ftran(basis, second) == lp_simplex_EXIT_FAILURE)
			return lp_simplex_EXIT_FAILURE;
		return lp_simplex_EXIT_SUCCESS;
	}
	lp_simplex_memcpy(mutable->impl->base_work, first,
		(size_t)basis->impl->rows * sizeof(double));
	lp_simplex_memcpy(mutable->impl->refine_work, second,
		(size_t)basis->impl->rows * sizeof(double));
	if (simplex_basis_base_ftran_pair_once(basis, mutable->impl->base_work,
			mutable->impl->refine_work, first, second) == lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	mutable->impl->compact_ftran_validation_countdown =
		__lp_simplex_MAX__(0,
			mutable->impl->compact_ftran_validation_countdown - 2);
	return lp_simplex_EXIT_SUCCESS;
}


static int simplex_basis_base_btran_once(
		const struct simplex_Basis *basis, const double *right, double *vector)
{
	struct simplex_Basis *mutable = (struct simplex_Basis *)basis;
	const struct simplex_CscMatrix *matrix = basis->impl->matrix;
	int i, k, row, position;
	lp_simplex_memset(vector, 0, (size_t)basis->impl->rows * sizeof(double));
	for (row = 0; row < basis->impl->rows; row++) {
		position = basis->impl->logical_position[row];
		if (position >= 0)
			vector[row] = -right[position];
	}
	for (i = 0; i < basis->impl->core_size; i++) {
		int column = basis->impl->core_basis[i];
		double value = right[basis->impl->core_position[i]];
		for (k = matrix->column_start[column];
		     k < matrix->column_start[column + 1]; k++) {
			row = matrix->row_index[k];
			if (basis->impl->logical_position[row] >= 0)
				value -= matrix->value[k] * vector[row];
		}
		mutable->impl->core_work[i] = value;
	}
	if (simplex_basis_core_solve(basis, mutable->impl->core_work, 'T') ==
	    lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	for (i = 0; i < basis->impl->core_size; i++)
		vector[basis->impl->core_row[i]] = mutable->impl->core_work[i];
	return lp_simplex_EXIT_SUCCESS;
}


static void simplex_basis_btran_residual(
		const struct simplex_Basis *basis, const double *right,
		const double *solution, double *residual, double *maximum)
{
	const struct simplex_CscMatrix *matrix = basis->impl->matrix;
	int i, k;
	double scale = 1.;
	for (i = 0; i < basis->impl->rows; i++) {
		int variable = basis->impl->base_index[i];
		double product = 0.;
		if (variable >= basis->impl->structural_columns)
			product = -solution[variable - basis->impl->structural_columns];
		else
			for (k = matrix->column_start[variable];
			     k < matrix->column_start[variable + 1]; k++)
				product += matrix->value[k] * solution[matrix->row_index[k]];
		residual[i] = right[i] - product;
		scale = __lp_simplex_MAX__(scale, __lp_simplex_ABS__(right[i]));
	}
	*maximum = 0.;
	for (i = 0; i < basis->impl->rows; i++)
		*maximum = __lp_simplex_MAX__(*maximum,
			__lp_simplex_ABS__(residual[i]) / scale);
}


static int simplex_basis_base_btran(
		const struct simplex_Basis *basis, double *vector)
{
	struct simplex_Basis *mutable = (struct simplex_Basis *)basis;
	int i, refinement, refined = 0;
	double residual;
	if (!basis->impl->compact_active)
		return simplex_basis_core_solve(basis, vector, 'T');
	lp_simplex_memcpy(mutable->impl->base_work, vector,
		(size_t)basis->impl->rows * sizeof(double));
	if (simplex_basis_base_btran_once(basis, mutable->impl->base_work, vector) ==
	    lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	if (mutable->impl->compact_btran_validation_countdown > 0) {
		mutable->impl->compact_btran_validation_countdown--;
		return lp_simplex_EXIT_SUCCESS;
	}
	mutable->impl->profile_compact_btran_validations++;
	for (refinement = 0; refinement < 2; refinement++) {
		simplex_basis_btran_residual(basis, mutable->impl->base_work, vector,
			mutable->impl->refine_work, &residual);
		if (residual <= 1e-11)
			break;
		refined = 1;
		mutable->impl->profile_compact_btran_refinements++;
		if (simplex_basis_base_btran_once(basis, mutable->impl->refine_work,
				mutable->impl->correction_work) == lp_simplex_EXIT_FAILURE)
			return lp_simplex_EXIT_FAILURE;
		for (i = 0; i < basis->impl->rows; i++)
			vector[i] += mutable->impl->correction_work[i];
	}
	mutable->impl->compact_btran_validation_countdown = refined ? 0 : 31;
	return lp_simplex_EXIT_SUCCESS;
}


static int simplex_basis_apply_eta(
		const struct simplex_Basis *basis, const int update, double *vector)
{
	int k;
	int dense_slot = basis->impl->eta_dense_slot[update];
	int p = basis->impl->eta_pivot[update];
	double pivot = basis->impl->eta_pivot_value[update];
	double value;
	if (__lp_simplex_ABS__(pivot) <= 1e-14)
		return lp_simplex_EXIT_FAILURE;
	if (vector[p] == 0.)
		return lp_simplex_EXIT_SUCCESS;
	value = vector[p] / pivot;
	if (dense_slot < 0)
		for (k = basis->impl->eta_start[update];
		     k < basis->impl->eta_start[update + 1]; k++)
			vector[basis->impl->eta_index[k]] -= value * basis->impl->eta_value[k];
	else
		lp_simplex_linalg_daxpy(basis->impl->rows, -value,
			basis->impl->dense_eta + (size_t)dense_slot * basis->impl->rows,
			1, vector, 1);
	vector[p] = value;
	return lp_simplex_EXIT_SUCCESS;
}


static int simplex_basis_apply_eta_pair(
		const struct simplex_Basis *basis, const int update,
		double *first, double *second)
{
	int k;
	int dense_slot = basis->impl->eta_dense_slot[update];
	int p = basis->impl->eta_pivot[update];
	double pivot = basis->impl->eta_pivot_value[update];
	double first_value, second_value;
	if (__lp_simplex_ABS__(pivot) <= 1e-14)
		return lp_simplex_EXIT_FAILURE;
	if (first[p] == 0. && second[p] == 0.)
		return lp_simplex_EXIT_SUCCESS;
	first_value = first[p] / pivot;
	second_value = second[p] / pivot;
	if (dense_slot < 0) {
		for (k = basis->impl->eta_start[update];
		     k < basis->impl->eta_start[update + 1]; k++) {
			int index = basis->impl->eta_index[k];
			double value = basis->impl->eta_value[k];
			first[index] -= first_value * value;
			second[index] -= second_value * value;
		}
	} else {
		const double *eta = basis->impl->dense_eta +
			(size_t)dense_slot * basis->impl->rows;
		for (k = 0; k < basis->impl->rows; k++) {
			first[k] -= first_value * eta[k];
			second[k] -= second_value * eta[k];
		}
	}
	first[p] = first_value;
	second[p] = second_value;
	return lp_simplex_EXIT_SUCCESS;
}


static int simplex_basis_apply_eta_transpose(
		const struct simplex_Basis *basis, const int update, double *vector)
{
	int k;
	int dense_slot = basis->impl->eta_dense_slot[update];
	int p = basis->impl->eta_pivot[update];
	double pivot = basis->impl->eta_pivot_value[update];
	double pivot_value = vector[p];
	double dot = 0.;
	if (__lp_simplex_ABS__(pivot) <= 1e-14)
		return lp_simplex_EXIT_FAILURE;
	if (dense_slot < 0)
		for (k = basis->impl->eta_start[update];
		     k < basis->impl->eta_start[update + 1]; k++)
			dot += basis->impl->eta_value[k] * vector[basis->impl->eta_index[k]];
	else
		dot = lp_simplex_linalg_ddot(basis->impl->rows,
			basis->impl->dense_eta +
			(size_t)dense_slot * basis->impl->rows, 1, vector, 1);
	vector[p] = (pivot_value - dot + pivot * pivot_value) / pivot;
	return lp_simplex_EXIT_SUCCESS;
}


static int simplex_basis_ftran_raw(
		const struct simplex_Basis *basis, double *vector)
{
	int update, first_update = 0;
	if (simplex_basis_base_ftran(basis, vector) == lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
#ifndef LP_SIMPLEX_HAVE_KLU
	if (!basis->impl->compact_active &&
	    simplex_sparse_lu_ft_active(&basis->impl->sparse))
		first_update = basis->impl->sparse.ft_update_count;
#endif
	for (update = first_update; update < basis->impl->update_count; update++) {
		((struct simplex_Basis *)basis)->impl->eta_apply_work +=
			basis->impl->eta_dense_slot[update] < 0
			? basis->impl->eta_start[update + 1] - basis->impl->eta_start[update]
			: basis->impl->rows;
		if (simplex_basis_apply_eta(basis, update, vector) ==
		    lp_simplex_EXIT_FAILURE)
			return lp_simplex_EXIT_FAILURE;
	}
	return lp_simplex_EXIT_SUCCESS;
}


static int simplex_basis_btran_raw(
		const struct simplex_Basis *basis, double *vector)
{
	int update, first_update = 0;
#ifndef LP_SIMPLEX_HAVE_KLU
	if (!basis->impl->compact_active &&
	    simplex_sparse_lu_ft_active(&basis->impl->sparse))
		first_update = basis->impl->sparse.ft_update_count;
#endif
	for (update = basis->impl->update_count - 1;
	     update >= first_update; update--) {
		((struct simplex_Basis *)basis)->impl->eta_apply_work +=
			basis->impl->eta_dense_slot[update] < 0
			? basis->impl->eta_start[update + 1] - basis->impl->eta_start[update]
			: basis->impl->rows;
		if (simplex_basis_apply_eta_transpose(basis, update, vector) ==
		    lp_simplex_EXIT_FAILURE)
			return lp_simplex_EXIT_FAILURE;
	}
	return simplex_basis_base_btran(basis, vector);
}


/* Componentwise backward error for B*x=right using the current basis, not
 * merely the basis captured by the last refactorization.  The denominator
 * |right|+|B||x| makes the certificate invariant to row/column magnitudes. */
static double simplex_basis_current_ftran_residual(
		const struct simplex_Basis *basis, const double *right,
		const double *solution, double *residual, double *denominator)
{
	const struct simplex_CscMatrix *matrix = basis->impl->matrix;
	double maximum = 0.;
	int i, k;
	for (i = 0; i < basis->impl->rows; i++) {
		residual[i] = right[i];
		denominator[i] = __lp_simplex_ABS__(right[i]);
	}
	for (i = 0; i < basis->impl->rows; i++) {
		int variable = basis->impl->index[i];
		double value = solution[i];
		if (variable >= basis->impl->structural_columns) {
			residual[variable - basis->impl->structural_columns] += value;
			denominator[variable - basis->impl->structural_columns] +=
				__lp_simplex_ABS__(value);
		} else
			for (k = matrix->column_start[variable];
			     k < matrix->column_start[variable + 1]; k++) {
				double term = matrix->value[k] * value;
				int row = matrix->row_index[k];
				residual[row] -= term;
				denominator[row] += __lp_simplex_ABS__(term);
			}
	}
	for (i = 0; i < basis->impl->rows; i++)
		maximum = __lp_simplex_MAX__(maximum,
			__lp_simplex_ABS__(residual[i]) /
			__lp_simplex_MAX__(denominator[i], DBL_MIN));
	return maximum;
}


static int simplex_basis_refine(
		const struct simplex_Basis *basis, double *vector)
{
	struct simplex_Basis *mutable = (struct simplex_Basis *)basis;
	double error, improved;
	double threshold = 16. * DBL_EPSILON *
		(double)(basis->impl->rows + 1);
	int i, refinement;
	lp_simplex_memcpy(mutable->impl->certify_right, vector,
		(size_t)basis->impl->rows * sizeof(double));
	if (simplex_basis_ftran_raw(basis, vector) == lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	error = simplex_basis_current_ftran_residual(basis,
			mutable->impl->certify_right, vector,
			mutable->impl->certify_residual, mutable->impl->core_work);
	for (refinement = 0; refinement < 3 && error > threshold; refinement++) {
		lp_simplex_memcpy(mutable->impl->certify_correction,
			mutable->impl->certify_residual,
			(size_t)basis->impl->rows * sizeof(double));
		if (simplex_basis_ftran_raw(
				basis, mutable->impl->certify_correction) ==
		    lp_simplex_EXIT_FAILURE)
			return lp_simplex_EXIT_FAILURE;
		for (i = 0; i < basis->impl->rows; i++)
			vector[i] += mutable->impl->certify_correction[i];
		improved = simplex_basis_current_ftran_residual(basis,
				mutable->impl->certify_right, vector,
				mutable->impl->certify_residual, mutable->impl->core_work);
		if (improved >= error) {
			for (i = 0; i < basis->impl->rows; i++)
				vector[i] -= mutable->impl->certify_correction[i];
			break;
		}
		error = improved;
	}
	return lp_simplex_EXIT_SUCCESS;
}


int simplex_basis_ftran_refined(
		const struct simplex_Basis *basis, double *vector)
{
	return simplex_basis_refine(basis, vector);
}


int simplex_basis_ftran_pair(
		const struct simplex_Basis *basis, double *first, double *second)
{
	struct simplex_Basis *mutable = (struct simplex_Basis *)basis;
	clock_t started = 0;
	int result = lp_simplex_EXIT_SUCCESS;
	int update, first_update = 0;
	if (basis->impl->profile_enabled)
		started = clock();
	if (simplex_basis_base_ftran_pair(basis, first, second) ==
		lp_simplex_EXIT_FAILURE) {
		result = lp_simplex_EXIT_FAILURE;
	}
	if (result == lp_simplex_EXIT_SUCCESS) {
#ifndef LP_SIMPLEX_HAVE_KLU
		if (!basis->impl->compact_active &&
		    simplex_sparse_lu_ft_active(&basis->impl->sparse)) {
			first_update = basis->impl->sparse.ft_update_count;
			for (update = first_update;
			     update < basis->impl->update_count; update++) {
				mutable->impl->eta_apply_work += 2L *
					(basis->impl->eta_dense_slot[update] < 0
					? basis->impl->eta_start[update + 1] -
						basis->impl->eta_start[update]
					: basis->impl->rows);
				if (simplex_basis_apply_eta_pair(
						basis, update, first, second) ==
				    lp_simplex_EXIT_FAILURE) {
					result = lp_simplex_EXIT_FAILURE;
					break;
				}
			}
		} else {
#endif
			for (update = 0; update < basis->impl->update_count; update++) {
				mutable->impl->eta_apply_work += 2L *
					(basis->impl->eta_dense_slot[update] < 0
					? basis->impl->eta_start[update + 1] -
						basis->impl->eta_start[update]
					: basis->impl->rows);
				if (simplex_basis_apply_eta_pair(
						basis, update, first, second) ==
				    lp_simplex_EXIT_FAILURE) {
					result = lp_simplex_EXIT_FAILURE;
					break;
				}
			}
#ifndef LP_SIMPLEX_HAVE_KLU
		}
#endif
	}
	if (basis->impl->profile_enabled) {
		mutable->impl->profile_ftran_seconds +=
			(double)(clock() - started) / (double)CLOCKS_PER_SEC;
		mutable->impl->profile_ftran_calls += 2;
	}
	return result;
}


int simplex_basis_ftran(const struct simplex_Basis *basis, double *vector)
{
	struct simplex_Basis *mutable = (struct simplex_Basis *)basis;
	clock_t started = 0;
	int result;
	if (basis->impl->profile_enabled)
		started = clock();
	result = simplex_basis_ftran_raw(basis, vector);
	if (basis->impl->profile_enabled) {
		mutable->impl->profile_ftran_seconds +=
			(double)(clock() - started) / (double)CLOCKS_PER_SEC;
		mutable->impl->profile_ftran_calls++;
	}
	return result;
}


int simplex_basis_btran(const struct simplex_Basis *basis, double *vector)
{
	struct simplex_Basis *mutable = (struct simplex_Basis *)basis;
	clock_t started = 0;
	int result;
	if (basis->impl->profile_enabled)
		started = clock();
	result = simplex_basis_btran_raw(basis, vector);
	if (basis->impl->profile_enabled) {
		mutable->impl->profile_btran_seconds +=
			(double)(clock() - started) / (double)CLOCKS_PER_SEC;
		mutable->impl->profile_btran_calls++;
	}
	return result;
}


static int simplex_basis_reserve_eta(
		struct simplex_Basis *basis, const int required)
{
	int capacity;
	double *value;
	int *index;
	if (required <= basis->impl->eta_capacity)
		return lp_simplex_EXIT_SUCCESS;
	capacity = basis->impl->eta_capacity;
	while (capacity < required) {
		int grown = capacity + capacity / 2;
		capacity = grown > capacity ? grown : required;
	}
	value = (double *)lp_simplex_realloc(basis->impl->eta_value,
		(size_t)capacity * sizeof(double));
	if (value == NULL)
		return lp_simplex_EXIT_FAILURE;
	basis->impl->eta_value = value;
	index = (int *)lp_simplex_realloc(basis->impl->eta_index,
		(size_t)capacity * sizeof(int));
	if (index == NULL)
		return lp_simplex_EXIT_FAILURE;
	basis->impl->eta_index = index;
	basis->impl->eta_capacity = capacity;
	return lp_simplex_EXIT_SUCCESS;
}


static int simplex_basis_reserve_dense_eta(
		struct simplex_Basis *basis, const int required_columns)
{
	int capacity;
	double *eta;
	if (required_columns <= basis->impl->dense_eta_capacity)
		return lp_simplex_EXIT_SUCCESS;
	capacity = basis->impl->dense_eta_capacity > 0
		? basis->impl->dense_eta_capacity : 1;
	while (capacity < required_columns) {
		int grown = capacity + capacity / 2 + 1;
		capacity = grown > capacity ? grown : required_columns;
	}
	eta = (double *)lp_simplex_realloc(basis->impl->dense_eta,
		(size_t)capacity * basis->impl->rows * sizeof(double));
	if (eta == NULL)
		return lp_simplex_EXIT_FAILURE;
	basis->impl->dense_eta = eta;
	basis->impl->dense_eta_capacity = capacity;
	return lp_simplex_EXIT_SUCCESS;
}


static int simplex_basis_update_impl(
		struct simplex_Basis *basis, const int leaving_position,
		const double *direction, const int *packed_index, int nonzeros)
{
	int i, count, first;
	int packed = packed_index != NULL;
	int sparse;
	double direction_maximum = 0.;
	double relative_pivot;
	/* On large bases an internal sparse refactorization is orders of magnitude
	 * dearer than another eta traversal.  Let the numerical hard limit govern
	 * those chains; smaller bases retain the deterministic ski-rental rule. */
	if (basis->impl->update_count >= basis->impl->update_limit)
		return 1;
	if (basis->impl->rows < 4096 && basis->impl->update_count > 0 &&
	    basis->impl->factor_work > 0 &&
	    basis->impl->eta_apply_work /
		(basis->impl->allow_sparse_eta ? 1 : 2) >=
		basis->impl->factor_work) {
		basis->impl->profile_fill_reinversions++;
		return 1;
	}
	if (__lp_simplex_ABS__(direction[leaving_position]) <= 1e-14)
		return lp_simplex_EXIT_FAILURE;
	first = basis->impl->eta_start[basis->impl->update_count];
	if (packed) {
		if (nonzeros <= 0 || nonzeros > basis->impl->rows)
			return lp_simplex_EXIT_FAILURE;
		for (i = 0; i < nonzeros; i++) {
			int row = packed_index[i];
			if (row < 0 || row >= basis->impl->rows || direction[row] == 0.)
				return lp_simplex_EXIT_FAILURE;
			direction_maximum = __lp_simplex_MAX__(direction_maximum,
				__lp_simplex_ABS__(direction[row]));
		}
	} else {
		nonzeros = 0;
		for (i = 0; i < basis->impl->rows; i++)
			if (direction[i] != 0.) {
				direction_maximum = __lp_simplex_MAX__(direction_maximum,
					__lp_simplex_ABS__(direction[i]));
				nonzeros++;
			}
	}
	relative_pivot = __lp_simplex_ABS__(direction[leaving_position]) /
		__lp_simplex_MAX__(direction_maximum, 1e-300);
	basis->impl->update_work += nonzeros;
	basis->impl->minimum_relative_pivot = __lp_simplex_MIN__(
		basis->impl->minimum_relative_pivot, relative_pivot);
	/* The basis mapping already describes the post-pivot basis when update() is
	 * called.  Refactorizing here therefore replaces, rather than repeats, a
	 * weak eta update and prevents its multiplier growth from contaminating all
	 * subsequent solves. */
	if (relative_pivot < SIMPLEX_BASIS_STABILITY_PIVOT) {
		basis->impl->profile_stability_reinversions++;
		return 1;
	}
#ifndef LP_SIMPLEX_HAVE_KLU
	if (!basis->impl->compact_active &&
	    simplex_sparse_lu_ft_active(&basis->impl->sparse) &&
	    basis->impl->update_count == basis->impl->sparse.ft_update_count &&
	    basis->impl->rows < 4096) {
		int result = simplex_sparse_lu_ft_update(
			&basis->impl->sparse, leaving_position, direction);
		if (result == lp_simplex_EXIT_SUCCESS) {
			basis->impl->update_count++;
			basis->impl->eta_start[basis->impl->update_count] =
				basis->impl->eta_start[basis->impl->update_count - 1];
			basis->impl->profile_eta_nonzeros += nonzeros;
			basis->impl->profile_eta_slots += basis->impl->rows;
		}
		return result;
	}
#endif
	sparse =
		basis->impl->compact_active ||
		(basis->impl->allow_sparse_eta && nonzeros * 8 <= basis->impl->rows);
	if (sparse) {
		if (simplex_basis_reserve_eta(basis, first + nonzeros) ==
		    lp_simplex_EXIT_FAILURE)
			return lp_simplex_EXIT_FAILURE;
		count = first;
		if (packed)
			for (i = 0; i < nonzeros; i++) {
				int row = packed_index[i];
				basis->impl->eta_index[count] = row;
				basis->impl->eta_value[count++] = direction[row];
			}
		else
			for (i = 0; i < basis->impl->rows; i++)
				if (direction[i] != 0.) {
					basis->impl->eta_index[count] = i;
					basis->impl->eta_value[count++] = direction[i];
				}
		basis->impl->eta_dense_slot[basis->impl->update_count] = -1;
		basis->impl->profile_eta_nonzeros += nonzeros;
		basis->impl->profile_eta_slots += basis->impl->rows;
	} else {
		int dense_slot = basis->impl->dense_eta_count;
		count = first;
		if (simplex_basis_reserve_dense_eta(basis, dense_slot + 1) ==
		    lp_simplex_EXIT_FAILURE)
			return lp_simplex_EXIT_FAILURE;
		lp_simplex_memcpy(basis->impl->dense_eta +
			(size_t)dense_slot * basis->impl->rows, direction,
			(size_t)basis->impl->rows * sizeof(double));
		basis->impl->eta_dense_slot[basis->impl->update_count] = dense_slot;
		basis->impl->dense_eta_count++;
	}
	basis->impl->eta_pivot[basis->impl->update_count] = leaving_position;
	basis->impl->eta_pivot_value[basis->impl->update_count] = direction[leaving_position];
	basis->impl->update_count++;
	basis->impl->eta_start[basis->impl->update_count] = count;
	return lp_simplex_EXIT_SUCCESS;
}


int simplex_basis_update(
		struct simplex_Basis *basis, const int leaving_position,
		const double *direction)
{
	return simplex_basis_update_impl(basis, leaving_position, direction,
		NULL, 0);
}


int simplex_basis_update_packed(
		struct simplex_Basis *basis, const int leaving_position,
		const double *direction, const int *index, const int nonzeros)
{
	return simplex_basis_update_impl(basis, leaving_position, direction,
		index, nonzeros);
}


void simplex_basis_set_sparse_eta(
		struct simplex_Basis *basis, const int allowed)
{
	basis->impl->allow_sparse_eta = allowed;
}


int simplex_basis_update_count(const struct simplex_Basis *basis)
{
	return basis->impl->update_count;
}


int simplex_basis_compact_active(const struct simplex_Basis *basis)
{
	return basis->impl->compact_active;
}


int simplex_basis_compact_ever_active(const struct simplex_Basis *basis)
{
	return basis->impl->compact_ever_active;
}


int simplex_basis_compact_requested(const struct simplex_Basis *basis)
{
	return basis->impl->compact_requested;
}


void simplex_basis_request_compact(struct simplex_Basis *basis)
{
	basis->impl->compact_requested = 1;
}


void simplex_basis_get_profile(
		const struct simplex_Basis *basis,
		struct simplex_BasisProfile *profile)
{
	lp_simplex_memset(profile, 0, sizeof(*profile));
	profile->factor_size = basis->impl->factor_size;
	profile->factor_seconds = basis->impl->profile_factor_seconds;
	profile->ftran_seconds = basis->impl->profile_ftran_seconds;
	profile->btran_seconds = basis->impl->profile_btran_seconds;
	profile->factor_calls = basis->impl->profile_factor_calls;
	profile->ftran_calls = basis->impl->profile_ftran_calls;
	profile->btran_calls = basis->impl->profile_btran_calls;
	profile->compact_calls = basis->impl->profile_compact_calls;
	profile->compact_min = basis->impl->profile_compact_min;
	profile->compact_max = basis->impl->profile_compact_max;
	profile->eta_nonzeros = basis->impl->profile_eta_nonzeros;
	profile->eta_slots = basis->impl->profile_eta_slots;
	profile->fill_reinversions = basis->impl->profile_fill_reinversions;
	profile->stability_reinversions =
		basis->impl->profile_stability_reinversions;
	profile->compact_ftran_validations =
		basis->impl->profile_compact_ftran_validations;
	profile->compact_btran_validations =
		basis->impl->profile_compact_btran_validations;
	profile->compact_ftran_refinements =
		basis->impl->profile_compact_ftran_refinements;
	profile->compact_btran_refinements =
		basis->impl->profile_compact_btran_refinements;
}
