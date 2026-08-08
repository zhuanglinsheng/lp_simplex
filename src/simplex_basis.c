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
	basis->core_size = 0;
	basis->factor_size = 0;
	basis->compact_active = 0;
	basis->compact_ever_active = 0;
	basis->compact_requested = 0;
	basis->core_basis = (int *)lp_simplex_malloc(
		(size_t)matrix->rows * sizeof(int));
	basis->core_position = (int *)lp_simplex_malloc(
		(size_t)matrix->rows * sizeof(int));
	basis->core_row = (int *)lp_simplex_malloc(
		(size_t)matrix->rows * sizeof(int));
	basis->row_to_core = (int *)lp_simplex_malloc(
		(size_t)matrix->rows * sizeof(int));
	basis->logical_position = (int *)lp_simplex_malloc(
		(size_t)matrix->rows * sizeof(int));
	basis->base_index = (int *)lp_simplex_malloc(
		(size_t)matrix->rows * sizeof(int));
	basis->base_work = (double *)lp_simplex_malloc(
		(size_t)matrix->rows * sizeof(double));
	basis->core_work = (double *)lp_simplex_malloc(
		(size_t)matrix->rows * sizeof(double));
	basis->refine_work = (double *)lp_simplex_malloc(
		(size_t)matrix->rows * sizeof(double));
	basis->correction_work = (double *)lp_simplex_malloc(
		(size_t)matrix->rows * sizeof(double));
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
	basis->eta_value = (double *)lp_simplex_malloc(
		(size_t)SIMPLEX_BASIS_UPDATE_LIMIT * matrix->rows * sizeof(double));
	basis->eta = (double *)lp_simplex_malloc(
		(size_t)SIMPLEX_BASIS_UPDATE_LIMIT * matrix->rows * sizeof(double));
	basis->eta_index = (int *)lp_simplex_malloc(
		(size_t)SIMPLEX_BASIS_UPDATE_LIMIT * matrix->rows * sizeof(int));
	basis->eta_start = (int *)lp_simplex_malloc(
		(size_t)(SIMPLEX_BASIS_UPDATE_LIMIT + 1) * sizeof(int));
	basis->eta_pivot = (int *)lp_simplex_malloc(
		(size_t)SIMPLEX_BASIS_UPDATE_LIMIT * sizeof(int));
	basis->eta_pivot_value = (double *)lp_simplex_malloc(
		(size_t)SIMPLEX_BASIS_UPDATE_LIMIT * sizeof(double));
	basis->eta_sparse = (unsigned char *)lp_simplex_malloc(
		(size_t)SIMPLEX_BASIS_UPDATE_LIMIT * sizeof(unsigned char));
	basis->update_count = 0;
	basis->compact_ftran_validation_countdown = 0;
	basis->compact_btran_validation_countdown = 0;
	if (basis->eta_start != NULL)
		basis->eta_start[0] = 0;
	basis->update_limit = SIMPLEX_BASIS_UPDATE_LIMIT;
	basis->profile_enabled = getenv("LP_SIMPLEX_PROFILE") != NULL;
	basis->profile_factor_seconds = 0.;
	basis->profile_ftran_seconds = 0.;
	basis->profile_btran_seconds = 0.;
	basis->profile_factor_calls = 0;
	basis->profile_ftran_calls = 0;
	basis->profile_btran_calls = 0;
	basis->profile_compact_calls = 0;
	basis->profile_compact_min = matrix->rows;
	basis->profile_compact_max = 0;
	basis->profile_eta_nonzeros = 0;
	basis->profile_eta_slots = 0;
	basis->profile_compact_ftran_validations = 0;
	basis->profile_compact_btran_validations = 0;
	basis->profile_compact_ftran_refinements = 0;
	basis->profile_compact_btran_refinements = 0;
	if (
#ifdef LP_SIMPLEX_HAVE_KLU
	    basis->base_column_start == NULL ||
#endif
	    basis->core_basis == NULL || basis->core_position == NULL ||
	    basis->core_row == NULL || basis->row_to_core == NULL ||
	    basis->logical_position == NULL || basis->base_work == NULL ||
	    basis->base_index == NULL || basis->core_work == NULL ||
	    basis->refine_work == NULL ||
	    basis->correction_work == NULL || basis->eta_value == NULL ||
	    basis->eta == NULL ||
	    basis->eta_index == NULL || basis->eta_start == NULL ||
	    basis->eta_pivot == NULL || basis->eta_pivot_value == NULL ||
	    basis->eta_sparse == NULL) {
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
	lp_simplex_free(basis->eta_value);
	lp_simplex_free(basis->eta);
	lp_simplex_free(basis->eta_index);
	lp_simplex_free(basis->eta_start);
	lp_simplex_free(basis->eta_pivot);
	lp_simplex_free(basis->eta_pivot_value);
	lp_simplex_free(basis->eta_sparse);
	lp_simplex_free(basis->core_basis);
	lp_simplex_free(basis->core_position);
	lp_simplex_free(basis->core_row);
	lp_simplex_free(basis->row_to_core);
	lp_simplex_free(basis->logical_position);
	lp_simplex_free(basis->base_index);
	lp_simplex_free(basis->base_work);
	lp_simplex_free(basis->core_work);
	lp_simplex_free(basis->refine_work);
	lp_simplex_free(basis->correction_work);
	basis->eta_value = NULL;
	basis->eta = NULL;
	basis->eta_index = NULL;
	basis->eta_start = NULL;
	basis->eta_pivot = NULL;
	basis->eta_pivot_value = NULL;
	basis->eta_sparse = NULL;
	basis->core_basis = NULL;
	basis->core_position = NULL;
	basis->core_row = NULL;
	basis->row_to_core = NULL;
	basis->logical_position = NULL;
	basis->base_index = NULL;
	basis->base_work = NULL;
	basis->core_work = NULL;
	basis->refine_work = NULL;
	basis->correction_work = NULL;
	basis->core_size = 0;
	basis->factor_size = 0;
	basis->compact_active = 0;
	basis->compact_ever_active = 0;
	basis->compact_requested = 0;
	basis->update_count = 0;
	basis->compact_ftran_validation_countdown = 0;
	basis->compact_btran_validation_countdown = 0;
}


static int simplex_basis_build_core(struct simplex_Basis *basis)
{
	int position, row;
	int structural_count = 0;
	int core_row_count = 0;
	int m = basis->rows;
	for (row = 0; row < m; row++) {
		basis->row_to_core[row] = -1;
		basis->logical_position[row] = -1;
	}
	for (position = 0; position < m; position++) {
		int variable = basis->index[position];
		if (variable < basis->structural_columns) {
			basis->core_basis[structural_count] = variable;
			basis->core_position[structural_count] = position;
			structural_count++;
		} else {
			row = variable - basis->structural_columns;
			if (row < 0 || row >= m || basis->logical_position[row] >= 0)
				return lp_simplex_EXIT_FAILURE;
			basis->logical_position[row] = position;
		}
	}
	for (row = 0; row < m; row++) {
		if (basis->logical_position[row] < 0) {
			basis->core_row[core_row_count] = row;
			basis->row_to_core[row] = core_row_count;
			core_row_count++;
		}
	}
	if (core_row_count != structural_count)
		return lp_simplex_EXIT_FAILURE;
	lp_simplex_memcpy(basis->base_index, basis->index,
		(size_t)m * sizeof(int));
	basis->core_size = structural_count;
	return lp_simplex_EXIT_SUCCESS;
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
	core_size = basis->core_size;
	basis->compact_active = basis->compact_requested &&
		basis->rows >= 4096 &&
		core_size < basis->rows &&
		core_size * 3 <= basis->rows;
	basis->factor_size = basis->compact_active ? core_size : basis->rows;
	basis->compact_ftran_validation_countdown = 0;
	basis->compact_btran_validation_countdown = 0;
	basis->update_limit = SIMPLEX_BASIS_UPDATE_LIMIT;
	if (basis->compact_active) {
		basis->compact_ever_active = 1;
		basis->profile_compact_calls++;
		basis->profile_compact_min = __lp_simplex_MIN__(
			basis->profile_compact_min, core_size);
		basis->profile_compact_max = __lp_simplex_MAX__(
			basis->profile_compact_max, core_size);
	}
#ifdef LP_SIMPLEX_HAVE_KLU
	factor_size = basis->factor_size;
	if (basis->numeric != NULL)
		klu_free_numeric(&basis->numeric, &basis->common);
	if (basis->symbolic != NULL)
		klu_free_symbolic(&basis->symbolic, &basis->common);
	lp_simplex_free(basis->base_row_index);
	lp_simplex_free(basis->base_value);
	basis->base_row_index = NULL;
	basis->base_value = NULL;
	if (factor_size == 0) {
		basis->update_count = 0;
		basis->eta_start[0] = 0;
		return lp_simplex_EXIT_SUCCESS;
	}
	for (column = 0; column < factor_size; column++) {
		int variable = basis->compact_active ? basis->core_basis[column]
			: basis->base_index[column];
		if (!basis->compact_active && variable >= basis->structural_columns)
			nonzeros++;
		else
			for (k = basis->matrix->column_start[variable];
			     k < basis->matrix->column_start[variable + 1]; k++)
				if (!basis->compact_active ||
				    basis->row_to_core[basis->matrix->row_index[k]] >= 0)
					nonzeros++;
	}
	basis->base_row_index = (int *)lp_simplex_malloc(
		(size_t)nonzeros * sizeof(int));
	basis->base_value = (double *)lp_simplex_malloc(
		(size_t)nonzeros * sizeof(double));
	if (basis->base_row_index == NULL || basis->base_value == NULL)
		return lp_simplex_EXIT_FAILURE;
	nonzeros = 0;
	for (column = 0; column < factor_size; column++) {
		int variable = basis->compact_active ? basis->core_basis[column]
			: basis->base_index[column];
		basis->base_column_start[column] = nonzeros;
		if (!basis->compact_active && variable >= basis->structural_columns) {
			basis->base_row_index[nonzeros] =
				variable - basis->structural_columns;
			basis->base_value[nonzeros++] = -1.;
		} else for (k = basis->matrix->column_start[variable];
		     k < basis->matrix->column_start[variable + 1]; k++) {
			int row = basis->compact_active
				? basis->row_to_core[basis->matrix->row_index[k]]
				: basis->matrix->row_index[k];
			if (row >= 0) {
				basis->base_row_index[nonzeros] =
					row;
				basis->base_value[nonzeros] = basis->matrix->value[k];
				nonzeros++;
			}
		}
	}
	basis->base_column_start[factor_size] = nonzeros;
	basis->symbolic = klu_analyze(factor_size, basis->base_column_start,
		basis->base_row_index, &basis->common);
	if (basis->symbolic == NULL)
		return lp_simplex_EXIT_FAILURE;
	basis->numeric = klu_factor(basis->base_column_start,
		basis->base_row_index, basis->base_value,
		basis->symbolic, &basis->common);
	if (basis->numeric == NULL)
		return lp_simplex_EXIT_FAILURE;
	basis->update_count = 0;
	basis->eta_start[0] = 0;
	return lp_simplex_EXIT_SUCCESS;
#else
	if (basis->sparse.dimension != basis->factor_size) {
		simplex_sparse_lu_destroy(&basis->sparse);
		if (basis->factor_size > 0 && simplex_sparse_lu_create(
				&basis->sparse, basis->factor_size) == lp_simplex_EXIT_FAILURE)
			return lp_simplex_EXIT_FAILURE;
	}
	if (basis->factor_size == 0 ||
	    (basis->compact_active && simplex_sparse_lu_factorize_submatrix(
		&basis->sparse, basis->matrix, basis->core_basis,
		basis->row_to_core) == lp_simplex_EXIT_SUCCESS) ||
	    (!basis->compact_active && simplex_sparse_lu_factorize(
		&basis->sparse, basis->matrix, basis->structural_columns,
		basis->base_index) == lp_simplex_EXIT_SUCCESS)) {
		basis->update_count = 0;
		basis->eta_start[0] = 0;
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


static int simplex_basis_core_solve(
		const struct simplex_Basis *basis, double *vector, const char trans)
{
	if (basis->factor_size == 0)
		return lp_simplex_EXIT_SUCCESS;
#ifdef LP_SIMPLEX_HAVE_KLU
	int success;
	if (trans == 'N')
		success = klu_solve(basis->symbolic, basis->numeric,
			basis->factor_size, 1, vector, (klu_common *)&basis->common);
	else
		success = klu_tsolve(basis->symbolic, basis->numeric,
			basis->factor_size, 1, vector, (klu_common *)&basis->common);
	return success ? lp_simplex_EXIT_SUCCESS : lp_simplex_EXIT_FAILURE;
#else
	return simplex_sparse_lu_solve(
		&((struct simplex_Basis *)basis)->sparse, vector, trans == 'T');
#endif
}


static int simplex_basis_core_solve_pair(
		const struct simplex_Basis *basis, double *first, double *second)
{
	if (basis->factor_size == 0)
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
		&((struct simplex_Basis *)basis)->sparse, first, second);
#endif
}


static int simplex_basis_base_ftran_once(
		const struct simplex_Basis *basis, const double *right, double *vector)
{
	struct simplex_Basis *mutable = (struct simplex_Basis *)basis;
	const struct simplex_CscMatrix *matrix = basis->matrix;
	int i, k, row, position;
	lp_simplex_memset(vector, 0, (size_t)basis->rows * sizeof(double));
	for (i = 0; i < basis->core_size; i++)
		mutable->core_work[i] = right[basis->core_row[i]];
	if (simplex_basis_core_solve(basis, mutable->core_work, 'N') ==
	    lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	for (i = 0; i < basis->core_size; i++)
		vector[basis->core_position[i]] = mutable->core_work[i];
	for (row = 0; row < basis->rows; row++) {
		position = basis->logical_position[row];
		if (position >= 0)
			vector[position] = -right[row];
	}
	for (i = 0; i < basis->core_size; i++) {
		int column = basis->core_basis[i];
		double value = mutable->core_work[i];
		for (k = matrix->column_start[column];
		     k < matrix->column_start[column + 1]; k++) {
			row = matrix->row_index[k];
			position = basis->logical_position[row];
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
	const struct simplex_CscMatrix *matrix = basis->matrix;
	int i, k, row, position;
	lp_simplex_memset(first, 0, (size_t)basis->rows * sizeof(double));
	lp_simplex_memset(second, 0, (size_t)basis->rows * sizeof(double));
	for (i = 0; i < basis->core_size; i++) {
		row = basis->core_row[i];
		mutable->core_work[i] = first_right[row];
		mutable->correction_work[i] = second_right[row];
	}
	if (simplex_basis_core_solve_pair(basis, mutable->core_work,
			mutable->correction_work) == lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	for (i = 0; i < basis->core_size; i++) {
		position = basis->core_position[i];
		first[position] = mutable->core_work[i];
		second[position] = mutable->correction_work[i];
	}
	for (row = 0; row < basis->rows; row++) {
		position = basis->logical_position[row];
		if (position >= 0) {
			first[position] = -first_right[row];
			second[position] = -second_right[row];
		}
	}
	for (i = 0; i < basis->core_size; i++) {
		int column = basis->core_basis[i];
		double first_value = mutable->core_work[i];
		double second_value = mutable->correction_work[i];
		for (k = matrix->column_start[column];
		     k < matrix->column_start[column + 1]; k++) {
			row = matrix->row_index[k];
			position = basis->logical_position[row];
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
	const struct simplex_CscMatrix *matrix = basis->matrix;
	int i, k;
	double scale = 1.;
	for (i = 0; i < basis->rows; i++) {
		residual[i] = right[i];
		scale = __lp_simplex_MAX__(scale, __lp_simplex_ABS__(right[i]));
	}
	for (i = 0; i < basis->rows; i++) {
		int variable = basis->base_index[i];
		double value = solution[i];
		if (variable >= basis->structural_columns)
			residual[variable - basis->structural_columns] += value;
		else
			for (k = matrix->column_start[variable];
			     k < matrix->column_start[variable + 1]; k++)
				residual[matrix->row_index[k]] -= matrix->value[k] * value;
	}
	*maximum = 0.;
	for (i = 0; i < basis->rows; i++)
		*maximum = __lp_simplex_MAX__(*maximum,
			__lp_simplex_ABS__(residual[i]) / scale);
}


static int simplex_basis_base_ftran(
		const struct simplex_Basis *basis, double *vector)
{
	struct simplex_Basis *mutable = (struct simplex_Basis *)basis;
	int i, refinement, refined = 0;
	double residual;
	if (!basis->compact_active)
		return simplex_basis_core_solve(basis, vector, 'N');
	lp_simplex_memcpy(mutable->base_work, vector,
		(size_t)basis->rows * sizeof(double));
	if (simplex_basis_base_ftran_once(basis, mutable->base_work, vector) ==
	    lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	if (mutable->compact_ftran_validation_countdown > 0) {
		mutable->compact_ftran_validation_countdown--;
		return lp_simplex_EXIT_SUCCESS;
	}
	mutable->profile_compact_ftran_validations++;
	for (refinement = 0; refinement < 2; refinement++) {
		simplex_basis_ftran_residual(basis, mutable->base_work, vector,
			mutable->refine_work, &residual);
		if (residual <= 1e-11)
			break;
		refined = 1;
		mutable->profile_compact_ftran_refinements++;
		if (simplex_basis_base_ftran_once(basis, mutable->refine_work,
				mutable->correction_work) == lp_simplex_EXIT_FAILURE)
			return lp_simplex_EXIT_FAILURE;
		for (i = 0; i < basis->rows; i++)
			vector[i] += mutable->correction_work[i];
	}
	mutable->compact_ftran_validation_countdown = refined ? 0 : 31;
	return lp_simplex_EXIT_SUCCESS;
}


static int simplex_basis_base_ftran_pair(
		const struct simplex_Basis *basis, double *first, double *second)
{
	struct simplex_Basis *mutable = (struct simplex_Basis *)basis;
	if (!basis->compact_active)
		return simplex_basis_core_solve_pair(basis, first, second);
	/* Keep the certified path unchanged; the other 31 calls share mapping,
	 * sparse-LU traversal, and recovery for both right-hand sides. */
	if (mutable->compact_ftran_validation_countdown <= 0) {
		if (simplex_basis_base_ftran(basis, first) == lp_simplex_EXIT_FAILURE ||
		    simplex_basis_base_ftran(basis, second) == lp_simplex_EXIT_FAILURE)
			return lp_simplex_EXIT_FAILURE;
		return lp_simplex_EXIT_SUCCESS;
	}
	lp_simplex_memcpy(mutable->base_work, first,
		(size_t)basis->rows * sizeof(double));
	lp_simplex_memcpy(mutable->refine_work, second,
		(size_t)basis->rows * sizeof(double));
	if (simplex_basis_base_ftran_pair_once(basis, mutable->base_work,
			mutable->refine_work, first, second) == lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	mutable->compact_ftran_validation_countdown =
		__lp_simplex_MAX__(0,
			mutable->compact_ftran_validation_countdown - 2);
	return lp_simplex_EXIT_SUCCESS;
}


static int simplex_basis_base_btran_once(
		const struct simplex_Basis *basis, const double *right, double *vector)
{
	struct simplex_Basis *mutable = (struct simplex_Basis *)basis;
	const struct simplex_CscMatrix *matrix = basis->matrix;
	int i, k, row, position;
	lp_simplex_memset(vector, 0, (size_t)basis->rows * sizeof(double));
	for (row = 0; row < basis->rows; row++) {
		position = basis->logical_position[row];
		if (position >= 0)
			vector[row] = -right[position];
	}
	for (i = 0; i < basis->core_size; i++) {
		int column = basis->core_basis[i];
		double value = right[basis->core_position[i]];
		for (k = matrix->column_start[column];
		     k < matrix->column_start[column + 1]; k++) {
			row = matrix->row_index[k];
			if (basis->logical_position[row] >= 0)
				value -= matrix->value[k] * vector[row];
		}
		mutable->core_work[i] = value;
	}
	if (simplex_basis_core_solve(basis, mutable->core_work, 'T') ==
	    lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	for (i = 0; i < basis->core_size; i++)
		vector[basis->core_row[i]] = mutable->core_work[i];
	return lp_simplex_EXIT_SUCCESS;
}


static void simplex_basis_btran_residual(
		const struct simplex_Basis *basis, const double *right,
		const double *solution, double *residual, double *maximum)
{
	const struct simplex_CscMatrix *matrix = basis->matrix;
	int i, k;
	double scale = 1.;
	for (i = 0; i < basis->rows; i++) {
		int variable = basis->base_index[i];
		double product = 0.;
		if (variable >= basis->structural_columns)
			product = -solution[variable - basis->structural_columns];
		else
			for (k = matrix->column_start[variable];
			     k < matrix->column_start[variable + 1]; k++)
				product += matrix->value[k] * solution[matrix->row_index[k]];
		residual[i] = right[i] - product;
		scale = __lp_simplex_MAX__(scale, __lp_simplex_ABS__(right[i]));
	}
	*maximum = 0.;
	for (i = 0; i < basis->rows; i++)
		*maximum = __lp_simplex_MAX__(*maximum,
			__lp_simplex_ABS__(residual[i]) / scale);
}


static int simplex_basis_base_btran(
		const struct simplex_Basis *basis, double *vector)
{
	struct simplex_Basis *mutable = (struct simplex_Basis *)basis;
	int i, refinement, refined = 0;
	double residual;
	if (!basis->compact_active)
		return simplex_basis_core_solve(basis, vector, 'T');
	lp_simplex_memcpy(mutable->base_work, vector,
		(size_t)basis->rows * sizeof(double));
	if (simplex_basis_base_btran_once(basis, mutable->base_work, vector) ==
	    lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	if (mutable->compact_btran_validation_countdown > 0) {
		mutable->compact_btran_validation_countdown--;
		return lp_simplex_EXIT_SUCCESS;
	}
	mutable->profile_compact_btran_validations++;
	for (refinement = 0; refinement < 2; refinement++) {
		simplex_basis_btran_residual(basis, mutable->base_work, vector,
			mutable->refine_work, &residual);
		if (residual <= 1e-11)
			break;
		refined = 1;
		mutable->profile_compact_btran_refinements++;
		if (simplex_basis_base_btran_once(basis, mutable->refine_work,
				mutable->correction_work) == lp_simplex_EXIT_FAILURE)
			return lp_simplex_EXIT_FAILURE;
		for (i = 0; i < basis->rows; i++)
			vector[i] += mutable->correction_work[i];
	}
	mutable->compact_btran_validation_countdown = refined ? 0 : 31;
	return lp_simplex_EXIT_SUCCESS;
}


static int simplex_basis_apply_eta(
		const struct simplex_Basis *basis, const int update, double *vector)
{
	int k;
	int p = basis->eta_pivot[update];
	double pivot = basis->eta_pivot_value[update];
	double value;
	if (__lp_simplex_ABS__(pivot) <= 1e-14)
		return lp_simplex_EXIT_FAILURE;
	value = vector[p] / pivot;
	if (basis->eta_sparse[update])
		for (k = basis->eta_start[update];
		     k < basis->eta_start[update + 1]; k++)
			vector[basis->eta_index[k]] -= value * basis->eta_value[k];
	else
		lp_simplex_linalg_daxpy(basis->rows, -value,
			basis->eta + (size_t)update * basis->rows, 1, vector, 1);
	vector[p] = value;
	return lp_simplex_EXIT_SUCCESS;
}


static int simplex_basis_apply_eta_pair(
		const struct simplex_Basis *basis, const int update,
		double *first, double *second)
{
	int k;
	int p = basis->eta_pivot[update];
	double pivot = basis->eta_pivot_value[update];
	double first_value, second_value;
	if (__lp_simplex_ABS__(pivot) <= 1e-14)
		return lp_simplex_EXIT_FAILURE;
	first_value = first[p] / pivot;
	second_value = second[p] / pivot;
	if (basis->eta_sparse[update]) {
		for (k = basis->eta_start[update];
		     k < basis->eta_start[update + 1]; k++) {
			int index = basis->eta_index[k];
			double value = basis->eta_value[k];
			first[index] -= first_value * value;
			second[index] -= second_value * value;
		}
	} else {
		const double *eta = basis->eta + (size_t)update * basis->rows;
		for (k = 0; k < basis->rows; k++) {
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
	int p = basis->eta_pivot[update];
	double pivot = basis->eta_pivot_value[update];
	double pivot_value = vector[p];
	double dot = 0.;
	if (__lp_simplex_ABS__(pivot) <= 1e-14)
		return lp_simplex_EXIT_FAILURE;
	if (basis->eta_sparse[update])
		for (k = basis->eta_start[update];
		     k < basis->eta_start[update + 1]; k++)
			dot += basis->eta_value[k] * vector[basis->eta_index[k]];
	else
		dot = lp_simplex_linalg_ddot(basis->rows,
			basis->eta + (size_t)update * basis->rows, 1, vector, 1);
	vector[p] = (pivot_value - dot + pivot * pivot_value) / pivot;
	return lp_simplex_EXIT_SUCCESS;
}


static int simplex_basis_ftran_raw(
		const struct simplex_Basis *basis, double *vector)
{
	int update;
	if (simplex_basis_base_ftran(basis, vector) == lp_simplex_EXIT_FAILURE)
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
	return simplex_basis_base_btran(basis, vector);
}


int simplex_basis_ftran_pair(
		const struct simplex_Basis *basis, double *first, double *second)
{
	struct simplex_Basis *mutable = (struct simplex_Basis *)basis;
	clock_t started = 0;
	int result = lp_simplex_EXIT_SUCCESS;
	int update;
	if (basis->profile_enabled)
		started = clock();
	if (simplex_basis_base_ftran_pair(basis, first, second) ==
		lp_simplex_EXIT_FAILURE) {
		result = lp_simplex_EXIT_FAILURE;
	}
	if (result == lp_simplex_EXIT_SUCCESS)
		for (update = 0; update < basis->update_count; update++)
			if (simplex_basis_apply_eta_pair(basis, update, first, second) ==
			    lp_simplex_EXIT_FAILURE) {
				result = lp_simplex_EXIT_FAILURE;
				break;
			}
	if (basis->profile_enabled) {
		mutable->profile_ftran_seconds +=
			(double)(clock() - started) / (double)CLOCKS_PER_SEC;
		mutable->profile_ftran_calls += 2;
	}
	return result;
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


int simplex_basis_update(
		struct simplex_Basis *basis, const int leaving_position,
		const double *direction)
{
	int i, count, first, nonzeros;
	if (basis->update_count >= basis->update_limit)
		return 1;
	if (__lp_simplex_ABS__(direction[leaving_position]) <= 1e-14)
		return lp_simplex_EXIT_FAILURE;
	first = basis->eta_start[basis->update_count];
	count = first;
	for (i = 0; i < basis->rows; i++) {
		if (direction[i] != 0.) {
			basis->eta_index[count] = i;
			basis->eta_value[count] = direction[i];
			count++;
		}
	}
	nonzeros = count - first;
	basis->eta_sparse[basis->update_count] = (unsigned char)(
		basis->compact_active ||
		(basis->rows >= 4096 && nonzeros * 8 <= basis->rows));
	if (basis->eta_sparse[basis->update_count]) {
		basis->profile_eta_nonzeros += nonzeros;
		basis->profile_eta_slots += basis->rows;
	} else {
		count = first;
		lp_simplex_memcpy(basis->eta +
			(size_t)basis->update_count * basis->rows, direction,
			(size_t)basis->rows * sizeof(double));
	}
	basis->eta_pivot[basis->update_count] = leaving_position;
	basis->eta_pivot_value[basis->update_count] = direction[leaving_position];
	basis->update_count++;
	basis->eta_start[basis->update_count] = count;
	return lp_simplex_EXIT_SUCCESS;
}
