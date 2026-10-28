/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_basis.c
 * @brief Basis factorization backend.
 *
 * Revised-simplex code depends only on this API.
 */

#include "policy/factorization.h"
#include "policy/controls.h"
#include "policy/resources.h"
#include "simplex_basis.h"
#include "linalg.h"
#include "simplex_sparse_lu.h"
#include "simplex_klu_sparse.h"
#include "simplex_sparse_vector.h"
#include "utils.h"
#include "simplex_storage.h"

#ifdef LP_SIMPLEX_HAVE_KLU
#include <klu.h>
#endif

#include <lp_simplex/status.h>

#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>

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
	int allow_ft;
	int *core_basis;
	int *core_position;
	int *core_row;
	int *row_to_core;
	int *logical_position;
	int *base_index;
	int *column_to_core;
	int core_nonzeros;
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
	struct simplex_KluSparse sparse_klu;
	struct simplex_SparseVector sparse_core, sparse_lift;
#endif
	struct simplex_SparseLu sparse;
	/* Backend actually in use.  When KLU is available both factorizations are
	 * compiled, so a run that breaks down numerically can be retried with the
	 * other one instead of being reported as a solve failure. */
	int backend;
	double *eta_value;
	double *dense_eta;
	double *eta_pivot_value;
	int *eta_start;
	int *eta_index;
	int *eta_pivot;
	int *eta_dense_slot;
	unsigned char eta_ordered[SIMPLEX_BASIS_UPDATE_LIMIT];
	int eta_capacity;
	int dense_eta_capacity;
	int dense_eta_count;
	int update_count;
	int storage_failed;
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
	long profile_klu_sparse_success, profile_klu_sparse_declined;
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

static void simplex_basis_reset_updates(struct simplex_Basis *basis);

/* True when the internal sparse LU backs this basis.  Both backends are
 * compiled when KLU is present, so every backend-specific branch is a runtime
 * decision; DEFAULT resolves to KLU in that case. */
static int simplex_basis_uses_internal(const struct simplex_Basis *basis)
{
#ifdef LP_SIMPLEX_HAVE_KLU
	return basis->impl->backend == SIMPLEX_BASIS_BACKEND_INTERNAL;
#else
	(void)basis;
	return 1;
#endif
}

void simplex_basis_set_backend(struct simplex_Basis *basis, const int backend)
{
	if (basis == NULL || basis->impl == NULL)
		return;
	basis->impl->backend = backend;
#ifdef LP_SIMPLEX_HAVE_KLU
	if (backend == SIMPLEX_BASIS_BACKEND_INTERNAL) {
		simplex_klu_sparse_destroy(&basis->impl->sparse_klu);
		if (basis->impl->numeric != NULL)
			klu_free_numeric(
				&basis->impl->numeric, &basis->impl->common);
		if (basis->impl->symbolic != NULL)
			klu_free_symbolic(
				&basis->impl->symbolic, &basis->impl->common);
	}
#endif
	simplex_basis_reset_updates(basis);
}

int simplex_basis_resolved_backend(const struct simplex_Basis *basis)
{
	if (basis == NULL || basis->impl == NULL)
		return SIMPLEX_BASIS_BACKEND_INTERNAL;
#ifdef LP_SIMPLEX_HAVE_KLU
	if (basis->impl->backend != SIMPLEX_BASIS_BACKEND_INTERNAL)
		return SIMPLEX_BASIS_BACKEND_KLU;
#endif
	return SIMPLEX_BASIS_BACKEND_INTERNAL;
}

int simplex_basis_create(struct simplex_Basis *basis,
	const struct simplex_CscMatrix *matrix, const int structural_columns,
	int *index)
{
	return simplex_basis_create_configured(basis, matrix,
		structural_columns, index,
		getenv(SIMPLEX_ENV_STABLE_KLU) != NULL);
}

/* Constructor-owned arrays only; numerical updates/factors have their own
 * storage lifetimes and admission. No global memory quota is implied. */
struct basis_StartupStorage {
	int eta_capacity;
	size_t core, work, starts, eta_values, eta_indices, eta_starts,
		eta_pivots;
};

static int basis_startup_storage(int rows, int structural, size_t limit,
	struct basis_StartupStorage *storage)
{
	struct basis_StartupStorage next = {0};
	*storage = next;
	if (rows < 0 || structural < 0 || structural > INT_MAX - rows ||
		SIMPLEX_BASIS_UPDATE_LIMIT <= 0 ||
		SIMPLEX_BASIS_UPDATE_LIMIT > (INT_MAX - 1) / 3)
		return 0;
	next.eta_capacity =
		rows > INT_MAX / 4 ? INT_MAX : __lp_simplex_MAX__(rows * 4, 16);
	if (!simplex_storage_bytes(
		    (size_t)rows, 6 * sizeof(int), limit, &next.core) ||
		!simplex_storage_bytes(
			(size_t)rows, 7 * sizeof(double), limit, &next.work) ||
		!simplex_storage_bytes(
			(size_t)rows + 1, sizeof(int), limit, &next.starts) ||
		!simplex_storage_bytes((size_t)next.eta_capacity,
			sizeof(double), limit, &next.eta_values) ||
		!simplex_storage_bytes((size_t)next.eta_capacity, sizeof(int),
			limit, &next.eta_indices) ||
		!simplex_storage_bytes(
			(size_t)3 * SIMPLEX_BASIS_UPDATE_LIMIT + 1, sizeof(int),
			limit, &next.eta_starts) ||
		!simplex_storage_bytes((size_t)SIMPLEX_BASIS_UPDATE_LIMIT,
			sizeof(double), limit, &next.eta_pivots))
		return 0;
	*storage = next;
	return 1;
}

int simplex_basis_create_configured(struct simplex_Basis *basis,
	const struct simplex_CscMatrix *matrix, const int structural_columns,
	int *index, int stable_klu)
{
	struct basis_StartupStorage storage;
	(void)stable_klu;
	if (basis == NULL || matrix == NULL)
		return lp_simplex_EXIT_FAILURE;

	lp_simplex_memset(basis, 0, sizeof(*basis));
	if (structural_columns > matrix->columns ||
		!basis_startup_storage(
			matrix->rows, structural_columns, (size_t)-1, &storage))
		return lp_simplex_EXIT_FAILURE;
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
	basis->impl->allow_ft = 1;
	basis->impl->core_basis = (int *)lp_simplex_malloc(storage.core);
	basis->impl->core_position =
		basis->impl->core_basis != NULL
			? basis->impl->core_basis + matrix->rows
			: NULL;
	basis->impl->core_row =
		basis->impl->core_basis != NULL
			? basis->impl->core_basis + (size_t)2 * matrix->rows
			: NULL;
	basis->impl->row_to_core =
		basis->impl->core_basis != NULL
			? basis->impl->core_basis + (size_t)3 * matrix->rows
			: NULL;
	basis->impl->logical_position =
		basis->impl->core_basis != NULL
			? basis->impl->core_basis + (size_t)4 * matrix->rows
			: NULL;
	basis->impl->base_index =
		basis->impl->core_basis != NULL
			? basis->impl->core_basis + (size_t)5 * matrix->rows
			: NULL;
	basis->impl->base_work = (double *)lp_simplex_malloc(storage.work);
	basis->impl->core_work = basis->impl->base_work != NULL
					 ? basis->impl->base_work + matrix->rows
					 : NULL;
	basis->impl->refine_work =
		basis->impl->base_work != NULL
			? basis->impl->base_work + (size_t)2 * matrix->rows
			: NULL;
	basis->impl->correction_work =
		basis->impl->base_work != NULL
			? basis->impl->base_work + (size_t)3 * matrix->rows
			: NULL;
	basis->impl->certify_right =
		basis->impl->base_work != NULL
			? basis->impl->base_work + (size_t)4 * matrix->rows
			: NULL;
	basis->impl->certify_residual =
		basis->impl->base_work != NULL
			? basis->impl->base_work + (size_t)5 * matrix->rows
			: NULL;
	basis->impl->certify_correction =
		basis->impl->base_work != NULL
			? basis->impl->base_work + (size_t)6 * matrix->rows
			: NULL;
#ifdef LP_SIMPLEX_HAVE_KLU
	basis->impl->base_column_start =
		(int *)lp_simplex_malloc(storage.starts);
	basis->impl->base_row_index = NULL;
	basis->impl->base_value = NULL;
	basis->impl->symbolic = NULL;
	basis->impl->numeric = NULL;
	klu_defaults(&basis->impl->common);
	if (stable_klu)
		basis->impl->common.tol = SIMPLEX_KLU_PIVOT_PREFERENCE;
#endif
	lp_simplex_memset(&basis->impl->sparse, 0, sizeof(basis->impl->sparse));
	basis->impl->backend = SIMPLEX_BASIS_BACKEND_DEFAULT;
	basis->impl->eta_capacity = storage.eta_capacity;
	basis->impl->eta_value =
		(double *)lp_simplex_malloc(storage.eta_values);
	/* Dense eta columns are grown by actual use.  Reserving update_limit * rows
	 * here makes the first dense update consume hundreds of megabytes on a
	 * large model even if reinversion follows immediately. */
	basis->impl->dense_eta = NULL;
	basis->impl->dense_eta_capacity = 0;
	basis->impl->dense_eta_count = 0;
	basis->impl->eta_index = (int *)lp_simplex_malloc(storage.eta_indices);
	basis->impl->eta_start = (int *)lp_simplex_malloc(storage.eta_starts);
	basis->impl->eta_pivot = basis->impl->eta_start != NULL
					 ? basis->impl->eta_start +
						   SIMPLEX_BASIS_UPDATE_LIMIT +
						   1
					 : NULL;
	basis->impl->eta_pivot_value =
		(double *)lp_simplex_malloc(storage.eta_pivots);
	basis->impl->eta_dense_slot =
		basis->impl->eta_pivot != NULL
			? basis->impl->eta_pivot + SIMPLEX_BASIS_UPDATE_LIMIT
			: NULL;
	basis->impl->update_count = 0;
	basis->impl->storage_failed = 0;
	basis->impl->update_work = 0;
	basis->impl->factor_work = 0;
	basis->impl->eta_apply_work = 0;
	basis->impl->minimum_relative_pivot = 1.;
	basis->impl->compact_ftran_validation_countdown = 0;
	basis->impl->compact_btran_validation_countdown = 0;
	if (basis->impl->eta_start != NULL)
		basis->impl->eta_start[0] = 0;
	basis->impl->update_limit = SIMPLEX_BASIS_UPDATE_LIMIT;
	basis->impl->profile_enabled = getenv(SIMPLEX_ENV_PROFILE) != NULL;
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
		basis->impl->core_basis == NULL ||
		basis->impl->core_position == NULL ||
		basis->impl->core_row == NULL ||
		basis->impl->row_to_core == NULL ||
		basis->impl->logical_position == NULL ||
		basis->impl->base_work == NULL ||
		basis->impl->base_index == NULL ||
		basis->impl->core_work == NULL ||
		basis->impl->refine_work == NULL ||
		basis->impl->correction_work == NULL ||
		basis->impl->certify_right == NULL ||
		basis->impl->certify_residual == NULL ||
		basis->impl->certify_correction == NULL ||
		basis->impl->eta_value == NULL ||
		basis->impl->eta_index == NULL ||
		basis->impl->eta_start == NULL ||
		basis->impl->eta_pivot == NULL ||
		basis->impl->eta_pivot_value == NULL ||
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
	simplex_klu_sparse_destroy(&basis->impl->sparse_klu);
	simplex_sparse_vector_destroy(&basis->impl->sparse_core);
	simplex_sparse_vector_destroy(&basis->impl->sparse_lift);
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
#endif
	simplex_sparse_lu_destroy(&basis->impl->sparse);
	lp_simplex_free(basis->impl->eta_value);
	lp_simplex_free(basis->impl->dense_eta);
	lp_simplex_free(basis->impl->eta_index);
	lp_simplex_free(basis->impl->eta_start);
	lp_simplex_free(basis->impl->eta_pivot_value);
	lp_simplex_free(basis->impl->core_basis);
	lp_simplex_free(basis->impl->base_work);
	lp_simplex_free(basis->impl->column_to_core);
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
			if (row < 0 || row >= m ||
				basis->impl->logical_position[row] >= 0)
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
	basis->impl->storage_failed = 0;
	basis->impl->update_count = 0;
	basis->impl->update_work = 0;
	basis->impl->minimum_relative_pivot = 1.;
	basis->impl->eta_start[0] = 0;
	basis->impl->dense_eta_count = 0;
}

#ifdef LP_SIMPLEX_HAVE_KLU
static void simplex_basis_prepare_packed(struct simplex_Basis *basis)
{
	struct simplex_BasisImpl *b = basis->impl;
	size_t bytes;
	if (!b->compact_active || b->column_to_core == NULL ||
		!simplex_storage_bytes((size_t)b->rows + b->core_size,
			2 * sizeof(int) + 2 * sizeof(double) +
				sizeof(unsigned int),
			SIMPLEX_BASIS_PACKED_BYTES, &bytes))
		return;
	if (b->sparse_core.dimension != b->core_size) {
		simplex_sparse_vector_destroy(&b->sparse_core);
		if (simplex_sparse_vector_create(&b->sparse_core, b->core_size))
			goto REFUSED;
	}
	if (!b->sparse_lift.index &&
		simplex_sparse_vector_create(&b->sparse_lift, b->rows))
		goto REFUSED;
	return;
REFUSED:
	simplex_sparse_vector_destroy(&b->sparse_core);
	simplex_sparse_vector_destroy(&b->sparse_lift);
}
#endif

static int simplex_basis_factorize_impl(struct simplex_Basis *basis)
{
	int core_size;
#ifdef LP_SIMPLEX_HAVE_KLU
	int column, k;
	int nonzeros = 0;
	int factor_size;
	size_t index_bytes, value_bytes;
#endif
#ifdef LP_SIMPLEX_HAVE_KLU
	simplex_klu_sparse_destroy(&basis->impl->sparse_klu);
#endif
	if (simplex_basis_build_core(basis) == lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	core_size = basis->impl->core_size;
	basis->impl->compact_active =
		basis->impl->compact_requested &&
		basis->impl->rows >= SIMPLEX_BASIS_COMPACT_MIN_ROWS &&
		core_size < basis->impl->rows &&
		core_size * SIMPLEX_BASIS_COMPACT_ROW_RATIO <=
			basis->impl->rows;
	basis->impl->factor_size =
		basis->impl->compact_active ? core_size : basis->impl->rows;
	basis->impl->core_nonzeros = 0;
	{
		int j;
		for (j = 0; j < core_size; j++) {
			int column = basis->impl->core_basis[j];
			int count =
				basis->impl->matrix->column_start[column + 1] -
				basis->impl->matrix->column_start[column];
			basis->impl->core_nonzeros =
				count > INT_MAX - basis->impl->core_nonzeros
					? INT_MAX
					: basis->impl->core_nonzeros + count;
		}
	}
	/* Optional inverse map: the immutable CSR view can then apply D^T using
	 * only nonzero logical RHS entries in B = [C 0; D -I]. Allocation refusal
	 * retains the CSC traversal; it is not a factorization failure. */
	if (basis->impl->compact_active &&
		basis->impl->matrix->row_start != NULL &&
		basis->impl->matrix->column_index != NULL &&
		basis->impl->matrix->row_value != NULL) {
		int j;
		size_t bytes;
		if (basis->impl->column_to_core == NULL &&
			simplex_storage_bytes(
				(size_t)basis->impl->structural_columns,
				sizeof(int), SIMPLEX_BASIS_CORE_MAP_BYTES,
				&bytes))
			basis->impl->column_to_core =
				(int *)lp_simplex_malloc(bytes);
		if (basis->impl->column_to_core != NULL) {
			for (j = 0; j < basis->impl->structural_columns; j++)
				basis->impl->column_to_core[j] = -1;
			for (j = 0; j < core_size; j++)
				basis->impl->column_to_core[basis->impl
						->core_basis[j]] = j;
		}
	}
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
	if (basis->impl->backend != SIMPLEX_BASIS_BACKEND_INTERNAL) {
		factor_size = basis->impl->factor_size;
		if (basis->impl->numeric != NULL)
			klu_free_numeric(
				&basis->impl->numeric, &basis->impl->common);
		if (basis->impl->symbolic != NULL)
			klu_free_symbolic(
				&basis->impl->symbolic, &basis->impl->common);
		lp_simplex_free(basis->impl->base_row_index);
		lp_simplex_free(basis->impl->base_value);
		basis->impl->base_row_index = NULL;
		basis->impl->base_value = NULL;
		if (factor_size == 0) {
			simplex_basis_prepare_packed(basis);
			simplex_basis_reset_updates(basis);
			return lp_simplex_EXIT_SUCCESS;
		}
		for (column = 0; column < factor_size; column++) {
			int variable =
				basis->impl->compact_active
					? basis->impl->core_basis[column]
					: basis->impl->base_index[column];

			if (!basis->impl->compact_active &&
				variable >= basis->impl->structural_columns) {
				if (nonzeros == INT_MAX) {
					basis->impl->storage_failed = 1;
					return lp_simplex_EXIT_FAILURE;
				}
				nonzeros++;
			} else
				for (k = basis->impl->matrix
						->column_start[variable];
					k <
					basis->impl->matrix
						->column_start[variable + 1];
					k++)
					if (!basis->impl->compact_active ||
						basis->impl->row_to_core[basis
								->impl->matrix
								->row_index
									[k]] >=
							0) {
						if (nonzeros == INT_MAX) {
							basis->impl
								->storage_failed =
								1;
							return lp_simplex_EXIT_FAILURE;
						}
						nonzeros++;
					}
		}
		if (!simplex_storage_bytes((size_t)nonzeros, sizeof(int),
			    (size_t)-1, &index_bytes) ||
			!simplex_storage_bytes((size_t)nonzeros, sizeof(double),
				(size_t)-1, &value_bytes)) {
			basis->impl->storage_failed = 1;
			return lp_simplex_EXIT_FAILURE;
		}
		basis->impl->base_row_index =
			(int *)lp_simplex_malloc(index_bytes);
		basis->impl->base_value =
			(double *)lp_simplex_malloc(value_bytes);
		if (basis->impl->base_row_index == NULL ||
			basis->impl->base_value == NULL) {
			basis->impl->storage_failed = 1;
			return lp_simplex_EXIT_FAILURE;
		}
		nonzeros = 0;
		for (column = 0; column < factor_size; column++) {
			int variable =
				basis->impl->compact_active
					? basis->impl->core_basis[column]
					: basis->impl->base_index[column];

			basis->impl->base_column_start[column] = nonzeros;
			if (!basis->impl->compact_active &&
				variable >= basis->impl->structural_columns) {
				basis->impl->base_row_index[nonzeros] =
					variable -
					basis->impl->structural_columns;
				basis->impl->base_value[nonzeros++] = -1.;
			} else
				for (k = basis->impl->matrix
						->column_start[variable];
					k <
					basis->impl->matrix
						->column_start[variable + 1];
					k++) {
					int row =
						basis->impl->compact_active
							? basis->impl->row_to_core[basis
									  ->impl
									  ->matrix
									  ->row_index
										  [k]]
							: basis->impl->matrix
								  ->row_index
									  [k];

					if (row >= 0) {
						basis->impl->base_row_index
							[nonzeros] = row;
						basis->impl
							->base_value[nonzeros] =
							basis->impl->matrix
								->value[k];
						nonzeros++;
					}
				}
		}
		basis->impl->base_column_start[factor_size] = nonzeros;
		basis->impl->symbolic = klu_analyze(factor_size,
			basis->impl->base_column_start,
			basis->impl->base_row_index, &basis->impl->common);
		if (basis->impl->symbolic == NULL) {
			basis->impl->storage_failed =
				basis->impl->common.status ==
					KLU_OUT_OF_MEMORY ||
				basis->impl->common.status == KLU_TOO_LARGE;
			return lp_simplex_EXIT_FAILURE;
		}
		basis->impl->numeric = klu_factor(
			basis->impl->base_column_start,
			basis->impl->base_row_index, basis->impl->base_value,
			basis->impl->symbolic, &basis->impl->common);
		if (basis->impl->numeric == NULL) {
			basis->impl->storage_failed =
				basis->impl->common.status ==
					KLU_OUT_OF_MEMORY ||
				basis->impl->common.status == KLU_TOO_LARGE;
			return lp_simplex_EXIT_FAILURE;
		}
		klu_flops(basis->impl->symbolic, basis->impl->numeric,
			&basis->impl->common);
		if (basis->impl->profile_enabled) {
			klu_common diagnostic = basis->impl->common;
			if (klu_rcond(basis->impl->symbolic,
				    basis->impl->numeric, &diagnostic))
				fprintf(stderr,
					"basis profile: KLU dimension=%d compact=%d reciprocal-diagonal=%.17g\n",
					factor_size,
					basis->impl->compact_active,
					diagnostic.rcond);
		}
		(void)simplex_klu_sparse_build(&basis->impl->sparse_klu,
			basis->impl->symbolic, basis->impl->numeric,
			&basis->impl->common);
		simplex_basis_prepare_packed(basis);
		simplex_basis_reset_updates(basis);
		return lp_simplex_EXIT_SUCCESS;
	}
#endif
	if (basis->impl->sparse.dimension != basis->impl->factor_size) {
		simplex_sparse_lu_destroy(&basis->impl->sparse);
		if (basis->impl->factor_size > 0 &&
			simplex_sparse_lu_create(&basis->impl->sparse,
				basis->impl->factor_size) ==
				lp_simplex_EXIT_FAILURE)
			return lp_simplex_EXIT_FAILURE;
	}
	if (basis->impl->factor_size == 0 ||
		(basis->impl->compact_active &&
			simplex_sparse_lu_factorize_submatrix(
				&basis->impl->sparse, basis->impl->matrix,
				basis->impl->core_basis,
				basis->impl->row_to_core) ==
				lp_simplex_EXIT_SUCCESS) ||
		(!basis->impl->compact_active &&
			simplex_sparse_lu_factorize(&basis->impl->sparse,
				basis->impl->matrix,
				basis->impl->structural_columns,
				basis->impl->base_index) ==
				lp_simplex_EXIT_SUCCESS)) {
		simplex_basis_reset_updates(basis);
		return lp_simplex_EXIT_SUCCESS;
	}
	return lp_simplex_EXIT_FAILURE;
}

int simplex_basis_factorize(struct simplex_Basis *basis)
{
	clock_t started = 0;
	double elapsed = 0.;
	int result;

	basis->impl->storage_failed = 0;
	basis->impl->sparse.storage_failed = 0;
	if (basis->impl->profile_enabled)
		started = clock();
	result = simplex_basis_factorize_impl(basis);
	if (result == lp_simplex_EXIT_FAILURE &&
		simplex_basis_uses_internal(basis))
		basis->impl->storage_failed =
			basis->impl->sparse.storage_failed;
	if (result == lp_simplex_EXIT_SUCCESS) {
		if (basis->impl->profile_enabled)
			elapsed = (double)(clock() - started) /
				  (double)CLOCKS_PER_SEC;
#ifdef LP_SIMPLEX_HAVE_KLU
		if (basis->impl->backend == SIMPLEX_BASIS_BACKEND_INTERNAL)
			basis->impl->factor_work =
				basis->impl->sparse.factor_work;
		else
			basis->impl->factor_work =
				(long)(basis->impl->common.flops +
					basis->impl->common.work +
					basis->impl->factor_size +
					(basis->impl->base_column_start != NULL
							? basis->impl->base_column_start
								  [basis->impl->factor_size]
							: 0) +
					(basis->impl->numeric != NULL
							? basis->impl->numeric
									  ->lnz +
								  basis->impl
									  ->numeric
									  ->unz +
								  basis->impl
									  ->numeric
									  ->nzoff
							: 0));
#else
		basis->impl->factor_work = basis->impl->sparse.factor_work;
#endif
		basis->impl->eta_apply_work = 0;
	}
	if (basis->impl->profile_enabled) {
		basis->impl->profile_factor_seconds += elapsed;
	}
	basis->impl->profile_factor_calls++;
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
	if (basis->impl->backend == SIMPLEX_BASIS_BACKEND_INTERNAL)
		return simplex_sparse_lu_solve(
			&((struct simplex_Basis *)basis)->impl->sparse, vector,
			trans == 'T');
	success = simplex_klu_sparse_solve(
		&((struct simplex_Basis *)basis)->impl->sparse_klu, vector,
		trans == 'T');
	if (success == 0)
		((struct simplex_Basis *)basis)
			->impl->profile_klu_sparse_success++;
	else if (success == 1)
		((struct simplex_Basis *)basis)
			->impl->profile_klu_sparse_declined++;
	if (success != 1)
		return success == 0 ? lp_simplex_EXIT_SUCCESS
				    : lp_simplex_EXIT_FAILURE;
	if (trans == 'N')
		success = klu_solve(basis->impl->symbolic, basis->impl->numeric,
			basis->impl->factor_size, 1, vector,
			(klu_common *)&basis->impl->common);
	else
		success = klu_tsolve(basis->impl->symbolic,
			basis->impl->numeric, basis->impl->factor_size, 1,
			vector, (klu_common *)&basis->impl->common);
	return success ? lp_simplex_EXIT_SUCCESS : lp_simplex_EXIT_FAILURE;
#else
	return simplex_sparse_lu_solve(
		&((struct simplex_Basis *)basis)->impl->sparse, vector,
		trans == 'T');
#endif
}

static int simplex_basis_core_solve_pair(
	const struct simplex_Basis *basis, double *first, double *second)
{
	if (basis->impl->factor_size == 0)
		return lp_simplex_EXIT_SUCCESS;
#ifdef LP_SIMPLEX_HAVE_KLU
	if (basis->impl->backend == SIMPLEX_BASIS_BACKEND_INTERNAL)
		return simplex_sparse_lu_solve_pair(
			&((struct simplex_Basis *)basis)->impl->sparse, first,
			second);
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

	/* Core and logical positions partition the complete basis. Each output
	 * is assigned below before use; a separate dense zero pass is redundant. */
	for (i = 0; i < basis->impl->core_size; i++)
		mutable->impl->core_work[i] = right[basis->impl->core_row[i]];
	if (simplex_basis_core_solve(basis, mutable->impl->core_work, 'N') ==
		lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	for (i = 0; i < basis->impl->core_size; i++)
		vector[basis->impl->core_position[i]] =
			mutable->impl->core_work[i];
	for (row = 0; row < basis->impl->rows; row++) {
		position = basis->impl->logical_position[row];
		if (position >= 0)
			vector[position] = -right[row];
	}
	for (i = 0; i < basis->impl->core_size; i++) {
		int column = basis->impl->core_basis[i];
		double value = mutable->impl->core_work[i];

		if (value == 0.)
			continue;
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

static int simplex_basis_base_ftran_pair_once(const struct simplex_Basis *basis,
	const double *first_right, const double *second_right, double *first,
	double *second)
{
	struct simplex_Basis *mutable = (struct simplex_Basis *)basis;
	const struct simplex_CscMatrix *matrix = basis->impl->matrix;
	int i, k, row, position;

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

		if (first_value == 0. && second_value == 0.)
			continue;
		for (k = matrix->column_start[column];
			k < matrix->column_start[column + 1]; k++) {
			row = matrix->row_index[k];
			position = basis->impl->logical_position[row];
			if (position >= 0) {
				first[position] +=
					matrix->value[k] * first_value;
				second[position] +=
					matrix->value[k] * second_value;
			}
		}
	}
	return lp_simplex_EXIT_SUCCESS;
}

static void simplex_basis_ftran_residual(const struct simplex_Basis *basis,
	const double *right, const double *solution, double *residual,
	double *maximum)
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
			residual[variable - basis->impl->structural_columns] +=
				value;
		else
			for (k = matrix->column_start[variable];
				k < matrix->column_start[variable + 1]; k++)
				residual[matrix->row_index[k]] -=
					matrix->value[k] * value;
	}
	*maximum = 0.;
	for (i = 0; i < basis->impl->rows; i++)
		*maximum = __lp_simplex_MAX__(
			*maximum, __lp_simplex_ABS__(residual[i]) / scale);
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
	if (simplex_basis_base_ftran_once(basis, mutable->impl->base_work,
		    vector) == lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	if (mutable->impl->compact_ftran_validation_countdown > 0) {
		mutable->impl->compact_ftran_validation_countdown--;
		return lp_simplex_EXIT_SUCCESS;
	}
	mutable->impl->profile_compact_ftran_validations++;
	for (refinement = 0; refinement < SIMPLEX_BASIS_COMPACT_REFINEMENTS;
		refinement++) {
		simplex_basis_ftran_residual(basis, mutable->impl->base_work,
			vector, mutable->impl->refine_work, &residual);
		if (residual <= SIMPLEX_BASIS_COMPACT_RESIDUAL_TARGET)
			break;
		refined = 1;
		mutable->impl->profile_compact_ftran_refinements++;
		if (simplex_basis_base_ftran_once(basis,
			    mutable->impl->refine_work,
			    mutable->impl->correction_work) ==
			lp_simplex_EXIT_FAILURE)
			return lp_simplex_EXIT_FAILURE;
		for (i = 0; i < basis->impl->rows; i++)
			vector[i] += mutable->impl->correction_work[i];
	}
	mutable->impl->compact_ftran_validation_countdown =
		refined ? 0 : SIMPLEX_BASIS_VALIDATION_PERIOD - 1;
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
		if (simplex_basis_base_ftran(basis, first) ==
				lp_simplex_EXIT_FAILURE ||
			simplex_basis_base_ftran(basis, second) ==
				lp_simplex_EXIT_FAILURE)
			return lp_simplex_EXIT_FAILURE;
		return lp_simplex_EXIT_SUCCESS;
	}
	lp_simplex_memcpy(mutable->impl->base_work, first,
		(size_t)basis->impl->rows * sizeof(double));
	lp_simplex_memcpy(mutable->impl->refine_work, second,
		(size_t)basis->impl->rows * sizeof(double));
	if (simplex_basis_base_ftran_pair_once(basis, mutable->impl->base_work,
		    mutable->impl->refine_work, first,
		    second) == lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	mutable->impl->compact_ftran_validation_countdown = __lp_simplex_MAX__(
		0, mutable->impl->compact_ftran_validation_countdown - 2);
	return lp_simplex_EXIT_SUCCESS;
}

static int simplex_basis_base_btran_once(
	const struct simplex_Basis *basis, const double *right, double *vector)
{
	struct simplex_Basis *mutable = (struct simplex_Basis *)basis;
	const struct simplex_CscMatrix *matrix = basis->impl->matrix;
	int i, k, row, position;
	int row_view = basis->impl->column_to_core != NULL &&
		       matrix->row_start != NULL &&
		       matrix->column_index != NULL &&
		       matrix->row_value != NULL;
	size_t row_work = 0;

	for (row = 0; row < basis->impl->rows; row++) {
		position = basis->impl->logical_position[row];
		vector[row] = position >= 0 ? -right[position] : 0.;
		if (row_view && vector[row] != 0.)
			row_work += matrix->row_start[row + 1] -
				    matrix->row_start[row];
	}
	if (row_view && row_work <= (size_t)basis->impl->core_nonzeros) {
		for (i = 0; i < basis->impl->core_size; i++)
			mutable->impl->core_work[i] =
				right[basis->impl->core_position[i]];
		/* Visit the cheaper support using the existing CSR view. No small
		 * coefficient or RHS is dropped, and the base residual audit and
		 * refinement path remain unchanged. */
		for (row = 0; row < basis->impl->rows; row++) {
			double value = vector[row];
			if (value == 0.)
				continue;
			for (k = matrix->row_start[row];
				k < matrix->row_start[row + 1]; k++) {
				int column = matrix->column_index[k];
				int core =
					column < basis->impl->structural_columns
						? basis->impl->column_to_core
							  [column]
						: -1;
				if (core >= 0)
					mutable->impl->core_work[core] -=
						matrix->row_value[k] * value;
			}
		}
	} else
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

static void simplex_basis_btran_residual(const struct simplex_Basis *basis,
	const double *right, const double *solution, double *residual,
	double *maximum)
{
	const struct simplex_CscMatrix *matrix = basis->impl->matrix;
	int i, k;
	double scale = 1.;

	for (i = 0; i < basis->impl->rows; i++) {
		int variable = basis->impl->base_index[i];
		double product = 0.;

		if (variable >= basis->impl->structural_columns)
			product = -solution[variable -
					    basis->impl->structural_columns];
		else
			for (k = matrix->column_start[variable];
				k < matrix->column_start[variable + 1]; k++)
				product += matrix->value[k] *
					   solution[matrix->row_index[k]];
		residual[i] = right[i] - product;
		scale = __lp_simplex_MAX__(scale, __lp_simplex_ABS__(right[i]));
	}
	*maximum = 0.;
	for (i = 0; i < basis->impl->rows; i++)
		*maximum = __lp_simplex_MAX__(
			*maximum, __lp_simplex_ABS__(residual[i]) / scale);
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
	if (simplex_basis_base_btran_once(basis, mutable->impl->base_work,
		    vector) == lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	if (mutable->impl->compact_btran_validation_countdown > 0) {
		mutable->impl->compact_btran_validation_countdown--;
		return lp_simplex_EXIT_SUCCESS;
	}
	mutable->impl->profile_compact_btran_validations++;
	for (refinement = 0; refinement < SIMPLEX_BASIS_COMPACT_REFINEMENTS;
		refinement++) {
		simplex_basis_btran_residual(basis, mutable->impl->base_work,
			vector, mutable->impl->refine_work, &residual);
		if (residual <= SIMPLEX_BASIS_COMPACT_RESIDUAL_TARGET)
			break;
		refined = 1;
		mutable->impl->profile_compact_btran_refinements++;
		if (simplex_basis_base_btran_once(basis,
			    mutable->impl->refine_work,
			    mutable->impl->correction_work) ==
			lp_simplex_EXIT_FAILURE)
			return lp_simplex_EXIT_FAILURE;
		for (i = 0; i < basis->impl->rows; i++)
			vector[i] += mutable->impl->correction_work[i];
	}
	mutable->impl->compact_btran_validation_countdown =
		refined ? 0 : SIMPLEX_BASIS_VALIDATION_PERIOD - 1;
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

	if (__lp_simplex_ABS__(pivot) <= SIMPLEX_BASIS_ABSOLUTE_PIVOT_FLOOR)
		return lp_simplex_EXIT_FAILURE;
	if (vector[p] == 0.)
		return lp_simplex_EXIT_SUCCESS;
	value = vector[p] / pivot;
	if (dense_slot < 0)
		for (k = basis->impl->eta_start[update];
			k < basis->impl->eta_start[update + 1]; k++)
			vector[basis->impl->eta_index[k]] -=
				value * basis->impl->eta_value[k];
	else
		lp_simplex_linalg_daxpy(basis->impl->rows, -value,
			basis->impl->dense_eta +
				(size_t)dense_slot * basis->impl->rows,
			1, vector, 1);
	vector[p] = value;
	return lp_simplex_EXIT_SUCCESS;
}

static int simplex_basis_apply_eta_pair(const struct simplex_Basis *basis,
	const int update, double *first, double *second)
{
	int k;
	int dense_slot = basis->impl->eta_dense_slot[update];
	int p = basis->impl->eta_pivot[update];
	double pivot = basis->impl->eta_pivot_value[update];
	double first_value, second_value;

	if (__lp_simplex_ABS__(pivot) <= SIMPLEX_BASIS_ABSOLUTE_PIVOT_FLOOR)
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

	if (__lp_simplex_ABS__(pivot) <= SIMPLEX_BASIS_ABSOLUTE_PIVOT_FLOOR)
		return lp_simplex_EXIT_FAILURE;
	if (dense_slot < 0)
		for (k = basis->impl->eta_start[update];
			k < basis->impl->eta_start[update + 1]; k++)
			dot += basis->impl->eta_value[k] *
			       vector[basis->impl->eta_index[k]];
	else
		dot = lp_simplex_linalg_ddot(basis->impl->rows,
			basis->impl->dense_eta +
				(size_t)dense_slot * basis->impl->rows,
			1, vector, 1);
	vector[p] = (pivot_value - dot + pivot * pivot_value) / pivot;
	return lp_simplex_EXIT_SUCCESS;
}

static int simplex_basis_ftran_raw(
	const struct simplex_Basis *basis, double *vector)
{
	int update, first_update = 0;

	if (simplex_basis_base_ftran(basis, vector) == lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	if (simplex_basis_uses_internal(basis) &&
		!basis->impl->compact_active &&
		simplex_sparse_lu_ft_active(&basis->impl->sparse))
		first_update = basis->impl->sparse.ft_update_count;
	for (update = first_update; update < basis->impl->update_count;
		update++) {
		((struct simplex_Basis *)basis)->impl->eta_apply_work +=
			basis->impl->eta_dense_slot[update] < 0
				? basis->impl->eta_start[update + 1] -
					  basis->impl->eta_start[update]
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

	if (simplex_basis_uses_internal(basis) &&
		!basis->impl->compact_active &&
		simplex_sparse_lu_ft_active(&basis->impl->sparse))
		first_update = basis->impl->sparse.ft_update_count;
	for (update = basis->impl->update_count - 1; update >= first_update;
		update--) {
		((struct simplex_Basis *)basis)->impl->eta_apply_work +=
			basis->impl->eta_dense_slot[update] < 0
				? basis->impl->eta_start[update + 1] -
					  basis->impl->eta_start[update]
				: basis->impl->rows;
		if (simplex_basis_apply_eta_transpose(basis, update, vector) ==
			lp_simplex_EXIT_FAILURE)
			return lp_simplex_EXIT_FAILURE;
	}
	return simplex_basis_base_btran(basis, vector);
}

#ifdef LP_SIMPLEX_HAVE_KLU
static int simplex_basis_base_ftran_packed_impl(
	struct simplex_Basis *basis, struct simplex_SparseVector *vector)
{
	struct simplex_BasisImpl *b = basis->impl;
	const struct simplex_CscMatrix *matrix = b->matrix;
	struct simplex_SparseVector *core = &b->sparse_core,
				    *lift = &b->sparse_lift;
	int i, k, result;
	if (simplex_basis_uses_internal(basis))
		return 1;
	if (!b->compact_active) {
		result = simplex_klu_sparse_solve_packed(
			&b->sparse_klu, vector, 0);
		if (!result)
			b->profile_klu_sparse_success++;
		else if (result == 1)
			b->profile_klu_sparse_declined++;
		return result;
	}
	if (b->compact_ftran_validation_countdown <= 0 || !lift->index ||
		core->dimension != b->core_size)
		return 1;
	simplex_sparse_vector_clear(core);
	simplex_sparse_vector_clear(lift);
	for (i = 0; i < vector->count; i++) {
		int row = vector->index[i], position = b->row_to_core[row];
		double value = vector->value[i];
		if (value == 0.)
			continue;
		if (position >= 0)
			(void)simplex_sparse_vector_set(core, position, value);
		else
			(void)simplex_sparse_vector_set(
				lift, b->logical_position[row], -value);
	}
	if (b->core_size > 0) {
		result = simplex_klu_sparse_solve_packed(
			&b->sparse_klu, core, 0);
		if (!result)
			b->profile_klu_sparse_success++;
		else if (result == 1) {
			b->profile_klu_sparse_declined++;
			lp_simplex_memcpy(b->core_work, core->dense,
				(size_t)b->core_size * sizeof(double));
			if (!klu_solve(b->symbolic, b->numeric, b->factor_size,
				    1, b->core_work, &b->common))
				return -1;
			(void)simplex_sparse_vector_pack(
				core, b->core_work, 0.);
		} else
			return -1;
	}
	/* Preserve the original ascending-core-column accumulation order. */
	simplex_sparse_vector_sort(core);
	for (i = 0; i < core->count; i++) {
		int column = b->core_basis[core->index[i]];
		double value = core->value[i];
		(void)simplex_sparse_vector_set(
			lift, b->core_position[core->index[i]], value);
		if (value == 0.)
			continue;
		for (k = matrix->column_start[column];
			k < matrix->column_start[column + 1]; k++) {
			int position =
				b->logical_position[matrix->row_index[k]];
			if (position < 0)
				continue;
			if (lift->mark[position] == lift->generation)
				lift->dense[position] +=
					matrix->value[k] * value;
			else
				(void)simplex_sparse_vector_set(lift, position,
					matrix->value[k] * value);
		}
	}
	simplex_sparse_vector_compact(lift);
	simplex_sparse_vector_clear(vector);
	for (i = 0; i < lift->count; i++)
		(void)simplex_sparse_vector_set(
			vector, lift->index[i], lift->value[i]);
	b->compact_ftran_validation_countdown--;
	return 0;
}

static int simplex_basis_base_btran_packed_impl(
	struct simplex_Basis *basis, struct simplex_SparseVector *vector)
{
	struct simplex_BasisImpl *b = basis->impl;
	const struct simplex_CscMatrix *matrix = b->matrix;
	struct simplex_SparseVector *core = &b->sparse_core,
				    *lift = &b->sparse_lift;
	size_t row_work = 0;
	int i, k, result;
	if (simplex_basis_uses_internal(basis))
		return 1;
	if (!b->compact_active) {
		result = simplex_klu_sparse_solve_packed(
			&b->sparse_klu, vector, 1);
		if (!result)
			b->profile_klu_sparse_success++;
		else if (result == 1)
			b->profile_klu_sparse_declined++;
		return result;
	}
	/* Preserve the existing periodic full residual audit and refinement.
	 * Only the intervening unvalidated lifts use this packed representation. */
	if (b->compact_btran_validation_countdown <= 0 || !lift->index ||
		core->dimension != b->core_size || b->column_to_core == NULL)
		return 1;
	simplex_sparse_vector_clear(core);
	simplex_sparse_vector_clear(lift);
	for (i = 0; i < vector->count; i++) {
		int column = b->base_index[vector->index[i]];
		double value = vector->value[i];
		if (value == 0.)
			continue;
		if (column >= b->structural_columns) {
			int row = column - b->structural_columns;
			row_work += matrix->row_start[row + 1] -
				    matrix->row_start[row];
			(void)simplex_sparse_vector_set(lift, row, -value);
		} else {
			int position = b->column_to_core[column];
			if (position < 0)
				return 1;
			(void)simplex_sparse_vector_set(core, position, value);
		}
	}
	if (row_work > (size_t)b->core_nonzeros)
		return 1;
	simplex_sparse_vector_sort(lift);
	for (i = 0; i < lift->count; i++) {
		int row = lift->index[i];
		double value = lift->value[i];
		for (k = matrix->row_start[row]; k < matrix->row_start[row + 1];
			k++) {
			int column = matrix->column_index[k];
			int position = column < b->structural_columns
					       ? b->column_to_core[column]
					       : -1;
			if (position < 0)
				continue;
			if (core->mark[position] == core->generation)
				core->dense[position] -=
					matrix->row_value[k] * value;
			else
				(void)simplex_sparse_vector_set(core, position,
					-matrix->row_value[k] * value);
		}
	}
	for (i = 0; i < core->count; i++)
		core->value[i] = core->dense[core->index[i]];
	if (b->core_size > 0) {
		result = simplex_klu_sparse_solve_packed(
			&b->sparse_klu, core, 1);
		if (!result)
			b->profile_klu_sparse_success++;
		else if (result == 1) {
			b->profile_klu_sparse_declined++;
			lp_simplex_memcpy(b->core_work, core->dense,
				(size_t)b->core_size * sizeof(double));
			if (!klu_tsolve(b->symbolic, b->numeric, b->factor_size,
				    1, b->core_work, &b->common))
				return -1;
			(void)simplex_sparse_vector_pack(
				core, b->core_work, 0.);
		} else
			return -1;
	}
	for (i = 0; i < core->count; i++)
		if (core->value[i] != 0.)
			(void)simplex_sparse_vector_set(lift,
				b->core_row[core->index[i]], core->value[i]);
	simplex_sparse_vector_clear(vector);
	for (i = 0; i < lift->count; i++)
		(void)simplex_sparse_vector_set(
			vector, lift->index[i], lift->value[i]);
	b->compact_btran_validation_countdown--;
	return 0;
}
#endif

/* Eta applications consume dense values and membership only. Materialize
 * packed values once at the end instead of duplicating every floating store. */
static void basis_packed_assign(
	struct simplex_SparseVector *vector, int row, double value)
{
	if (vector->mark[row] != vector->generation) {
		if (value == 0.)
			return;
		vector->mark[row] = vector->generation;
		vector->slot[row] = vector->count;
		vector->index[vector->count++] = row;
	}
	vector->dense[row] = value;
}

int simplex_basis_ftran_packed(
	const struct simplex_Basis *basis, struct simplex_SparseVector *vector)
{
	struct simplex_Basis *mutable = (struct simplex_Basis *)basis;
	struct simplex_BasisImpl *b = mutable->impl;
	int update, i, result = 1;
	clock_t started = 0;
	if (b->profile_enabled)
		started = clock();
	/* Internal FT retains its established raw kernel. */
	if (simplex_basis_uses_internal(basis))
		goto DENSE;
#ifdef LP_SIMPLEX_HAVE_KLU
	result = simplex_basis_base_ftran_packed_impl(mutable, vector);
#endif
	if (result == 1) {
		lp_simplex_memcpy(b->certify_correction, vector->dense,
			(size_t)b->rows * sizeof(double));
		result = simplex_basis_base_ftran(basis, b->certify_correction);
		if (!result)
			(void)simplex_sparse_vector_pack(
				vector, b->certify_correction, 0.);
	}
	if (result)
		goto FINISH;
	for (update = 0; update < b->update_count; update++) {
		int p = b->eta_pivot[update];
		double value, pivot = b->eta_pivot_value[update];
		b->eta_apply_work += b->eta_dense_slot[update] < 0
					     ? b->eta_start[update + 1] -
						       b->eta_start[update]
					     : b->rows;
		if (fabs(pivot) <= SIMPLEX_BASIS_ABSOLUTE_PIVOT_FLOOR) {
			result = -1;
			goto FINISH;
		}
		if (vector->dense[p] == 0.)
			continue;
		value = vector->dense[p] / pivot;
		if (b->eta_dense_slot[update] < 0) {
			for (i = b->eta_start[update];
				i < b->eta_start[update + 1]; i++) {
				int row = b->eta_index[i];
				double next = vector->dense[row] -
					      value * b->eta_value[i];
				basis_packed_assign(vector, row, next);
			}
		} else {
			const double *eta =
				b->dense_eta +
				(size_t)b->eta_dense_slot[update] * b->rows;
			for (i = 0; i < b->rows; i++)
				if (eta[i] != 0.)
					basis_packed_assign(vector, i,
						vector->dense[i] -
							value * eta[i]);
		}
		(void)simplex_sparse_vector_set(vector, p, value);
	}
	simplex_sparse_vector_compact(vector);
	simplex_sparse_vector_sort(vector);
	goto FINISH;
DENSE:
	lp_simplex_memcpy(b->certify_correction, vector->dense,
		(size_t)b->rows * sizeof(double));
	result = simplex_basis_ftran_raw(basis, b->certify_correction);
	if (!result)
		(void)simplex_sparse_vector_pack(
			vector, b->certify_correction, 0.);
FINISH:
	if (b->profile_enabled)
		b->profile_ftran_seconds +=
			(double)(clock() - started) / CLOCKS_PER_SEC;
	b->profile_ftran_calls++;
	return result;
}

/* Eta transpose changes only its pivot coordinate. Keep the existing ordered
 * support ordered with one insertion, so dense eta dots can visit only live
 * RHS coordinates in row order. No second dense index owner is needed. */
static void basis_packed_pivot(
	struct simplex_SparseVector *vector, int position, double value)
{
	int previous = vector->count, i;
	(void)simplex_sparse_vector_set(vector, position, value);
	if (vector->count == previous)
		return;
	i = previous;
	while (i > 0 && vector->index[i - 1] > position) {
		int row = vector->index[i - 1];
		vector->index[i] = row;
		vector->value[i] = vector->value[i - 1];
		vector->slot[row] = i;
		i--;
	}
	vector->index[i] = position;
	vector->value[i] = value;
	vector->slot[position] = i;
}

int simplex_basis_btran_packed(
	const struct simplex_Basis *basis, struct simplex_SparseVector *vector)
{
	struct simplex_Basis *mutable = (struct simplex_Basis *)basis;
	int update, first_update = 0, result = 1;
	clock_t started = 0;
	if (basis->impl->profile_enabled)
		started = clock();
	simplex_sparse_vector_sort(vector);
	if (simplex_basis_uses_internal(basis) &&
		!basis->impl->compact_active &&
		simplex_sparse_lu_ft_active(&basis->impl->sparse))
		first_update = basis->impl->sparse.ft_update_count;
	for (update = basis->impl->update_count - 1; update >= first_update;
		update--) {
		int position = basis->impl->eta_pivot[update];
		mutable->impl->eta_apply_work +=
			basis->impl->eta_dense_slot[update] < 0
				? basis->impl->eta_start[update + 1] -
					  basis->impl->eta_start[update]
				: basis->impl->rows;
		if (basis->impl->eta_dense_slot[update] >= 0 &&
			vector->count < basis->impl->rows / 4) {
			const double *eta =
				basis->impl->dense_eta +
				(size_t)basis->impl->eta_dense_slot[update] *
					basis->impl->rows;
			double pivot = basis->impl->eta_pivot_value[update];
			double old = vector->dense[position], dot = 0.;
			int i;
			if (fabs(pivot) <= SIMPLEX_BASIS_ABSOLUTE_PIVOT_FLOOR) {
				result = -1;
				goto FINISH;
			}
			for (i = 0; i < vector->count; i++) {
				int row = vector->index[i];
				dot += eta[row] * vector->dense[row];
			}
			basis_packed_pivot(vector, position,
				(old - dot + pivot * old) / pivot);
		} else if (basis->impl->eta_dense_slot[update] < 0 &&
			   basis->impl->eta_ordered[update]) {
			int first = basis->impl->eta_start[update];
			int end = basis->impl->eta_start[update + 1];
			int length = end - first, depth = 1, scale = length;
			double pivot = basis->impl->eta_pivot_value[update];
			double old = vector->dense[position], dot = 0.;
			int i;
			while (scale > 1) {
				depth++;
				scale /= 2;
			}
			if (vector->count >= length / depth) {
				if (simplex_basis_apply_eta_transpose(
					    basis, update, vector->dense)) {
					result = -1;
					goto FINISH;
				}
				basis_packed_pivot(vector, position,
					vector->dense[position]);
				continue;
			}
			if (fabs(pivot) <= SIMPLEX_BASIS_ABSOLUTE_PIVOT_FLOOR) {
				result = -1;
				goto FINISH;
			}
			/* Both supports are ordered. Search the larger eta for each RHS
			 * coordinate, then accumulate matches in the original eta order. */
			for (i = 0; i < vector->count; i++) {
				int row = vector->index[i], lo = first,
				    hi = end;
				if (vector->dense[row] == 0.)
					continue;
				while (lo < hi) {
					int middle = lo + (hi - lo) / 2;
					if (basis->impl->eta_index[middle] <
						row)
						lo = middle + 1;
					else
						hi = middle;
				}
				if (lo < end &&
					basis->impl->eta_index[lo] == row)
					dot += basis->impl->eta_value[lo] *
					       vector->dense[row];
			}
			basis_packed_pivot(vector, position,
				(old - dot + pivot * old) / pivot);
		} else {
			if (simplex_basis_apply_eta_transpose(
				    basis, update, vector->dense)) {
				result = -1;
				goto FINISH;
			}
			basis_packed_pivot(
				vector, position, vector->dense[position]);
		}
	}
#ifdef LP_SIMPLEX_HAVE_KLU
	result = simplex_basis_base_btran_packed_impl(mutable, vector);
#endif
	if (result == 1) {
		lp_simplex_memcpy(mutable->impl->certify_correction,
			vector->dense,
			(size_t)basis->impl->rows * sizeof(double));
		result = simplex_basis_base_btran(
			basis, mutable->impl->certify_correction);
		if (!result)
			(void)simplex_sparse_vector_pack(
				vector, mutable->impl->certify_correction, 0.);
	}
	if (!result)
		simplex_sparse_vector_sort(vector);
FINISH:
	if (basis->impl->profile_enabled)
		mutable->impl->profile_btran_seconds +=
			(double)(clock() - started) / CLOCKS_PER_SEC;
	mutable->impl->profile_btran_calls++;
	return result;
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
			residual[variable - basis->impl->structural_columns] +=
				value;
			denominator[variable -
				    basis->impl->structural_columns] +=
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
	for (i = 0; i < basis->impl->rows; i++) {
		if (!isfinite(residual[i]) || !isfinite(denominator[i]))
			return HUGE_VAL;
		maximum = __lp_simplex_MAX__(maximum,
			__lp_simplex_ABS__(residual[i]) /
				__lp_simplex_MAX__(denominator[i], DBL_MIN));
	}
	return maximum;
}

static double simplex_basis_current_btran_residual(
	const struct simplex_Basis *basis, const double *right,
	const double *solution, double *residual, double *denominator)
{
	const struct simplex_CscMatrix *matrix = basis->impl->matrix;
	double maximum = 0.;
	int i, k;

	for (i = 0; i < basis->impl->rows; i++) {
		int variable = basis->impl->index[i];
		long double sum = right[i],
			    scale = fabsl((long double)right[i]);

		if (variable >= basis->impl->structural_columns) {
			double value =
				solution[variable -
					 basis->impl->structural_columns];

			sum += value;
			scale += __lp_simplex_ABS__(value);
		} else
			for (k = matrix->column_start[variable];
				k < matrix->column_start[variable + 1]; k++) {
				long double term =
					(long double)matrix->value[k] *
					solution[matrix->row_index[k]];

				sum -= term;
				scale += fabsl(term);
			}
		residual[i] = (double)sum;
		denominator[i] = (double)scale;
		if (!isfinite(residual[i]) || !isfinite(denominator[i]))
			return HUGE_VAL;
		maximum = __lp_simplex_MAX__(
			maximum, (double)(fabsl(sum) /
					  (scale > DBL_MIN ? scale : DBL_MIN)));
	}
	return maximum;
}

static int simplex_basis_refine(
	const struct simplex_Basis *basis, double *vector, const int transpose)
{
	struct simplex_Basis *mutable = (struct simplex_Basis *)basis;
	int (*solve)(const struct simplex_Basis *, double *) =
		transpose ? simplex_basis_btran_raw : simplex_basis_ftran_raw;
	double (*residual)(const struct simplex_Basis *, const double *,
		const double *, double *, double *) =
		transpose ? simplex_basis_current_btran_residual
			  : simplex_basis_current_ftran_residual;
	double error, improved;
	double threshold = 16. * DBL_EPSILON * (double)(basis->impl->rows + 1);
	int i, refinement;

	lp_simplex_memcpy(mutable->impl->certify_right, vector,
		(size_t)basis->impl->rows * sizeof(double));
	if (solve(basis, vector) == lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	error = residual(basis, mutable->impl->certify_right, vector,
		mutable->impl->certify_residual, mutable->impl->core_work);
	if (!isfinite(error))
		return lp_simplex_EXIT_FAILURE;
	for (refinement = 0; refinement < SIMPLEX_BASIS_TERMINAL_REFINEMENTS &&
			     error > threshold;
		refinement++) {
		lp_simplex_memcpy(mutable->impl->certify_correction,
			mutable->impl->certify_residual,
			(size_t)basis->impl->rows * sizeof(double));
		if (solve(basis, mutable->impl->certify_correction) ==
			lp_simplex_EXIT_FAILURE)
			return lp_simplex_EXIT_FAILURE;
		for (i = 0; i < basis->impl->rows; i++)
			vector[i] += mutable->impl->certify_correction[i];
		improved = residual(basis, mutable->impl->certify_right, vector,
			mutable->impl->certify_residual,
			mutable->impl->core_work);
		if (!isfinite(improved))
			return lp_simplex_EXIT_FAILURE;
		if (improved >= error) {
			for (i = 0; i < basis->impl->rows; i++)
				vector[i] -=
					mutable->impl->certify_correction[i];
			break;
		}
		error = improved;
	}
	return lp_simplex_EXIT_SUCCESS;
}

int simplex_basis_ftran_refined(
	const struct simplex_Basis *basis, double *vector)
{
	return simplex_basis_refine(basis, vector, 0);
}

int simplex_basis_btran_refined(
	const struct simplex_Basis *basis, double *vector)
{
	return simplex_basis_refine(basis, vector, 1);
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
		if (simplex_basis_uses_internal(basis) &&
			!basis->impl->compact_active &&
			simplex_sparse_lu_ft_active(&basis->impl->sparse))
			first_update = basis->impl->sparse.ft_update_count;
		for (update = first_update; update < basis->impl->update_count;
			update++) {
			mutable->impl->eta_apply_work +=
				2L *
				(basis->impl->eta_dense_slot[update] < 0
						? basis->impl->eta_start
								  [update + 1] -
							  basis->impl->eta_start
								  [update]
						: basis->impl->rows);
			if (simplex_basis_apply_eta_pair(basis, update, first,
				    second) == lp_simplex_EXIT_FAILURE) {
				result = lp_simplex_EXIT_FAILURE;
				break;
			}
		}
	}
	if (basis->impl->profile_enabled) {
		mutable->impl->profile_ftran_seconds +=
			(double)(clock() - started) / (double)CLOCKS_PER_SEC;
	}
	mutable->impl->profile_ftran_calls += 2;
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
	}
	mutable->impl->profile_ftran_calls++;
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
	}
	mutable->impl->profile_btran_calls++;
	return result;
}

/* Keep the ordinary growth sequence, but never evaluate an overflowing int.
 * The exact request is sufficient when geometric growth is unrepresentable. */
static int basis_eta_capacity(int capacity, int required, int dense)
{
	if (capacity < 0 || required < 0)
		return -1;
	while (capacity < required) {
		int increment = capacity / 2 + dense;

		if (!increment || capacity > INT_MAX - increment)
			return required;
		capacity += increment;
	}
	return capacity;
}

static int basis_eta_bytes(int capacity, int rows, size_t limit, size_t *bytes)
{
	size_t column;
	*bytes = 0;
	return capacity >= 0 && rows >= 0 &&
	       simplex_storage_bytes(
		       (size_t)rows, sizeof(double), limit, &column) &&
	       simplex_storage_bytes((size_t)capacity, column, limit, bytes);
}

static int simplex_basis_reserve_eta(
	struct simplex_Basis *basis, const int required)
{
	int capacity;
	double *value;
	int *index;
	size_t value_bytes, index_bytes;

	if (required < 0)
		return lp_simplex_EXIT_FAILURE;
	if (required <= basis->impl->eta_capacity)
		return lp_simplex_EXIT_SUCCESS;
	capacity = basis_eta_capacity(basis->impl->eta_capacity, required, 0);
	if (capacity < 0 ||
		!simplex_storage_bytes((size_t)capacity, sizeof(double),
			(size_t)-1, &value_bytes) ||
		!simplex_storage_bytes((size_t)capacity, sizeof(int),
			(size_t)-1, &index_bytes))
		return lp_simplex_EXIT_FAILURE;
	value = (double *)lp_simplex_realloc(
		basis->impl->eta_value, value_bytes);
	if (value == NULL)
		return lp_simplex_EXIT_FAILURE;
	basis->impl->eta_value = value;
	/* A second allocation failure keeps the old usable capacity and both
	 * prefixes. The enlarged value allocation is owned and safely reusable. */
	index = (int *)lp_simplex_realloc(basis->impl->eta_index, index_bytes);
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
	size_t bytes;

	if (required_columns < 0)
		return lp_simplex_EXIT_FAILURE;
	if (required_columns <= basis->impl->dense_eta_capacity)
		return lp_simplex_EXIT_SUCCESS;
	capacity = basis->impl->dense_eta_capacity > 0
			   ? basis->impl->dense_eta_capacity
			   : 1;
	capacity = basis_eta_capacity(capacity, required_columns, 1);
	if (!basis_eta_bytes(capacity, basis->impl->rows, (size_t)-1, &bytes))
		return lp_simplex_EXIT_FAILURE;
	eta = (double *)lp_simplex_realloc(basis->impl->dense_eta, bytes);
	if (eta == NULL)
		return lp_simplex_EXIT_FAILURE;
	basis->impl->dense_eta = eta;
	basis->impl->dense_eta_capacity = capacity;
	return lp_simplex_EXIT_SUCCESS;
}

static int simplex_basis_update_impl(struct simplex_Basis *basis,
	const int leaving_position, const double *direction,
	const int *packed_index, int nonzeros)
{
	int i, count, first;
	int packed = packed_index != NULL;
	int sparse;
	double direction_maximum = 0.;
	double relative_pivot;
	basis->impl->storage_failed = 0;
	/* Large internal sparse refactorizations retain their hard-limit policy.
	 * KLU has its own factor-work estimate: do not suppress its deterministic
	 * amortization rule merely because the original basis has many rows. */
	if (basis->impl->update_count >= basis->impl->update_limit)
		return 1;
	if ((basis->impl->rows < SIMPLEX_BASIS_LARGE_ROWS ||
		    !simplex_basis_uses_internal(basis)) &&
		basis->impl->update_count > 0 && basis->impl->factor_work > 0 &&
		basis->impl->eta_apply_work /
				(basis->impl->allow_sparse_eta ? 1 : 2) >=
			basis->impl->factor_work) {
		basis->impl->profile_fill_reinversions++;
		return 1;
	}
	if (__lp_simplex_ABS__(direction[leaving_position]) <=
		SIMPLEX_BASIS_ABSOLUTE_PIVOT_FLOOR)
		return lp_simplex_EXIT_FAILURE;
	first = basis->impl->eta_start[basis->impl->update_count];
	if (packed) {
		if (nonzeros <= 0 || nonzeros > basis->impl->rows)
			return lp_simplex_EXIT_FAILURE;
		for (i = 0; i < nonzeros; i++) {
			int row = packed_index[i];

			if (row < 0 || row >= basis->impl->rows ||
				direction[row] == 0.)
				return lp_simplex_EXIT_FAILURE;
			direction_maximum =
				__lp_simplex_MAX__(direction_maximum,
					__lp_simplex_ABS__(direction[row]));
		}
	} else {
		nonzeros = 0;
		for (i = 0; i < basis->impl->rows; i++)
			if (direction[i] != 0.) {
				direction_maximum = __lp_simplex_MAX__(
					direction_maximum,
					__lp_simplex_ABS__(direction[i]));
				nonzeros++;
			}
	}
	/* The accepted nonzero leaving pivot proves direction_maximum > 0. */
	relative_pivot = __lp_simplex_ABS__(direction[leaving_position]) /
			 direction_maximum;
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
	if (basis->impl->allow_ft && simplex_basis_uses_internal(basis) &&
		!basis->impl->compact_active &&
		simplex_sparse_lu_ft_active(&basis->impl->sparse) &&
		basis->impl->update_count ==
			basis->impl->sparse.ft_update_count &&
		basis->impl->rows < SIMPLEX_BASIS_LARGE_ROWS) {
		int result = simplex_sparse_lu_ft_update(
			&basis->impl->sparse, leaving_position, direction);

		if (result == lp_simplex_EXIT_FAILURE &&
			basis->impl->sparse.storage_failed)
			basis->impl->storage_failed = 1;
		if (result == lp_simplex_EXIT_SUCCESS) {
			basis->impl->update_count++;
			basis->impl->eta_start[basis->impl->update_count] =
				basis->impl
					->eta_start[basis->impl->update_count -
						    1];
			basis->impl->profile_eta_nonzeros += nonzeros;
			basis->impl->profile_eta_slots += basis->impl->rows;
		}
		return result;
	}
	sparse = basis->impl->compact_active ||
		 (basis->impl->allow_sparse_eta &&
			 nonzeros <= basis->impl->rows / 8);
	if (sparse) {
		if (first < 0 || first > INT_MAX - nonzeros ||
			simplex_basis_reserve_eta(basis, first + nonzeros) ==
				lp_simplex_EXIT_FAILURE) {
			basis->impl->storage_failed = 1;
			return lp_simplex_EXIT_FAILURE;
		}
		count = first;
		if (packed)
			for (i = 0; i < nonzeros; i++) {
				int row = packed_index[i];

				basis->impl->eta_index[count] = row;
				basis->impl->eta_value[count++] =
					direction[row];
			}
		else
			for (i = 0; i < basis->impl->rows; i++)
				if (direction[i] != 0.) {
					basis->impl->eta_index[count] = i;
					basis->impl->eta_value[count++] =
						direction[i];
				}
		basis->impl->eta_dense_slot[basis->impl->update_count] = -1;
		basis->impl->eta_ordered[basis->impl->update_count] = 1;
		if (packed)
			for (i = first + 1; i < count; i++)
				if (basis->impl->eta_index[i - 1] >=
					basis->impl->eta_index[i]) {
					basis->impl->eta_ordered[basis->impl
							->update_count] = 0;
					break;
				}
		basis->impl->profile_eta_nonzeros += nonzeros;
		basis->impl->profile_eta_slots += basis->impl->rows;
	} else {
		int dense_slot = basis->impl->dense_eta_count;

		count = first;
		if (dense_slot < 0 || dense_slot == INT_MAX ||
			simplex_basis_reserve_dense_eta(basis,
				dense_slot + 1) == lp_simplex_EXIT_FAILURE) {
			basis->impl->storage_failed = 1;
			return lp_simplex_EXIT_FAILURE;
		}
		lp_simplex_memcpy(
			basis->impl->dense_eta +
				(size_t)dense_slot * basis->impl->rows,
			direction, (size_t)basis->impl->rows * sizeof(double));
		basis->impl->eta_dense_slot[basis->impl->update_count] =
			dense_slot;
		basis->impl->dense_eta_count++;
	}
	basis->impl->eta_pivot[basis->impl->update_count] = leaving_position;
	basis->impl->eta_pivot_value[basis->impl->update_count] =
		direction[leaving_position];
	basis->impl->update_count++;
	basis->impl->eta_start[basis->impl->update_count] = count;
	return lp_simplex_EXIT_SUCCESS;
}

int simplex_basis_storage_failed(const struct simplex_Basis *basis)
{
	return basis->impl->storage_failed;
}

int simplex_basis_update(struct simplex_Basis *basis,
	const int leaving_position, const double *direction)
{
	return simplex_basis_update_impl(
		basis, leaving_position, direction, NULL, 0);
}

int simplex_basis_update_packed(struct simplex_Basis *basis,
	const int leaving_position, const double *direction, const int *index,
	const int nonzeros)
{
	return simplex_basis_update_impl(
		basis, leaving_position, direction, index, nonzeros);
}

void simplex_basis_set_sparse_eta(
	struct simplex_Basis *basis, const int allowed)
{
	basis->impl->allow_sparse_eta = allowed;
}

void simplex_basis_use_product_updates(struct simplex_Basis *basis)
{
	/* Existing FT updates remain part of the base solve.  Subsequent pivots
	 * are kept in the independent eta file, including after reinversion. */
	basis->impl->allow_ft = 0;
}

void simplex_basis_use_default_updates(struct simplex_Basis *basis)
{
	basis->impl->allow_ft = 1;
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
	const struct simplex_Basis *basis, struct simplex_BasisProfile *profile)
{
	lp_simplex_memset(profile, 0, sizeof(*profile));

	profile->factor_size = basis->impl->factor_size;
	profile->factor_seconds = basis->impl->profile_factor_seconds;
	profile->ftran_seconds = basis->impl->profile_ftran_seconds;
	profile->btran_seconds = basis->impl->profile_btran_seconds;
	profile->factor_calls = basis->impl->profile_factor_calls;
	profile->ftran_calls = basis->impl->profile_ftran_calls;
	profile->btran_calls = basis->impl->profile_btran_calls;
	profile->klu_sparse_success = basis->impl->profile_klu_sparse_success;
	profile->klu_sparse_declined = basis->impl->profile_klu_sparse_declined;
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
