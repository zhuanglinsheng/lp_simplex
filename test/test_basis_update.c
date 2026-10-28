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
#include "simplex_sparse_vector.h"

#include <lp_simplex/status.h>

#include <stdio.h>
#include <math.h>
#include <string.h>

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

/* Compact bases retain a permuted logical identity and a small structural
 * core. Check both RHS directions against independently formed B products,
 * including dense/tiny supports, eta updates, and a changed refactorization. */
static void compact_product(const struct simplex_CscMatrix *matrix,
	const int *index, const double *x, double *right, int transpose)
{
	int i, j, k;
	memset(right, 0, 4096 * sizeof(double));
	for (i = 0; i < 4096; i++) {
		j = index[i];
		if (j >= 4) {
			if (transpose)
				right[i] = -x[j - 4];
			else
				right[j - 4] -= x[i];
		} else
			for (k = matrix->column_start[j];
				k < matrix->column_start[j + 1]; k++) {
				int row = matrix->row_index[k];
				if (transpose)
					right[i] += matrix->value[k] * x[row];
				else
					right[row] += matrix->value[k] * x[i];
			}
	}
}

static int test_compact_lift(int backend, int with_rows)
{
	int starts[5] = {0, 4, 8, 11, 15};
	int indices[15] = {1, 17, 23, 4000, 5, 23, 3000, 4095, 17, 3000, 4000,
		1, 17, 23, 3000};
	double values[15] = {.25, 2., .5, -1., .75, 3., -.25, .5, 1., 4., .125,
		.5, 3., -.5, .25};
	int rows[4097], columns[15], cursor[4096], index[4096];
	double row_values[15], right[4096], second[4096];
	double expected[4096], direction[4096];
	struct simplex_CscMatrix matrix;
	struct simplex_Basis basis;
	int i, j, k, pass, epoch, transpose, valid = 1;
	struct simplex_SparseVector packed;
	memset(&matrix, 0, sizeof(matrix));
	memset(rows, 0, sizeof(rows));
	for (k = 0; k < 15; k++)
		rows[indices[k] + 1]++;
	for (i = 0; i < 4096; i++) {
		rows[i + 1] += rows[i];
		cursor[i] = rows[i];
		index[i] = 4 + (i + 7) % 4096;
	}
	for (j = 0; j < 4; j++)
		for (k = starts[j]; k < starts[j + 1]; k++) {
			i = cursor[indices[k]]++;
			columns[i] = j;
			row_values[i] = values[k];
		}
	index[10] = 0;
	index[16] = 1;
	index[2993] = 2;
	matrix.rows = 4096;
	matrix.columns = 4;
	matrix.nonzeros = 15;
	matrix.column_start = starts;
	matrix.row_index = indices;
	matrix.value = values;
	matrix.row_start = with_rows ? rows : NULL;
	matrix.column_index = with_rows ? columns : NULL;
	matrix.row_value = with_rows ? row_values : NULL;
	if (simplex_basis_create(&basis, &matrix, 4, index))
		return 0;
	if (simplex_sparse_vector_create(&packed, 4096)) {
		simplex_basis_destroy(&basis);
		return 0;
	}
	simplex_basis_set_backend(&basis, backend);
	if (with_rows != 2)
		simplex_basis_request_compact(&basis);
	if (simplex_basis_factorize(&basis) ||
		simplex_basis_compact_active(&basis) != (with_rows != 2))
		valid = 0;
	for (epoch = 0; valid && epoch < 3; epoch++) {
		for (pass = 0; valid && pass < 40; pass++) {
			double scale = pass % 3 == 2 ? 1e-200 : 1.;
			for (transpose = 0; valid && transpose < 2;
				transpose++) {
				for (i = 0; i < 4096; i++) {
					double value =
						pass % 3 == 0
							? (double)((i % 7) - 3)
							: (i == 1 || i == 17 ||
								  i == 3000);
					expected[i] = scale * value;
				}
				compact_product(&matrix, index, expected, right,
					transpose);
				memcpy(second, right, sizeof(right));
				if (pass % 4 == 0)
					memset(second, 0, sizeof(second));
				if (!transpose && pass % 2 == 0) {
					(void)simplex_sparse_vector_pack(
						&packed, right, 0.);
					if (simplex_basis_ftran_packed(
						    &basis, &packed) ||
						simplex_basis_ftran(
							&basis, second))
						valid = 0;
					memcpy(right, packed.dense,
						sizeof(right));
				} else if (transpose && pass % 2 == 0) {
					(void)simplex_sparse_vector_pack(
						&packed, right, 0.);
					if (simplex_basis_btran_packed(
						    &basis, &packed))
						valid = 0;
					memcpy(right, packed.dense,
						sizeof(right));
					for (i = 1; i < packed.count; i++)
						if (packed.index[i - 1] >=
							packed.index[i])
							valid = 0;
				} else if (transpose ? simplex_basis_btran(
							       &basis, right)
						     : simplex_basis_ftran_pair(
							       &basis, right,
							       second))
					valid = 0;
				for (i = 0; valid && i < 4096; i++) {
					double target = pass % 4 == 0
								? 0.
								: expected[i];
					if (!isfinite(right[i]) ||
						!isfinite(second[i]) ||
						fabs((right[i] - expected[i]) /
							scale) > 1e-12 ||
						(!transpose &&
							fabs((second[i] -
								     target) /
								scale) > 1e-12))
						valid = 0;
				}
			}
		}
		if (valid && epoch == 0) {
			simplex_csc_column_to_dense(&matrix, 4, 3, direction);
			if (simplex_basis_ftran(&basis, direction))
				valid = 0;
			index[10] = 3;
			i = simplex_basis_update(&basis, 10, direction);
			if (i == lp_simplex_EXIT_FAILURE ||
				(i == 1 && simplex_basis_factorize(&basis)))
				valid = 0;
		} else if (valid && epoch == 1 &&
			   simplex_basis_factorize(&basis))
			valid = 0;
	}
	if (!valid)
		fprintf(stderr, "compact backend=%d CSR=%d epoch=%d\n", backend,
			with_rows, epoch);
	simplex_basis_destroy(&basis);
	simplex_sparse_vector_destroy(&packed);
	return valid;
}

static int test_mixed_eta_packed(int dense_update)
{
	enum {
		N = 4096
	};

	int start[3] = {0, N, N + 1}, row[N + 1], index[N];
	double values[N + 1], expected[N], right[N], dense[N], direction[N];
	struct simplex_CscMatrix matrix = {0};
	struct simplex_Basis basis;
	struct simplex_SparseVector packed;
	int i, j, k, epoch, pass, transpose, valid = 1;
	for (i = 0; i < N; i++) {
		row[i] = i;
		values[i] = i == 17			? 2.
			    : (dense_update || i < 256) ? (i % 7 - 3) * .125
							: 0.;
		index[i] = i + 2;
	}
	row[N] = 23;
	values[N] = 1.;
	matrix.rows = N;
	matrix.columns = 2;
	matrix.nonzeros = N + 1;
	matrix.column_start = start;
	matrix.row_index = row;
	matrix.value = values;
	if (simplex_basis_create(&basis, &matrix, 2, index))
		return 0;
	if (simplex_sparse_vector_create(&packed, N)) {
		simplex_basis_destroy(&basis);
		return 0;
	}
	simplex_basis_set_sparse_eta(&basis, 1);
	valid = !simplex_basis_factorize(&basis);
	for (epoch = 0; valid && epoch < 3; epoch++) {
		if (epoch < 2) {
			int pivot = epoch == 0 ? 17 : 23, result;
			simplex_csc_column_to_dense(
				&matrix, 2, epoch, direction);
			valid = !simplex_basis_ftran(&basis, direction);
			index[pivot] = epoch;
			result = simplex_basis_update(&basis, pivot, direction);
			if (result == 1)
				result = simplex_basis_factorize(&basis);
			if (result != 0)
				fprintf(stderr,
					"mixed eta update epoch=%d result=%d\n",
					epoch, result);
			valid = valid && result == 0;
		} else
			valid = !simplex_basis_factorize(&basis);
		for (pass = 0; valid && pass < 6; pass++)
			for (transpose = 0; valid && transpose < 2;
				transpose++) {
				double scale = pass >= 3 ? 1e-200 : 1.;
				memset(right, 0, sizeof(right));
				for (i = 0; i < N; i++)
					expected[i] =
						scale *
						(pass % 3 == 0 ? (i % 5 - 2)
							       : i == 17 ||
									 i == 23 ||
									 i == 251);
				for (i = 0; i < N; i++) {
					j = index[i];
					if (j >= 2) {
						if (transpose)
							right[i] = -expected[j -
									     2];
						else
							right[j - 2] -=
								expected[i];
					} else
						for (k = start[j];
							k < start[j + 1]; k++)
							if (transpose)
								right[i] +=
									values[k] *
									expected[row[k]];
							else
								right[row[k]] +=
									values[k] *
									expected[i];
				}
				memcpy(dense, right, sizeof(right));
				(void)simplex_sparse_vector_pack(
					&packed, right, 0.);
				valid = !(transpose ? simplex_basis_btran_packed(
							      &basis, &packed)
						    : simplex_basis_ftran_packed(
							      &basis,
							      &packed)) &&
					!(transpose ? simplex_basis_btran(
							      &basis, dense)
						    : simplex_basis_ftran(
							      &basis, dense));
				for (i = 0; valid && i < N; i++)
					valid = isfinite(packed.dense[i]) &&
						fabs((packed.dense[i] -
							     expected[i]) /
							scale) < 1e-12 &&
						fabs((dense[i] - expected[i]) /
							scale) < 1e-12;
				if (!valid)
					fprintf(stderr,
						"mixed eta epoch=%d pass=%d transpose=%d row=%d\n",
						epoch, pass, transpose, i);
			}
	}
	simplex_sparse_vector_destroy(&packed);
	simplex_basis_destroy(&basis);
	return valid;
}

/* Unit RHS at the last row propagates through a long reverse eta chain. Alternate packed eta order to exercise both ordered and raw dots. */
static int test_eta_chain(int amortized)
{
	enum {
		M = 4096,
		K = 129
	};

	int starts[K + 1], rows[2 * K], index[M], support[2];
	double values[2 * K], direction[M], expected[M], rhs[M];
	struct simplex_CscMatrix matrix = {0};
	struct simplex_Basis basis;
	struct simplex_SparseVector packed;
	int i, update, pass, valid = 1;
	for (i = 0; i < K; i++) {
		starts[i] = 2 * i;
		rows[2 * i] = i;
		rows[2 * i + 1] = i + 1;
		values[2 * i] = -1.;
		values[2 * i + 1] = -.25;
	}
	starts[K] = 2 * K;
	for (i = 0; i < M; i++)
		index[i] = K + i;
	matrix.rows = M;
	matrix.columns = K;
	matrix.nonzeros = 2 * K;
	matrix.column_start = starts;
	matrix.row_index = rows;
	matrix.value = values;
	if (simplex_basis_create(&basis, &matrix, K, index))
		return 0;
	/* The internal fixture deliberately keeps a long chain; the KLU fixture
	 * must honor cost-driven rebuilds even beyond the large-row boundary. */
	simplex_basis_set_backend(
		&basis, amortized ? SIMPLEX_BASIS_BACKEND_KLU
				  : SIMPLEX_BASIS_BACKEND_INTERNAL);
	if (amortized && simplex_basis_resolved_backend(&basis) !=
				 SIMPLEX_BASIS_BACKEND_KLU) {
		simplex_basis_destroy(&basis);
		return 1;
	}
	if (simplex_sparse_vector_create(&packed, M)) {
		simplex_basis_destroy(&basis);
		return 0;
	}
	simplex_basis_set_sparse_eta(&basis, 1);
	valid = !simplex_basis_factorize(&basis);
	for (update = 0; valid && update < K; update++) {
		int result, solve;
		for (solve = 0; amortized && valid && solve < 64; solve++) {
			memset(rhs, 0, sizeof(rhs));
			rhs[0] = 1.;
			valid = !simplex_basis_ftran(&basis, rhs);
		}
		simplex_csc_column_to_dense(&matrix, K, update, direction);
		valid = !simplex_basis_ftran(&basis, direction);
		index[update] = update;
		support[update % 2] = update;
		support[1 - update % 2] = update + 1;
		result = simplex_basis_update_packed(
			&basis, update, direction, support, 2);
		if (amortized && result == 1)
			result = simplex_basis_factorize(&basis);
		valid = valid && !result;
	}
	if (amortized) {
		struct simplex_BasisProfile profile;
		simplex_basis_get_profile(&basis, &profile);
		valid = valid && profile.fill_reinversions > 0;
	} else
		valid = valid && simplex_basis_update_count(&basis) == K;
	for (pass = 0; valid && pass < 2; pass++) {
		double scale = pass ? 1e-200 : 1.;
		memset(expected, 0, sizeof(expected));
		memset(rhs, 0, sizeof(rhs));
		expected[K] = -scale;
		rhs[K] = scale;
		for (i = K - 1; i >= 0; i--)
			expected[i] = -.25 * expected[i + 1];
		(void)simplex_sparse_vector_pack(&packed, rhs, 0.);
		valid = !simplex_basis_btran_packed(&basis, &packed);
		for (i = 0; valid && i < M; i++)
			valid = isfinite(packed.dense[i]) &&
				fabs((packed.dense[i] - expected[i]) / scale) <
					1e-12;
		/* Check the same RHS against the independent original bidiagonal basis. */
		for (i = 0; valid && i < K; i++)
			valid = fabs((-packed.dense[i] -
					     .25 * packed.dense[i + 1]) /
					scale) < 1e-12;
	}
	simplex_sparse_vector_destroy(&packed);
	simplex_basis_destroy(&basis);
	return valid;
}

int main(void)
{
	if (!test_eta_chain(0) || !test_eta_chain(1) ||
		!test_mixed_eta_packed(0) || !test_mixed_eta_packed(1) ||
		!test_compact_lift(SIMPLEX_BASIS_BACKEND_DEFAULT, 0) ||
		!test_compact_lift(SIMPLEX_BASIS_BACKEND_DEFAULT, 1) ||
		!test_compact_lift(SIMPLEX_BASIS_BACKEND_DEFAULT, 2) ||
		!test_compact_lift(SIMPLEX_BASIS_BACKEND_INTERNAL, 0) ||
		!test_compact_lift(SIMPLEX_BASIS_BACKEND_INTERNAL, 1) ||
		!test_compact_lift(SIMPLEX_BASIS_BACKEND_INTERNAL, 2)) {
		printf("compact basis lift regression failed\n");
		return 1;
	}
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
