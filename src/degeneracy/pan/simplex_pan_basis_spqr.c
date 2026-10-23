/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_pan_basis_spqr.c
 * @brief Sparse linear-algebra backend for Pan's dynamically deficient basis.
 *
 * The hot path maintains a sparse dynamic R factor and applies Q implicitly
 * through the normalized active columns.  SuiteSparseQR is reserved for rank-
 * revealing recovery/certification.
 */

#include "simplex_pan_basis.h"
#include "simplex_pan_config.h"

#ifdef LP_SIMPLEX_HAVE_SPQR

#include "pan_dynamic_sparse_qr.h"
#include "pan_multifrontal_qr.h"
#include "utils.h"

#include <lp_simplex/status.h>

#include <SuiteSparseQR_C.h>

#include <math.h>
#include <stdlib.h>
#include <float.h>
#include <time.h>
#include <stdlib.h>


struct pan_SparseBackend {
	cholmod_common common;
	cholmod_sparse *matrix;
	struct pan_DynamicSparseQR factor;
	struct pan_MultifrontalQR orthogonal;
	int matrix_dirty;
	double *condition_work;
	unsigned char *stability_mask;
	double *pricing_work;
	double *product_correction;
};

int simplex_pan_dual_newton_qr(const struct simplex_PanStandard *a,
	const double *point, const unsigned char *free_mask, const double *map,
	const double *multiplier, double scale, double *correction)
{
	cholmod_common common;
	cholmod_sparse *matrix = NULL, *transpose = NULL;
	cholmod_dense *rhs = NULL, *solution = NULL, *coefficient = NULL;
	int *index = NULL, *cursor = NULL, *p, *row;
	double *value, *right;
	int i, j, k, count = 0, nonzeros = 0, result = lp_simplex_EXIT_FAILURE;

	int32_t rank = -1;
	clock_t start = clock();
	cholmod_start(&common);

	common.print = 0;
	index = (int *)lp_simplex_malloc((size_t)a->columns * sizeof(int));
	cursor = (int *)lp_simplex_malloc((size_t)a->rows * sizeof(int));
	if (!index || !cursor)
		goto done;
	for (j = 0; j < a->columns; j++) {
		double norm = a->column_norm[j], r;

		index[j] = -1;
		if (norm == 0.)
			continue;
		r = (a->objective[j] / scale) / norm;
		for (k = a->column_start[j]; k < a->column_start[j + 1]; k++)
			r += (a->value[k] / norm) * map[a->row_index[k]] *
			     multiplier[a->row_index[k]];
		if ((free_mask ? free_mask[j] : point[j] > 0.) || r < 0.) {
			index[j] = count++;
			nonzeros += a->column_start[j + 1] - a->column_start[j];
		}
	}
	if (count == 0)
		goto done;
	matrix = cholmod_allocate_sparse((size_t)count, (size_t)a->rows,
		(size_t)nonzeros, 1, 1, 0, CHOLMOD_REAL, &common);
	rhs = cholmod_allocate_dense(
		(size_t)count, 1, (size_t)count, CHOLMOD_REAL, &common);
	if (!matrix || !rhs)
		goto done;
	p = (int *)matrix->p;
	row = (int *)matrix->i;
	value = (double *)matrix->x;
	right = (double *)rhs->x;
	lp_simplex_memset(p, 0, ((size_t)a->rows + 1) * sizeof(int));
	for (j = 0; j < a->columns; j++)
		if (index[j] >= 0)
			for (k = a->column_start[j]; k < a->column_start[j + 1];
				k++)
				p[a->row_index[k] + 1]++;
	for (i = 0; i < a->rows; i++) {
		p[i + 1] += p[i];
		cursor[i] = p[i];
	}
	for (j = 0; j < a->columns; j++)
		if (index[j] >= 0) {
			double norm = a->column_norm[j],
			       r = (a->objective[j] / scale) / norm;

			for (k = a->column_start[j]; k < a->column_start[j + 1];
				k++) {
				int position = cursor[a->row_index[k]]++;
				double v = (a->value[k] / norm) *
					   map[a->row_index[k]];

				row[position] = index[j];
				value[position] = v;
				r += v * multiplier[a->row_index[k]];
			}
			right[index[j]] = -r;
		}
	rank = SuiteSparseQR_i_C(SPQR_ORDERING_COLAMD, 1e-12, (int32_t)a->rows,
		2, matrix, NULL, rhs, NULL, &solution, NULL, NULL, NULL, NULL,
		NULL, &common);
	if (rank < 0 || !solution || solution->nrow < (size_t)a->rows)
		goto done;
	/* The rank-deficient backslash returns a basic solution. Its arbitrary
	 * null-space component can unnecessarily violate inactive inequalities.
	 * Project it onto range(U) with a second rectangular QR to obtain the
	 * minimum-norm correction, without forming U U^T. */
	transpose = cholmod_transpose(matrix, 1, &common);
	if (!transpose)
		goto done;
	if (SuiteSparseQR_i_C(SPQR_ORDERING_COLAMD, 1e-12, (int32_t)count, 2,
		    transpose, NULL, solution, NULL, &coefficient, NULL, NULL,
		    NULL, NULL, NULL, &common) < 0 ||
		!coefficient || coefficient->nrow < (size_t)count)
		goto done;
	lp_simplex_memset(correction, 0, (size_t)a->rows * sizeof(double));
	p = (int *)transpose->p;
	row = (int *)transpose->i;
	value = (double *)transpose->x;
	for (j = 0; j < count; j++)
		for (k = p[j]; k < p[j + 1]; k++)
			correction[row[k]] +=
				value[k] * ((double *)coefficient->x)[j];
	for (i = 0; i < a->rows; i++) {
		if (!isfinite(correction[i]))
			goto done;
	}
	result = lp_simplex_EXIT_SUCCESS;
done:
	if (getenv("LP_SIMPLEX_PROFILE"))
		fprintf(stderr,
			"pan-bda dual newton QR: active=%d rows=%d rank=%d seconds=%.6f result=%d\n",
			count, a->rows, (int)rank,
			(double)(clock() - start) / CLOCKS_PER_SEC, result);
	cholmod_free_dense(&coefficient, &common);
	cholmod_free_dense(&solution, &common);
	cholmod_free_dense(&rhs, &common);
	cholmod_free_sparse(&transpose, &common);
	cholmod_free_sparse(&matrix, &common);
	cholmod_finish(&common);

	lp_simplex_free(index);
	lp_simplex_free(cursor);
	return result;
}


static cholmod_sparse *pan_sparse_basis_matrix(struct simplex_PanBasis *basis);


/* A fresh multifrontal factor may reorder columns to reduce sparse fill.
 * Make that order explicit in every active-set map, never only inside R. */
static void pan_accept_factor_order(struct simplex_PanBasis *basis)
{
	struct pan_SparseBackend *backend =
		(struct pan_SparseBackend *)basis->sparse_backend;
	int i;

	for (i = 0; i < basis->count; i++) {
		int column = backend->factor.tag[i];

		basis->column[i] = column;
		basis->position[column] = i;
		basis->column_norm[i] =
			simplex_pan_column_norm(basis->standard, column);
	}
}


static void pan_sync_pack_statistics(struct simplex_PanBasis *basis)
{
	struct pan_SparseBackend *backend =
		(struct pan_SparseBackend *)basis->sparse_backend;

	basis->numeric_packs = backend->factor.numeric_packs;
	basis->numeric_packed_nodes = backend->factor.numeric_packed_nodes;
	basis->border_solves = backend->factor.border_solves;
}


static void pan_sync_dynamic_statistics(struct simplex_PanBasis *basis)
{
	struct pan_SparseBackend *backend =
		(struct pan_SparseBackend *)basis->sparse_backend;

	pan_multifrontal_qr_refresh_structure(
		&backend->orthogonal, &backend->factor);
	basis->extensions = backend->factor.extensions;
	basis->downdates = backend->factor.downdates;
	basis->rotations = backend->factor.rotations;
	basis->symbolic_updates = backend->factor.symbolic_updates;
	basis->anchor_rebuilds = backend->orthogonal.anchor_rebuilds;
	basis->local_retriangularizations =
		backend->orthogonal.local_retriangularizations;
	basis->block_householders = backend->orthogonal.block_householders;
	basis->maximum_front = backend->orthogonal.maximum_front;
	basis->factor_nonzeros = backend->orthogonal.factor_nonzeros;
	basis->orthogonal_error = backend->orthogonal.orthogonal_error;
	basis->backward_error = backend->orthogonal.backward_error;
	pan_sync_pack_statistics(basis);
}


static void pan_record_backward_error(struct simplex_PanBasis *basis,
	const double residual, const double right_norm)
{
	struct pan_SparseBackend *backend =
		(struct pan_SparseBackend *)basis->sparse_backend;
	double error = residual / (1. + right_norm);

	backend->orthogonal.backward_error = error;
	/* A solve does not change R's sparsity, parent tree or front sizes.
	 * Refreshing all of them here added an O(nnz(R)) scan to every solve. */
	basis->backward_error = error;
}


static int pan_dynamic_factor_matches(const struct simplex_PanBasis *basis)
{
	const struct pan_SparseBackend *backend =
		(const struct pan_SparseBackend *)basis->sparse_backend;
	int i;

	if (backend->factor.count != basis->count)
		return 0;
	for (i = 0; i < basis->count; i++)
		if (backend->factor.tag[i] != basis->column[i])
			return 0;
	return 1;
}


static double pan_dense_norm(const double *vector, const int count)
{
	double scale = 0.;
	double sum = 1.;
	int i;

	for (i = 0; i < count; i++) {
		double value = vector[i] < 0. ? -vector[i] : vector[i];

		if (value == 0.)
			continue;
		if (scale < value) {
			double ratio = scale / value;

			sum = 1. + sum * ratio * ratio;
			scale = value;
		} else {
			double ratio = value / scale;

			sum += ratio * ratio;
		}
	}
	return scale == 0. ? 0. : scale * sqrt(sum);
}


double simplex_pan_column_norm(
	const struct simplex_PanStandard *standard, const int column)
{
	return standard->column_norm[column];
}


static void pan_column_axpy(const struct simplex_PanStandard *standard,
	const int column, const double scale, double *dense)
{
	int k;

	for (k = standard->column_start[column];
		k < standard->column_start[column + 1]; k++)
		dense[standard->row_index[k]] += scale * standard->value[k];
}


static double pan_column_dense_dot(const struct simplex_PanStandard *standard,
	const int column, const double *dense)
{
	return simplex_pan_column_dot(standard, column, dense);
}


static void pan_normalized_column_axpy(
	const struct simplex_PanStandard *standard, const int column,
	const double scale, double *dense)
{
	int k;

	for (k = standard->column_start[column];
		k < standard->column_start[column + 1]; k++)
		dense[standard->row_index[k]] +=
			scale * standard->normalized_value[k];
}


static void pan_basis_multiply(const struct simplex_PanBasis *basis,
	const double *coefficient, double *result)
{
	struct pan_SparseBackend *backend =
		(struct pan_SparseBackend *)basis->sparse_backend;
	double *correction = backend->product_correction;
	int i, k;

	lp_simplex_memset(
		result, 0, (size_t)basis->standard->rows * sizeof(double));
	lp_simplex_memset(
		correction, 0, (size_t)basis->standard->rows * sizeof(double));
	for (i = 0; i < basis->count; i++) {
		int column = basis->column[i];

		for (k = basis->standard->column_start[column];
			k < basis->standard->column_start[column + 1]; k++) {
			int row = basis->standard->row_index[k];
			double a = basis->standard->value[k],
			       p = a * coefficient[i], next = result[row] + p;

			correction[row] +=
				(fabs(result[row]) >= fabs(p)
						? (result[row] - next) + p
						: (p - next) + result[row]) +
				fma(a, coefficient[i], -p);
			result[row] = next;
		}
	}
	for (i = 0; i < basis->standard->rows; i++)
		result[i] += correction[i];
}


static int pan_ls_stationary(const struct simplex_PanBasis *basis,
	const double *residual, double right_norm)
{
	int i;

	if (basis->full_refinement)
		return 0;
	for (i = 0; i < basis->count; i++) {
		double gradient = simplex_pan_column_dot(basis->standard,
					  basis->column[i], residual) /
				  basis->column_norm[i];

		if (!isfinite(gradient) ||
			fabs(gradient) > 8. * DBL_EPSILON * right_norm)
			return 0;
	}
	return 1;
}


static int pan_implicit_q_transpose_multiply(
	const struct simplex_PanBasis *basis, const double *vector,
	double *result)
{
	struct pan_SparseBackend *backend =
		(struct pan_SparseBackend *)basis->sparse_backend;
	int i;

	if (pan_multifrontal_qr_apply_transpose(&backend->orthogonal, vector,
		    basis->right_work) == lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	for (i = 0; i < backend->factor.count; i++)
		result[i] = basis->right_work[i];
	return lp_simplex_EXIT_SUCCESS;
}


static int pan_implicit_q_multiply(const struct simplex_PanBasis *basis,
	const double *coefficient, double *solve_work, double *result)
{
	struct pan_SparseBackend *backend =
		(struct pan_SparseBackend *)basis->sparse_backend;

	lp_simplex_memset(
		solve_work, 0, (size_t)basis->standard->rows * sizeof(double));
	if (backend->factor.count > 0)
		lp_simplex_memcpy(solve_work, coefficient,
			(size_t)backend->factor.count * sizeof(double));
	return pan_multifrontal_qr_apply(
		&backend->orthogonal, solve_work, result);
}


/* Apply the accumulated sparse Householder/Givens product to the normalized
 * entering column.  The tail Householder is retained as part of Q, so the
 * bordered R column and its rank certificate come from a true QR update. */
static int pan_dynamic_append_column(
	struct simplex_PanBasis *basis, const int column)
{
	struct pan_SparseBackend *backend =
		(struct pan_SparseBackend *)basis->sparse_backend;
	const struct simplex_PanStandard *standard = basis->standard;
	int count = backend->factor.count;
	double diagonal;
	int transform_count = backend->orthogonal.transform_count;
	int value_count = backend->orthogonal.value_count;
	long long local_update_work = backend->orthogonal.local_update_work;
	long local_retriangularizations =
		backend->orthogonal.local_retriangularizations;
	long block_householders = backend->orthogonal.block_householders;

	lp_simplex_memset(
		basis->row_work, 0, (size_t)standard->rows * sizeof(double));
	pan_normalized_column_axpy(standard, column, 1., basis->row_work);
	if (pan_multifrontal_qr_append(&backend->orthogonal, basis->row_work,
		    count, basis->basis_work,
		    &diagonal) == lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	if (pan_dynamic_sparse_qr_append(&backend->factor, column,
		    basis->basis_work, diagonal,
		    basis->rank_tolerance) == lp_simplex_EXIT_FAILURE) {
		backend->orthogonal.transform_count = transform_count;
		backend->orthogonal.value_count = value_count;
		backend->orthogonal.local_update_work = local_update_work;
		backend->orthogonal.local_retriangularizations =
			local_retriangularizations;
		backend->orthogonal.block_householders = block_householders;
		return lp_simplex_EXIT_FAILURE;
	}
	return lp_simplex_EXIT_SUCCESS;
}


static int pan_basis_reserve(struct simplex_PanBasis *basis, const int capacity)
{
	double *norm;
	double *work;
	int *columns;

	if (capacity <= basis->capacity)
		return lp_simplex_EXIT_SUCCESS;
	columns = (int *)lp_simplex_malloc((size_t)capacity * sizeof(int));
	norm = (double *)lp_simplex_malloc((size_t)capacity * sizeof(double));
	work = (double *)lp_simplex_malloc((size_t)capacity * sizeof(double));
	if (columns == NULL || norm == NULL || work == NULL) {
		lp_simplex_free(columns);
		lp_simplex_free(norm);
		lp_simplex_free(work);
		return lp_simplex_EXIT_FAILURE;
	}
	if (basis->count > 0) {
		lp_simplex_memcpy(columns, basis->column,
			(size_t)basis->count * sizeof(int));
		lp_simplex_memcpy(norm, basis->column_norm,
			(size_t)basis->count * sizeof(double));
	}
	lp_simplex_free(basis->column);
	lp_simplex_free(basis->column_norm);
	lp_simplex_free(basis->basis_work);
	basis->column = columns;
	basis->column_norm = norm;
	basis->basis_work = work;
	basis->capacity = capacity;
	return lp_simplex_EXIT_SUCCESS;
}


int simplex_pan_basis_create(struct simplex_PanBasis *basis,
	const struct simplex_PanStandard *standard, const double rank_tolerance)
{
	struct simplex_PanConfig config;
	struct pan_SparseBackend *backend;

	lp_simplex_memset(basis, 0, sizeof(*basis));
	basis->standard = standard;
	simplex_pan_config_init(&config, NULL);
	simplex_pan_basis_configure(basis, &config);
	basis->rank_tolerance = rank_tolerance;
	basis->position = (int *)lp_simplex_malloc(
		(size_t)(standard->columns > 0 ? standard->columns : 1) *
		sizeof(int));
	basis->row_work = (double *)lp_simplex_malloc(
		(size_t)standard->rows * sizeof(double));
	basis->right_work = (double *)lp_simplex_malloc(
		(size_t)standard->rows * sizeof(double));
	backend =
		(struct pan_SparseBackend *)lp_simplex_malloc(sizeof(*backend));
	if (basis->position == NULL || basis->row_work == NULL ||
		basis->right_work == NULL || backend == NULL) {
		lp_simplex_free(backend);
		simplex_pan_basis_destroy(basis);
		return lp_simplex_EXIT_FAILURE;
	}
	lp_simplex_memset(backend, 0, sizeof(*backend));
	if (!cholmod_start(&backend->common)) {
		lp_simplex_free(backend);
		simplex_pan_basis_destroy(basis);
		return lp_simplex_EXIT_FAILURE;
	}
	backend->common.supernodal = CHOLMOD_AUTO;
	backend->common.quick_return_if_not_posdef = 1;
	basis->sparse_backend = backend;
	if (pan_dynamic_sparse_qr_initialize(&backend->factor) ==
			lp_simplex_EXIT_FAILURE ||
		pan_multifrontal_qr_initialize(&backend->orthogonal,
			standard->rows,
			&backend->common) == lp_simplex_EXIT_FAILURE) {
		simplex_pan_basis_destroy(basis);
		return lp_simplex_EXIT_FAILURE;
	}
	backend->product_correction = (double *)lp_simplex_malloc(
		(size_t)standard->rows * sizeof(double));
	if (!backend->product_correction) {
		simplex_pan_basis_destroy(basis);
		return lp_simplex_EXIT_FAILURE;
	}
	if (standard->columns > 0)
		lp_simplex_memset(basis->position, 0xff,
			(size_t)standard->columns * sizeof(int));
	if (pan_basis_reserve(basis, standard->rows > 8 ? 8 : standard->rows) ==
		lp_simplex_EXIT_FAILURE) {
		simplex_pan_basis_destroy(basis);
		return lp_simplex_EXIT_FAILURE;
	}
	return lp_simplex_EXIT_SUCCESS;
}


void simplex_pan_basis_destroy(struct simplex_PanBasis *basis)
{
	struct pan_SparseBackend *backend;

	if (basis == NULL)
		return;
	backend = (struct pan_SparseBackend *)basis->sparse_backend;
	if (backend != NULL) {
		pan_multifrontal_qr_destroy(&backend->orthogonal);
		pan_dynamic_sparse_qr_destroy(&backend->factor);
		cholmod_free_sparse(&backend->matrix, &backend->common);
		cholmod_finish(&backend->common);

		lp_simplex_free(backend->condition_work);
		lp_simplex_free(backend->stability_mask);
		lp_simplex_free(backend->pricing_work);
		lp_simplex_free(backend->product_correction);
		lp_simplex_free(backend);
	}
	lp_simplex_free(basis->column);
	lp_simplex_free(basis->position);
	lp_simplex_free(basis->column_norm);
	lp_simplex_free(basis->basis_work);
	lp_simplex_free(basis->row_work);
	lp_simplex_free(basis->right_work);
	lp_simplex_memset(basis, 0, sizeof(*basis));
}


int simplex_pan_basis_add(struct simplex_PanBasis *basis, const int column)
{
	int capacity;

	if (column < 0 || column >= basis->standard->columns ||
		basis->position[column] >= 0 ||
		basis->count >= basis->standard->rows)
		return lp_simplex_EXIT_FAILURE;
	if (basis->count == basis->capacity) {
		capacity = basis->capacity > 0 ? 2 * basis->capacity : 8;
		if (capacity > basis->standard->rows)
			capacity = basis->standard->rows;
		if (pan_basis_reserve(basis, capacity) ==
			lp_simplex_EXIT_FAILURE)
			return lp_simplex_EXIT_FAILURE;
	}
	basis->position[column] = basis->count;
	basis->column[basis->count++] = column;
	return lp_simplex_EXIT_SUCCESS;
}


void simplex_pan_basis_remove_marked(
	struct simplex_PanBasis *basis, const unsigned char *remove)
{
	struct pan_SparseBackend *backend =
		(struct pan_SparseBackend *)basis->sparse_backend;
	int synchronized = pan_dynamic_factor_matches(basis);
	int destination = 0;
	int i;

	if (synchronized)
		for (i = basis->count - 1; i >= 0; i--)
			if (remove[i]) {
				pan_multifrontal_qr_mark_path(
					&backend->orthogonal, i);
				if (pan_dynamic_sparse_qr_delete(
					    &backend->factor, i,
					    basis->rank_tolerance) ==
						lp_simplex_EXIT_FAILURE ||
					pan_multifrontal_qr_record_givens(
						&backend->orthogonal,
						&backend->factor) ==
						lp_simplex_EXIT_FAILURE) {
					pan_dynamic_sparse_qr_clear(
						&backend->factor);
					pan_multifrontal_qr_reset(
						&backend->orthogonal);
					synchronized = 0;
					break;
				}
			}
	for (i = 0; i < basis->count; i++) {
		int column = basis->column[i];

		if (remove[i]) {
			basis->position[column] = -1;
			continue;
		}
		basis->column[destination] = column;
		basis->column_norm[destination] = basis->column_norm[i];
		basis->position[column] = destination++;
	}
	basis->count = destination;
	backend->matrix_dirty = 1;
	pan_sync_dynamic_statistics(basis);
}


static cholmod_sparse *pan_sparse_basis_matrix(struct simplex_PanBasis *basis)
{
	struct pan_SparseBackend *backend =
		(struct pan_SparseBackend *)basis->sparse_backend;
	cholmod_sparse *matrix;
	double *values;
	int *column_start;
	int *row_index;
	size_t nonzeros = 0;
	int i, k, position = 0;

	for (i = 0; i < basis->count; i++)
		nonzeros +=
			(size_t)(basis->standard
					 ->column_start[basis->column[i] + 1] -
				 basis->standard
					 ->column_start[basis->column[i]]);
	matrix = cholmod_allocate_sparse((size_t)basis->standard->rows,
		(size_t)basis->count, nonzeros, 1, 1, 0, CHOLMOD_REAL,
		&backend->common);
	if (matrix == NULL)
		return NULL;
	column_start = (int *)matrix->p;
	row_index = (int *)matrix->i;
	values = (double *)matrix->x;
	for (i = 0; i < basis->count; i++) {
		int column = basis->column[i];
		double norm = simplex_pan_column_norm(basis->standard, column);

		column_start[i] = position;
		basis->column_norm[i] = norm;
		if (norm == 0.) {
			cholmod_free_sparse(&matrix, &backend->common);

			return NULL;
		}
		for (k = basis->standard->column_start[column];
			k < basis->standard->column_start[column + 1]; k++) {
			row_index[position] = basis->standard->row_index[k];
			values[position++] =
				basis->standard->normalized_value[k];
		}
	}
	column_start[basis->count] = position;
	return matrix;
}


static int pan_maybe_rebuild_orthogonal_anchor(struct simplex_PanBasis *basis)
{
	struct pan_SparseBackend *backend =
		(struct pan_SparseBackend *)basis->sparse_backend;
	cholmod_sparse *matrix;

	if (!pan_multifrontal_qr_should_rebuild(&backend->orthogonal))
		return lp_simplex_EXIT_SUCCESS;
	if (!pan_dynamic_factor_matches(basis))
		return lp_simplex_EXIT_FAILURE;
	matrix = pan_sparse_basis_matrix(basis);
	if (matrix == NULL)
		return lp_simplex_EXIT_FAILURE;
	if (pan_multifrontal_qr_rebuild(&backend->orthogonal, matrix,
		    &backend->factor, basis->column,
		    basis->rank_tolerance) == lp_simplex_EXIT_FAILURE) {
		cholmod_free_sparse(&matrix, &backend->common);

		return lp_simplex_EXIT_FAILURE;
	}
	cholmod_free_sparse(&backend->matrix, &backend->common);

	backend->matrix = matrix;
	backend->matrix_dirty = 0;
	pan_accept_factor_order(basis);
	/* matrix has the old column order; rebuild its cached view on demand. */
	backend->matrix_dirty = 1;
	basis->factorizations++;
	basis->factor_flops += backend->orthogonal.anchor_factor_work;
	pan_sync_dynamic_statistics(basis);
	return lp_simplex_EXIT_SUCCESS;
}


int simplex_pan_basis_rebuild(struct simplex_PanBasis *basis)
{
	struct pan_SparseBackend *backend =
		(struct pan_SparseBackend *)basis->sparse_backend;
	cholmod_sparse *matrix;

	if (basis->count == 0) {
		pan_dynamic_sparse_qr_clear(&backend->factor);
		pan_multifrontal_qr_reset(&backend->orthogonal);
		backend->matrix_dirty = 1;
		return lp_simplex_EXIT_SUCCESS;
	}
	matrix = pan_sparse_basis_matrix(basis);
	if (matrix == NULL)
		return lp_simplex_EXIT_FAILURE;
	if (pan_multifrontal_qr_rebuild(&backend->orthogonal, matrix,
		    &backend->factor, basis->column,
		    basis->rank_tolerance) == lp_simplex_EXIT_FAILURE) {
		cholmod_free_sparse(&matrix, &backend->common);

		return lp_simplex_EXIT_FAILURE;
	}
	cholmod_free_sparse(&matrix, &backend->common);

	backend->matrix_dirty = 1;
	basis->factorizations++;
	basis->factor_flops += backend->orthogonal.anchor_factor_work;
	pan_accept_factor_order(basis);
	pan_sync_dynamic_statistics(basis);
	return lp_simplex_EXIT_SUCCESS;
}


int simplex_pan_basis_repair_rank(
	struct simplex_PanBasis *basis, const double *value)
{
	struct pan_SparseBackend *backend =
		(struct pan_SparseBackend *)basis->sparse_backend;
	cholmod_sparse *matrix;
	int *order;
	double *coefficient;
	int i, positive = 0, count = basis->count;
	int result = lp_simplex_EXIT_FAILURE;

	if (count == 0)
		return simplex_pan_basis_rebuild(basis);
	order = (int *)lp_simplex_malloc((size_t)count * sizeof(int));
	coefficient =
		(double *)lp_simplex_malloc((size_t)count * sizeof(double));
	if (order == NULL || coefficient == NULL) {
		lp_simplex_free(order);
		lp_simplex_free(coefficient);
		return lp_simplex_EXIT_FAILURE;
	}
	for (i = 0; i < count; i++)
		if (value[basis->column[i]] != 0.)
			order[positive++] = basis->column[i];
	count = positive;
	for (i = 0; i < basis->count; i++)
		if (value[basis->column[i]] == 0.)
			order[count++] = basis->column[i];
	lp_simplex_memcpy(basis->column, order, (size_t)count * sizeof(int));
	for (i = 0; i < count; i++)
		basis->position[order[i]] = i;
	matrix = pan_sparse_basis_matrix(basis);
	if (matrix == NULL)
		goto cleanup;
	/* First obtain a well-conditioned independent subset. Its numerical
	 * ordering may omit required columns; the exchange loop restores them. */
	if (pan_multifrontal_qr_repair_rank(&backend->orthogonal, matrix,
		    &backend->factor, basis->column,
		    __lp_simplex_MAX__(
			    basis->rank_tolerance, sqrt(DBL_EPSILON)),
		    0) == lp_simplex_EXIT_FAILURE) {
		cholmod_free_sparse(&matrix, &backend->common);

		goto cleanup;
	}
	cholmod_free_sparse(&matrix, &backend->common);

	for (i = 0; i < count; i++)
		basis->position[order[i]] = -1;
	basis->count = backend->factor.count;
	pan_accept_factor_order(basis);
	basis->factorizations++;
	basis->factor_flops += backend->orthogonal.anchor_factor_work;
	backend->matrix_dirty = 1;
	if (basis->count < positive)
		goto cleanup;
	for (i = 0; i < positive; i++) {
		int j, entering = order[i], leaving = -1;
		double residual, best = 0.;

		if (basis->position[entering] >= 0)
			continue;
		if (simplex_pan_basis_column_projection(basis, entering,
			    coefficient,
			    &residual) == lp_simplex_EXIT_FAILURE ||
			residual / simplex_pan_column_norm(
					   basis->standard, entering) >
				sqrt(DBL_EPSILON))
			goto cleanup;
		for (j = 0; j < basis->count; j++) {
			double strength =
				fabs(coefficient[j]) * basis->column_norm[j];

			if (value[basis->column[j]] == 0. && strength > best) {
				best = strength;
				leaving = j;
			}
		}
		if (leaving < 0 || !isfinite(best) ||
			best / simplex_pan_column_norm(
				       basis->standard, entering) <=
				basis->rank_tolerance ||
			simplex_pan_basis_replace_factorized(basis, leaving,
				entering) == lp_simplex_EXIT_FAILURE)
			goto cleanup;
	}
	basis->rank_repairs++;
	pan_sync_dynamic_statistics(basis);
	result = lp_simplex_EXIT_SUCCESS;
cleanup:
	lp_simplex_free(order);
	lp_simplex_free(coefficient);
	return result;
}


/* Hager-style sparse estimate of ||R^-1||_1. Its large response vector also
 * identifies optional columns involved in a weak direction. */
static int pan_estimate_condition(struct simplex_PanBasis *basis)
{
	struct pan_SparseBackend *backend =
		(struct pan_SparseBackend *)basis->sparse_backend;
	struct pan_DynamicSparseQR *r = &backend->factor;
	int n = basis->count, i, j, pass, previous = -1;
	double norm = 0., inverse = 0.;
	double *x = backend->condition_work;
	double *y = x + basis->standard->rows;
	double *z = y + basis->standard->rows;
	double *weak = z + basis->standard->rows;

	if (n == 0) {
		basis->reciprocal_condition = 1.;
		return lp_simplex_EXIT_SUCCESS;
	}
	for (j = 0; j < n; j++) {
		int node = r->column_head[j];
		double sum = 0.;

		while (node >= 0) {
			sum += fabs(r->node[node].value);
			node = r->node[node].column_next;
		}
		norm = fmax(norm, sum);
		x[j] = 1. / n;
		weak[j] = 0.;
	}
	for (pass = 0; pass < 5; pass++) {
		double estimate = 0., maximum = 0., pairing = 0.;
		int index = -1;

		lp_simplex_memcpy(y, x, (size_t)n * sizeof(double));
		if (pan_dynamic_sparse_qr_solve_upper(r, y))
			return lp_simplex_EXIT_FAILURE;
		for (i = 0; i < n; i++) {
			estimate += fabs(y[i]);
			z[i] = y[i] >= 0. ? 1. : -1.;
		}
		if (estimate > inverse) {
			inverse = estimate;
			lp_simplex_memcpy(weak, y, (size_t)n * sizeof(double));
		}
		if (pan_dynamic_sparse_qr_solve_upper_transpose(r, z))
			return lp_simplex_EXIT_FAILURE;
		for (i = 0; i < n; i++) {
			pairing += z[i] * x[i];
			if (fabs(z[i]) > maximum) {
				maximum = fabs(z[i]);
				index = i;
			}
		}
		if (index < 0 || index == previous || maximum <= pairing)
			break;
		lp_simplex_memset(x, 0, (size_t)n * sizeof(double));
		x[index] = 1.;
		previous = index;
	}
	basis->reciprocal_condition =
		norm > 0. && inverse > 0. ? 1. / (norm * inverse) : 0.;
	return isfinite(basis->reciprocal_condition) ? lp_simplex_EXIT_SUCCESS
						     : lp_simplex_EXIT_FAILURE;
}


int simplex_pan_basis_stabilize(struct simplex_PanBasis *basis,
	const double *value, const int protected_column)
{
	struct pan_SparseBackend *backend =
		(struct pan_SparseBackend *)basis->sparse_backend;
	int attempt;

	if (backend->condition_work == NULL) {
		backend->condition_work = (double *)lp_simplex_malloc(
			(size_t)basis->standard->rows * 4 * sizeof(double));
		backend->stability_mask = (unsigned char *)lp_simplex_malloc(
			(size_t)basis->standard->rows);
		if (backend->condition_work == NULL ||
			backend->stability_mask == NULL) {
			lp_simplex_free(backend->condition_work);
			lp_simplex_free(backend->stability_mask);
			backend->condition_work = NULL;
			backend->stability_mask = NULL;
			return lp_simplex_EXIT_FAILURE;
		}
	}
	for (attempt = 0; attempt <= basis->standard->rows; attempt++) {
		int i, leaving = -1;
		double best = 0.;
		double *weak = backend->condition_work +
			       (size_t)basis->standard->rows * 3;

		if (pan_estimate_condition(basis))
			return lp_simplex_EXIT_FAILURE;
		if (basis->reciprocal_condition >= sqrt(DBL_EPSILON))
			return lp_simplex_EXIT_SUCCESS;
		for (i = 0; i < basis->count; i++) {
			int column = basis->column[i];

			if (value[column] == 0. && column != protected_column &&
				fabs(weak[i]) > best) {
				best = fabs(weak[i]);
				leaving = i;
			}
		}
		if (leaving < 0)
			return lp_simplex_EXIT_FAILURE;
		lp_simplex_memset(
			backend->stability_mask, 0, (size_t)basis->count);
		backend->stability_mask[leaving] = 1;
		simplex_pan_basis_remove_marked(basis, backend->stability_mask);
		basis->stability_drops++;
		if (simplex_pan_basis_factorize(basis))
			return lp_simplex_EXIT_FAILURE;
	}
	return lp_simplex_EXIT_FAILURE;
}


int simplex_pan_basis_pricing_probes(struct simplex_PanBasis *basis,
	const double *value, double *span_probe, double *edge_probe)
{
	struct pan_SparseBackend *backend =
		(struct pan_SparseBackend *)basis->sparse_backend;
	const int rows = basis->standard->rows;
	int i, zeros = 0;
	double norm, residual, saved_error = backend->orthogonal.backward_error;
	double scale = rows > basis->count
			       ? 1. / sqrt((double)(rows - basis->count))
			       : 0.;

	if (backend->pricing_work == NULL) {
		backend->pricing_work = (double *)lp_simplex_malloc(
			(size_t)rows * sizeof(double));
		if (backend->pricing_work == NULL)
			return lp_simplex_EXIT_FAILURE;
	}
	for (i = 0; i < rows; i++) {
		unsigned int hash = ((unsigned int)i + 1U) * 2654435761U;

		hash ^= hash >> 16;
		basis->row_work[i] =
			i < basis->count ? 0. : (hash & 1U ? scale : -scale);
	}
	if (pan_multifrontal_qr_apply(
		    &backend->orthogonal, basis->row_work, span_probe))
		return lp_simplex_EXIT_FAILURE;
	for (i = 0; i < basis->count; i++) {
		int column = basis->column[i];
		unsigned int hash = ((unsigned int)column + 1U) * 2654435761U;

		backend->pricing_work[i] =
			value[column] == 0.
				? .5 + (double)(hash & 1023U) / 1024.
				: 0.;
		zeros += value[column] == 0.;
	}
	if (zeros == 0) {
		lp_simplex_memset(edge_probe, 0, (size_t)rows * sizeof(double));
		return lp_simplex_EXIT_SUCCESS;
	}
	if (simplex_pan_basis_minimum_norm(
		    basis, backend->pricing_work, edge_probe, &residual))
		return lp_simplex_EXIT_FAILURE;
	backend->orthogonal.backward_error =
		fmax(saved_error, backend->orthogonal.backward_error);
	norm = pan_dense_norm(edge_probe, rows);
	if (!isfinite(norm) || norm == 0. ||
		residual >
			sqrt(DBL_EPSILON) * (1. + sqrt((double)basis->count)))
		return lp_simplex_EXIT_FAILURE;
	for (i = 0; i < rows; i++)
		edge_probe[i] /= norm;
	return lp_simplex_EXIT_SUCCESS;
}


int simplex_pan_basis_select_columns(struct simplex_PanBasis *basis,
	const int *column, int count, const double *required)
{
	int i;

	if (basis->count != 0 || count < 0 ||
		count > basis->standard->columns ||
		pan_basis_reserve(basis, count))
		return lp_simplex_EXIT_FAILURE;
	for (i = 0; i < count; i++) {
		int j = column[i];

		if (j < 0 || j >= basis->standard->columns ||
			basis->position[j] >= 0)
			return lp_simplex_EXIT_FAILURE;
		basis->column[basis->count] = j;
		basis->position[j] = basis->count++;
	}
	return simplex_pan_basis_repair_rank(basis, required);
}


int simplex_pan_basis_complete(
	struct simplex_PanBasis *basis, const double *value)
{
	int *saved;
	int count = basis->count, j, result;

	saved = (int *)lp_simplex_malloc(
		(size_t)(count > 0 ? count : 1) * sizeof(int));
	if (saved == NULL)
		return lp_simplex_EXIT_FAILURE;
	lp_simplex_memcpy(saved, basis->column, (size_t)count * sizeof(int));
	if (pan_basis_reserve(basis, basis->standard->columns) ==
		lp_simplex_EXIT_FAILURE) {
		lp_simplex_free(saved);
		return lp_simplex_EXIT_FAILURE;
	}
	/* Only the batch rank-revealing call sees the overcomplete set. No
	 * triangular solve or dynamic update is performed until it is reduced. */
	basis->count = 0;
	for (j = 0; j < basis->standard->columns; j++) {
		basis->position[j] = -1;
		if (simplex_pan_column_norm(basis->standard, j) == 0.)
			continue;
		basis->column[basis->count] = j;
		basis->position[j] = basis->count++;
	}
	result = simplex_pan_basis_repair_rank(basis, value);
	if (result == lp_simplex_EXIT_FAILURE) {
		for (j = 0; j < basis->standard->columns; j++)
			basis->position[j] = -1;
		basis->count = count;
		lp_simplex_memcpy(
			basis->column, saved, (size_t)count * sizeof(int));
		for (j = 0; j < count; j++)
			basis->position[saved[j]] = j;
		simplex_pan_basis_rebuild(basis);
	}
	lp_simplex_free(saved);
	return result;
}


static int pan_timed_factorize_impl(struct simplex_PanBasis *basis)
{
	struct pan_SparseBackend *backend =
		(struct pan_SparseBackend *)basis->sparse_backend;
	int append = backend->factor.count + 1 == basis->count;
	int i;
	int j;

	if (basis->count == 0) {
		pan_dynamic_sparse_qr_clear(&backend->factor);
		pan_multifrontal_qr_reset(&backend->orthogonal);
		backend->matrix_dirty = 1;
		return lp_simplex_EXIT_SUCCESS;
	}
	if (pan_dynamic_factor_matches(basis))
		return pan_maybe_rebuild_orthogonal_anchor(basis);
	if (append)
		for (i = 0; i < backend->factor.count; i++)
			if (backend->factor.tag[i] != basis->column[i]) {
				append = 0;
				break;
			}
	if (append) {
		j = basis->count - 1;
		basis->column_norm[j] = simplex_pan_column_norm(
			basis->standard, basis->column[j]);
		if (basis->column_norm[j] == 0.)
			return lp_simplex_EXIT_FAILURE;
		if (pan_dynamic_append_column(basis, basis->column[j]) ==
			lp_simplex_EXIT_FAILURE)
			return lp_simplex_EXIT_FAILURE;
		basis->updates++;
		backend->matrix_dirty = 1;
		pan_sync_dynamic_statistics(basis);
		return pan_maybe_rebuild_orthogonal_anchor(basis);
	}
	return simplex_pan_basis_rebuild(basis);
}


int simplex_pan_basis_replace_factorized(
	struct simplex_PanBasis *basis, const int position, const int column)
{
	struct pan_SparseBackend *backend =
		(struct pan_SparseBackend *)basis->sparse_backend;
	double entering_norm;
	int leaving;

	if (position < 0 || position >= basis->count || column < 0 ||
		column >= basis->standard->columns ||
		basis->position[column] >= 0)
		return lp_simplex_EXIT_FAILURE;
	entering_norm = simplex_pan_column_norm(basis->standard, column);
	if (entering_norm == 0.)
		return lp_simplex_EXIT_FAILURE;
	leaving = basis->column[position];
	if (pan_dynamic_factor_matches(basis) &&
		pan_dynamic_sparse_qr_delete(&backend->factor, position,
			basis->rank_tolerance) == lp_simplex_EXIT_SUCCESS) {
		pan_multifrontal_qr_mark_path(&backend->orthogonal, position);
		if (pan_multifrontal_qr_record_givens(&backend->orthogonal,
			    &backend->factor) == lp_simplex_EXIT_SUCCESS &&
			pan_dynamic_append_column(basis, column) ==
				lp_simplex_EXIT_SUCCESS &&
			pan_dynamic_sparse_qr_move_last(&backend->factor,
				position, basis->rank_tolerance) ==
				lp_simplex_EXIT_SUCCESS &&
			pan_multifrontal_qr_record_givens(&backend->orthogonal,
				&backend->factor) == lp_simplex_EXIT_SUCCESS) {
			basis->position[leaving] = -1;
			basis->column[position] = column;
			basis->position[column] = position;
			basis->column_norm[position] = entering_norm;
			basis->updates++;
			backend->matrix_dirty = 1;
			pan_sync_dynamic_statistics(basis);
			return pan_maybe_rebuild_orthogonal_anchor(basis);
		}
	}
	pan_dynamic_sparse_qr_clear(&backend->factor);
	pan_multifrontal_qr_reset(&backend->orthogonal);
	basis->position[leaving] = -1;
	basis->column[position] = column;
	basis->position[column] = position;
	basis->column_norm[position] = entering_norm;
	backend->matrix_dirty = 1;
	return simplex_pan_basis_factorize(basis);
}


static cholmod_dense *pan_dense_create(
	struct pan_SparseBackend *backend, const double *value, const int rows)
{
	cholmod_dense *dense = cholmod_allocate_dense(
		(size_t)rows, 1, (size_t)rows, CHOLMOD_REAL, &backend->common);

	if (dense != NULL && rows > 0)
		lp_simplex_memcpy(
			dense->x, value, (size_t)rows * sizeof(double));
	return dense;
}


static int pan_timed_least_squares_impl(struct simplex_PanBasis *basis,
	const double *right, double *solution, double *residual_norm)
{
	struct pan_SparseBackend *backend =
		(struct pan_SparseBackend *)basis->sparse_backend;
	int refinement, i;
	double right_norm = pan_dense_norm(right, basis->standard->rows);

	if (basis->count == 0) {
		lp_simplex_memcpy(basis->row_work, right,
			(size_t)basis->standard->rows * sizeof(double));
		*residual_norm =
			pan_dense_norm(basis->row_work, basis->standard->rows);
		return lp_simplex_EXIT_SUCCESS;
	}
	if (!pan_dynamic_factor_matches(basis))
		return lp_simplex_EXIT_FAILURE;
	if (pan_implicit_q_transpose_multiply(basis, right, solution) ==
			lp_simplex_EXIT_FAILURE ||
		pan_dynamic_sparse_qr_solve_upper(&backend->factor, solution) ==
			lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	for (i = 0; i < basis->count; i++)
		solution[i] /= basis->column_norm[i];
	for (refinement = 0; refinement < 2; refinement++) {
		pan_basis_multiply(basis, solution, basis->row_work);
		for (i = 0; i < basis->standard->rows; i++)
			basis->row_work[i] = right[i] - basis->row_work[i];
		if (pan_ls_stationary(basis, basis->row_work, right_norm))
			break;
		if (pan_implicit_q_transpose_multiply(basis, basis->row_work,
			    basis->basis_work) == lp_simplex_EXIT_FAILURE ||
			pan_dynamic_sparse_qr_solve_upper(&backend->factor,
				basis->basis_work) == lp_simplex_EXIT_FAILURE)
			return lp_simplex_EXIT_FAILURE;
		for (i = 0; i < basis->count; i++)
			solution[i] +=
				basis->basis_work[i] / basis->column_norm[i];
		basis->refinements++;
	}
	pan_basis_multiply(basis, solution, basis->row_work);
	for (i = 0; i < basis->standard->rows; i++)
		basis->row_work[i] = right[i] - basis->row_work[i];
	*residual_norm = pan_dense_norm(basis->row_work, basis->standard->rows);
	return lp_simplex_EXIT_SUCCESS;
}


int simplex_pan_basis_least_squares_qr(struct simplex_PanBasis *basis,
	const double *right, double *solution, double *residual_norm)
{
	struct pan_SparseBackend *backend =
		(struct pan_SparseBackend *)basis->sparse_backend;
	cholmod_dense *right_dense = NULL;
	cholmod_dense *scaled_solution = NULL;
	double *values;
	int i;

	if (basis->count == 0)
		return simplex_pan_basis_least_squares(
			basis, right, solution, residual_norm);
	if (backend->matrix_dirty) {
		cholmod_sparse *matrix = pan_sparse_basis_matrix(basis);

		if (matrix == NULL)
			goto failure;
		cholmod_free_sparse(&backend->matrix, &backend->common);

		backend->matrix = matrix;
		backend->matrix_dirty = 0;
	}
	right_dense = pan_dense_create(backend, right, basis->standard->rows);
	if (right_dense == NULL)
		goto failure;
	{
		int32_t rank = SuiteSparseQR_i_C(SPQR_ORDERING_COLAMD,
			basis->rank_tolerance, (int32_t)basis->count, 2,
			backend->matrix, NULL, right_dense, NULL,
			&scaled_solution, NULL, NULL, NULL, NULL, NULL,
			&backend->common);
		basis->factorizations++;
		basis->factor_flops += backend->common.SPQR_flopcount;
		if (rank < 0)
			goto failure;
	}
	if (scaled_solution == NULL ||
		scaled_solution->nrow < (size_t)basis->count)
		goto failure;
	values = (double *)scaled_solution->x;
	for (i = 0; i < basis->count; i++)
		solution[i] = values[i] / basis->column_norm[i];
	pan_basis_multiply(basis, solution, basis->row_work);
	for (i = 0; i < basis->standard->rows; i++)
		basis->row_work[i] = right[i] - basis->row_work[i];
	*residual_norm = pan_dense_norm(basis->row_work, basis->standard->rows);
	cholmod_free_dense(&scaled_solution, &backend->common);
	cholmod_free_dense(&right_dense, &backend->common);

	return lp_simplex_EXIT_SUCCESS;
failure:
	cholmod_free_dense(&scaled_solution, &backend->common);
	cholmod_free_dense(&right_dense, &backend->common);

	return lp_simplex_EXIT_FAILURE;
}


static int pan_timed_minimum_norm_impl(struct simplex_PanBasis *basis,
	const double *right, double *solution, double *residual_norm)
{
	struct pan_SparseBackend *backend =
		(struct pan_SparseBackend *)basis->sparse_backend;
	int refinement, i;

	if (basis->count == 0) {
		lp_simplex_memset(solution, 0,
			(size_t)basis->standard->rows * sizeof(double));
		*residual_norm = 0.;
		return lp_simplex_EXIT_SUCCESS;
	}
	for (i = 0; i < basis->count; i++)
		basis->basis_work[i] = right[i] / basis->column_norm[i];
	if (!pan_dynamic_factor_matches(basis) ||
		pan_dynamic_sparse_qr_solve_upper_transpose(&backend->factor,
			basis->basis_work) == lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	if (pan_implicit_q_multiply(basis, basis->basis_work, basis->right_work,
		    solution) == lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	for (refinement = 0; refinement < 2; refinement++) {
		double solution_norm =
			pan_dense_norm(solution, basis->standard->rows);
		int accurate = !basis->full_refinement;

		for (i = 0; i < basis->count; i++)
			basis->basis_work[i] =
				right[i] - pan_column_dense_dot(basis->standard,
						   basis->column[i], solution);
		for (i = 0; i < basis->count; i++) {
			double normalized =
				basis->basis_work[i] / basis->column_norm[i];

			if (!isfinite(normalized) ||
				fabs(normalized) >
					4. * DBL_EPSILON *
						(fabs(right[i] /
							 basis->column_norm
								 [i]) +
							solution_norm))
				accurate = 0;
		}
		if (accurate)
			break;
		for (i = 0; i < basis->count; i++)
			basis->basis_work[i] /= basis->column_norm[i];
		if (pan_dynamic_sparse_qr_solve_upper_transpose(
			    &backend->factor, basis->basis_work) ==
			lp_simplex_EXIT_FAILURE)
			return lp_simplex_EXIT_FAILURE;
		if (pan_implicit_q_multiply(basis, basis->basis_work,
			    basis->row_work,
			    basis->right_work) == lp_simplex_EXIT_FAILURE)
			return lp_simplex_EXIT_FAILURE;
		for (i = 0; i < basis->standard->rows; i++)
			solution[i] += basis->right_work[i];
		basis->refinements++;
	}
	for (i = 0; i < basis->count; i++)
		basis->basis_work[i] =
			right[i] - pan_column_dense_dot(basis->standard,
					   basis->column[i], solution);
	*residual_norm = pan_dense_norm(basis->basis_work, basis->count);
	pan_record_backward_error(
		basis, *residual_norm, pan_dense_norm(right, basis->count));
	return lp_simplex_EXIT_SUCCESS;
}


static int pan_timed_column_projection_impl(struct simplex_PanBasis *basis,
	const int column, double *coefficient, double *residual_norm)
{
	struct pan_SparseBackend *backend =
		(struct pan_SparseBackend *)basis->sparse_backend;
	int refinement, i, k;
	double orthogonal_residual;

	if (basis->count == 0) {
		*residual_norm =
			simplex_pan_column_norm(basis->standard, column);
		basis->projection_solve_error = *residual_norm;
		return lp_simplex_EXIT_SUCCESS;
	}
	if (!pan_dynamic_factor_matches(basis))
		return lp_simplex_EXIT_FAILURE;
	lp_simplex_memset(basis->row_work, 0,
		(size_t)basis->standard->rows * sizeof(double));
	pan_column_axpy(basis->standard, column, 1., basis->row_work);
	if (pan_implicit_q_transpose_multiply(basis, basis->row_work,
		    coefficient) == lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	/* The orthogonal tail measures rank without amplifying the triangular
	 * condition number. ||a-B*delta|| also includes solve error and must not
	 * be mistaken for a new independent direction. */
	orthogonal_residual = pan_dense_norm(basis->right_work + basis->count,
		basis->standard->rows - basis->count);
	if (pan_dynamic_sparse_qr_solve_upper(&backend->factor, coefficient) ==
		lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	for (i = 0; i < basis->count; i++)
		coefficient[i] /= basis->column_norm[i];
	for (refinement = 0; refinement < 2; refinement++) {
		pan_basis_multiply(basis, coefficient, basis->row_work);
		for (i = 0; i < basis->standard->rows; i++)
			basis->row_work[i] = -basis->row_work[i];
		for (k = basis->standard->column_start[column];
			k < basis->standard->column_start[column + 1]; k++)
			basis->row_work[basis->standard->row_index[k]] +=
				basis->standard->value[k];
		if (pan_ls_stationary(basis, basis->row_work,
			    simplex_pan_column_norm(basis->standard, column)))
			break;
		if (pan_implicit_q_transpose_multiply(basis, basis->row_work,
			    basis->basis_work) == lp_simplex_EXIT_FAILURE ||
			pan_dynamic_sparse_qr_solve_upper(&backend->factor,
				basis->basis_work) == lp_simplex_EXIT_FAILURE)
			return lp_simplex_EXIT_FAILURE;
		for (i = 0; i < basis->count; i++)
			coefficient[i] +=
				basis->basis_work[i] / basis->column_norm[i];
		basis->refinements++;
	}
	pan_basis_multiply(basis, coefficient, basis->row_work);
	for (i = 0; i < basis->standard->rows; i++)
		basis->row_work[i] = -basis->row_work[i];
	for (k = basis->standard->column_start[column];
		k < basis->standard->column_start[column + 1]; k++)
		basis->row_work[basis->standard->row_index[k]] +=
			basis->standard->value[k];
	basis->projection_solve_error =
		pan_dense_norm(basis->row_work, basis->standard->rows);
	*residual_norm = orthogonal_residual;
	if (orthogonal_residual / simplex_pan_column_norm(
					  basis->standard, column) <=
			basis->rank_tolerance &&
		basis->projection_solve_error /
				simplex_pan_column_norm(
					basis->standard, column) >
			basis->rank_tolerance)
		basis->rank_disagreements++;
	return lp_simplex_EXIT_SUCCESS;
}


int simplex_pan_basis_factorize(struct simplex_PanBasis *basis)
{
	int profile = getenv("LP_SIMPLEX_PROFILE") != NULL;
	clock_t start = profile ? clock() : 0;
	int result = pan_timed_factorize_impl(basis);

	if (profile)
		basis->factor_seconds +=
			(double)(clock() - start) / CLOCKS_PER_SEC;
	return result;
}


int simplex_pan_basis_least_squares(struct simplex_PanBasis *basis,
	const double *right, double *solution, double *residual_norm)
{
	int profile = getenv("LP_SIMPLEX_PROFILE") != NULL;
	clock_t start = profile ? clock() : 0;
	int result = pan_timed_least_squares_impl(
		basis, right, solution, residual_norm);

	pan_sync_pack_statistics(basis);
	if (profile)
		basis->least_squares_seconds +=
			(double)(clock() - start) / CLOCKS_PER_SEC;
	return result;
}


int simplex_pan_basis_minimum_norm(struct simplex_PanBasis *basis,
	const double *right, double *solution, double *residual_norm)
{
	int profile = getenv("LP_SIMPLEX_PROFILE") != NULL;
	clock_t start = profile ? clock() : 0;
	int result = pan_timed_minimum_norm_impl(
		basis, right, solution, residual_norm);

	pan_sync_pack_statistics(basis);
	if (profile)
		basis->minimum_norm_seconds +=
			(double)(clock() - start) / CLOCKS_PER_SEC;
	return result;
}


int simplex_pan_basis_column_projection(struct simplex_PanBasis *basis,
	int column, double *coefficient, double *residual_norm)
{
	int profile = getenv("LP_SIMPLEX_PROFILE") != NULL;
	clock_t start = profile ? clock() : 0;
	int result = pan_timed_column_projection_impl(
		basis, column, coefficient, residual_norm);

	pan_sync_pack_statistics(basis);
	if (profile)
		basis->projection_seconds +=
			(double)(clock() - start) / CLOCKS_PER_SEC;
	return result;
}
#endif
