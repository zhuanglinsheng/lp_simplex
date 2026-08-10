/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 *
 * Sparse linear-algebra backend for Pan's dynamically deficient basis.
 * The hot path maintains a sparse dynamic R factor and applies Q implicitly
 * through the normalized active columns.  SuiteSparseQR is reserved for
 * rank-revealing recovery/certification.
 */
#include "simplex_pan_basis.h"

#ifdef LP_SIMPLEX_HAVE_SPQR

#include "pan_dynamic_sparse_qr.h"
#include "pan_multifrontal_qr.h"
#include "utils.h"

#include <lp_simplex/status.h>

#include <SuiteSparseQR_C.h>

#include <math.h>


struct pan_SparseBackend {
	cholmod_common common;
	cholmod_sparse *matrix;
	struct pan_DynamicSparseQR factor;
	struct pan_MultifrontalQR orthogonal;
	int matrix_dirty;
};


static cholmod_sparse *pan_sparse_basis_matrix(
		struct simplex_PanBasis *basis);


static void pan_sync_dynamic_statistics(struct simplex_PanBasis *basis)
{
	struct pan_SparseBackend *backend =
		(struct pan_SparseBackend *)basis->sparse_backend;
	pan_multifrontal_qr_refresh_structure(&backend->orthogonal,
		&backend->factor);
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
}


static void pan_record_backward_error(
		struct simplex_PanBasis *basis, const double residual,
		const double right_norm)
{
	struct pan_SparseBackend *backend =
		(struct pan_SparseBackend *)basis->sparse_backend;
	double error = residual / (1. + right_norm);
	backend->orthogonal.backward_error = error;
	pan_sync_dynamic_statistics(basis);
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


static void pan_column_axpy(
		const struct simplex_PanStandard *standard,
		const int column, const double scale, double *dense)
{
	int k;
	for (k = standard->column_start[column];
	     k < standard->column_start[column + 1]; k++)
		dense[standard->row_index[k]] += scale * standard->value[k];
}


static double pan_column_dense_dot(
		const struct simplex_PanStandard *standard,
		const int column, const double *dense)
{
	double dot = 0.;
	int k;
	for (k = standard->column_start[column];
	     k < standard->column_start[column + 1]; k++)
		dot += standard->value[k] * dense[standard->row_index[k]];
	return dot;
}


static void pan_normalized_column_axpy(
		const struct simplex_PanStandard *standard,
		const int column, const double scale, double *dense)
{
	int k;
	for (k = standard->column_start[column];
	     k < standard->column_start[column + 1]; k++)
		dense[standard->row_index[k]] +=
			scale * standard->normalized_value[k];
}


static void pan_basis_multiply(
		const struct simplex_PanBasis *basis,
		const double *coefficient, double *result)
{
	int i;
	lp_simplex_memset(result, 0,
		(size_t)basis->standard->rows * sizeof(double));
	for (i = 0; i < basis->count; i++)
		pan_column_axpy(basis->standard, basis->column[i],
			coefficient[i], result);
}


static int pan_implicit_q_transpose_multiply(
		const struct simplex_PanBasis *basis,
		const double *vector, double *result)
{
	struct pan_SparseBackend *backend =
		(struct pan_SparseBackend *)basis->sparse_backend;
	int i;
	if (pan_multifrontal_qr_apply_transpose(&backend->orthogonal,
		vector, basis->right_work) == lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	for (i = 0; i < backend->factor.count; i++)
		result[i] = basis->right_work[i];
	return lp_simplex_EXIT_SUCCESS;
}


static int pan_implicit_q_multiply(
		const struct simplex_PanBasis *basis,
		const double *coefficient, double *solve_work, double *result)
{
	struct pan_SparseBackend *backend =
		(struct pan_SparseBackend *)basis->sparse_backend;
	lp_simplex_memset(solve_work, 0,
		(size_t)basis->standard->rows * sizeof(double));
	if (backend->factor.count > 0)
		lp_simplex_memcpy(solve_work, coefficient,
			(size_t)backend->factor.count * sizeof(double));
	return pan_multifrontal_qr_apply(&backend->orthogonal,
		solve_work, result);
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
	lp_simplex_memset(basis->row_work, 0,
		(size_t)standard->rows * sizeof(double));
	pan_normalized_column_axpy(standard, column, 1., basis->row_work);
	if (pan_multifrontal_qr_append(&backend->orthogonal,
		basis->row_work, count, basis->basis_work, &diagonal) ==
	    lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	if (pan_dynamic_sparse_qr_append(&backend->factor, column,
		basis->basis_work, diagonal, basis->rank_tolerance) ==
	    lp_simplex_EXIT_FAILURE) {
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


static int pan_basis_reserve(
		struct simplex_PanBasis *basis, const int capacity)
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


int simplex_pan_basis_create(
		struct simplex_PanBasis *basis,
		const struct simplex_PanStandard *standard,
		const double rank_tolerance)
{
	struct pan_SparseBackend *backend;
	lp_simplex_memset(basis, 0, sizeof(*basis));
	basis->standard = standard;
	basis->rank_tolerance = rank_tolerance;
	basis->position = (int *)lp_simplex_malloc(
		(size_t)(standard->columns > 0 ? standard->columns : 1) *
		sizeof(int));
	basis->row_work = (double *)lp_simplex_malloc(
		(size_t)standard->rows * sizeof(double));
	basis->right_work = (double *)lp_simplex_malloc(
		(size_t)standard->rows * sizeof(double));
	backend = (struct pan_SparseBackend *)lp_simplex_malloc(sizeof(*backend));
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
		standard->rows, &backend->common) ==
	    lp_simplex_EXIT_FAILURE) {
		simplex_pan_basis_destroy(basis);
		return lp_simplex_EXIT_FAILURE;
	}
	if (standard->columns > 0)
		lp_simplex_memset(basis->position, 0xff,
			(size_t)standard->columns * sizeof(int));
	if (pan_basis_reserve(basis,
		standard->rows > 8 ? 8 : standard->rows) ==
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
	    basis->position[column] >= 0 || basis->count >= basis->standard->rows)
		return lp_simplex_EXIT_FAILURE;
	if (basis->count == basis->capacity) {
		capacity = basis->capacity > 0 ? 2 * basis->capacity : 8;
		if (capacity > basis->standard->rows)
			capacity = basis->standard->rows;
		if (pan_basis_reserve(basis, capacity) == lp_simplex_EXIT_FAILURE)
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
				pan_multifrontal_qr_mark_path(&backend->orthogonal, i);
				if (pan_dynamic_sparse_qr_delete(&backend->factor, i,
					basis->rank_tolerance) == lp_simplex_EXIT_FAILURE ||
				    pan_multifrontal_qr_record_givens(
					&backend->orthogonal, &backend->factor) ==
				    lp_simplex_EXIT_FAILURE) {
					pan_dynamic_sparse_qr_clear(&backend->factor);
					pan_multifrontal_qr_reset(&backend->orthogonal);
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


static cholmod_sparse *pan_sparse_basis_matrix(
		struct simplex_PanBasis *basis)
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
		nonzeros += (size_t)(basis->standard->column_start[
			basis->column[i] + 1] - basis->standard->column_start[
			basis->column[i]]);
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
			values[position++] = basis->standard->normalized_value[k];
		}
	}
	column_start[basis->count] = position;
	return matrix;
}


static int pan_maybe_rebuild_orthogonal_anchor(
		struct simplex_PanBasis *basis)
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
		&backend->factor, basis->column, basis->rank_tolerance) ==
	    lp_simplex_EXIT_FAILURE) {
		cholmod_free_sparse(&matrix, &backend->common);
		return lp_simplex_EXIT_FAILURE;
	}
	cholmod_free_sparse(&backend->matrix, &backend->common);
	backend->matrix = matrix;
	backend->matrix_dirty = 0;
	basis->factorizations++;
	basis->factor_flops += backend->orthogonal.anchor_factor_work;
	return lp_simplex_EXIT_SUCCESS;
}


int simplex_pan_basis_factorize(struct simplex_PanBasis *basis)
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
	/* Recovery path: rebuild by a sequence of bordered extensions.  This is
	 * intentionally not a symbolic/numeric monolithic refactorization. */
	pan_dynamic_sparse_qr_clear(&backend->factor);
	pan_multifrontal_qr_reset(&backend->orthogonal);
	basis->factorizations++;
	for (j = 0; j < basis->count; j++) {
		basis->column_norm[j] = simplex_pan_column_norm(
			basis->standard, basis->column[j]);
		if (basis->column_norm[j] == 0.)
			goto failure;
		if (pan_dynamic_append_column(basis, basis->column[j]) ==
		    lp_simplex_EXIT_FAILURE)
			goto failure;
	}
	backend->matrix_dirty = 1;
	pan_sync_dynamic_statistics(basis);
	return pan_maybe_rebuild_orthogonal_anchor(basis);
failure:
	pan_dynamic_sparse_qr_clear(&backend->factor);
	pan_multifrontal_qr_reset(&backend->orthogonal);
	pan_sync_dynamic_statistics(basis);
	return lp_simplex_EXIT_FAILURE;
}


int simplex_pan_basis_replace_factorized(
		struct simplex_PanBasis *basis, const int position, const int column)
{
	struct pan_SparseBackend *backend =
		(struct pan_SparseBackend *)basis->sparse_backend;
	double entering_norm;
	int leaving;
	if (position < 0 || position >= basis->count || column < 0 ||
	    column >= basis->standard->columns || basis->position[column] >= 0)
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
		    pan_dynamic_sparse_qr_move_last(&backend->factor, position,
			basis->rank_tolerance) == lp_simplex_EXIT_SUCCESS &&
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
	cholmod_dense *dense = cholmod_allocate_dense((size_t)rows, 1,
		(size_t)rows, CHOLMOD_REAL, &backend->common);
	if (dense != NULL && rows > 0)
		lp_simplex_memcpy(dense->x, value, (size_t)rows * sizeof(double));
	return dense;
}


int simplex_pan_basis_least_squares(
		struct simplex_PanBasis *basis, const double *right,
		double *solution, double *residual_norm)
{
	struct pan_SparseBackend *backend =
		(struct pan_SparseBackend *)basis->sparse_backend;
	int refinement, i;
	if (basis->count == 0) {
		lp_simplex_memcpy(basis->row_work, right,
			(size_t)basis->standard->rows * sizeof(double));
		*residual_norm = pan_dense_norm(basis->row_work,
			basis->standard->rows);
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
		if (pan_implicit_q_transpose_multiply(basis,
			basis->row_work, basis->basis_work) ==
		    lp_simplex_EXIT_FAILURE ||
		    pan_dynamic_sparse_qr_solve_upper(&backend->factor,
			basis->basis_work) ==
		    lp_simplex_EXIT_FAILURE)
			return lp_simplex_EXIT_FAILURE;
		for (i = 0; i < basis->count; i++)
			solution[i] += basis->basis_work[i] /
				basis->column_norm[i];
		basis->refinements++;
	}
	pan_basis_multiply(basis, solution, basis->row_work);
	for (i = 0; i < basis->standard->rows; i++)
		basis->row_work[i] = right[i] - basis->row_work[i];
	*residual_norm = pan_dense_norm(basis->row_work, basis->standard->rows);
	return lp_simplex_EXIT_SUCCESS;
}


int simplex_pan_basis_least_squares_qr(
		struct simplex_PanBasis *basis, const double *right,
		double *solution, double *residual_norm)
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
	scaled_solution = SuiteSparseQR_C_backslash(SPQR_ORDERING_FIXED,
		basis->rank_tolerance, backend->matrix, right_dense,
		&backend->common);
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


int simplex_pan_basis_minimum_norm(
		struct simplex_PanBasis *basis, const double *right,
		double *solution, double *residual_norm)
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
		basis->basis_work) ==
	    lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	if (pan_implicit_q_multiply(basis, basis->basis_work,
		basis->right_work, solution) == lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	for (refinement = 0; refinement < 2; refinement++) {
		for (i = 0; i < basis->count; i++)
			basis->basis_work[i] = right[i] -
				pan_column_dense_dot(basis->standard,
					basis->column[i], solution);
		for (i = 0; i < basis->count; i++)
			basis->basis_work[i] /= basis->column_norm[i];
		if (pan_dynamic_sparse_qr_solve_upper_transpose(&backend->factor,
			basis->basis_work) ==
		    lp_simplex_EXIT_FAILURE)
			return lp_simplex_EXIT_FAILURE;
		if (pan_implicit_q_multiply(basis, basis->basis_work,
			basis->row_work, basis->right_work) ==
		    lp_simplex_EXIT_FAILURE)
			return lp_simplex_EXIT_FAILURE;
		for (i = 0; i < basis->standard->rows; i++)
			solution[i] += basis->right_work[i];
		basis->refinements++;
	}
	for (i = 0; i < basis->count; i++)
		basis->basis_work[i] = right[i] -
			pan_column_dense_dot(basis->standard,
				basis->column[i], solution);
	*residual_norm = pan_dense_norm(basis->basis_work, basis->count);
	pan_record_backward_error(basis, *residual_norm,
		pan_dense_norm(right, basis->count));
	return lp_simplex_EXIT_SUCCESS;
}


int simplex_pan_basis_column_projection(
		struct simplex_PanBasis *basis, const int column,
		double *coefficient, double *residual_norm)
{
	struct pan_SparseBackend *backend =
		(struct pan_SparseBackend *)basis->sparse_backend;
	int refinement, i, k;
	if (basis->count == 0) {
		*residual_norm = simplex_pan_column_norm(basis->standard, column);
		return lp_simplex_EXIT_SUCCESS;
	}
	if (!pan_dynamic_factor_matches(basis))
		return lp_simplex_EXIT_FAILURE;
	lp_simplex_memset(basis->row_work, 0,
		(size_t)basis->standard->rows * sizeof(double));
	pan_column_axpy(basis->standard, column, 1., basis->row_work);
	if (pan_implicit_q_transpose_multiply(basis,
		basis->row_work, coefficient) == lp_simplex_EXIT_FAILURE ||
	    pan_dynamic_sparse_qr_solve_upper(&backend->factor, coefficient) ==
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
		if (pan_implicit_q_transpose_multiply(basis,
			basis->row_work, basis->basis_work) ==
		    lp_simplex_EXIT_FAILURE ||
		    pan_dynamic_sparse_qr_solve_upper(&backend->factor,
			basis->basis_work) ==
		    lp_simplex_EXIT_FAILURE)
			return lp_simplex_EXIT_FAILURE;
		for (i = 0; i < basis->count; i++)
			coefficient[i] += basis->basis_work[i] /
				basis->column_norm[i];
		basis->refinements++;
	}
	pan_basis_multiply(basis, coefficient, basis->row_work);
	for (i = 0; i < basis->standard->rows; i++)
		basis->row_work[i] = -basis->row_work[i];
	for (k = basis->standard->column_start[column];
	     k < basis->standard->column_start[column + 1]; k++)
		basis->row_work[basis->standard->row_index[k]] +=
			basis->standard->value[k];
	*residual_norm = pan_dense_norm(basis->row_work, basis->standard->rows);
	return lp_simplex_EXIT_SUCCESS;
}

#endif
