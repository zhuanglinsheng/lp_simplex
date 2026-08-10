/*
 * Orthogonal part of Pan's dynamic sparse QR.
 *
 * SuiteSparseQR supplies a multifrontal sparse-Householder anchor.  Structural
 * modifications after that anchor are represented exactly by right-side
 * Givens and tail Householder transformations.  Applying Q or Q' therefore
 * never forms R^-1 R^-T.  A measured work crossover, not an iteration count,
 * decides when the accumulated local transformations should be folded into a
 * new multifrontal anchor.
 */
#include "pan_multifrontal_qr.h"
#include "utils.h"

#include <lp_simplex/status.h>

#include <float.h>
#include <math.h>


static int pan_transform_reserve(
		struct pan_MultifrontalQR *qr, const int count)
{
	struct pan_QRTransform *replacement;
	int capacity;
	if (count <= qr->transform_capacity)
		return lp_simplex_EXIT_SUCCESS;
	capacity = qr->transform_capacity > 0 ? 2 * qr->transform_capacity :
		(qr->dimension > 0 ? qr->dimension : 8);
	while (capacity < count)
		capacity *= 2;
	replacement = (struct pan_QRTransform *)lp_simplex_malloc(
		(size_t)capacity * sizeof(*replacement));
	if (replacement == NULL)
		return lp_simplex_EXIT_FAILURE;
	if (qr->transform_count > 0)
		lp_simplex_memcpy(replacement, qr->transform,
			(size_t)qr->transform_count * sizeof(*replacement));
	lp_simplex_free(qr->transform);
	qr->transform = replacement;
	qr->transform_capacity = capacity;
	return lp_simplex_EXIT_SUCCESS;
}


static int pan_value_reserve(struct pan_MultifrontalQR *qr, const int count)
{
	double *replacement;
	int capacity;
	if (count <= qr->value_capacity)
		return lp_simplex_EXIT_SUCCESS;
	capacity = qr->value_capacity > 0 ? 2 * qr->value_capacity :
		(qr->dimension > 0 ? qr->dimension : 8);
	while (capacity < count)
		capacity *= 2;
	replacement = (double *)lp_simplex_malloc(
		(size_t)capacity * sizeof(double));
	if (replacement == NULL)
		return lp_simplex_EXIT_FAILURE;
	if (qr->value_count > 0)
		lp_simplex_memcpy(replacement, qr->transform_value,
			(size_t)qr->value_count * sizeof(double));
	lp_simplex_free(qr->transform_value);
	qr->transform_value = replacement;
	qr->value_capacity = capacity;
	return lp_simplex_EXIT_SUCCESS;
}


static void pan_apply_householder(
		const double *v, const int first, const int length,
		const double tau, double *x)
{
	double dot = 0.;
	int i;
	if (tau == 0.)
		return;
	for (i = 0; i < length; i++)
		dot += v[i] * x[first + i];
	dot *= tau;
	for (i = 0; i < length; i++)
		x[first + i] -= v[i] * dot;
}


static void pan_apply_transform(
		const struct pan_MultifrontalQR *qr,
		const struct pan_QRTransform *transform,
		const int transpose, double *x)
{
	if (transform->kind == pan_QR_HOUSEHOLDER) {
		pan_apply_householder(qr->transform_value +
			transform->value_offset, transform->first,
			transform->length, transform->tau, x);
	} else {
		int first = transform->first;
		double a = x[first];
		double b = x[first + 1];
		double c = transform->cosine;
		double s = transform->sine;
		if (transpose) {
			x[first] = c * a + s * b;
			x[first + 1] = -s * a + c * b;
		} else {
			x[first] = c * a - s * b;
			x[first + 1] = s * a + c * b;
		}
	}
}


static int pan_apply_anchor_transpose(
		const struct pan_MultifrontalQR *qr,
		const double *vector, double *result)
{
	const cholmod_sparse *h = qr->householder;
	int i;
	if (h == NULL) {
		lp_simplex_memcpy(result, vector,
			(size_t)qr->dimension * sizeof(double));
		return lp_simplex_EXIT_SUCCESS;
	}
	{
		const int *start = (const int *)h->p;
		const int *row = (const int *)h->i;
		const double *value = (const double *)h->x;
		const double *tau = (const double *)qr->householder_tau->x;
		int column;
		for (i = 0; i < qr->dimension; i++)
			result[qr->row_permutation[i]] = vector[i];
		for (column = 0; column < (int)h->ncol; column++) {
			double dot = 0.;
			int k;
			for (k = start[column]; k < start[column + 1]; k++)
				dot += value[k] * result[row[k]];
			dot *= tau[column];
			for (k = start[column]; k < start[column + 1]; k++)
				result[row[k]] -= value[k] * dot;
		}
	}
	return lp_simplex_EXIT_SUCCESS;
}


static int pan_apply_anchor(
		const struct pan_MultifrontalQR *qr,
		const double *vector, double *result)
{
	const cholmod_sparse *h = qr->householder;
	int i;
	if (h == NULL) {
		lp_simplex_memcpy(result, vector,
			(size_t)qr->dimension * sizeof(double));
		return lp_simplex_EXIT_SUCCESS;
	}
	{
		const int *start = (const int *)h->p;
		const int *row = (const int *)h->i;
		const double *value = (const double *)h->x;
		const double *tau = (const double *)qr->householder_tau->x;
		int column;
		lp_simplex_memcpy(qr->work, vector,
			(size_t)qr->dimension * sizeof(double));
		for (column = (int)h->ncol - 1; column >= 0; column--) {
			double dot = 0.;
			int k;
			for (k = start[column]; k < start[column + 1]; k++)
				dot += value[k] * qr->work[row[k]];
			dot *= tau[column];
			for (k = start[column]; k < start[column + 1]; k++)
				qr->work[row[k]] -= value[k] * dot;
		}
		for (i = 0; i < qr->dimension; i++)
			result[i] = qr->work[qr->row_permutation[i]];
	}
	return lp_simplex_EXIT_SUCCESS;
}


int pan_multifrontal_qr_initialize(
		struct pan_MultifrontalQR *qr, const int dimension,
		cholmod_common *common)
{
	lp_simplex_memset(qr, 0, sizeof(*qr));
	if (dimension <= 0 || common == NULL)
		return lp_simplex_EXIT_FAILURE;
	qr->dimension = dimension;
	qr->common = common;
	qr->work = (double *)lp_simplex_malloc(
		(size_t)dimension * sizeof(double));
	qr->parent = (int *)lp_simplex_malloc((size_t)dimension * sizeof(int));
	qr->front_size = (int *)lp_simplex_malloc(
		(size_t)dimension * sizeof(int));
	qr->dirty_front = (unsigned char *)lp_simplex_malloc(
		(size_t)dimension);
	if (qr->work == NULL || qr->parent == NULL || qr->front_size == NULL ||
	    qr->dirty_front == NULL) {
		pan_multifrontal_qr_destroy(qr);
		return lp_simplex_EXIT_FAILURE;
	}
	lp_simplex_memset(qr->dirty_front, 0, (size_t)dimension);
	return lp_simplex_EXIT_SUCCESS;
}


void pan_multifrontal_qr_clear_updates(struct pan_MultifrontalQR *qr)
{
	qr->transform_count = 0;
	qr->value_count = 0;
	qr->pending_apply_work = 0;
	qr->local_update_work = 0;
	qr->backward_error = 0.;
	lp_simplex_memset(qr->dirty_front, 0, (size_t)qr->dimension);
}


void pan_multifrontal_qr_reset(struct pan_MultifrontalQR *qr)
{
	cholmod_free_sparse(&qr->householder, qr->common);
	cholmod_free_dense(&qr->householder_tau, qr->common);
	if (qr->row_permutation != NULL) {
		SuiteSparse_free(qr->row_permutation);
		qr->row_permutation = NULL;
	}
	qr->anchor_apply_work = 0;
	qr->anchor_factor_work = 0.;
	pan_multifrontal_qr_clear_updates(qr);
}


void pan_multifrontal_qr_destroy(struct pan_MultifrontalQR *qr)
{
	if (qr == NULL)
		return;
	if (qr->common != NULL) {
		cholmod_free_sparse(&qr->householder, qr->common);
		cholmod_free_dense(&qr->householder_tau, qr->common);
		if (qr->row_permutation != NULL)
			SuiteSparse_free(qr->row_permutation);
	}
	lp_simplex_free(qr->transform);
	lp_simplex_free(qr->transform_value);
	lp_simplex_free(qr->parent);
	lp_simplex_free(qr->front_size);
	lp_simplex_free(qr->dirty_front);
	lp_simplex_free(qr->work);
	lp_simplex_memset(qr, 0, sizeof(*qr));
}


void pan_multifrontal_qr_refresh_structure(
		struct pan_MultifrontalQR *qr,
		const struct pan_DynamicSparseQR *factor)
{
	int i;
	qr->maximum_front = 0;
	qr->factor_nonzeros = 0;
	for (i = 0; i < factor->count; i++) {
		const struct pan_DynamicSparseRow *row = factor->row + i;
		int parent = -1;
		int node = row->head;
		while (node >= 0) {
			if (factor->node[node].column > i) {
				parent = factor->node[node].column;
				break;
			}
			node = factor->node[node].row_next;
		}
		qr->parent[i] = parent;
		qr->front_size[i] = row->count;
		qr->factor_nonzeros += row->count;
		if (row->count > qr->maximum_front)
			qr->maximum_front = row->count;
	}
}


int pan_multifrontal_qr_rebuild(
		struct pan_MultifrontalQR *qr, cholmod_sparse *matrix,
		struct pan_DynamicSparseQR *factor, const int *tag,
		const double rank_tolerance)
{
	cholmod_sparse *r = NULL;
	cholmod_sparse *h = NULL;
	cholmod_dense *tau = NULL;
	int32_t *permutation = NULL;
	int32_t *column_permutation = NULL;
	int32_t rank;
	int i;
	if (matrix == NULL || matrix->nrow != (size_t)qr->dimension)
		return lp_simplex_EXIT_FAILURE;
	rank = SuiteSparseQR_i_C(SPQR_ORDERING_FIXED, rank_tolerance,
		(int32_t)matrix->ncol, 0, matrix, NULL, NULL, NULL, NULL,
		&r, &column_permutation, &h, &permutation, &tau, qr->common);
	if (rank != (int32_t)matrix->ncol || r == NULL || h == NULL ||
	    tau == NULL || permutation == NULL)
		goto failure;
	if (column_permutation != NULL)
		for (i = 0; i < (int)matrix->ncol; i++)
			if (column_permutation[i] != i)
				goto failure;
	if (r->itype != CHOLMOD_INT ||
	    pan_dynamic_sparse_qr_load_csc(factor, tag, (int)matrix->ncol,
		(const int *)r->p, (const int *)r->i, (const double *)r->x,
		rank_tolerance) == lp_simplex_EXIT_FAILURE)
		goto failure;
	cholmod_free_sparse(&qr->householder, qr->common);
	cholmod_free_dense(&qr->householder_tau, qr->common);
	if (qr->row_permutation != NULL)
		SuiteSparse_free(qr->row_permutation);
	qr->householder = h;
	qr->householder_tau = tau;
	qr->row_permutation = permutation;
	h = NULL;
	tau = NULL;
	permutation = NULL;
	qr->anchor_factor_work = qr->common->SPQR_flopcount;
	qr->anchor_apply_work = 4 * (long long)qr->householder->nzmax;
	qr->anchor_rebuilds++;
	pan_multifrontal_qr_clear_updates(qr);
	pan_multifrontal_qr_refresh_structure(qr, factor);
	cholmod_free_sparse(&r, qr->common);
	if (column_permutation != NULL)
		SuiteSparse_free(column_permutation);
	return lp_simplex_EXIT_SUCCESS;
failure:
	cholmod_free_sparse(&r, qr->common);
	cholmod_free_sparse(&h, qr->common);
	cholmod_free_dense(&tau, qr->common);
	if (permutation != NULL)
		SuiteSparse_free(permutation);
	if (column_permutation != NULL)
		SuiteSparse_free(column_permutation);
	return lp_simplex_EXIT_FAILURE;
}


int pan_multifrontal_qr_apply_transpose(
		struct pan_MultifrontalQR *qr, const double *vector,
		double *result)
{
	int i;
	if (pan_apply_anchor_transpose(qr, vector, result) ==
	    lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	for (i = 0; i < qr->transform_count; i++)
		pan_apply_transform(qr, qr->transform + i, 1, result);
	for (i = 0; i < qr->transform_count; i++)
		qr->pending_apply_work += qr->transform[i].kind == pan_QR_GIVENS
			? 6 : 4 * qr->transform[i].length;
	return lp_simplex_EXIT_SUCCESS;
}


int pan_multifrontal_qr_apply(
		struct pan_MultifrontalQR *qr, const double *vector,
		double *result)
{
	int i;
	lp_simplex_memcpy(qr->work, vector,
		(size_t)qr->dimension * sizeof(double));
	for (i = qr->transform_count - 1; i >= 0; i--)
		pan_apply_transform(qr, qr->transform + i, 0, qr->work);
	if (pan_apply_anchor(qr, qr->work, result) == lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	for (i = 0; i < qr->transform_count; i++)
		qr->pending_apply_work += qr->transform[i].kind == pan_QR_GIVENS
			? 6 : 4 * qr->transform[i].length;
	return lp_simplex_EXIT_SUCCESS;
}


int pan_multifrontal_qr_append(
		struct pan_MultifrontalQR *qr, const double *column,
		const int active_count, double *border, double *diagonal)
{
	struct pan_QRTransform *transform;
	double alpha;
	double tail_norm = 0.;
	double beta;
	double tau;
	double norm_error;
	int length;
	int offset;
	int i;
	if (active_count < 0 || active_count >= qr->dimension ||
	    pan_multifrontal_qr_apply_transpose(qr, column, qr->work) ==
	    lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	for (i = 0; i < active_count; i++)
		border[i] = qr->work[i];
	length = qr->dimension - active_count;
	alpha = qr->work[active_count];
	for (i = active_count + 1; i < qr->dimension; i++)
		tail_norm = hypot(tail_norm, qr->work[i]);
	if (tail_norm == 0.) {
		beta = alpha;
		tau = 0.;
	} else {
		beta = -copysign(hypot(alpha, tail_norm), alpha);
		tau = (beta - alpha) / beta;
	}
	if (pan_transform_reserve(qr, qr->transform_count + 1) ==
	    lp_simplex_EXIT_FAILURE ||
	    pan_value_reserve(qr, qr->value_count + length) ==
	    lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	offset = qr->value_count;
	qr->transform_value[offset] = 1.;
	if (tau == 0.)
		for (i = 1; i < length; i++)
			qr->transform_value[offset + i] = 0.;
	else
		for (i = 1; i < length; i++)
			qr->transform_value[offset + i] =
				qr->work[active_count + i] / (alpha - beta);
	transform = qr->transform + qr->transform_count++;
	transform->kind = pan_QR_HOUSEHOLDER;
	transform->first = active_count;
	transform->length = length;
	transform->value_offset = offset;
	transform->cosine = 0.;
	transform->sine = 0.;
	transform->tau = tau;
	qr->value_count += length;
	qr->local_update_work += 4 * (long long)length;
	qr->local_retriangularizations++;
	qr->block_householders++;
	norm_error = 0.;
	for (i = 0; i < length; i++)
		norm_error += qr->transform_value[offset + i] *
			qr->transform_value[offset + i];
	norm_error = tau == 0. ? 0. : fabs(tau * norm_error - 2.);
	if (norm_error > qr->orthogonal_error)
		qr->orthogonal_error = norm_error;
	*diagonal = beta;
	return lp_simplex_EXIT_SUCCESS;
}


int pan_multifrontal_qr_record_givens(
		struct pan_MultifrontalQR *qr,
		const struct pan_DynamicSparseQR *factor)
{
	int i;
	if (pan_transform_reserve(qr,
		qr->transform_count + factor->event_count) ==
	    lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	for (i = 0; i < factor->event_count; i++) {
		struct pan_QRTransform *transform =
			qr->transform + qr->transform_count++;
		transform->kind = pan_QR_GIVENS;
		transform->first = factor->event_row[i];
		transform->length = 2;
		transform->value_offset = 0;
		transform->cosine = factor->event_cosine[i];
		transform->sine = factor->event_sine[i];
		transform->tau = 0.;
		qr->local_update_work += 6;
	}
	return lp_simplex_EXIT_SUCCESS;
}


void pan_multifrontal_qr_mark_path(
		struct pan_MultifrontalQR *qr, const int first)
{
	int front = first;
	while (front >= 0 && front < qr->dimension &&
	       !qr->dirty_front[front]) {
		qr->dirty_front[front] = 1;
		front = qr->parent[front];
	}
}


int pan_multifrontal_qr_should_rebuild(
		const struct pan_MultifrontalQR *qr)
{
	/* Exact work crossover: another use of the accumulated transforms costs
	 * more arithmetic than the last measured multifrontal factorization. */
	long long next_apply = qr->anchor_apply_work + qr->local_update_work;
	if (qr->transform_count == 0)
		return 0;
	if (qr->backward_error > DBL_EPSILON * (double)(qr->dimension + 1))
		return 1;
	if (qr->anchor_factor_work <= 0.)
		return qr->pending_apply_work >= qr->local_update_work;
	next_apply = qr->local_update_work;
	return (double)(qr->pending_apply_work + next_apply) >=
		qr->anchor_factor_work;
}
