/* True orthogonal-transform backend layered over the dynamic sparse R. */
#ifndef LP_SIMPLEX_PAN_MULTIFRONTAL_QR_H
#define LP_SIMPLEX_PAN_MULTIFRONTAL_QR_H

#include "pan_dynamic_sparse_qr.h"

#include <SuiteSparseQR_C.h>


enum pan_QRTransformKind {
	pan_QR_GIVENS = 1,
	pan_QR_HOUSEHOLDER = 2
};


struct pan_QRTransform {
	int kind;
	int first;
	int length;
	int value_offset;
	double cosine;
	double sine;
	double tau;
};


struct pan_MultifrontalQR {
	cholmod_common *common;
	cholmod_sparse *householder;
	cholmod_dense *householder_tau;
	int32_t *row_permutation;
	struct pan_QRTransform *transform;
	double *transform_value;
	int *parent;
	int *front_size;
	unsigned char *dirty_front;
	double *work;
	int dimension;
	int transform_count;
	int transform_capacity;
	int value_count;
	int value_capacity;
	long long anchor_apply_work;
	long long pending_apply_work;
	long long local_update_work;
	double anchor_factor_work;
	double orthogonal_error;
	double backward_error;
	long anchor_rebuilds;
	long local_retriangularizations;
	long block_householders;
	long maximum_front;
	long factor_nonzeros;
};


int pan_multifrontal_qr_initialize(
		struct pan_MultifrontalQR *qr, int dimension,
		cholmod_common *common);

void pan_multifrontal_qr_destroy(struct pan_MultifrontalQR *qr);

void pan_multifrontal_qr_clear_updates(struct pan_MultifrontalQR *qr);

void pan_multifrontal_qr_reset(struct pan_MultifrontalQR *qr);

int pan_multifrontal_qr_rebuild(
		struct pan_MultifrontalQR *qr, cholmod_sparse *matrix,
		struct pan_DynamicSparseQR *factor, const int *tag,
		double rank_tolerance);

int pan_multifrontal_qr_apply_transpose(
		struct pan_MultifrontalQR *qr, const double *vector,
		double *result);

int pan_multifrontal_qr_apply(
		struct pan_MultifrontalQR *qr, const double *vector,
		double *result);

int pan_multifrontal_qr_append(
		struct pan_MultifrontalQR *qr, const double *column,
		int active_count, double *border, double *diagonal);

int pan_multifrontal_qr_record_givens(
		struct pan_MultifrontalQR *qr,
		const struct pan_DynamicSparseQR *factor);

void pan_multifrontal_qr_mark_path(
		struct pan_MultifrontalQR *qr, int first);

void pan_multifrontal_qr_refresh_structure(
		struct pan_MultifrontalQR *qr,
		const struct pan_DynamicSparseQR *factor);

int pan_multifrontal_qr_should_rebuild(
		const struct pan_MultifrontalQR *qr);

#endif
