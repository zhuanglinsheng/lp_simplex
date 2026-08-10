/* Unit tests for Pan's dynamic sparse QR update engine. */
#include "pan_dynamic_sparse_qr.h"

#include <lp_simplex/status.h>

#include <math.h>
#include <stdio.h>


static int pan_check_solution(
		struct pan_DynamicSparseQR *factor, const double *gram,
		const double *expected, const int count)
{
	double right[4];
	int i;
	int j;
	for (i = 0; i < count; i++) {
		right[i] = 0.;
		for (j = 0; j < count; j++)
			right[i] += gram[i * count + j] * expected[j];
	}
	if (pan_dynamic_sparse_qr_solve(factor, right) ==
	    lp_simplex_EXIT_FAILURE)
		return 0;
	for (i = 0; i < count; i++)
		if (fabs(right[i] - expected[i]) > 1e-10)
			return 0;
	return 1;
}


int main(void)
{
	struct pan_DynamicSparseQR factor;
	double border1[1] = {0.6};
	double border2[2] = {0.3, 0.4};
	double border3[2] = {0., 0.9078412990032035};
	double duplicate[3] = {1., 0., 0.};
	double gram3[9] = {
		1., 0.6, 0.3,
		0.6, 1., 0.5,
		0.3, 0.5, 1.
	};
	double after_delete[4] = {1., 0.3, 0.3, 1.};
	double after_insert[9] = {
		1., 0., 0.3,
		0., 1., 0.8660254037844386,
		0.3, 0.8660254037844386, 1.
	};
	double expected3[3] = {1., -2., 0.5};
	double expected2[2] = {-0.25, 2.};
	double expected_insert[3] = {0.5, -1., 2.};
	int failed = 0;
	if (pan_dynamic_sparse_qr_initialize(&factor) ||
	    pan_dynamic_sparse_qr_append(&factor, 10, NULL, 1., 1e-12) ||
	    pan_dynamic_sparse_qr_append(&factor, 11, border1, 0.8, 1e-12) ||
	    pan_dynamic_sparse_qr_append(&factor, 12, border2,
		0.8660254037844386, 1e-12) ||
	    pan_dynamic_sparse_qr_validate(&factor) ||
	    !pan_check_solution(&factor, gram3, expected3, 3)) {
		fprintf(stderr, "dynamic QR bordered extension failed\n");
		failed = 1;
	}
	if (!failed &&
	    (pan_dynamic_sparse_qr_delete(&factor, 1, 1e-12) ||
	     pan_dynamic_sparse_qr_validate(&factor) ||
	     factor.count != 2 || factor.tag[0] != 10 || factor.tag[1] != 12 ||
	     !pan_check_solution(&factor, after_delete, expected2, 2))) {
		fprintf(stderr, "dynamic QR downdate failed\n");
		failed = 1;
	}
	if (!failed &&
	    (pan_dynamic_sparse_qr_append(&factor, 13, border3,
		0.4193139346887673, 1e-12) ||
	     pan_dynamic_sparse_qr_move_last(&factor, 1, 1e-12) ||
	     pan_dynamic_sparse_qr_validate(&factor) ||
	     factor.tag[0] != 10 || factor.tag[1] != 13 || factor.tag[2] != 12 ||
	     !pan_check_solution(&factor, after_insert, expected_insert, 3))) {
		fprintf(stderr, "dynamic QR insertion/retriangularization failed\n");
		failed = 1;
	}
	if (!failed &&
	    (pan_dynamic_sparse_qr_append(&factor, 14, duplicate, 0.,
		1e-12) != lp_simplex_EXIT_FAILURE || factor.count != 3)) {
		fprintf(stderr, "dynamic QR rank rejection was not transactional\n");
		failed = 1;
	}
	if (!failed && (factor.extensions < 4 || factor.downdates != 1 ||
	    factor.rotations == 0)) {
		fprintf(stderr, "dynamic QR update statistics missing\n");
		failed = 1;
	}
	pan_dynamic_sparse_qr_destroy(&factor);
	if (!failed)
		printf("dynamic sparse QR tests passed\n");
	return failed;
}
