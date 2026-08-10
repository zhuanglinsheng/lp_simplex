#include "pan_dynamic_sparse_qr.h"
#include "pan_multifrontal_qr.h"

#include <lp_simplex/status.h>

#include <math.h>
#include <stdio.h>


int main(void)
{
	cholmod_common common;
	struct pan_DynamicSparseQR factor;
	struct pan_MultifrontalQR orthogonal;
	double a0[3] = {0.6, 0.8, 0.};
	double a1[3] = {0.3, 0.4, 0.8660254037844386};
	double border[3];
	double diagonal;
	double vector[3] = {0.25, -0.5, 0.75};
	double transformed[3];
	double restored[3];
	cholmod_sparse *matrix = NULL;
	int tag[2] = {0, 1};
	int failed = 0;
	cholmod_start(&common);
	if (pan_dynamic_sparse_qr_initialize(&factor) ||
	    pan_multifrontal_qr_initialize(&orthogonal, 3, &common) ||
	    pan_multifrontal_qr_append(&orthogonal, a0, 0, border, &diagonal) ||
	    pan_dynamic_sparse_qr_append(&factor, 0, border, diagonal, 1e-12) ||
	    pan_multifrontal_qr_append(&orthogonal, a1, 1, border, &diagonal) ||
	    pan_dynamic_sparse_qr_append(&factor, 1, border, diagonal, 1e-12) ||
	    pan_multifrontal_qr_apply_transpose(&orthogonal, vector,
		transformed) ||
	    pan_multifrontal_qr_apply(&orthogonal, transformed, restored))
		failed = 1;
	if (!failed && (fabs(restored[0] - vector[0]) > 1e-12 ||
	    fabs(restored[1] - vector[1]) > 1e-12 ||
	    fabs(restored[2] - vector[2]) > 1e-12))
		{
			fprintf(stderr, "restore %.17g %.17g %.17g\n",
				restored[0], restored[1], restored[2]);
			failed = 1;
		}
	if (!failed && (fabs(fabs(border[0]) - 0.5) > 1e-12 ||
	    fabs(fabs(diagonal) - 0.8660254037844386) > 1e-12))
		{
			fprintf(stderr, "border %.17g diagonal %.17g\n",
				border[0], diagonal);
			failed = 1;
		}
	if (!failed && (pan_dynamic_sparse_qr_delete(&factor, 0, 1e-12) ||
	    pan_multifrontal_qr_record_givens(&orthogonal, &factor) ||
	    pan_multifrontal_qr_apply_transpose(&orthogonal, a1, transformed) ||
	    fabs(fabs(transformed[0]) - 1.) > 1e-12 ||
	    fabs(transformed[1]) > 1e-12 || fabs(transformed[2]) > 1e-12))
		{
			fprintf(stderr, "delete qt %.17g %.17g %.17g\n",
				transformed[0], transformed[1], transformed[2]);
			failed = 1;
		}
	if (!failed && (pan_multifrontal_qr_append(&orthogonal, a0, 1,
	    border, &diagonal) || pan_dynamic_sparse_qr_append(&factor, 0,
	    border, diagonal, 1e-12) || pan_dynamic_sparse_qr_move_last(
	    &factor, 0, 1e-12) || pan_multifrontal_qr_record_givens(
	    &orthogonal, &factor) || pan_multifrontal_qr_apply_transpose(
	    &orthogonal, a0, transformed) ||
	    fabs(fabs(transformed[0]) - 1.) > 1e-12 ||
	    fabs(transformed[1]) > 1e-12))
		{
			fprintf(stderr, "move qt %.17g %.17g %.17g\n",
				transformed[0], transformed[1], transformed[2]);
			failed = 1;
		}
	if (!failed) {
		int *start;
		int *row;
		double *value;
		matrix = cholmod_allocate_sparse(3, 2, 5, 1, 1, 0,
			CHOLMOD_REAL, &common);
		if (matrix == NULL)
			failed = 1;
		else {
			start = (int *)matrix->p;
			row = (int *)matrix->i;
			value = (double *)matrix->x;
			start[0] = 0; start[1] = 2; start[2] = 5;
			row[0] = 0; row[1] = 1; row[2] = 0; row[3] = 1;
			row[4] = 2;
			value[0] = 0.6; value[1] = 0.8;
			value[2] = 0.3; value[3] = 0.4;
			value[4] = 0.8660254037844386;
			if (pan_multifrontal_qr_rebuild(&orthogonal, matrix,
				&factor, tag, 1e-12) ||
			    pan_multifrontal_qr_apply_transpose(&orthogonal,
				a0, transformed) || fabs(transformed[1]) > 1e-12 ||
			    fabs(transformed[2]) > 1e-12 ||
			    pan_multifrontal_qr_apply(&orthogonal, transformed,
				restored) || fabs(restored[0] - a0[0]) > 1e-12 ||
			    fabs(restored[1] - a0[1]) > 1e-12)
				failed = 1;
		}
	}
	cholmod_free_sparse(&matrix, &common);
	pan_multifrontal_qr_destroy(&orthogonal);
	pan_dynamic_sparse_qr_destroy(&factor);
	cholmod_finish(&common);
	if (failed)
		fprintf(stderr, "multifrontal transform invariant failed\n");
	return failed;
}
