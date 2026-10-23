/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file test_pan_qr_stress.c
 * @brief Audit QR against original CSC columns, including nontrivial
 * permutations.
 *
 * Orthogonality alone would not catch pairing a correct Q with the wrong R.
 */

#include "pan_dynamic_sparse_qr.h"
#include "pan_multifrontal_qr.h"

#include <lp_simplex/status.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>


static int run_stress(const int m, const int n)
{
	cholmod_common common;
	cholmod_sparse *matrix;
	struct pan_DynamicSparseQR factor;
	struct pan_MultifrontalQR qr;
	double *column, *qt, *restored;
	int *tag, *start, *row;
	double *value;
	int i, j, k, valid = 1;
	double worst = 0.;
	cholmod_start(&common);

	matrix = cholmod_allocate_sparse((size_t)m, (size_t)n, (size_t)m * n, 1,
		1, 0, CHOLMOD_REAL, &common);
	column = (double *)calloc((size_t)m * 3, sizeof(double));
	tag = (int *)malloc((size_t)n * sizeof(int));
	if (!matrix || !column || !tag)
		return 0;
	qt = column + m;
	restored = qt + m;
	start = (int *)matrix->p;
	row = (int *)matrix->i;
	value = (double *)matrix->x;
	k = 0;
	for (j = 0; j < n; j++) {
		start[j] = k;
		tag[j] = j;
		for (i = 0; i < m; i++) {
			double a = sin((i + 3.) * (j + 2.)) * .07;

			if (i == (j * 17 + 11) % m)
				a += 2.;
			if (i == (j * 7 + 5) % m)
				a += .3;
			if ((i * 7 + j * 3) % 5 == 0 || fabs(a) > .1) {
				row[k] = i;
				value[k++] = a;
			}
		}
	}
	start[n] = k;
	if (pan_dynamic_sparse_qr_initialize(&factor) ||
		pan_multifrontal_qr_initialize(&qr, m, &common) ||
		pan_multifrontal_qr_rebuild(&qr, matrix, &factor, tag, 1e-12))
		return 0;
	for (j = 0; j < n; j++) {
		int original = factor.tag[j];

		for (i = 0; i < m; i++)
			column[i] = 0.;
		for (k = start[original]; k < start[original + 1]; k++)
			column[row[k]] = value[k];
		if (pan_multifrontal_qr_apply_transpose(&qr, column, qt) ||
			pan_multifrontal_qr_apply(&qr, qt, restored))
			valid = 0;
		for (i = 0; i < m; i++) {
			double r = 0., error;
			int node = factor.column_head[j];

			while (node >= 0) {
				if (factor.node[node].owner.row == i)
					r = factor.node[node].value;
				node = factor.node[node].column_next;
			}
			error = fmax(
				fabs(qt[i] - r), fabs(restored[i] - column[i]));
			worst = fmax(worst, error);
		}
	}
	valid = valid && worst < 1e-11;
	printf("sparse QR audit m=%d n=%d error=%.6g %s\n", m, n, worst,
		valid ? "passed" : "failed");
	pan_multifrontal_qr_destroy(&qr);
	pan_dynamic_sparse_qr_destroy(&factor);
	cholmod_free_sparse(&matrix, &common);
	cholmod_finish(&common);

	free(column);
	free(tag);
	return valid;
}


int main(void)
{
	return run_stress(53, 17) && run_stress(127, 61) && run_stress(263, 113)
		       ? 0
		       : 1;
}
