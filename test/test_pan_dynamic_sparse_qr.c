/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file test_pan_dynamic_sparse_qr.c
 * @brief Unit tests for Pan's dynamic sparse QR update engine.
 */

#include "pan_dynamic_sparse_qr.h"

#include <lp_simplex/status.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>


static int border_cache_stress(int count)
{
	struct pan_DynamicSparseQR factor;
	double *storage, *border, *expected, *right, *up, *down;
	double worst = 0.;
	clock_t start = clock();
	int i, j, good = 1;

	storage = (double *)calloc((size_t)count * 5, sizeof(double));
	if (!storage)
		return 0;
	border = storage;
	expected = border + count;
	right = expected + count;
	up = right + count;
	down = up + count;
	pan_dynamic_sparse_qr_initialize(&factor);
	/* Also cover an empty factor and capacity growth with a live cache. */
	if (pan_dynamic_sparse_qr_solve_upper(&factor, right) ||
		pan_dynamic_sparse_qr_solve_upper_transpose(&factor, right))
		good = 0;
	for (j = 0; j < count && good; j++) {
		expected[j] = sin(j + 1.);
		down[j] = 2. * expected[j];
		up[j] = 2. * expected[j];
		for (i = 0; i < j; i++) {
			border[i] = (i + j) % 3 == 0
					    ? 0.
					    : .01 * sin((i + 1.) * (j + 3.));
			up[i] += border[i] * expected[j];
			down[j] += border[i] * expected[i];
		}
		if (pan_dynamic_sparse_qr_append(
			    &factor, j, border, 2., 1e-12)) {
			good = 0;
			break;
		}
		for (i = 0; i <= j; i++)
			right[i] = up[i];
		if (pan_dynamic_sparse_qr_solve_upper(&factor, right))
			good = 0;
		for (i = 0; i <= j; i++)
			worst = fmax(worst, fabs(right[i] - expected[i]));
		for (i = 0; i <= j; i++)
			right[i] = down[i];
		if (pan_dynamic_sparse_qr_solve_upper_transpose(&factor, right))
			good = 0;
		for (i = 0; i <= j; i++)
			worst = fmax(worst, fabs(right[i] - expected[i]));
	}
	good = good && worst < 1e-11 && factor.numeric_packs > 0 &&
	       factor.numeric_packs <= 1 + count / 32 &&
	       factor.border_solves > count &&
	       pan_dynamic_sparse_qr_validate(&factor) == 0;
	printf("QR append cache n=%d packs=%ld copied=%ld border-solves=%ld error=%.6g seconds=%.6f %s\n",
		count, factor.numeric_packs, factor.numeric_packed_nodes,
		factor.border_solves, worst,
		(double)(clock() - start) / CLOCKS_PER_SEC,
		good ? "passed" : "failed");
	pan_dynamic_sparse_qr_clear(&factor);
	if (pan_dynamic_sparse_qr_solve_upper(&factor, right))
		good = 0;
	pan_dynamic_sparse_qr_destroy(&factor);
	free(storage);
	return good;
}


static int pan_check_solution(struct pan_DynamicSparseQR *factor,
	const double *gram, const double *expected, const int count)
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


int main(int argc, char **argv)
{
	struct pan_DynamicSparseQR factor;
	double border1[1] = {0.6};
	double border2[2] = {0.3, 0.4};
	double border3[2] = {0., 0.9078412990032035};
	double duplicate[3] = {1., 0., 0.};
	double gram3[9] = {1., 0.6, 0.3, 0.6, 1., 0.5, 0.3, 0.5, 1.};
	double after_delete[4] = {1., 0.3, 0.3, 1.};
	double after_insert[9] = {1., 0., 0.3, 0., 1., 0.8660254037844386, 0.3,
		0.8660254037844386, 1.};
	double expected3[3] = {1., -2., 0.5};
	double expected2[2] = {-0.25, 2.};
	double expected_insert[3] = {0.5, -1., 2.};
	int failed = 0;

	if (argc == 2)
		return border_cache_stress(
			       atoi(argv[1]) > 0 ? atoi(argv[1]) : 129)
			       ? 0
			       : 1;
	if (!border_cache_stress(129))
		return 1;
	if (pan_dynamic_sparse_qr_initialize(&factor) ||
		pan_dynamic_sparse_qr_append(&factor, 10, NULL, 1., 1e-12) ||
		pan_dynamic_sparse_qr_append(
			&factor, 11, border1, 0.8, 1e-12) ||
		pan_dynamic_sparse_qr_append(
			&factor, 12, border2, 0.8660254037844386, 1e-12) ||
		pan_dynamic_sparse_qr_validate(&factor) ||
		!pan_check_solution(&factor, gram3, expected3, 3)) {
		fprintf(stderr, "dynamic QR bordered extension failed\n");
		failed = 1;
	}
	if (!failed && (pan_dynamic_sparse_qr_delete(&factor, 1, 1e-12) ||
			       pan_dynamic_sparse_qr_validate(&factor) ||
			       factor.count != 2 || factor.tag[0] != 10 ||
			       factor.tag[1] != 12 ||
			       !pan_check_solution(
				       &factor, after_delete, expected2, 2))) {
		fprintf(stderr, "dynamic QR downdate failed\n");
		failed = 1;
	}
	if (!failed &&
		(pan_dynamic_sparse_qr_append(
			 &factor, 13, border3, 0.4193139346887673, 1e-12) ||
			pan_dynamic_sparse_qr_move_last(&factor, 1, 1e-12) ||
			pan_dynamic_sparse_qr_validate(&factor) ||
			factor.tag[0] != 10 || factor.tag[1] != 13 ||
			factor.tag[2] != 12 ||
			!pan_check_solution(
				&factor, after_insert, expected_insert, 3))) {
		fprintf(stderr,
			"dynamic QR insertion/retriangularization failed\n");
		failed = 1;
	}
	if (!failed && (pan_dynamic_sparse_qr_append(&factor, 14, duplicate, 0.,
				1e-12) != lp_simplex_EXIT_FAILURE ||
			       factor.count != 3)) {
		fprintf(stderr,
			"dynamic QR rank rejection was not transactional\n");
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
