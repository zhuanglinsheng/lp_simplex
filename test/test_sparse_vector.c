/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file test_sparse_vector.c
 * @brief Packing, generation reuse and fused maximum validation.
 */

#include "simplex_sparse_vector.h"

#include <lp_simplex/status.h>

#include <math.h>
#include <stdio.h>

static int packing(void)
{
	struct simplex_SparseVector vector;
	double dense[5] = {0., -3., .5, 0., 2.}, maximum = -1.;
	int good;

	if (simplex_sparse_vector_create(&vector, 5) != lp_simplex_EXIT_SUCCESS)
		return 0;
	good = simplex_sparse_vector_pack_maximum(&vector, dense, .5,
		       &maximum) == lp_simplex_EXIT_SUCCESS &&
	       maximum == 3. && vector.count == 2 && vector.index[0] == 1 &&
	       vector.index[1] == 4 && vector.value[0] == -3.;
	dense[1] = 0.;
	dense[3] = -4.;
	good = good &&
	       simplex_sparse_vector_pack_maximum(&vector, dense, 0.,
		       &maximum) == lp_simplex_EXIT_SUCCESS &&
	       maximum == 4. && vector.count == 3 &&
	       simplex_sparse_vector_get(&vector, 1) == 0. &&
	       vector.dense[1] == 0. && vector.slot[3] == 1;
	dense[2] = NAN;
	good = good &&
	       simplex_sparse_vector_pack_maximum(&vector, dense, 0.,
		       &maximum) == lp_simplex_EXIT_SUCCESS &&
	       isinf(maximum) && isnan(simplex_sparse_vector_get(&vector, 2));
	dense[2] = HUGE_VAL;
	good = good &&
	       simplex_sparse_vector_pack_maximum(&vector, dense, 0.,
		       &maximum) == lp_simplex_EXIT_SUCCESS &&
	       isinf(maximum);
	dense[2] = dense[3] = dense[4] = 0.;
	good = good &&
	       simplex_sparse_vector_pack_maximum(&vector, dense, 0.,
		       &maximum) == lp_simplex_EXIT_SUCCESS &&
	       maximum == 0. && vector.count == 0;
	simplex_sparse_vector_destroy(&vector);
	return good;
}

int main(void)
{
	if (!packing()) {
		puts("Sparse vector packing regression failed");
		return 1;
	}
	puts("Sparse vector packing handles reuse, tolerance and nonfinite values");
	return 0;
}
