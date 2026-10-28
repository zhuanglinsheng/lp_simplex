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

static int compact_reuse(void)
{
	struct simplex_SparseVector vector;
	double tiny = nextafter(0., 1.);
	int valid;
	if (simplex_sparse_vector_create(&vector, 5))
		return 0;
	(void)simplex_sparse_vector_set(&vector, 4, 3.);
	(void)simplex_sparse_vector_set(&vector, 1, tiny);
	(void)simplex_sparse_vector_set(&vector, 3, -2.);
	(void)simplex_sparse_vector_set(&vector, 4, 0.);
	simplex_sparse_vector_compact(&vector);
	simplex_sparse_vector_sort(&vector);
	valid = vector.count == 2 && vector.index[0] == 1 &&
		vector.value[0] == tiny && vector.index[1] == 3 &&
		vector.value[1] == -2.;
	(void)simplex_sparse_vector_set(&vector, 4, 7.);
	valid = valid && vector.count == 3 && vector.slot[4] == 2 &&
		vector.dense[4] == 7.;
	simplex_sparse_vector_clear(&vector);
	valid = valid && vector.dense[1] == 0. && vector.dense[3] == 0. &&
		vector.dense[4] == 0.;
	simplex_sparse_vector_destroy(&vector);
	return valid;
}

static int integer_ordering(void)
{
	struct simplex_SparseVector vector;
	const int counts[] = {1, 16, 32, 33, 257, 4097};
	int test, i, repeat, valid = 1;
	if (simplex_sparse_vector_create(&vector, 262147))
		return 0;
	for (test = 0; test < 6; test++) {
		simplex_sparse_vector_clear(&vector);
		for (i = counts[test] - 1; i >= 0; i--) {
			int row = (i * 7919) % vector.dimension;
			(void)simplex_sparse_vector_set(&vector, row, i + 1.);
		}
		for (repeat = 0; repeat < 2; repeat++) {
			simplex_sparse_vector_sort(&vector);
			valid = valid && vector.count == counts[test];
			for (i = 0; i < vector.count; i++) {
				int row = vector.index[i];
				valid = valid &&
					(i == 0 || vector.index[i - 1] < row) &&
					vector.slot[row] == i &&
					vector.value[i] == vector.dense[row] &&
					simplex_sparse_vector_get(&vector,
						row) == vector.value[i];
			}
		}
		for (i = 0; i < counts[test]; i++) {
			int row = (i * 7919) % vector.dimension;
			(void)simplex_sparse_vector_add(&vector, row, 1.);
			valid = valid && simplex_sparse_vector_get(
						 &vector, row) == i + 2.;
		}
	}
	simplex_sparse_vector_destroy(&vector);
	return valid;
}

int main(void)
{
	if (!packing() || !compact_reuse() || !integer_ordering()) {
		puts("Sparse vector packing regression failed");
		return 1;
	}
	puts("Sparse vector packing handles reuse, tolerance and nonfinite values");
	return 0;
}
