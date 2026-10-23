/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file test_pan_phase1.c
 * @brief Regression tests for Pan Phase I feasibility recovery.
 */

#include "simplex_pan_phase1.h"

#include <lp_simplex/status.h>

#include <math.h>
#include <stdio.h>
#include <string.h>


/* Orthogonal singleton minimizers are an exact NNLS warm start.  Opposite
 * signs and duplicate columns must not seed a negative or dependent value. */
int main(void)
{
	struct simplex_PanStandard standard;
	struct simplex_PanBasis basis;
	int starts[5] = {0, 1, 2, 3, 4};
	int rows[4] = {0, 0, 1, 2};
	double entries[4] = {2., 1., -3., 1.};
	double normalized[4] = {1., 1., -1., 1.};
	double norms[4] = {2., 1., 3., 1.};
	double rhs[3] = {4., -6., 0.};
	double value[4];
	int iterations = 0, status = -1, result, valid;

	memset(&standard, 0, sizeof(standard));
	standard.rows = 3;
	standard.columns = 4;
	standard.nonzeros = 4;
	standard.column_start = starts;
	standard.row_index = rows;
	standard.value = entries;
	standard.normalized_value = normalized;
	standard.column_norm = norms;
	standard.rhs = rhs;
	if (simplex_pan_basis_create(&basis, &standard, 1e-10) ==
		lp_simplex_EXIT_FAILURE)
		return 1;
	result = simplex_pan_phase1_run(
		&basis, value, 1e-8, 1e-8, 10, &iterations, &status);
	valid = result == lp_simplex_EXIT_SUCCESS && basis.count == 2 &&
		iterations == 0 && fabs(value[0] - 2.) < 1e-12 &&
		fabs(value[2] - 2.) < 1e-12 && value[1] == 0. && value[3] == 0.;
	simplex_pan_basis_destroy(&basis);
	printf("Pan singleton crash %s\n", valid ? "passed" : "failed");
	return valid ? 0 : 1;
}
