/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file test_pan_mixed_nnls.c
 * @brief Regression tests for mixed free and nonnegative tangent-cone
 * coordinates.
 */

#include "simplex_pan_phase1.h"

#include <lp_simplex/status.h>

#include <math.h>
#include <stdio.h>
#include <string.h>


static int check(int feasible)
{
	struct simplex_PanStandard standard;
	struct simplex_PanBasis basis;
	int start[3] = {0, 1, 3}, row[3] = {0, 0, 1};
	double a[3] = {1., 1., 1.},
	       normalized[3] = {1., .7071067811865475, .7071067811865475};
	double norm[2] = {1., 1.4142135623730951}, rhs[2] = {0., 1.}, x[2],
	       residual[2];
	unsigned char free_mask[2] = {1, 0};
	int steps = 0, status = -1, state, valid;

	memset(&standard, 0, sizeof(standard));
	standard.rows = 2;
	standard.columns = 2;
	standard.nonzeros = 3;
	standard.column_start = start;
	standard.row_index = row;
	standard.value = a;
	standard.normalized_value = normalized;
	standard.column_norm = norm;
	standard.rhs = rhs;
	if (!feasible)
		rhs[1] = -1.;
	if (simplex_pan_basis_create(&basis, &standard, 1e-12) ||
		simplex_pan_basis_add(&basis, 0) ||
		simplex_pan_basis_rebuild(&basis))
		return 0;
	state = simplex_pan_phase1_mixed(&basis, x, free_mask, residual, 1e-10,
		1e-10, 100, &steps, &status);
	valid = feasible ? state == 0 && fabs(x[0] + 1.) < 1e-12 &&
				   fabs(x[1] - 1.) < 1e-12
			 : state == 1 && fabs(residual[0]) < 1e-12 &&
				   fabs(residual[1] + 1.) < 1e-12;
	printf("mixed NNLS fit=%d preserves signed free coordinates %s\n",
		feasible, valid ? "passed" : "failed");
	simplex_pan_basis_destroy(&basis);
	return valid;
}


int main(void)
{
	return check(1) && check(0) ? 0 : 1;
}
