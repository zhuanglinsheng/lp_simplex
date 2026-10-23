/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file test_pan_phase1_accuracy.c
 * @brief Regression tests for Pan Phase I feasibility accuracy.
 */

#include "simplex_pan_phase1.h"

#include <lp_simplex/status.h>

#include <math.h>
#include <stdio.h>
#include <string.h>


static int small_positive(void)
{
	struct simplex_PanStandard standard;
	struct simplex_PanBasis basis;
	int start[2] = {0, 2}, row[2] = {0, 1};
	double a[2] = {1e8, 1e8},
	       normalized[2] = {0.7071067811865475244, 0.7071067811865475244};
	double norm[1] = {1.4142135623730950488e8}, rhs[2] = {1., 1.}, x[1];
	int steps = 0, status = -1, state, valid;

	memset(&standard, 0, sizeof(standard));
	standard.rows = 2;
	standard.columns = 1;
	standard.nonzeros = 2;
	standard.column_start = start;
	standard.row_index = row;
	standard.value = a;
	standard.normalized_value = normalized;
	standard.column_norm = norm;
	standard.rhs = rhs;
	if (simplex_pan_basis_create(&basis, &standard, 1e-10))
		return 0;
	state = simplex_pan_phase1_run(
		&basis, x, 1e-7, 1e-7, 10, &steps, &status);
	valid = state == lp_simplex_EXIT_SUCCESS && basis.count == 1 &&
		steps == 1 && fabs(x[0] - 1e-8) < 1e-20 &&
		fabs(a[0] * x[0] - rhs[0]) <= 1e-7;
	simplex_pan_basis_destroy(&basis);
	printf("NNLS retains small positive values %s\n",
		valid ? "passed" : "failed");
	return valid;
}


static int infeasible_witness(void)
{
	struct simplex_PanStandard standard;
	struct simplex_PanBasis basis;
	int start[2] = {0, 1}, row[1] = {0}, steps = 0, status = -1, state;
	double a[1] = {-1.}, norm[1] = {1.}, rhs[1] = {1e-4}, x[1];

	memset(&standard, 0, sizeof(standard));
	standard.rows = 1;
	standard.columns = 1;
	standard.nonzeros = 1;
	standard.column_start = start;
	standard.row_index = row;
	standard.value = a;
	standard.normalized_value = a;
	standard.column_norm = norm;
	standard.rhs = rhs;
	if (simplex_pan_basis_create(&basis, &standard, 1e-10))
		return 0;
	state = simplex_pan_phase1_run(
		&basis, x, 1e-7, 1e-7, 10, &steps, &status);
	simplex_pan_basis_destroy(&basis);
	return state == lp_simplex_EXIT_FAILURE &&
	       status == lp_simplex_Infeasibility;
}


int main(void)
{
	return small_positive() && infeasible_witness() ? 0 : 1;
}
