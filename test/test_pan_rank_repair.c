/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file test_pan_rank_repair.c
 * @brief Regression tests for dependent support and Pan basis rank repair.
 */

#include "simplex_pan_config.h"
#include "simplex_pan_basis.h"

#include <lp_simplex/status.h>

#include <math.h>
#include <stdio.h>
#include <string.h>

static int check_repair(const int required_dependent)
{
	struct simplex_PanConfig basis_controls;
	struct simplex_PanStandard standard;
	struct simplex_PanBasis basis;
	int start[4] = {0, 1, 2, 4}, row[4] = {0, 1, 0, 1};
	double entry[4] = {1., 1., 1., 1.};
	double normalized[4] = {
		1., 1., 0.7071067811865475244, 0.7071067811865475244};
	double norm[3] = {1., 1., 1.4142135623730950488};
	double value[3] = {2., 3., 0.}, rhs[3] = {2., 3., 0.}, solution[3],
	       residual;
	int i, result, valid;

	simplex_pan_config_init(&basis_controls, NULL);
	memset(&standard, 0, sizeof(standard));
	standard.rows = 3;
	standard.columns = 3;
	standard.nonzeros = 4;
	standard.column_start = start;
	standard.row_index = row;
	standard.value = entry;
	standard.normalized_value = normalized;
	standard.column_norm = norm;
	if (simplex_pan_basis_create_configured(
		    &basis, &standard, 1e-10, &basis_controls))
		return 0;
	for (i = 0; i < 3; i++)
		if (simplex_pan_basis_add(&basis, i))
			return 0;
	if (required_dependent)
		value[2] = 1.;
	result = simplex_pan_basis_repair_rank(&basis, value);
	valid = required_dependent
			? result == lp_simplex_EXIT_FAILURE
			: result == lp_simplex_EXIT_SUCCESS &&
				  basis.count == 2 && basis.position[0] >= 0 &&
				  basis.position[1] >= 0 &&
				  basis.position[2] < 0 &&
				  !simplex_pan_basis_least_squares(
					  &basis, rhs, solution, &residual) &&
				  residual < 1e-12;
	printf("Pan rank repair protected=%d %s\n", required_dependent,
		valid ? "passed" : "failed");
	simplex_pan_basis_destroy(&basis);
	return valid;
}

static int check_completion(int scaled)
{
	struct simplex_PanConfig basis_controls;
	struct simplex_PanStandard standard;
	struct simplex_PanBasis basis;
	int start[6] = {0, 1, 2, 3, 5, 7}, row[7] = {0, 1, 2, 0, 1, 1, 2};
	double a[7] = {1., 1., 1., 1., 1., 1., 1.};
	double normalized[7] = {1., 1., 1., .7071067811865475,
		.7071067811865475, .7071067811865475, .7071067811865475};
	double norm[5] = {1., 1., 1., 1.4142135623730951, 1.4142135623730951};
	double value[5] = {0., 0., 0., 1., 2.}, rhs[3] = {1., 3., 2.},
	       solution[3], residual;
	int i, valid;

	simplex_pan_config_init(&basis_controls, NULL);
	memset(&standard, 0, sizeof(standard));
	if (scaled) {
		a[3] *= 1e-12;
		a[4] *= 1e-12;
		norm[3] *= 1e-12;
		value[3] *= 1e12;
	}
	standard.rows = 3;
	standard.columns = 5;
	standard.nonzeros = 7;
	standard.column_start = start;
	standard.row_index = row;
	standard.value = a;
	standard.normalized_value = normalized;
	standard.column_norm = norm;
	if (simplex_pan_basis_create_configured(
		    &basis, &standard, 1e-10, &basis_controls))
		return 0;
	if (simplex_pan_basis_add(&basis, 3) ||
		simplex_pan_basis_add(&basis, 4) ||
		simplex_pan_basis_rebuild(&basis))
		return 0;
	valid = !simplex_pan_basis_complete(&basis, value) &&
		basis.count == 3 && basis.position[3] >= 0 &&
		basis.position[4] >= 0 &&
		!simplex_pan_basis_least_squares(
			&basis, rhs, solution, &residual) &&
		residual < 1e-12;
	for (i = 0; i < basis.count; i++)
		valid = valid &&
			fabs(solution[i] - value[basis.column[i]]) <
				1e-12 * (1. + fabs(value[basis.column[i]]));
	printf("Pan overcomplete bootstrap scaled=%d preserves positive support %s\n",
		scaled, valid ? "passed" : "failed");
	simplex_pan_basis_destroy(&basis);
	return valid;
}

static int check_stability(void)
{
	struct simplex_PanConfig basis_controls;
	struct simplex_PanStandard standard;
	struct simplex_PanBasis basis;
	int start[3] = {0, 1, 3}, row[3] = {0, 0, 1};
	double a[3] = {1., 1., 1e-10}, norm[2] = {1., 1.}, value[2] = {1., 0.};
	int valid;

	simplex_pan_config_init(&basis_controls, NULL);
	memset(&standard, 0, sizeof(standard));
	standard.rows = 2;
	standard.columns = 2;
	standard.nonzeros = 3;
	standard.column_start = start;
	standard.row_index = row;
	standard.value = a;
	standard.normalized_value = a;
	standard.column_norm = norm;
	if (simplex_pan_basis_create_configured(
		    &basis, &standard, 1e-14, &basis_controls))
		return 0;
	if (simplex_pan_basis_add(&basis, 0) ||
		simplex_pan_basis_add(&basis, 1) ||
		simplex_pan_basis_rebuild(&basis))
		return 0;
	valid = !simplex_pan_basis_stabilize(&basis, value, -1) &&
		basis.count == 1 && basis.position[0] >= 0 &&
		basis.position[1] < 0 && basis.stability_drops == 1 &&
		basis.reciprocal_condition > .9;
	printf("Pan weak zero-column removal %s\n",
		valid ? "passed" : "failed");
	simplex_pan_basis_destroy(&basis);
	return valid;
}

int main(void)
{
	return check_repair(0) && check_repair(1) && check_completion(0) &&
			       check_completion(1) && check_stability()
		       ? 0
		       : 1;
}
