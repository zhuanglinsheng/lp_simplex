/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file test_pan_cone_oracle.c
 * @brief Regression tests for tangent-cone directions on degenerate models.
 */

#include "simplex_pan_config.h"
#include "simplex_pan_cone.h"
#include "simplex_pan_basis.h"

#include <lp_simplex/status.h>

#include <math.h>
#include <stdio.h>
#include <string.h>

/* Analytic tangent-cone problems: a feasible decreasing edge and a stationary
 * constrained fit. The latter has an improving cost column but no feasible
 * decreasing direction, so it must not be returned as an exact fit. */
static int run(int stationary)
{
	struct simplex_PanConfig basis_controls;
	struct simplex_PanStandard a;
	struct simplex_PanBasis b;
	int start[4] = {0, 1, 2, 4}, rows[4] = {0, 1, 0, 1};
	double values[4] = {1., 1., 1., 1.};
	double norms[3] = {1., 1., 1.4142135623730951};
	double normalized[4] = {
		1., 1., .7071067811865475244, .7071067811865475244};
	double cost[3] = {0., 0., -1.}, rhs[2] = {0., 0.}, d[3], r[3];
	unsigned char free_mask[3] = {1, 0, 0};
	int steps = 0, status = lp_simplex_CondUnsatisfied, state, valid;

	simplex_pan_config_init(&basis_controls, NULL);
	memset(&a, 0, sizeof(a));
	a.rows = stationary ? 2 : 1;
	a.columns = stationary ? 3 : 2;
	a.nonzeros = stationary ? 4 : 2;
	if (!stationary) {
		rows[1] = 0;
		cost[0] = 1.;
		cost[1] = 0.;
	}
	a.column_start = start;
	a.row_index = rows;
	a.value = values;
	a.column_norm = norms;
	a.normalized_value = normalized;
	a.objective = cost;
	a.rhs = rhs;
	if (simplex_pan_basis_create_configured(&b, &a, 1e-10, &basis_controls))
		return 0;
	if (simplex_pan_basis_add(&b, 0) || simplex_pan_basis_rebuild(&b)) {
		simplex_pan_basis_destroy(&b);
		return 0;
	}
	state = simplex_pan_cone_implicit_audited(
		&b, d, free_mask, r, 1e-10, 100, &steps, &status, NULL, NULL);
	valid = stationary
			? state == 1 && fabs(d[0] + .5) < 1e-9 &&
				  fabs(d[1]) < 1e-9 && fabs(d[2] - .5) < 1e-9 &&
				  fabs(r[0]) < 1e-9 && fabs(r[1] + .5) < 1e-9 &&
				  fabs(r[2] + .5) < 1e-9
			: state == 0 && fabs(d[0] + 1.) < 1e-9 &&
				  fabs(d[1] - 1.) < 1e-9 && fabs(r[0]) < 1e-9 &&
				  fabs(r[1]) < 1e-9;
	simplex_pan_basis_destroy(&b);
	printf("implicit cone %s %s\n",
		stationary ? "stationary witness" : "decreasing edge",
		valid ? "passed" : "failed");
	return valid;
}

int main(void)
{
	return run(0) && run(1) ? 0 : 1;
}
