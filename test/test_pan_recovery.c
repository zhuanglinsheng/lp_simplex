/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file test_pan_recovery.c
 * @brief Regression tests for Pan feasible-point recovery and acceptance
 * checks.
 */

#include "simplex_pan_recovery.h"
#include "simplex_pan_cone.h"
#include "simplex_pan_basis.h"

#include <lp_simplex/lp_simplex.h>

#include <math.h>
#include <stdio.h>
#include <string.h>


static int restore_cases(void)
{
	struct simplex_PanStandard a;
	struct lp_simplex_Options o;
	int starts[3] = {0, 1, 2}, rows[2] = {0, 1}, used, state;
	double values[2] = {1., 1e12}, norm[2] = {1., 1e12}, unit[2] = {1., 1.};
	double rhs[2] = {1., 1.}, cost[2] = {1., 0.}, point[2], saved[2];

	memset(&a, 0, sizeof(a));
	a.rows = 2;
	a.columns = 2;
	a.nonzeros = 2;
	a.column_start = starts;
	a.row_index = rows;
	a.value = values;
	a.normalized_value = unit;
	a.column_norm = norm;
	a.rhs = rhs;
	a.objective = cost;
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_PAN_BDA);
	o.primal_tolerance = 1e-7;
	o.iteration_limit = 100;
	point[0] = 1. - 9e-8;
	point[1] = (1. + 8e-8) * 1e-12;
	used = 0;
	state = simplex_pan_restore(&a, point, &o, 1.L, &used);
	if (state != 1 || fabs(point[0] - 1.) > 1e-14 ||
		fabs(point[1] * 1e12 - 1.) > 1e-14 || used <= 0 || used > 100)
		return 0;
	/* Repair must use original costs even when working costs are shifted. */
	a.audit_objective = cost;
	{
		double shifted[2] = {-1e18, -1e18};

		a.objective = shifted;
		point[0] = 1. - 9e-8;
		point[1] = 1e-12;
		memcpy(saved, point, sizeof(point));
		used = 0;
		if (simplex_pan_restore(&a, point, &o, .5L, &used) != -1 ||
			memcmp(point, saved, sizeof(point)))
			return 0;
	}
	a.objective = cost;
	a.audit_objective = NULL;
	/* No remaining budget and an invalid negative coordinate are failures
	 * with exact rollback; neither is silently clipped. */
	point[0] = 1. - 9e-8;
	point[1] = 1e-12;
	memcpy(saved, point, sizeof(point));
	used = 100;
	if (simplex_pan_restore(&a, point, &o, 1.L, &used) != -1 ||
		used != 100 || memcmp(point, saved, sizeof(point)))
		return 0;
	point[1] = -1e-20;
	memcpy(saved, point, sizeof(point));
	used = 0;
	if (simplex_pan_restore(&a, point, &o, 1.L, &used) != -1 ||
		memcmp(point, saved, sizeof(point)))
		return 0;
	/* An empty positive face cannot repair an inconsistent equality. */
	point[0] = point[1] = 0.;
	used = 0;
	if (simplex_pan_restore(&a, point, &o, 1.L, &used) != -1 ||
		point[0] != 0. || point[1] != 0.)
		return 0;
	/* A correction crossing the nonnegative boundary is shortened. It must
	 * still reduce the full original residual enough, or roll back. */
	rhs[0] = -1e-6;
	point[0] = 1e-8;
	point[1] = 1e-12;
	memcpy(saved, point, sizeof(point));
	used = 0;
	return simplex_pan_restore(&a, point, &o, 1.L, &used) == -1 &&
	       !memcmp(point, saved, sizeof(point));
}


static int expanding_face(void)
{
	struct simplex_PanStandard a;
	struct lp_simplex_Options o;
	int starts[71], rows[70], j, used = 0;
	double values[70], norm[70], rhs[70], cost[70], point[70];

	memset(&a, 0, sizeof(a));
	for (j = 0; j < 70; j++) {
		starts[j] = rows[j] = j;
		values[j] = norm[j] = rhs[j] = point[j] = 1.;
		cost[j] = 0.;
	}
	starts[70] = 70;
	point[69] -= 9e-8;
	a.rows = a.columns = a.nonzeros = 70;
	a.column_start = starts;
	a.row_index = rows;
	a.value = a.normalized_value = values;
	a.column_norm = norm;
	a.rhs = rhs;
	a.objective = cost;
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_PAN_BDA);
	o.primal_tolerance = 1e-7;
	o.iteration_limit = 200;
	return simplex_pan_restore(&a, point, &o, 0.L, &used) == 1 &&
	       used > 65 && used <= 200 && fabs(point[69] - 1.) < 1e-14;
}


/* Exercise the real trial -> residual repair -> commit -> dual audit path.
 * The tangent direction preserves an initial residual near the tolerance;
 * repair must regain headroom without increasing the original objective. */
static int repair_in_solver(void)
{
	struct simplex_PanStandard a;
	struct simplex_PanBasis b;
	struct lp_simplex_Options o;
	int starts[4] = {0, 1, 2, 4}, rows[4] = {0, 1, 0, 1};
	double values[4] = {1., 1., 1., 1.},
	       normalized[4] = {1., 1., .7071067811865475, .7071067811865475};
	double norms[3] = {1., 1., 1.4142135623730951}, rhs[2] = {1., 1.},
	       cost[3] = {1., 0., 0.};
	double point[3] = {1. - 9e-8, 1., 0.}, pi[2], dual, residual[2],
	       work[2];
	int used = 0, status = -1, result, good;

	memset(&a, 0, sizeof(a));
	a.rows = 2;
	a.columns = 3;
	a.nonzeros = 4;
	a.column_start = starts;
	a.row_index = rows;
	a.value = values;
	a.normalized_value = normalized;
	a.column_norm = norms;
	a.rhs = rhs;
	a.objective = cost;
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_PAN_BDA);
	o.primal_tolerance = 1e-7;
	o.dual_tolerance = 1e-9;
	o.iteration_limit = 1000;
	if (simplex_pan_basis_create(&b, &a, o.pivot_tolerance))
		return 0;
	result = simplex_pan_cone_phase2(
		&a, &b, point, pi, &o, &used, &status, &dual);
	simplex_pan_basis_destroy(&b);
	good = result == 0 && status == lp_simplex_Success && used <= 1000 &&
	       point[0] >= 0. && point[1] >= 0. && point[2] >= 0. &&
	       point[0] < 1e-9 &&
	       simplex_pan_residual(&a, point, residual, work) < 1e-12;
	/* Exhaustion during a trial/repair cannot commit an extra LP step. */
	point[0] = 1. - 9e-8;
	point[1] = 1.;
	point[2] = 0.;
	used = 0;
	o.iteration_limit = 2;
	if (simplex_pan_basis_create(&b, &a, o.pivot_tolerance))
		return 0;
	result = simplex_pan_cone_phase2(
		&a, &b, point, pi, &o, &used, &status, &dual);
	simplex_pan_basis_destroy(&b);
	return good && result != 0 && status == lp_simplex_ExceedIterLimit &&
	       used <= 2 && point[0] == 1. - 9e-8 && point[1] == 1. &&
	       point[2] == 0.;
}


int main(void)
{
	struct simplex_PanRecovery s;
	double point[2] = {1., 2.};
	int good;

	memset(&s, 0, sizeof(s));
	s.low_progress = 2;
	good = simplex_pan_recovery_compress(&s, point, 2) == 1 &&
	       simplex_pan_recovery_compress(&s, point, 2) == 0 &&
	       s.low_progress == 2;
	point[1] = 3.;
	good = good && simplex_pan_recovery_compress(&s, point, 2) == 1 &&
	       s.compressions == 2 && s.compression_skips == 1;
	good = good && restore_cases() && expanding_face() &&
	       repair_in_solver();
	printf("Pan bounded recovery and original-model repair %s\n",
		good ? "passed" : "failed");
	return good ? 0 : 1;
}
