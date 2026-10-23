/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file test_pan_trial.c
 * @brief Regression tests for finite tangent-cone trial-step acceptance.
 */

#include "simplex_pan_cone.h"
#include "simplex_pan_phase1.h"

#include <lp_simplex/lp_simplex.h>

#include <math.h>
#include <stdio.h>
#include <string.h>


static int actual_step(void)
{
	struct simplex_PanStandard a;
	struct simplex_PanDirectionTrial t;
	struct lp_simplex_Options o;
	int starts[3] = {0, 1, 2}, row[2] = {0, 0}, used = 0;
	double values[2] = {1., 1.}, norm[2] = {1., 1.}, rhs[1] = {1.};
	double cost[2] = {1., 0.}, shifted[2] = {1e18, 1e18},
	       point[2] = {1e-7, 1. - 1e-7};
	double trial[2], residual[1], work[1], direction[2] = {-1., .9};

	memset(&a, 0, sizeof(a));
	memset(&t, 0, sizeof(t));
	a.rows = 1;
	a.columns = 2;
	a.nonzeros = 2;
	a.column_start = starts;
	a.row_index = row;
	a.value = values;
	a.normalized_value = values;
	a.column_norm = norm;
	a.rhs = rhs;
	a.objective = shifted;
	a.audit_objective = cost;
	a.objective_offset = -1e18;
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_PAN_BDA);
	o.primal_tolerance = 1e-7;
	o.iteration_limit = 10;
	t.a = &a;
	t.options = &o;
	t.point = point;
	t.trial = trial;
	t.residual = residual;
	t.work = work;
	t.iterations = &used;
	t.calls = 31;
	if (!simplex_pan_cone_trial_step(&t, direction) || used != 1 ||
		t.accepted != 1 || trial[0] != 0. ||
		fabs(trial[1] - (1. - 1e-8)) > 1e-15 || point[0] != 1e-7 ||
		point[1] != 1. - 1e-7)
		return 0;
	/* Same objective benefit but too much original equality error. */
	direction[1] = -1.;
	t.calls = 31;
	if (simplex_pan_cone_trial_step(&t, direction) || used != 2)
		return 0;
	/* A cost increase hidden by a large centered objective is refused. */
	direction[0] = 1.;
	direction[1] = -1.;
	t.calls = 31;
	if (simplex_pan_cone_trial_step(&t, direction) || used != 3)
		return 0;
	/* Do not audit rays using a finite-step relaxation. */
	cost[0] = cost[1] = -1.;
	direction[0] = direction[1] = 1.;
	t.calls = 31;
	if (simplex_pan_cone_trial_step(&t, direction))
		return 0;
	used = o.iteration_limit;
	t.calls = 31;
	return !simplex_pan_cone_trial_step(&t, direction) &&
	       used == o.iteration_limit;
}


static int exhaust(void *context, const double *direction)
{
	int *used = (int *)context;

	(void)direction;
	*used = 1;
	return 0;
}


static int callback_budget(int implicit)
{
	struct simplex_PanStandard a;
	struct simplex_PanBasis b;
	int starts[3] = {0, 1, 2}, rows[3] = {0, 0, 0}, used = 0, status = -1,
	    state;
	double values[3] = {1., 1., 1.}, normalized[3] = {1., 1., 1.},
	       norm[2] = {1., 1.};
	double rhs[2] = {0., -1.}, cost[2] = {1., 0.}, d[2], r[2];
	unsigned char free_mask[2] = {1, 0};

	memset(&a, 0, sizeof(a));
	a.rows = 1;
	a.columns = 2;
	a.nonzeros = 2;
	if (!implicit) {
		a.rows = 2;
		a.nonzeros = 3;
		starts[1] = 2;
		starts[2] = 3;
		rows[1] = 1;
		norm[0] = sqrt(2.);
		normalized[0] = normalized[1] = 1. / sqrt(2.);
	}
	a.column_start = starts;
	a.row_index = rows;
	a.value = values;
	a.normalized_value = normalized;
	a.column_norm = norm;
	a.rhs = rhs;
	a.objective = cost;
	if (simplex_pan_basis_create(&b, &a, 1e-10) ||
		simplex_pan_basis_add(&b, 0) || simplex_pan_basis_rebuild(&b))
		return 0;
	state = implicit ? simplex_pan_cone_implicit_audited(&b, d, free_mask,
				   r, 1e-10, 1, &used, &status, exhaust, &used)
			 : simplex_pan_phase1_mixed_audited(&b, d, free_mask, r,
				   1e-10, 1e-10, 1, &used, &status, exhaust,
				   &used);
	simplex_pan_basis_destroy(&b);
	return state < 0 && status == lp_simplex_ExceedIterLimit && used == 1;
}


int main(void)
{
	int good = actual_step() && callback_budget(0) && callback_budget(1);

	printf("Pan original-model NNLS trial and callback budget %s\n",
		good ? "passed" : "failed");
	return good ? 0 : 1;
}
