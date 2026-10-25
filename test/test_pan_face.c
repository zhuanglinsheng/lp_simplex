/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file test_pan_face.c
 * @brief Regression tests for Pan degenerate-face control and transition rules.
 */

#include "simplex_pan_config.h"
#include "simplex_pan_cone.h"
#include "simplex_pan_basis.h"
#include "simplex_pan.h"
#include "simplex_problem.h"

#include <lp_simplex/lp_simplex.h>
#include <lp_simplex/status.h>

#include <math.h>
#include <stdio.h>
#include <string.h>

static int compression(void)
{
	struct simplex_PanConfig controls;
	struct simplex_PanStandard a;
	struct lp_simplex_Options options;
	int start[5] = {0, 1, 2, 4, 5}, row[5] = {0, 1, 0, 1, 2};
	double entry[5] = {1., 1., 1., 1., 1e12},
	       norm[4] = {1., 1., 1.4142135623730951, 1e12};
	double normalized[5] = {
		1., 1., .7071067811865475, .7071067811865475, 1.};
	double rhs[3] = {1., 1., 1.}, cost[4] = {1., 1., 0., 0.},
	       point[4] = {.5, .5, .5, 1e-12};
	double r[3], work[3];
	int state, j, support = 0;

	memset(&a, 0, sizeof(a));
	a.rows = 3;
	a.columns = 4;
	a.nonzeros = 5;
	a.column_start = start;
	a.row_index = row;
	a.value = entry;
	a.normalized_value = normalized;
	a.column_norm = norm;
	a.rhs = rhs;
	a.objective = cost;
	lp_simplex_default_options(&options, lp_simplex_ALGORITHM_PAN_BDA);
	simplex_pan_config_init(&controls, &options);
	state = simplex_pan_support_compress_configured(
		&a, point, &options, &controls);
	for (j = 0; j < 4; j++)
		if (point[j] > 0.)
			support++;
	return state == 1 && support <= 3 && point[0] >= 0. && point[1] >= 0. &&
	       point[2] >= 0. && point[0] + point[1] <= 1. + 1e-12 &&
	       fabs(1e12 * point[3] - 1.) < 1e-12 &&
	       simplex_pan_residual(&a, point, r, work) < 1e-12;
}

/* A weak positive coordinate must decrease to reach the optimum. Restricting
 * its direction initially is safe only if failed original complementarity
 * reopens it. Separate rows keep this positive support independent. */
static int reopen(void)
{
	struct simplex_PanConfig basis_controls;
	struct simplex_PanConfig controls;
	struct simplex_PanStandard a;
	struct simplex_PanBasis basis;
	struct lp_simplex_Options options;
	int start[4] = {0, 1, 2, 4}, row[4] = {0, 1, 0, 1};
	double entry[4] = {1., 1., 1., 1.},
	       norm[3] = {1., 1., 1.4142135623730951};
	double normalized[4] = {1., 1., .7071067811865475, .7071067811865475};
	double rhs[2] = {1e-9, 1.}, cost[3] = {1., 0., 0.},
	       point[3] = {1e-9, 1., 0.}, pi[2], dual;
	int iterations = 0, status = -1, state;

	simplex_pan_config_init(&basis_controls, NULL);
	memset(&a, 0, sizeof(a));
	a.rows = 2;
	a.columns = 3;
	a.nonzeros = 4;
	a.column_start = start;
	a.row_index = row;
	a.value = entry;
	a.normalized_value = normalized;
	a.column_norm = norm;
	a.rhs = rhs;
	a.objective = cost;
	lp_simplex_default_options(&options, lp_simplex_ALGORITHM_PAN_BDA);
	simplex_pan_config_init(&controls, &options);
	options.primal_tolerance = 1e-12;
	options.dual_tolerance = 1e-12;
	options.iteration_limit = 1000;
	if (simplex_pan_basis_create_configured(
		    &basis, &a, options.pivot_tolerance, &basis_controls))
		return 0;
	state = simplex_pan_cone_phase2_configured(&a, &basis, point, pi,
		&options, &iterations, &status, &dual, &controls);
	simplex_pan_basis_destroy(&basis);
	return state == 0 && status == lp_simplex_Success &&
	       fabs(point[0]) < 1e-12 && fabs(point[2] - 1e-9) < 1e-12 &&
	       fabs(point[1] + point[2] - 1.) < 1e-12 && dual < 1e-12;
}

static int minimum_norm_newton(void)
{
	struct simplex_PanStandard a;
	int start[3] = {0, 2, 4}, row[4] = {0, 1, 0, 1};
	double entry[4] = {1., 1., 1., 1.},
	       norm[2] = {1.4142135623730951, 1.4142135623730951};
	double point[2] = {1., 1.}, cost[2] = {1., 1.}, map[2] = {1., 1.},
	       y[2] = {0., 0.}, delta[2];
	unsigned char free_mask[2] = {1, 1};

	memset(&a, 0, sizeof(a));
	a.rows = 2;
	a.columns = 2;
	a.nonzeros = 4;
	a.column_start = start;
	a.row_index = row;
	a.value = entry;
	a.column_norm = norm;
	a.objective = cost;
	return simplex_pan_dual_newton_qr(
		       &a, point, free_mask, map, y, 1., delta) == 0 &&
	       fabs(delta[0] + .5) < 1e-12 && fabs(delta[1] + .5) < 1e-12;
}

static int rollback(void)
{
	struct simplex_PanConfig controls;
	struct simplex_PanStandard a;
	struct lp_simplex_Options options;
	int start[3] = {0, 1, 2}, row[2] = {0, 0};
	double entry[2] = {1., -1.}, norm[2] = {1., 1.}, rhs[1] = {0.},
	       cost[2] = {-1., -1.}, point[2] = {1., 1.};

	memset(&a, 0, sizeof(a));
	a.rows = 1;
	a.columns = 2;
	a.nonzeros = 2;
	a.column_start = start;
	a.row_index = row;
	a.value = entry;
	a.normalized_value = entry;
	a.column_norm = norm;
	a.rhs = rhs;
	a.objective = cost;
	lp_simplex_default_options(&options, lp_simplex_ALGORITHM_PAN_BDA);
	simplex_pan_config_init(&controls, &options);
	/* Both finite support-reducing endpoints increase cost. Compression must
	 * restore the point and leave recession-ray handling to the LP solver. */
	return simplex_pan_support_compress_configured(
		       &a, point, &options, &controls) == -1 &&
	       point[0] == 1. && point[1] == 1.;
}

static int dependent_warm_point(void)
{
	struct lp_Model *model = lp_model_create(2, 3);
	struct simplex_Problem problem;
	struct lp_simplex_Options options;
	struct lp_simplex_Result result;
	double point[3] = {.5, .5, .5}, x[3];
	int state, good, i;

	if (!model)
		return 0;
	for (i = 0; i < 2; i++) {
		model->constraints[i].coef[i] = 1.;
		model->constraints[i].coef[2] = 1.;
		model->constraints[i].rhs = 1.;
		model->constraints[i].type = optm_CONS_T_EQ;
		model->objective[i] = 1.;
	}
	if (simplex_problem_from_model(&problem, model)) {
		lp_model_free(model);
		return 0;
	}
	lp_simplex_default_options(&options, lp_simplex_ALGORITHM_PAN_BDA);
	state = simplex_pan_solve_from_feasible(
		&problem, &options, point, x, &result);
	good = state == 0 && result.status == lp_simplex_Success &&
	       fabs(result.objective) < 1e-10 && fabs(x[2] - 1.) < 1e-10 &&
	       fabs(x[0]) < 1e-10 && fabs(x[1]) < 1e-10;
	simplex_problem_destroy(&problem);
	lp_model_free(model);
	return good;
}

/* An equivalent large dual shift can round away the original cost difference.
 * Both orientations must still compress toward the true lower-cost endpoint. */
static int original_objective_audit(void)
{
	struct simplex_PanConfig controls;
	struct simplex_PanStandard a;
	struct lp_simplex_Options options;
	int start[3] = {0, 1, 2}, row[2] = {0, 0}, orientation;
	double entry[2] = {1., 1.}, norm[2] = {1., 1.}, rhs[1] = {1.};
	double shifted[2] = {1e18, 1e18}, original[2], point[2];

	memset(&a, 0, sizeof(a));
	a.rows = 1;
	a.columns = 2;
	a.nonzeros = 2;
	a.column_start = start;
	a.row_index = row;
	a.value = entry;
	a.normalized_value = entry;
	a.column_norm = norm;
	a.rhs = rhs;
	a.objective = shifted;
	a.objective_offset = -1e18;
	a.audit_objective = original;
	lp_simplex_default_options(&options, lp_simplex_ALGORITHM_PAN_BDA);
	simplex_pan_config_init(&controls, &options);
	for (orientation = 0; orientation < 2; orientation++) {
		original[0] = orientation ? 0. : 1.;
		original[1] = orientation ? 1. : 0.;
		point[0] = point[1] = .5;
		if (simplex_pan_support_compress_configured(
			    &a, point, &options, &controls) != 1 ||
			fabs(point[0] + point[1] - 1.) > 1e-12 ||
			point[0] < 0. || point[1] < 0. ||
			original[0] * point[0] + original[1] * point[1] > 1e-12)
			return 0;
	}
	return 1;
}

static int select_required(void)
{
	struct simplex_PanConfig basis_controls;
	struct simplex_PanStandard a;
	struct simplex_PanBasis b;
	int starts[5] = {0, 1, 2, 4, 5}, rows[5] = {0, 1, 0, 1, 0},
	    columns[4] = {0, 1, 2, 3};
	double values[5] = {1., 1., 1., 1., 1e-12},
	       norm[4] = {1., 1., 1.4142135623730951, 1e-12};
	double normalized[5] = {
		1., 1., .7071067811865475, .7071067811865475, 1.};
	double required[4] = {0., 0., 0., 1.};
	int good;

	simplex_pan_config_init(&basis_controls, NULL);
	memset(&a, 0, sizeof(a));
	a.rows = 2;
	a.columns = 4;
	a.nonzeros = 5;
	a.column_start = starts;
	a.row_index = rows;
	a.value = values;
	a.column_norm = norm;
	a.normalized_value = normalized;
	if (simplex_pan_basis_create_configured(&b, &a, 1e-10, &basis_controls))
		return 0;
	good = simplex_pan_basis_select_columns(&b, columns, 4, required) ==
		       0 &&
	       b.count == 2 && b.position[3] >= 0;
	simplex_pan_basis_destroy(&b);
	return good;
}

int main(void)
{
	int good = select_required() && original_objective_audit() &&
		   minimum_norm_newton() && rollback() && compression() &&
		   reopen() && dependent_warm_point();

	printf("Pan feasible support compression and face reopening %s\n",
		good ? "passed" : "failed");
	return good ? 0 : 1;
}
