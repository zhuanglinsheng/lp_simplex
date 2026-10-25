/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file test_pan_phase1.c
 * @brief Regression tests for Pan Phase I feasibility recovery.
 */

#include "simplex_pan_config.h"
#include "simplex_pan_phase1.h"

#include <lp_simplex/lp_simplex.h>
#include <lp_simplex/status.h>

#include <math.h>
#include <stdio.h>
#include <string.h>

/* The search correlation can overflow/underflow even when all input data,
 * norms, and a feasible witness are finite. Such a stop is numerical, not
 * mathematical infeasibility. Use powers of two for exact witness checks. */
static int test_classification(void)
{
	struct simplex_PanConfig basis_controls;
	int variant;
	simplex_pan_config_init(&basis_controls, NULL);
	for (variant = 0; variant < 6; variant++) {
		struct simplex_PanStandard standard;
		struct simplex_PanBasis basis;
		int starts[2] = {0, 2}, rows[2] = {0, 1};
		int work = 0, status = -1, code, expected;
		double entries[2], normalized[2], norm[1], rhs[2], value[1];
		double tolerance = 1e-8, witness = 1.;
		unsigned char free_mask[1] = {1};

		entries[0] = entries[1] = 1.;
		rhs[0] = rhs[1] = 1.;
		expected = lp_simplex_Success;
		if (variant == 1) {
			entries[0] = entries[1] = ldexp(1., 1022);
			rhs[0] = rhs[1] = 8.;
			witness = ldexp(1., -1019);
			expected = lp_simplex_PrecisionError;
		} else if (variant == 2) {
			entries[0] = entries[1] = ldexp(1., -1000);
			rhs[0] = rhs[1] = ldexp(1., -100);
			tolerance = ldexp(1., -105);
			witness = ldexp(1., 900);
			expected = lp_simplex_PrecisionError;
		} else if (variant == 3) {
			entries[0] = entries[1] = -1.;
			expected = lp_simplex_Infeasibility;
		} else if (variant >= 4) {
			entries[1] = -1.;
			expected = variant == 4 ? lp_simplex_Infeasibility
						: lp_simplex_CondUnsatisfied;
		}
		if (variant <= 2 && (entries[0] * witness != rhs[0] ||
					    entries[1] * witness != rhs[1]))
			return 0;
		norm[0] = hypot(entries[0], entries[1]);
		normalized[0] = entries[0] / norm[0];
		normalized[1] = entries[1] / norm[0];
		memset(&standard, 0, sizeof(standard));
		standard.rows = 2;
		standard.columns = 1;
		standard.nonzeros = 2;
		standard.column_start = starts;
		standard.row_index = rows;
		standard.value = entries;
		standard.normalized_value = normalized;
		standard.column_norm = norm;
		standard.rhs = rhs;
		if (simplex_pan_basis_create_configured(
			    &basis, &standard, 1e-9, &basis_controls))
			return 0;
		if (variant == 5)
			code = simplex_pan_phase1_mixed(&basis, value,
				free_mask, NULL, tolerance, tolerance, 100,
				&work, &status);
		else
			code = simplex_pan_phase1_run(&basis, value, tolerance,
				tolerance, 100, &work, &status);
		simplex_pan_basis_destroy(&basis);
		/* Feasible NNLS reports a handoff, not LP optimality. */
		if ((variant == 0 &&
			    (code != 0 ||
				    status != lp_simplex_CondUnsatisfied)) ||
			(variant > 0 && status != expected) ||
			(variant > 0 && variant < 5 && code != -1) ||
			(variant == 5 && code != 1) || work < 0 || work > 100) {
			fprintf(stderr,
				"Phase I classification %d: code=%d status=%d work=%d\n",
				variant, code, status, work);
			return 0;
		}
	}
	return 1;
}

static int test_public_classification(void)
{
	struct lp_Model *model = lp_model_create(2, 1);
	struct lp_simplex_Options options;
	struct lp_simplex_Result result;
	double point[1], a = ldexp(1., 1022), witness = ldexp(1., -1019);
	int i, state, valid;

	if (!model || a * witness != 8.) {
		lp_model_free(model);
		return 0;
	}
	for (i = 0; i < 2; i++) {
		model->constraints[i].type = optm_CONS_T_EQ;
		model->constraints[i].rhs = 8.;
		model->constraints[i].coef[0] = a;
	}
	lp_simplex_default_options(&options, lp_simplex_ALGORITHM_PAN_BDA);
	options.presolve = 0;
	state = lp_simplex_solve(model, &options, point, &result);
	valid = state == lp_simplex_EXIT_FAILURE &&
		result.status == lp_simplex_PrecisionError &&
		result.iterations >= 0 &&
		result.iterations <= options.iteration_limit;
	lp_model_free(model);
	return valid;
}

/* Orthogonal singleton minimizers are an exact NNLS warm start.  Opposite
 * signs and duplicate columns must not seed a negative or dependent value. */
int main(void)
{
	struct simplex_PanConfig basis_controls;
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

	simplex_pan_config_init(&basis_controls, NULL);
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
	if (simplex_pan_basis_create_configured(&basis, &standard, 1e-10,
		    &basis_controls) == lp_simplex_EXIT_FAILURE)
		return 1;
	result = simplex_pan_phase1_run(
		&basis, value, 1e-8, 1e-8, 10, &iterations, &status);
	valid = result == lp_simplex_EXIT_SUCCESS && basis.count == 2 &&
		iterations == 0 && fabs(value[0] - 2.) < 1e-12 &&
		fabs(value[2] - 2.) < 1e-12 && value[1] == 0. && value[3] == 0.;
	simplex_pan_basis_destroy(&basis);
	printf("Pan singleton crash %s\n", valid ? "passed" : "failed");
	return valid && test_classification() && test_public_classification()
		       ? 0
		       : 1;
}
