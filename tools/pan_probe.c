/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file pan_probe.c
 * @brief Isolate Pan Phase II with an explicitly supplied feasible original
 * point.
 *
 * External reference points are diagnostics, never a production dependency.
 */

#include <lp_simplex/lp_simplex.h>

#include "simplex_pan.h"
#include "simplex_problem.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>


int main(int argc, char **argv)
{
	struct lp_Model *model;
	struct simplex_Problem problem;
	struct lp_simplex_Options options;
	struct lp_simplex_Result result;
	double *start, *x, violation = 0., objective = 0.;
	FILE *input;
	clock_t begun;
	int j, i, k, state, passed;

	if (argc < 4 || argc > 6) {
		fprintf(stderr,
			"usage: pan_probe model.mps point.txt iterations [objective-or-dash [save-point]]\n");
		return 2;
	}
	model = lp_read_mps(argv[1]);
	if (model == NULL)
		return 2;
	if (simplex_problem_from_model(&problem, model)) {
		lp_model_free(model);
		return 2;
	}
	start = (double *)calloc(
		(size_t)(model->n > 0 ? model->n : 1) * 2, sizeof(double));
	if (start == NULL) {
		simplex_problem_destroy(&problem);
		lp_model_free(model);
		return 2;
	}
	x = start + model->n;
	input = fopen(argv[2], "r");
	if (input == NULL) {
		free(start);
		simplex_problem_destroy(&problem);
		lp_model_free(model);
		return 2;
	}
	for (j = 0; j < model->n; j++)
		if (fscanf(input, "%lf", start + j) != 1 ||
			!isfinite(start[j])) {
			fprintf(stderr, "invalid point at column %d\n", j);
			fclose(input);
			free(start);
			simplex_problem_destroy(&problem);
			lp_model_free(model);
			return 2;
		}
	fclose(input);
	lp_simplex_default_options(&options, lp_simplex_ALGORITHM_PAN_BDA);
	options.iteration_limit = atoi(argv[3]);
	result.status = lp_simplex_CondUnsatisfied;
	result.iterations = 0;
	result.dual_infeasibility = HUGE_VAL;
	begun = clock();
	state = simplex_pan_solve_from_feasible(
		&problem, &options, start, x, &result);
	for (j = 0; j < problem.columns; j++) {
		const struct optm_VariableBound *bound = problem.bounds + j;
		double error = 0.;

		objective += problem.objective[j] * x[j];
		if (!isfinite(x[j]))
			violation = HUGE_VAL;
		if (bound->b_type == optm_BOUND_T_LO ||
			bound->b_type == optm_BOUND_T_BS)
			error = fmax(error, bound->lb - x[j]);
		if (bound->b_type == optm_BOUND_T_UP ||
			bound->b_type == optm_BOUND_T_BS)
			error = fmax(error, x[j] - bound->ub);
		violation = fmax(violation, error);
	}
	for (i = 0; i < problem.rows; i++) {
		long double activity = 0.;
		double error;

		for (k = problem.matrix.row_start[i];
			k < problem.matrix.row_start[i + 1]; k++)
			activity += (long double)problem.matrix.row_value[k] *
				    x[problem.matrix.column_index[k]];
		error = (double)(activity - problem.rhs[i]);
		if (problem.row_type[i] == optm_CONS_T_EQ)
			error = fabs(error);
		else if (problem.row_type[i] == optm_CONS_T_GE)
			error = -error;
		violation = fmax(violation, error);
	}
	objective += model->objective_offset;
	passed = state == lp_simplex_EXIT_SUCCESS &&
		 result.status == lp_simplex_Success && isfinite(objective) &&
		 violation <= options.primal_tolerance &&
		 result.dual_infeasibility <= options.dual_tolerance;
	if (argc >= 5 && strcmp(argv[4], "-") != 0) {
		double expected = atof(argv[4]);

		passed = passed &&
			 fabs(objective - expected) <=
				 options.dual_tolerance * (1. + fabs(expected));
	}
	if (argc == 6 && isfinite(objective) &&
		violation <= options.primal_tolerance) {
		int written = 1;

		input = fopen(argv[5], "w");
		if (input == NULL)
			written = 0;
		else {
			for (j = 0; j < model->n; j++)
				if (fprintf(input, "%.17g\n", x[j]) < 0)
					written = 0;
			if (fclose(input) != 0)
				written = 0;
		}
		if (written)
			printf("saved point: primal=%.6g status=%d path=%s\n",
				violation, result.status, argv[5]);
		else {
			fprintf(stderr, "could not save point: %s\n", argv[5]);
			passed = 0;
		}
	}
	printf("pan probe: status=%d iterations=%d objective=%.17g primal=%.6g dual=%.6g seconds=%.6f verified=%d\n",
		result.status, result.iterations, objective, violation,
		result.dual_infeasibility,
		(double)(clock() - begun) / CLOCKS_PER_SEC, passed);
	free(start);
	simplex_problem_destroy(&problem);
	lp_model_free(model);
	return passed ? 0 : 1;
}
