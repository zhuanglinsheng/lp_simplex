/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file test_regressions.c
 * @brief Regressions for input handling and for the success contract.
 *
 * Two invariants are pinned here:
 *
 *   1. `lp_simplex_solve` may only return success when the vector it reports
 *      is finite and feasible for the caller's model.  Before this test the
 *      tableau path could report success with a residual of 5e5 and the
 *      presolve path could report success with an infinite objective.
 *   2. The MPS reader accepts the dialect features that occur in real files
 *      (comments, blank lines, OBJSENSE MAX, integer markers) and refuses the
 *      ones it cannot honour -- loudly, and without reading past the data.
 */

#include <lp_simplex/lp_simplex.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef SIMPLE_EXAMPLES_DIR
#define SIMPLE_EXAMPLES_DIR "data/simple_examples"
#endif

#define TEMP_FILE "lp_simplex_regression.mps"


static int failures = 0;


static void fail(const char *label, const char *detail)
{
	printf("[FAIL] %s: %s\n", label, detail);
	failures++;
}


static void pass(const char *label)
{
	printf("[ OK ] %s\n", label);
}


/* Check the returned vector independently of the solver's residual field. */
static double primal_violation(const struct lp_Model *model, const double *x)
{
	double maximum = 0.;
	int i, j, k;

	for (j = 0; j < model->n; j++) {
		const struct optm_VariableBound *bound = model->bounds + j;

		if (!isfinite(x[j]))
			return HUGE_VAL;
		if ((bound->b_type == optm_BOUND_T_LO ||
			    bound->b_type == optm_BOUND_T_BS) &&
			bound->lb - x[j] > maximum)
			maximum = bound->lb - x[j];
		if ((bound->b_type == optm_BOUND_T_UP ||
			    bound->b_type == optm_BOUND_T_BS) &&
			x[j] - bound->ub > maximum)
			maximum = x[j] - bound->ub;
	}
	for (i = 0; i < model->m; i++) {
		double activity = 0., violation;

		for (k = model->row_start[i]; k < model->row_start[i + 1]; k++)
			activity +=
				model->row_value[k] * x[model->column_index[k]];
		if (!isfinite(activity))
			return HUGE_VAL;
		violation = activity - model->constraints[i].rhs;
		if (model->constraints[i].type == optm_CONS_T_EQ)
			violation = fabs(violation);
		else if (model->constraints[i].type == optm_CONS_T_GE)
			violation = -violation;
		if (violation > maximum)
			maximum = violation;
	}
	return maximum;
}


/**
 * Write `lines` (NULL terminated) to a scratch file and return its path.
 *
 * Fixtures are given line by line rather than as one literal because C90 only
 * requires compilers to support 509-character string literals, and these
 * models are longer than that.
 */
static const char *write_fixture(const char *const *lines)
{
	FILE *stream = fopen(TEMP_FILE, "w");
	int i;

	if (stream == NULL)
		return NULL;
	for (i = 0; lines[i] != NULL; i++)
		fputs(lines[i], stream);
	fclose(stream);
	return TEMP_FILE;
}


/**
 * Solve `file` and report the status.  Returns 0 when the call failed and 1
 * when it succeeded, so callers can assert on both `state` and `status`.
 */
static int solve_file(const char *file, int algorithm, double *x,
	struct lp_simplex_Result *result)
{
	struct lp_simplex_Options options;
	struct lp_Model *model;
	int state;

	model = lp_read_mps(file);
	if (model == NULL)
		return -1;
	lp_simplex_default_options(&options, algorithm);
	state = lp_simplex_solve(model, &options, x, result);
	lp_model_free(model);
	return state == lp_simplex_EXIT_SUCCESS;
}


/* Success must imply a finite objective and a feasible reported vector. */
static void check_success_is_honest(
	const char *file, int algorithm, const char *label)
{
	struct lp_simplex_Options options;
	struct lp_Model *model;
	struct lp_simplex_Result result;
	double *x;
	int state;

	model = lp_read_mps(file);
	if (model == NULL) {
		fail(label, "could not read model");
		return;
	}
	lp_simplex_default_options(&options, algorithm);
	x = (double *)malloc((size_t)model->n * sizeof(double));
	if (x == NULL) {
		lp_model_free(model);
		fail(label, "allocation failed");
		return;
	}
	memset(&result, 0, sizeof(result));
	state = lp_simplex_solve(model, &options, x, &result);
	if (state == lp_simplex_EXIT_SUCCESS) {
		if (result.status != lp_simplex_Success)
			fail(label, "success returned an inconsistent status");
		else if (!isfinite(result.objective))
			fail(label, "success returned a non-finite objective");
		else if (primal_violation(model, x) >
			 10. * options.primal_tolerance)
			fail(label, "success returned an infeasible x");
		else
			pass(label);
	} else {
		fail(label, "feasible example was not solved");
	}
	free(x);
	lp_model_free(model);
}


/* Comments and blank lines appear in real MPS files; they must not shift the
 * counting pass away from the fill pass. */
static void test_comments_and_blank_lines(void)
{
	static const char *const lines[] = {"NAME          CMTS\n",
		"* a leading comment\n", "\n", "ROWS\n", " N  COST\n",
		"* a comment inside ROWS\n", " L  R0\n", "\n", "COLUMNS\n",
		"    X0        COST                 1   R0                   1\n",
		"\n", "RHS\n", "    RHS1      R0                   5\n",
		"ENDATA\n", NULL};
	struct lp_simplex_Result result;
	const char *file;
	double x[1];
	int state;

	file = write_fixture(lines);
	if (file == NULL) {
		fail("comments", "could not write fixture");
		return;
	}
	state = solve_file(file, lp_simplex_ALGORITHM_DUAL_REVISED, x, &result);
	remove(TEMP_FILE);
	if (state != 1)
		fail("comments", "model rejected");
	else if (result.primal_infeasibility > 1e-7)
		fail("comments", "infeasible result");
	else
		pass("comments and blank lines tolerated");
}


/* Both OBJSENSE layouts must preserve the direction and constant term. */
static void test_objective_senses(void)
{
	static const char *const senses[] = {"OBJSENSE\n    MAXIMIZE\n",
		"OBJSENSE MAXIMIZE\n", "OBJSENSE max\n",
		"OBJSENSE\n* comment\n\n    MINIMIZE\n", "OBJSENSE MIN\n",
		"OBJSENSE minimize\n", "OBJSENSE MAXIMUM\n", "OBJSENSE\n",
		"OBJSENSE MAX EXTRA\n"};
	const char *lines[] = {"NAME          SENSE\n", NULL, "ROWS\n",
		" N  COST\n", " L  R0\n", "COLUMNS\n",
		"    X0        COST                 1   R0                   1\n",
		"RHS\n",
		"    RHS1      R0                   5   COST                 7\n",
		"BOUNDS\n", " LO BND1      X0                   0\n",
		" UP BND1      X0                  10\n", "ENDATA\n", NULL};
	int variant, algorithm, presolve;

	for (variant = 0; variant < 9; variant++) {
		struct lp_Model *model;
		const char *file;

		lines[1] = senses[variant];
		file = write_fixture(lines);
		if (file == NULL) {
			fail("objsense", "could not write fixture");
			return;
		}
		model = lp_read_mps(file);
		remove(TEMP_FILE);
		if (variant >= 6) {
			if (model != NULL)
				fail("objsense", "invalid sense was accepted");
			lp_model_free(model);
			continue;
		}
		if (model == NULL) {
			fail("objsense", "valid sense was rejected");
			continue;
		}
		for (algorithm = 0; algorithm < 3; algorithm++)
			for (presolve = 0; presolve < 2; presolve++) {
				struct lp_simplex_Options options;
				struct lp_simplex_Result result;
				double x[1], expected_x = variant < 3 ? 5. : 0.;
				int state;

				lp_simplex_default_options(&options, algorithm);
				options.presolve = presolve;
				state = lp_simplex_solve(
					model, &options, x, &result);
				if (state != lp_simplex_EXIT_SUCCESS ||
					result.status != lp_simplex_Success)
					fail("objsense", "solve failed");
				else if (!isfinite(result.objective) ||
					 fabs(result.objective -
						 (expected_x - 7.)) > 1e-7 ||
					 !isfinite(x[0]) ||
					 fabs(x[0] - expected_x) > 1e-7)
					fail("objsense",
						"wrong optimum or objective offset");
			}
		lp_model_free(model);
	}
}


/* A file without a ROWS section used to reach a negative-size allocation. */
static void test_malformed_without_rows(void)
{
	static const char *const lines[] = {"NAME          BAD\n", " N  COST\n",
		" L  R0\n", "COLUMNS\n",
		"    X0        COST                 1   R0                   1\n",
		"RHS\n", "    RHS1      R0                  10\n", "RANGES\n",
		"    RNG1      R0                   5\n", "ENDATA\n", NULL};
	struct lp_Model *model;
	const char *file;

	file = write_fixture(lines);
	if (file == NULL) {
		fail("malformed", "could not write fixture");
		return;
	}
	model = lp_read_mps(file);
	remove(TEMP_FILE);
	if (model != NULL) {
		lp_model_free(model);
		fail("malformed", "accepted a file without ROWS");
	} else {
		pass("malformed MPS refused");
	}
}


/* Integer models are out of scope; they must be refused rather than relaxed. */
static void test_integer_markers_rejected(void)
{
	static const char *const lines[] = {"NAME          INTT\n", "ROWS\n",
		" N  COST\n", " L  R0\n", "COLUMNS\n",
		"    MARKER                 'MARKER'                 'INTORG'\n",
		"    X0        COST                 1   R0                   1\n",
		"    MARKER                 'MARKER'                 'INTEND'\n",
		"RHS\n", "    RHS1      R0                   5\n", "ENDATA\n",
		NULL};
	struct lp_simplex_Result result;
	struct lp_Model *model;
	const char *file;
	double x[1];

	file = write_fixture(lines);
	if (file == NULL) {
		fail("integer", "could not write fixture");
		return;
	}
	model = lp_read_mps(file);
	remove(TEMP_FILE);
	if (model == NULL) {
		fail("integer", "model with integer markers could not be read");
		return;
	}
	{
		struct lp_simplex_Options options;
		int state;

		lp_simplex_default_options(
			&options, lp_simplex_ALGORITHM_DUAL_REVISED);
		state = lp_simplex_solve(model, &options, x, &result);
		lp_model_free(model);
		if (state == lp_simplex_EXIT_SUCCESS)
			fail("integer",
				"integer markers were silently relaxed");
		else
			pass("integer model refused");
	}
}


/*
 * An infeasible model whose data spans several orders of magnitude.  Bound
 * propagation used to drive its bounds towards overflow, the reduced solve
 * produced an infinite objective, and the result was reported as success.
 */
static void test_infeasible_scaled_model(void)
{
	static const char *const lines[] = {"NAME          HARD\n", "ROWS\n",
		" N  COST\n", " L  R0\n", " E  R1\n", " E  R2\n", " G  R3\n",
		" E  R4\n", " G  R5\n", "COLUMNS\n",
		"    X0        COST                 2   R0               0.001\n",
		"    X0        R1                   1   R2                   1\n",
		"    X0        R4                   5   R5                   3\n",
		"    X1        COST                 0   R0                   3\n",
		"    X1        R1               0.001   R2                   1\n",
		"    X1        R3              0.0001   R4                  -1\n",
		"    X2        COST                -3   R0               0.001\n",
		"    X2        R1                0.01   R2                0.01\n",
		"    X2        R5                   2\n",
		"    X3        COST                -3   R0                   2\n",
		"    X3        R1              0.0001   R2                0.01\n",
		"    X3        R4                0.01   R5                   2\n",
		"RHS\n",
		"    RHS1      R0                 100   R1                 100\n",
		"    RHS1      R2                 -10   R3                  30\n",
		"    RHS1      R5                 100\n", "BOUNDS\n",
		" LO BND1      X0                 -50\n", " FR BND1      X1\n",
		" LO BND1      X2                   0\n",
		" UP BND1      X2                   2\n", " FR BND1      X3\n",
		"ENDATA\n", NULL};
	struct lp_simplex_Result result;
	struct lp_Model *model;
	const char *file;
	double *x;
	int state;

	file = write_fixture(lines);
	if (file == NULL) {
		fail("scaled infeasible", "could not write fixture");
		return;
	}
	model = lp_read_mps(file);
	remove(TEMP_FILE);
	if (model == NULL) {
		fail("scaled infeasible", "could not read fixture");
		return;
	}
	x = (double *)malloc((size_t)model->n * sizeof(double));
	if (x == NULL) {
		lp_model_free(model);
		fail("scaled infeasible", "allocation failed");
		return;
	}
	{
		struct lp_simplex_Options options;

		lp_simplex_default_options(
			&options, lp_simplex_ALGORITHM_DUAL_REVISED);
		memset(&result, 0, sizeof(result));
		state = lp_simplex_solve(model, &options, x, &result);
		if (state == lp_simplex_EXIT_SUCCESS)
			fail("scaled infeasible",
				"infeasible model reported as success");
		else if (result.objective != result.objective ||
			 result.objective - result.objective != 0.)
			fail("scaled infeasible",
				"non-finite objective reported");
		else
			pass("scaled infeasible model refused");
	}
	free(x);
	lp_model_free(model);
}


/* The invariant above, applied to every bundled example and every algorithm. */
static void test_success_contract_on_examples(void)
{
	static const char *names[] = {"degenerate_pivot.mps",
		"energy_system.mps", "equality_system.mps",
		"free_and_lower_bounds.mps", "free_variables.mps",
		"mixed_inequalities.mps", "production_planning.mps",
		"ranged_rows.mps", "variable_bounds.mps"};
	static const int algorithms[] = {lp_simplex_ALGORITHM_TABLEAU,
		lp_simplex_ALGORITHM_DUAL_REVISED,
		lp_simplex_ALGORITHM_PAN_BDA};
	char path[1024];
	int i, k;
	int checked = 0;

	for (k = 0; k < 3; k++)
		for (i = 0; i < (int)(sizeof(names) / sizeof(names[0])); i++) {
			if (strlen(SIMPLE_EXAMPLES_DIR) + strlen(names[i]) + 2 >
				sizeof(path)) {
				fail("success contract",
					"path buffer too small");
				return;
			}
			strcpy(path, SIMPLE_EXAMPLES_DIR);
			strcat(path, "/");
			strcat(path, names[i]);
			check_success_is_honest(path, algorithms[k], names[i]);
			checked++;
		}
	printf("       success contract checked on %d solver/example pairs\n",
		checked);
}


/* Beale's cycling example, with repeated constraints and column permutations.
 * Verify the original objective and constraints with presolve both on and off. */
static void test_degenerate_cycling_model(void)
{
	static const double coefficients[3][4] = {
		{.5, -5.5, -2.5, 9.}, {.5, -1.5, -.5, 1.}, {1., 0., 0., 0.}};
	static const double costs[4] = {-10., 57., 9., 24.};
	int permutation, presolve, i, j;

	for (presolve = 0; presolve <= 1; presolve++)
		for (permutation = 0; permutation < 4; permutation++) {
			struct lp_Model *model = lp_model_create(5, 4);
			struct lp_simplex_Options options;
			struct lp_simplex_Result result;
			double x[4];
			int state;

			if (model == NULL) {
				fail("degenerate cycling", "allocation failed");
				return;
			}
			for (j = 0; j < 4; j++)
				model->objective[j] =
					costs[(j + permutation) % 4];
			for (i = 0; i < 5; i++) {
				int row = i < 3 ? i : i - 3;

				model->constraints[i].type = optm_CONS_T_LE;
				model->constraints[i].rhs = row == 2 ? 1. : 0.;
				for (j = 0; j < 4; j++)
					model->constraints[i].coef[j] =
						coefficients[row]
							    [(j + permutation) %
								    4];
			}
			if (lp_model_build_sparse(model) !=
				lp_simplex_EXIT_SUCCESS) {
				lp_model_free(model);
				fail("degenerate cycling",
					"sparse conversion failed");
				return;
			}
			lp_simplex_default_options(
				&options, lp_simplex_ALGORITHM_DUAL_REVISED);
			options.presolve = presolve;
			options.iteration_limit = 1000;
			state = lp_simplex_solve(model, &options, x, &result);
			if (state != lp_simplex_EXIT_SUCCESS ||
				result.status != lp_simplex_Success ||
				fabs(result.objective + 1.) > 1e-7 ||
				primal_violation(model, x) >
					options.primal_tolerance)
				fail("degenerate cycling",
					"wrong optimum or infeasible original vector");
			lp_model_free(model);
		}
}


/* Harris may accept a small local dual violation.  Across a wide variable
 * bound that can hide a material objective gap, even with exact feasibility. */
static void test_wide_bound_objective_accuracy(void)
{
	struct lp_Model *model = lp_model_create(1, 2);
	struct lp_simplex_Options options;
	struct lp_simplex_Result result;
	double x[2];
	int j, state;

	if (model == NULL) {
		fail("wide-bound optimality", "allocation failed");
		return;
	}
	lp_simplex_default_options(&options, lp_simplex_ALGORITHM_DUAL_REVISED);
	options.presolve = 0;
	model->objective[0] = .5 * options.dual_tolerance;
	model->constraints[0].type = optm_CONS_T_EQ;
	model->constraints[0].rhs = 1e8;
	for (j = 0; j < 2; j++) {
		model->constraints[0].coef[j] = 1.;
		model->bounds[j].b_type = optm_BOUND_T_BS;
		model->bounds[j].ub = 1e8;
	}
	if (lp_model_build_sparse(model) != lp_simplex_EXIT_SUCCESS) {
		fail("wide-bound optimality", "sparse conversion failed");
		lp_model_free(model);
		return;
	}
	state = lp_simplex_solve(model, &options, x, &result);
	if (state != lp_simplex_EXIT_SUCCESS ||
		result.status != lp_simplex_Success ||
		fabs(result.objective) > options.dual_tolerance ||
		primal_violation(model, x) > options.primal_tolerance)
		fail("wide-bound optimality",
			"local reduced costs hid a large objective gap");
	else
		pass("wide-bound optimality");
	lp_model_free(model);
}


/* Extended-precision MPS numbers must not be truncated to twelve columns.
 * Check both numeric positions and the RHS/bound paths used by Netlib files. */
static void test_long_mps_numbers(void)
{
	FILE *stream = fopen(TEMP_FILE, "w");
	struct lp_Model *model;

	if (stream == NULL) {
		fail("long MPS numbers", "could not write fixture");
		return;
	}
	fputs("NAME          LONGNUM\nROWS\n N  COST\n E  R0\n L  R1\nCOLUMNS\n", stream);
	fprintf(stream, "    %-8s  %-8s  %s\n", "X0", "COST",
		"1.2345678901234567e-20");
	fprintf(stream, "    %-8s  %-8s  %12s   %-8s  %s\n", "X0", "R0",
		"1.25", "R1", "-0.17999999999999999");
	fputs("RHS\n", stream);
	fprintf(stream, "    %-8s  %-8s  %s\n", "RHS1", "R0",
		"-2447.8299999999999");
	fprintf(stream, "    %-8s  %-8s  %s\n", "RHS1", "R1",
		"1.3999999999999999");
	fputs("BOUNDS\n", stream);
	fprintf(stream, " FX %-8s  %-8s  %s\n", "BND1", "X0",
		"12345.678901234567");
	fputs("ENDATA\n", stream);
	fclose(stream);
	model = lp_read_mps(TEMP_FILE);
	remove(TEMP_FILE);
	if (model == NULL) {
		fail("long MPS numbers", "could not read fixture");
		return;
	}
	if (model->m != 2 || model->n != 1 ||
		model->objective[0] != 1.2345678901234567e-20 ||
		model->row_value[model->row_start[0]] != 1.25 ||
		model->row_value[model->row_start[1]] != -0.17999999999999999 ||
		model->constraints[0].rhs != -2447.8299999999999 ||
		model->constraints[1].rhs != 1.3999999999999999 ||
		model->bounds[0].lb != 12345.678901234567 ||
		model->bounds[0].ub != 12345.678901234567)
		fail("long MPS numbers", "numeric tokens changed while reading");
	else
		pass("long MPS numbers");
	lp_model_free(model);
}


int main(void)
{
	test_long_mps_numbers();
	test_comments_and_blank_lines();
	test_objective_senses();
	test_malformed_without_rows();
	test_integer_markers_rejected();
	test_infeasible_scaled_model();
	test_success_contract_on_examples();
	test_degenerate_cycling_model();
	test_wide_bound_objective_accuracy();

	if (failures != 0) {
		printf("%d regression check(s) failed\n", failures);
		return 1;
	}
	printf("all regression checks passed\n");
	return 0;
}
