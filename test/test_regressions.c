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
static void check_success_is_honest(const char *file, int algorithm,
		const char *label)
{
	struct lp_simplex_Options options;
	struct lp_Model *model;
	struct lp_simplex_Result result;
	double *x;
	int state;

	(void)algorithm;
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
		int j;
		for (j = 0; j < model->n; j++) {
			if (x[j] != x[j] || x[j] - x[j] != 0.) {
				fail(label, "success returned a non-finite x");
				break;
			}
		}
		if (result.objective != result.objective ||
		    result.objective - result.objective != 0.)
			fail(label, "success returned a non-finite objective");
		else if (result.primal_infeasibility >
			 10. * options.primal_tolerance)
			fail(label, "success returned an infeasible x");
		else
			pass(label);
	} else {
		/* A failure is always acceptable here; only wrong success is not. */
		pass(label);
	}
	free(x);
	lp_model_free(model);
}


/* Comments and blank lines appear in real MPS files; they must not shift the
 * counting pass away from the fill pass. */
static void test_comments_and_blank_lines(void)
{
	static const char *const lines[] = {
		"NAME          CMTS\n",
		"* a leading comment\n",
		"\n",
		"ROWS\n",
		" N  COST\n",
		"* a comment inside ROWS\n",
		" L  R0\n",
		"\n",
		"COLUMNS\n",
		"    X0        COST                 1   R0                   1\n",
		"\n",
		"RHS\n",
		"    RHS1      R0                   5\n",
		"ENDATA\n",
		NULL
	};
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


/* OBJSENSE MAXIMIZE must be honoured, not silently solved as a minimisation. */
static void test_objsense_maximize(void)
{
	static const char *const lines[] = {
		"NAME          MAXT\n",
		"OBJSENSE\n",
		"    MAXIMIZE\n",
		"ROWS\n",
		" N  COST\n",
		" L  R0\n",
		"COLUMNS\n",
		"    X0        COST                 1   R0                   1\n",
		"RHS\n",
		"    RHS1      R0                   5\n",
		"BOUNDS\n",
		" LO BND1      X0                   0\n",
		" UP BND1      X0                  10\n",
		"ENDATA\n",
		NULL
	};
	struct lp_simplex_Result result;
	const char *file;
	double x[1];
	int state;

	file = write_fixture(lines);
	if (file == NULL) {
		fail("objsense", "could not write fixture");
		return;
	}
	state = solve_file(file, lp_simplex_ALGORITHM_DUAL_REVISED, x, &result);
	remove(TEMP_FILE);
	if (state != 1)
		fail("objsense", "maximisation model rejected");
	else if (result.objective < 5. - 1e-7)
		fail("objsense", "MAX was ignored (objective is not 5)");
	else
		pass("OBJSENSE MAXIMIZE honoured");
}


/* A file without a ROWS section used to reach a negative-size allocation. */
static void test_malformed_without_rows(void)
{
	static const char *const lines[] = {
		"NAME          BAD\n",
		" N  COST\n",
		" L  R0\n",
		"COLUMNS\n",
		"    X0        COST                 1   R0                   1\n",
		"RHS\n",
		"    RHS1      R0                  10\n",
		"RANGES\n",
		"    RNG1      R0                   5\n",
		"ENDATA\n",
		NULL
	};
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
	static const char *const lines[] = {
		"NAME          INTT\n",
		"ROWS\n",
		" N  COST\n",
		" L  R0\n",
		"COLUMNS\n",
		"    MARKER                 'MARKER'                 'INTORG'\n",
		"    X0        COST                 1   R0                   1\n",
		"    MARKER                 'MARKER'                 'INTEND'\n",
		"RHS\n",
		"    RHS1      R0                   5\n",
		"ENDATA\n",
		NULL
	};
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
		lp_simplex_default_options(&options,
			lp_simplex_ALGORITHM_DUAL_REVISED);
		state = lp_simplex_solve(model, &options, x, &result);
		lp_model_free(model);
		if (state == lp_simplex_EXIT_SUCCESS)
			fail("integer", "integer markers were silently relaxed");
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
	static const char *const lines[] = {
		"NAME          HARD\n",
		"ROWS\n",
		" N  COST\n",
		" L  R0\n",
		" E  R1\n",
		" E  R2\n",
		" G  R3\n",
		" E  R4\n",
		" G  R5\n",
		"COLUMNS\n",
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
		"    RHS1      R5                 100\n",
		"BOUNDS\n",
		" LO BND1      X0                 -50\n",
		" FR BND1      X1\n",
		" LO BND1      X2                   0\n",
		" UP BND1      X2                   2\n",
		" FR BND1      X3\n",
		"ENDATA\n",
		NULL
	};
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
		lp_simplex_default_options(&options,
			lp_simplex_ALGORITHM_DUAL_REVISED);
		memset(&result, 0, sizeof(result));
		state = lp_simplex_solve(model, &options, x, &result);
		if (state == lp_simplex_EXIT_SUCCESS)
			fail("scaled infeasible",
				"infeasible model reported as success");
		else if (result.objective != result.objective ||
			 result.objective - result.objective != 0.)
			fail("scaled infeasible", "non-finite objective reported");
		else
			pass("scaled infeasible model refused");
	}
	free(x);
	lp_model_free(model);
}


/* The invariant above, applied to every bundled example and every algorithm. */
static void test_success_contract_on_examples(void)
{
	static const char *names[] = {
		"degenerate_pivot.mps", "energy_system.mps",
		"equality_system.mps", "free_and_lower_bounds.mps",
		"free_variables.mps", "mixed_inequalities.mps",
		"production_planning.mps", "ranged_rows.mps",
		"variable_bounds.mps"
	};
	static const int algorithms[] = {
		lp_simplex_ALGORITHM_TABLEAU,
		lp_simplex_ALGORITHM_DUAL_REVISED,
		lp_simplex_ALGORITHM_PAN_BDA
	};
	char path[1024];
	int i, k;
	int checked = 0;

	for (k = 0; k < 3; k++)
		for (i = 0; i < (int)(sizeof(names) / sizeof(names[0])); i++) {
			if (strlen(SIMPLE_EXAMPLES_DIR) + strlen(names[i]) + 2 >
			    sizeof(path)) {
				fail("success contract", "path buffer too small");
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


int main(void)
{
	test_comments_and_blank_lines();
	test_objsense_maximize();
	test_malformed_without_rows();
	test_integer_markers_rejected();
	test_infeasible_scaled_model();
	test_success_contract_on_examples();

	if (failures != 0) {
		printf("%d regression check(s) failed\n", failures);
		return 1;
	}
	printf("all regression checks passed\n");
	return 0;
}
