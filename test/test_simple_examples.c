/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file test_simple_examples.c
 * @brief Table-driven regression test for the simple MPS example collection.
 *
 * Every example is loaded through lp_read_mps(), solved through the
 * public model wrapper, and compared with a documented prediction.  The test
 * checks model dimensions, successful termination, the optimal objective, and
 * where the reference solution is stable, the complete primal vector.
 *
 * Degenerate models can have several equally optimal basic feasible solutions.
 * Such examples deliberately set expected_x to NULL and test the invariant
 * optimal objective instead of requiring one arbitrary basis.
 */

#include <lp_simplex/lp_simplex.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef SIMPLE_EXAMPLES_DIR
/**
 * Fallback for editors or manual compilation that do not import CMake target
 * definitions.  CMake overrides this with an absolute source-tree path; the
 * fallback assumes the executable is launched from the repository root.
 */
#define SIMPLE_EXAMPLES_DIR "data/simple_examples"
#endif


/** Description and prediction for one MPS regression example. */
struct ExampleSpec {
	const char *label;          /**< Human-readable test name. */
	const char *filename;       /**< File relative to SIMPLE_EXAMPLES_DIR. */
	const char *criteria;       /**< Simplex pivot rule used by the example. */
	int iteration_limit;        /**< Maximum Phase I plus Phase II pivots. */
	int expected_m;             /**< Predicted number of constraints. */
	int expected_n;             /**< Predicted number of variables. */
	double expected_value;      /**< Predicted optimal objective value. */
	double tolerance;           /**< Absolute and relative comparison scale. */
	const double *expected_x;   /**< Predicted primal vector, or NULL. */
};

static const double production_solution[] = {4.5, 4.5};
static const double mixed_solution[] = {6., 4.};
static const double free_lower_solution[] = {10., -3.};
static const double free_solution[] = {2. / 3., 4. / 3.};
static const double equality_solution[] = {0.5, 1.25, 0., 1.};
static const double bounds_solution[] = {-4., 3., 4.};
static const double ranged_solution[] = {6.};

/**
 * Expected results for the complete simple example collection.
 *
 * The degenerate-pivot and energy-system examples may admit basis-dependent
 * primal representations, so their success criterion is the optimal value.
 */
static const struct ExampleSpec examples[] = {
	{
		"production planning", "production_planning.mps", "bland", 1000,
		4, 2, -22.5, 1e-8, production_solution
	},
	{
		"mixed inequalities", "mixed_inequalities.mps", "bland", 1000,
		3, 2, -34., 1e-8, mixed_solution
	},
	{
		"degenerate pivot", "degenerate_pivot.mps", "bland", 1000,
		3, 4, -1.25, 1e-8, NULL
	},
	{
		"free and lower bounds", "free_and_lower_bounds.mps", "bland", 1000,
		2, 2, -22., 1e-8, free_lower_solution
	},
	{
		"free variables", "free_variables.mps", "bland", 1000,
		6, 2, -10. / 9., 1e-7, free_solution
	},
	{
		"equality system", "equality_system.mps", "bland", 1000,
		4, 4, 7. / 4., 1e-8, equality_solution
	},
	{
		"energy system", "energy_system.mps", "dantzig", 10000,
		27, 32, -464.753142857143, 1e-7, NULL
	},
	{
		"variable bounds", "variable_bounds.mps", "bland", 1000,
		1, 3, -5., 1e-8, bounds_solution
	},
	{
		"ranged rows", "ranged_rows.mps", "bland", 1000,
		8, 1, -6., 1e-8, ranged_solution
	}
};

/** Return the absolute value without relying on non-C90 math helpers. */
static double absolute_value(double value)
{
	return value < 0. ? -value : value;
}

/** Compare floating-point values with a combined absolute/relative tolerance. */
static int nearly_equal(double actual, double expected, double tolerance)
{
	double scale = absolute_value(expected);
	if (scale < 1.)
		scale = 1.;
	return absolute_value(actual - expected) <= tolerance * scale;
}

/** Build a path into caller-provided storage, returning zero if it will not fit. */
static int build_example_path(char *path, size_t capacity, const char *filename)
{
	const char *directory = SIMPLE_EXAMPLES_DIR;
	size_t directory_length = strlen(directory);
	size_t filename_length = strlen(filename);

	if (directory_length + 1 + filename_length + 1 > capacity)
		return 0;
	memcpy(path, directory, directory_length);
	path[directory_length] = '/';
	memcpy(path + directory_length + 1, filename, filename_length + 1);
	return 1;
}

/**
 * Load, solve, and verify one example.
 *
 * @return 1 when every prediction is reached; otherwise 0 after printing the
 *         exact failed condition.
 */
static int run_example(const struct ExampleSpec *example, const int algorithm)
{
	char path[1024];
	struct lp_Model *model;
	struct lp_simplex_Options options;
	struct lp_simplex_Result result;
	double *x;
	double value = 0.;
	int state;
	int i;
	int passed = 1;

	if (!build_example_path(path, sizeof(path), example->filename)) {
		printf("[FAIL] %s: MPS path is too long\n", example->label);
		return 0;
	}
	model = lp_read_mps(path);
	if (model == NULL) {
		printf("[FAIL] %s: could not read %s\n", example->label, path);
		return 0;
	}
	if (model->m != example->expected_m || model->n != example->expected_n) {
		printf("[FAIL] %s: dimensions=(%d,%d), expected=(%d,%d)\n",
		       example->label, model->m, model->n,
		       example->expected_m, example->expected_n);
		lp_model_free(model);
		return 0;
	}

	x = (double *)malloc((size_t)model->n * sizeof(double));
	if (x == NULL) {
		printf("[FAIL] %s: solution allocation failed\n", example->label);
		lp_model_free(model);
		return 0;
	}
	lp_simplex_default_options(&options, algorithm);
	options.iteration_limit = example->iteration_limit;
	if (algorithm == lp_simplex_ALGORITHM_TABLEAU)
		options.pricing = strcmp(example->criteria, "dantzig") == 0
			? lp_simplex_PRICING_DANTZIG : lp_simplex_PRICING_BLAND;
	state = lp_simplex_solve(model, &options, x, &result);
	value = result.objective;
	if (state != lp_simplex_EXIT_SUCCESS || result.status != lp_simplex_Success) {
		printf("[FAIL] %s: state=%d, code=%d; expected successful optimum\n",
		       example->label, state, result.status);
		passed = 0;
	} else if (algorithm == lp_simplex_ALGORITHM_TABLEAU &&
		   result.iterations <= 0) {
		printf("[FAIL] %s: tableau did not report its pivot iterations\n",
		       example->label);
		passed = 0;
	} else if (result.primal_infeasibility > 10. *
		   (example->tolerance > options.primal_tolerance
		    ? example->tolerance : options.primal_tolerance)) {
		printf("[FAIL] %s: primal infeasibility=%.15g\n",
		       example->label, result.primal_infeasibility);
		passed = 0;
	} else if (!nearly_equal(value, example->expected_value, example->tolerance)) {
		printf("[FAIL] %s: objective=%.15g, expected=%.15g\n",
		       example->label, value, example->expected_value);
		passed = 0;
	} else if (example->expected_x != NULL) {
		double solution_tolerance = example->tolerance >
			options.primal_tolerance ? example->tolerance :
			options.primal_tolerance;
		for (i = 0; i < model->n; i++) {
			if (!nearly_equal(x[i], example->expected_x[i],
				10. * solution_tolerance)) {
				printf("[FAIL] %s: x[%d]=%.15g, expected=%.15g\n",
				       example->label, i, x[i], example->expected_x[i]);
				passed = 0;
				break;
			}
		}
	}
	if (passed) {
		const char *algorithm_name = algorithm == lp_simplex_ALGORITHM_TABLEAU
			? "tableau" : algorithm == lp_simplex_ALGORITHM_PAN_BDA
			? "pan-bda" : "dual-revised";
		printf("[PASS] %-12s %-23s objective=% .12g (expected % .12g)\n",
		       algorithm_name,
		       example->label, value, example->expected_value);
	}

	free(x);
	lp_model_free(model);
	return passed;
}

/** Run all MPS examples and fail the process if any prediction is missed. */
int main(void)
{
	const int count = (int)(sizeof(examples) / sizeof(examples[0]));
	int failures = 0;
	int algorithm, i;

	for (algorithm = lp_simplex_ALGORITHM_TABLEAU;
	     algorithm <= lp_simplex_ALGORITHM_PAN_BDA; algorithm++) {
		for (i = 0; i < count; i++) {
			if (!run_example(examples + i, algorithm))
				failures++;
		}
	}
	printf("%d/%d solver/example pairs reached their predicted results.\n",
	       3 * count - failures, 3 * count);
	return failures == 0 ? 0 : 1;
}
