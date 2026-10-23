/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file lp_solve_mps.c
 * @brief Solve one MPS model and optionally validate an expected result.
 *
 * Dataset preparation and reference lookup belong to benchmark.py. This
 * driver only uses explicit --expected or --infeasible predictions; without
 * either option, a successful solve is required.
 */

#include <lp_simplex/lp_simplex.h>

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#define DEFAULT_ITERATION_LIMIT LP_SIMPLEX_DEFAULT_AUTO_ITERATION_LIMIT
#define DEFAULT_TOLERANCE 1e-7
#define PATH_BUFFER_SIZE 1024

/** Parsed command-line configuration. */
struct Options {
	const char *mps_path;
	const char *criteria;
	int algorithm;
	int iteration_limit;
	int presolve;
	double tolerance;
	int print_solution;
	int progress;
	double progress_interval;
	const char *save_point;
	int expected_kind; /**< 0 = success, 1 = objective, 2 = infeasible. */
	double expected_value;
};

/** Print command-line documentation. */
static void print_usage(const char *program)
{
	printf("Usage: %s [options] MODEL.mps\n", program);
	printf("\nOptions:\n");
	printf("  --algorithm NAME      auto (default), dual-revised, pan-bda, or tableau\n");
	printf("  --criteria RULE       manual tableau/Pan pricing: dantzig, normalized, bland\n");
	printf("  --iterations N        total work limit (default: %d)\n",
		DEFAULT_ITERATION_LIMIT);
	printf("  --no-presolve         bypass presolve for differential diagnostics\n");
	printf("  --tolerance VALUE     objective relative tolerance (default: %.0e)\n",
		DEFAULT_TOLERANCE);
	printf("  --expected VALUE      check the objective against VALUE\n");
	printf("  --infeasible          expect an infeasible model\n");
	printf("  --print-solution      print every primal variable\n");
	printf("  --progress            report original-model bounds and gap to stderr\n");
	printf("  --progress-interval S positive CPU seconds between reports (default 5)\n");
	printf("  --save-point PATH     atomically checkpoint improved feasible points during solve\n");
	printf("  --help                show this help text\n");
}


/** Parse a complete decimal floating-point argument. */
static int parse_double(const char *text, double *value)
{
	char *end;

	errno = 0;
	*value = strtod(text, &end);
	return errno == 0 && end != text && *end == '\0';
}


/** Parse a positive decimal integer argument. */
static int parse_positive_int(const char *text, int *value)
{
	char *end;
	long parsed;

	errno = 0;
	parsed = strtol(text, &end, 10);
	if (errno != 0 || end == text || *end != '\0' || parsed <= 0 ||
		parsed > 2147483647L)
		return 0;
	*value = (int)parsed;
	return 1;
}


/** Consume the argument following an option, or report a missing value. */
static const char *option_value(int argc, char **argv, int *index)
{
	if (*index + 1 >= argc)
		return NULL;
	(*index)++;
	return argv[*index];
}


/** Parse command-line arguments; return 1 on success, 0 on a usage error. */
static int parse_options(int argc, char **argv, struct Options *options)
{
	int i;

	options->mps_path = NULL;
	options->criteria = "dantzig";
	options->algorithm = lp_simplex_ALGORITHM_AUTO;
	options->iteration_limit = DEFAULT_ITERATION_LIMIT;
	options->presolve = 1;
	options->tolerance = DEFAULT_TOLERANCE;
	options->print_solution = 0;
	options->progress = 0;
	options->progress_interval = 5.;
	options->save_point = NULL;
	options->expected_kind = 0;
	options->expected_value = 0.;

	for (i = 1; i < argc; i++) {
		const char *argument = argv[i];
		const char *value;

		if (strcmp(argument, "--help") == 0) {
			print_usage(argv[0]);
			exit(0);
		} else if (strcmp(argument, "--progress") == 0) {
			options->progress = 1;
		} else if (strcmp(argument, "--progress-interval") == 0) {
			value = option_value(argc, argv, &i);
			if (!value ||
				!parse_double(
					value, &options->progress_interval) ||
				!isfinite(options->progress_interval) ||
				options->progress_interval <= 0.)
				return 0;
			options->progress = 1;
		} else if (strcmp(argument, "--print-solution") == 0) {
			options->print_solution = 1;
		} else if (strcmp(argument, "--save-point") == 0) {
			options->save_point = option_value(argc, argv, &i);
			if (options->save_point == NULL)
				return 0;
		} else if (strcmp(argument, "--infeasible") == 0) {
			options->expected_kind = 2;
		} else if (strcmp(argument, "--no-presolve") == 0) {
			options->presolve = 0;
		} else if (strcmp(argument, "--algorithm") == 0) {
			value = option_value(argc, argv, &i);
			if (value == NULL)
				return 0;
			if (strcmp(value, "auto") == 0)
				options->algorithm = lp_simplex_ALGORITHM_AUTO;
			else if (strcmp(value, "tableau") == 0)
				options->algorithm =
					lp_simplex_ALGORITHM_TABLEAU;
			else if (strcmp(value, "dual-revised") == 0)
				options->algorithm =
					lp_simplex_ALGORITHM_DUAL_REVISED;
			else if (strcmp(value, "pan-bda") == 0)
				options->algorithm =
					lp_simplex_ALGORITHM_PAN_BDA;
			else
				return 0;
		} else if (strcmp(argument, "--criteria") == 0) {
			value = option_value(argc, argv, &i);
			if (value == NULL ||
				(strcmp(value, "bland") != 0 &&
					strcmp(value, "dantzig") != 0 &&
					strcmp(value, "normalized") != 0 &&
					strcmp(value, "pan97") != 0))
				return 0;
			options->criteria = value;
		} else if (strcmp(argument, "--iterations") == 0) {
			value = option_value(argc, argv, &i);
			if (value == NULL || !parse_positive_int(value,
						     &options->iteration_limit))
				return 0;
		} else if (strcmp(argument, "--tolerance") == 0) {
			value = option_value(argc, argv, &i);
			if (value == NULL ||
				!parse_double(value, &options->tolerance) ||
				options->tolerance <= 0.)
				return 0;
		} else if (strcmp(argument, "--expected") == 0) {
			value = option_value(argc, argv, &i);
			if (value == NULL ||
				!parse_double(value, &options->expected_value))
				return 0;
			options->expected_kind = 1;
		} else if (argument[0] == '-') {
			return 0;
		} else if (options->mps_path == NULL) {
			options->mps_path = argument;
		} else {
			return 0;
		}
	}
	return options->mps_path != NULL;
}


/** Return the filename component of a POSIX or Windows path. */
static const char *path_basename(const char *path)
{
	const char *name = path;
	const char *cursor;

	for (cursor = path; *cursor != '\0'; cursor++) {
		if (*cursor == '/' || *cursor == '\\')
			name = cursor + 1;
	}
	return name;
}


/** Copy the MPS basename without its extension into fixed storage. */
static int model_name_from_path(const char *path, char *name, size_t capacity)
{
	const char *base = path_basename(path);
	size_t length = strlen(base);

	if (length > 4 && strcmp(base + length - 4, ".mps") == 0)
		length -= 4;
	if (length + 1 > capacity)
		return 0;
	memcpy(name, base, length);
	name[length] = '\0';
	return 1;
}


/** Return a stable human-readable label for a public solver code. */
static const char *solver_code_name(int code)
{
	switch (code) {
	case lp_simplex_Success:
		return "success";
	case lp_simplex_MemoryAllocError:
		return "memory allocation error";
	case lp_simplex_CondUnsatisfied:
		return "invalid input";
	case lp_simplex_ExceedIterLimit:
		return "iteration limit";
	case lp_simplex_Singularity:
		return "singularity";
	case lp_simplex_OverDetermination:
		return "overdetermination";
	case lp_simplex_Unboundedness:
		return "unbounded";
	case lp_simplex_Infeasibility:
		return "infeasible";
	case lp_simplex_Degeneracy:
		return "degeneracy";
	case lp_simplex_PrecisionError:
		return "precision error";
	default:
		return "unknown";
	}
}


/** Compare objectives with a combined absolute and relative tolerance. */
static int objectives_match(double actual, double expected, double tolerance)
{
	double scale = expected < 0. ? -expected : expected;
	double difference = actual - expected;

	if (difference < 0.)
		difference = -difference;
	if (scale < 1.)
		scale = 1.;
	return difference <= tolerance * scale;
}


/** Print the primal vector when requested. */
static void print_solution(const struct lp_Model *model, const double *x)
{
	int i;

	printf("solution:\n");
	for (i = 0; i < model->n; i++)
		printf("  %-8.8s = %.15g\n", model->bounds[i].name, x[i]);
}


/** Measure the returned point against the public model, after postsolve. */
static double original_primal_infeasibility(const struct lp_Model *model,
	const double *x, int *worst_row, int *worst_column)
{
	double maximum = 0.;
	int i, j, k;
	*worst_row = -1;
	*worst_column = -1;
	for (j = 0; j < model->n; j++) {
		double violation = 0.;

		if (!isfinite(x[j]))
			return HUGE_VAL;
		if ((model->bounds[j].b_type == optm_BOUND_T_LO ||
			    model->bounds[j].b_type == optm_BOUND_T_BS) &&
			x[j] < model->bounds[j].lb)
			violation = model->bounds[j].lb - x[j];
		if ((model->bounds[j].b_type == optm_BOUND_T_UP ||
			    model->bounds[j].b_type == optm_BOUND_T_BS) &&
			x[j] > model->bounds[j].ub &&
			x[j] - model->bounds[j].ub > violation)
			violation = x[j] - model->bounds[j].ub;
		if (violation > maximum) {
			maximum = violation;
			*worst_column = j;
			*worst_row = -1;
		}
	}
	for (i = 0; i < model->m; i++) {
		double activity = 0., correction = 0., violation = 0.;

		if (model->row_start != NULL)
			for (k = model->row_start[i];
				k < model->row_start[i + 1]; k++) {
				double a = model->row_value[k],
				       v = x[model->column_index[k]];
				double product = a * v,
				       next = activity + product;

				if (!isfinite(product) || !isfinite(next))
					return HUGE_VAL;
				correction +=
					(fabs(activity) >= fabs(product)
							? (activity - next) +
								  product
							: (product - next) +
								  activity) +
					fma(a, v, -product);
				activity = next;
			}
		else
			for (j = 0; j < model->n; j++) {
				double a = model->constraints[i].coef[j],
				       product = a * x[j],
				       next = activity + product;

				if (!isfinite(product) || !isfinite(next))
					return HUGE_VAL;
				correction +=
					(fabs(activity) >= fabs(product)
							? (activity - next) +
								  product
							: (product - next) +
								  activity) +
					fma(a, x[j], -product);
				activity = next;
			}
		activity += correction;
		if (!isfinite(activity))
			return HUGE_VAL;
		if (model->constraints[i].type == optm_CONS_T_EQ)
			violation = activity - model->constraints[i].rhs;
		else if (model->constraints[i].type == optm_CONS_T_GE &&
			 activity < model->constraints[i].rhs)
			violation = model->constraints[i].rhs - activity;
		else if (model->constraints[i].type == optm_CONS_T_LE &&
			 activity > model->constraints[i].rhs)
			violation = activity - model->constraints[i].rhs;
		if (violation < 0.)
			violation = -violation;
		if (violation > maximum) {
			maximum = violation;
			*worst_row = i;
			*worst_column = -1;
		}
	}
	return maximum;
}


/** Load, solve, report, and validate one MPS model. */
static const char *dual_residual_name(int space)
{
	switch (space) {
	case lp_simplex_RESIDUAL_ORIGINAL:
		return "original-model";
	case lp_simplex_RESIDUAL_PRESOLVED:
		return "presolved-model";
	case lp_simplex_RESIDUAL_STANDARD:
		return "standard-form";
	case lp_simplex_RESIDUAL_SCALED:
		return "scaled-model";
	case lp_simplex_RESIDUAL_AUXILIARY:
		return "auxiliary-model";
	default:
		return "unavailable";
	}
}


struct DriverProgress {
	const struct Options *options;
	int columns;
	int checkpoint_failed;
};

/* Write a sibling file first: a killed solve cannot truncate the last good
 * point. rename() replaces it only after every coordinate has been closed. */
static int checkpoint_point(const char *path, const double *x, int n)
{
	char *temporary = malloc(strlen(path) + 5);
	FILE *stream;
	int j, written = 1;
	if (!temporary)
		return 0;
	sprintf(temporary, "%s.tmp", path);
	stream = fopen(temporary, "w");
	if (!stream) {
		free(temporary);
		return 0;
	}
	for (j = 0; j < n; j++)
		if (fprintf(stream, "%.17g\n", x[j]) < 0) {
			written = 0;
			break;
		}
	if (fclose(stream))
		written = 0;
	if (written && rename(temporary, path))
		written = 0;
	if (!written)
		remove(temporary);
	free(temporary);
	return written;
}
static void print_bound(const char *name, double bound)
{
	fprintf(stderr, " %s=", name);
	if (isfinite(bound))
		fprintf(stderr, "%.15g", bound);
	else
		fprintf(stderr, "unavailable");
}
static void driver_progress(void *context, const struct lp_simplex_Progress *p)
{
	struct DriverProgress *c = context;
	if (c->options->progress) {
		fprintf(stderr,
			"progress: cpu=%.3f stage=%s stage-iterations=%d",
			p->cpu_seconds, p->stage, p->stage_iterations);
		print_bound("lower", p->lower_bound);
		print_bound("upper", p->upper_bound);
		if (isfinite(p->relative_gap))
			fprintf(stderr, " gap=%.9g%% absolute-gap=%.9g",
				100. * p->relative_gap, p->absolute_gap);
		else
			fprintf(stderr, " gap=unavailable");
		fprintf(stderr, "\n");
	}
	if (p->point && c->options->save_point) {
		if (checkpoint_point(
			    c->options->save_point, p->point, c->columns)) {
			c->checkpoint_failed = 0;
			fprintf(stderr, "saved point: primal=%.6g path=%s\n",
				p->primal_infeasibility,
				c->options->save_point);
		} else {
			c->checkpoint_failed = 1;
			fprintf(stderr, "could not checkpoint point %s: %s\n",
				c->options->save_point, strerror(errno));
		}
	}
	fflush(stderr);
}


int main(int argc, char **argv)
{
	struct Options options;
	struct lp_Model *model;
	struct lp_simplex_Options solve_options;
	struct lp_simplex_Result result;
	struct DriverProgress progress;
	double *x;
	double objective = 0.;
	clock_t started, finished;
	double elapsed;
	double original_infeasibility = 0.;
	char model_name[128];
	int state, worst_row = -1, worst_column = -1;
	int passed;

	if (!parse_options(argc, argv, &options)) {
		print_usage(argv[0]);
		return 2;
	}
	if (!model_name_from_path(
		    options.mps_path, model_name, sizeof(model_name))) {
		fprintf(stderr, "error: model filename is too long\n");
		return 2;
	}
	model = lp_read_mps(options.mps_path);
	if (model == NULL) {
		fprintf(stderr, "[FAIL] %s: could not read MPS model\n",
			model_name);
		return 1;
	}
	x = (double *)malloc((size_t)model->n * sizeof(double));
	if (x == NULL) {
		fprintf(stderr, "[FAIL] %s: solution allocation failed\n",
			model_name);
		lp_model_free(model);
		return 1;
	}

	started = clock();
	lp_simplex_default_options(&solve_options, options.algorithm);
	solve_options.iteration_limit = options.iteration_limit;
	solve_options.presolve = options.presolve;
	progress.options = &options;
	progress.columns = model->n;
	progress.checkpoint_failed = 0;
	if (options.progress || options.save_point) {
		solve_options.progress_callback = driver_progress;
		solve_options.progress_context = &progress;
		solve_options.progress_interval = options.progress_interval;
	}
	if (options.algorithm == lp_simplex_ALGORITHM_TABLEAU)
		solve_options.pricing = strcmp(options.criteria, "dantzig") == 0
						? lp_simplex_PRICING_DANTZIG
						: lp_simplex_PRICING_BLAND;
	else if (options.algorithm == lp_simplex_ALGORITHM_PAN_BDA)
		solve_options.pricing =
			strcmp(options.criteria, "normalized") == 0
				? lp_simplex_PRICING_PAN_NORMALIZED
				: lp_simplex_PRICING_DANTZIG;
	state = lp_simplex_solve(model, &solve_options, x, &result);
	objective = result.objective;
	if (state == lp_simplex_EXIT_SUCCESS)
		original_infeasibility = original_primal_infeasibility(
			model, x, &worst_row, &worst_column);
	finished = clock();
	elapsed = (double)(finished - started) / (double)CLOCKS_PER_SEC;

	printf("model:      %s\n", model_name);
	printf("dimensions: %d constraints, %d variables\n", model->m,
		model->n);
	printf("algorithm:  %s\n",
		options.algorithm == lp_simplex_ALGORITHM_TABLEAU   ? "tableau"
		: options.algorithm == lp_simplex_ALGORITHM_PAN_BDA ? "pan-bda"
		: options.algorithm == lp_simplex_ALGORITHM_AUTO
			? "auto"
			: "dual-revised");
	printf("pricing:    %s\n",
		options.algorithm == lp_simplex_ALGORITHM_TABLEAU
			? options.criteria
		: options.algorithm == lp_simplex_ALGORITHM_PAN_BDA
			? (strcmp(options.criteria, "normalized") == 0
					  ? "normalized-violation"
					  : "dantzig-violation")
		: options.algorithm == lp_simplex_ALGORITHM_AUTO
			? "automatic"
			: "dual-steepest-edge");
	printf("result:     state=%d, code=%d (%s)\n", state, result.status,
		solver_code_name(result.status));
	printf("iterations: %d\n", result.iterations);
	printf("residuals:  primal=%.3g, dual=%.3g (%s)\n",
		result.primal_infeasibility, result.dual_infeasibility,
		dual_residual_name(result.dual_residual_space));
	if (state == lp_simplex_EXIT_SUCCESS)
		printf("postsolve:  primal=%.3g, worst-row=%d, worst-column=%d\n",
			original_infeasibility, worst_row, worst_column);
	printf("time:       %.6f seconds\n", elapsed);
	if (state == lp_simplex_EXIT_SUCCESS)
		printf("objective:  %.15g\n", objective);
	if (state == lp_simplex_EXIT_SUCCESS && options.print_solution)
		print_solution(model, x);

	if (options.expected_kind == 2) {
		passed = state == lp_simplex_EXIT_FAILURE &&
			 result.status == lp_simplex_Infeasibility;
		printf("expected:   infeasible\n");
	} else if (options.expected_kind == 1) {
		passed = state == lp_simplex_EXIT_SUCCESS &&
			 result.status == lp_simplex_Success &&
			 objectives_match(objective, options.expected_value,
				 options.tolerance);
		printf("expected:   objective %.15g (tolerance %.3g)\n",
			options.expected_value, options.tolerance);
	} else {
		passed = state == lp_simplex_EXIT_SUCCESS &&
			 result.status == lp_simplex_Success;
		printf("expected:   no reference; successful solve required\n");
	}
	printf("[%s] prediction %s\n", passed ? "PASS" : "FAIL",
		passed ? "reached" : "not reached");

	free(x);
	lp_model_free(model);
	return passed && !progress.checkpoint_failed ? 0 : 1;
}
