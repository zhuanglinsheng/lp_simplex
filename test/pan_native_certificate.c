/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file pan_native_certificate.c
 * @brief Original-model native solve with independently readable primal/dual
 * export.
 *
 * No presolve, reference objective or external solver is used.
 */

#include "simplex_pan_config.h"
#include <lp_simplex/lp_simplex.h>

#include "simplex_pan_bounded.h"
#include "simplex_problem.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static int write_vector(
	const char *prefix, const char *suffix, const double *x, int count)
{
	char *name = (char *)malloc(strlen(prefix) + strlen(suffix) + 1);
	FILE *f;
	int i, good = 1;

	if (!name)
		return 0;
	strcpy(name, prefix);
	strcat(name, suffix);
	f = fopen(name, "w");
	free(name);
	if (!f)
		return 0;
	for (i = 0; i < count; i++)
		if (fprintf(f, "%.17g\n", x[i]) < 0)
			good = 0;
	if (fclose(f))
		good = 0;
	return good;
}

int main(int argc, char **argv)
{
	struct simplex_PanConfig controls;
	struct lp_Model *model;
	struct simplex_Problem p;
	struct lp_simplex_Options o;
	struct lp_simplex_Result result;
	double *storage, *x, *y, *start;
	clock_t begun;
	int state, j, good = 0;
	FILE *f;

	if (argc < 4 || argc > 5) {
		fprintf(stderr,
			"usage: pan_native_certificate model.mps iterations output-prefix [start-point]\n");
		return 2;
	}
	model = lp_read_mps(argv[1]);
	if (!model)
		return 2;
	if (simplex_problem_from_model(&p, model)) {
		lp_model_free(model);
		return 2;
	}
	storage = (double *)calloc(
		2 * ((size_t)p.columns + 1) + (size_t)p.rows + 1,
		sizeof(double));
	if (!storage) {
		simplex_problem_destroy(&p);
		lp_model_free(model);
		return 2;
	}
	x = storage;
	start = x + p.columns + 1;
	y = start + p.columns + 1;
	if (argc == 5) {
		f = fopen(argv[4], "r");
		if (!f)
			goto done;
		for (j = 0; j < p.columns; j++)
			if (fscanf(f, "%lf", start + j) != 1 ||
				!isfinite(start[j])) {
				fclose(f);
				goto done;
			}
		fclose(f);
	}
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_PAN_BDA);
	o.presolve = 0;
	o.iteration_limit = atoi(argv[2]);
	if (o.iteration_limit < 0)
		goto done;
	begun = clock();
	simplex_pan_config_resolve(&controls, &p, &o);
	controls.model_form = SIMPLEX_PAN_NATIVE_BOUNDS;
	controls.direction_method = SIMPLEX_PAN_CONE_DIRECTION;
	if (controls.oracle == SIMPLEX_PAN_NO_ORACLE)
		controls.oracle = SIMPLEX_PAN_EXPLICIT_ORACLE;
	state = simplex_pan_bounded_solve_certified_configured(
		&p, &o, argc == 5 ? start : NULL, x, y, &result, &controls);
	good = state == 0 && result.status == lp_simplex_Success;
	if (good)
		good = write_vector(argv[3], "-primal.txt", x, p.columns) &&
		       write_vector(argv[3], "-dual.txt", y, p.rows);
	printf("native certificate: status=%d iterations=%d objective=%.17g primal=%.9g dual=%.9g cpu=%.6f exported=%d\n",
		result.status, result.iterations,
		result.objective + model->objective_offset,
		result.primal_infeasibility, result.dual_infeasibility,
		(double)(clock() - begun) / CLOCKS_PER_SEC, good);
done:
	free(storage);
	simplex_problem_destroy(&p);
	lp_model_free(model);
	return good ? 0 : 1;
}
