/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file test_pan_basis_stress.c
 * @brief Regression tests for Pan basis operations under repeated numerical
 * updates.
 */

#include "simplex_pan_config.h"
#include "simplex_pan_basis.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Long mixed update/rebuild sequences; equations are checked from original
 * columns, independently of the factor's own orthogonality statistics. */
int main(int argc, char **argv)
{
	struct simplex_PanConfig basis_controls;
	struct simplex_PanStandard standard;
	struct simplex_PanBasis basis;
	int large = argc > 1 && strcmp(argv[1], "--large") == 0;
	int m = large ? 10007 : 127, n = large ? 2048 : 62,
	    size = large ? 1024 : 31;
	int width = large ? 5 : m;
	int *start, *row;
	double *entry, *normalized, *norm, *rhs, *solution, *right, *dual;
	int i, j, k, t, failed = 0;
	double residual, worst = 0.;

	simplex_pan_config_init(&basis_controls, NULL);
	memset(&standard, 0, sizeof(standard));
	start = (int *)malloc((size_t)(n + 1) * sizeof(int));
	row = (int *)malloc((size_t)width * n * sizeof(int));
	entry = (double *)malloc((size_t)width * n * sizeof(double));
	normalized = (double *)malloc((size_t)width * n * sizeof(double));
	norm = (double *)calloc((size_t)n, sizeof(double));
	rhs = (double *)calloc((size_t)m * 4, sizeof(double));
	if (!start || !row || !entry || !normalized || !norm || !rhs)
		return 2;
	solution = rhs + m;
	right = solution + m;
	dual = right + m;
	for (j = 0; j < n; j++) {
		start[j] = j * width;
		for (i = 0; i < width; i++) {
			int actual =
				large ? (i == 0 ? j
						: n + (j * 37 + i * 997) %
									(m - n))
				      : i;

			k = j * width + i;
			row[k] = actual;
			entry[k] =
				large ? ldexp(i == 0 ? 2.
						     : .03 * sin((actual + 1.) *
								     (j + 3.)),
						j % 41 - 20)
				      : .03 * sin((i + 1.) * (j + 3.)) +
						(i == (j * 17 + 11) % m ? 2.
									: 0.);
			norm[j] = hypot(norm[j], entry[k]);
		}
		for (i = 0; i < width; i++)
			normalized[j * width + i] =
				entry[j * width + i] / norm[j];
	}
	start[n] = width * n;
	standard.rows = m;
	standard.columns = n;
	standard.nonzeros = width * n;
	standard.column_start = start;
	standard.row_index = row;
	standard.value = entry;
	standard.normalized_value = normalized;
	standard.column_norm = norm;
	if (simplex_pan_basis_create_configured(
		    &basis, &standard, 1e-10, &basis_controls))
		return 2;
	for (j = 0; j < size; j++)
		if (simplex_pan_basis_add(&basis, j))
			return 2;
	if (simplex_pan_basis_rebuild(&basis))
		return 2;
	for (t = 0; t < 1000 && !failed; t++) {
		int entering = (t * 7 + 13) % n;

		while (basis.position[entering] >= 0)
			entering = (entering + 1) % n;
		if (simplex_pan_basis_replace_factorized(
			    &basis, t % size, entering))
			failed = 1;
		if (t % 137 == 0 && simplex_pan_basis_rebuild(&basis))
			failed = 1;
		if (t % 17 != 0)
			continue;
		memset(rhs, 0, (size_t)m * sizeof(double));
		for (j = 0; j < size; j++) {
			int col = basis.column[j];
			double x = sin(col + 1.) / (large ? norm[col] : 1.);

			if (basis.position[col] != j)
				failed = 1;
			for (k = start[col]; k < start[col + 1]; k++)
				rhs[row[k]] += entry[k] * x;
		}
		if (simplex_pan_basis_least_squares(
			    &basis, rhs, solution, &residual))
			failed = 1;
		worst = fmax(worst, residual);
		for (j = 0; j < size; j++)
			worst = fmax(worst,
				fabs(solution[j] *
						(large ? norm[basis.column[j]]
						       : 1.) -
					sin(basis.column[j] + 1.)));
		for (i = 0; i < m; i++)
			dual[i] = cos(i + 2.);
		for (j = 0; j < size; j++)
			right[j] = simplex_pan_column_dot(
				&standard, basis.column[j], dual);
		if (simplex_pan_basis_minimum_norm(
			    &basis, right, solution, &residual))
			failed = 1;
		if (!large)
			worst = fmax(worst, residual);
		for (j = 0; j < size; j++)
			worst = fmax(worst,
				fabs(simplex_pan_column_dot(&standard,
					     basis.column[j], solution) -
					right[j]) /
					(large ? 1. + fabs(right[j]) : 1.));
	}
	failed = failed || worst > 1e-10;
	printf("Pan 1000 mixed QR updates rows=%d columns=%d active=%d error=%.6g rebuilds=%ld %s\n",
		m, n, size, worst, basis.factorizations,
		failed ? "failed" : "passed");
	simplex_pan_basis_destroy(&basis);
	free(start);
	free(row);
	free(entry);
	free(normalized);
	free(norm);
	free(rhs);
	return failed ? 1 : 0;
}
