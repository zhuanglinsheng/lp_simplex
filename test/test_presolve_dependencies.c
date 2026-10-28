/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file test_presolve_dependencies.c
 * @brief Exact conservation dependencies and conservative refusal.
 */
#include <lp_simplex/lp_simplex.h>
#include "simplex_presolve_dependencies.h"
#include "simplex_presolve.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Independent dense elimination of the augmented equations checks that
 * publication preserves their affine row space, including partial searches. */
static int affine_rank(const struct lp_Model *m, const unsigned char *active)
{
	int i, j, k, row = 0, pivot, width = m->n + 1;
	double *a =
		(double *)calloc((size_t)m->m * (size_t)width, sizeof(double));
	if (!a)
		return -1;
	for (i = 0; i < m->m; i++)
		if (!active || active[i]) {
			for (j = 0; j < m->n; j++)
				a[i * width + j] = m->constraints[i].coef[j];
			a[i * width + m->n] = m->constraints[i].rhs;
		}
	for (j = 0; j < width && row < m->m; j++) {
		pivot = row;
		for (i = row + 1; i < m->m; i++)
			if (fabs(a[i * width + j]) > fabs(a[pivot * width + j]))
				pivot = i;
		if (fabs(a[pivot * width + j]) < 1e-9)
			continue;
		for (k = j; k < width; k++) {
			double v = a[row * width + k];
			a[row * width + k] = a[pivot * width + k];
			a[pivot * width + k] = v;
		}
		for (i = row + 1; i < m->m; i++) {
			double v = a[i * width + j] / a[row * width + j];
			for (k = j; k < width; k++)
				a[i * width + k] -= v * a[row * width + k];
		}
		row++;
	}
	free(a);
	return row;
}

static int marginal(int n, int i, int k, int side, int target)
{
	int own = side ? k : i;
	return 2 * n + (i * n + k) * 2 * (n - 1) + side * (n - 1) + target -
	       (target > own);
}

static struct lp_Model *conservation_model(int n, int permute)
{
	int rows = 2 * n + 2 * n * n * (n - 1);
	int cols = n * n + n * n * (n - 1) * (n - 1) / 2;
	struct lp_Model *m = lp_model_create(rows, cols);
	int i, j, k, l, a, b, q = n * n, r, c;
	double *matrix;
	if (!m)
		return NULL;
	for (i = 0; i < n; i++)
		for (k = 0; k < n; k++) {
			c = i * n + k;
			m->constraints[i].coef[c] = 1.;
			m->constraints[n + k].coef[c] = 1.;
			for (j = 0; j < n; j++) {
				if (j != i)
					m->constraints[marginal(n, i, k, 0, j)]
						.coef[c] = -1.;
				if (j != k)
					m->constraints[marginal(n, i, k, 1, j)]
						.coef[c] = -1.;
			}
		}
	for (i = 0; i < n; i++)
		for (j = i + 1; j < n; j++)
			for (k = 0; k < n; k++)
				for (l = 0; l < n; l++)
					if (k != l) {
						m->constraints[marginal(n, i, k,
								       0, j)]
							.coef[q] = 1.;
						m->constraints[marginal(n, i, k,
								       1, l)]
							.coef[q] = 1.;
						m->constraints[marginal(n, j, l,
								       0, i)]
							.coef[q] = 1.;
						m->constraints[marginal(n, j, l,
								       1, k)]
							.coef[q++] = 1.;
					}
	for (r = 0; r < rows; r++) {
		m->constraints[r].type = optm_CONS_T_EQ;
		m->constraints[r].rhs = r < 2 * n ? 1. : 0.;
	}
	if (permute) {
		/* Reverse rows/columns and negate alternating equations, with RHS.
		 * No names are supplied; the rule must follow coefficients only. */
		for (r = 0; r < rows / 2; r++) {
			a = rows - 1 - r;
			for (c = 0; c < cols; c++) {
				double v = m->constraints[r].coef[c];
				m->constraints[r].coef[c] =
					m->constraints[a].coef[c];
				m->constraints[a].coef[c] = v;
			}
			{
				double v = m->constraints[r].rhs;
				m->constraints[r].rhs = m->constraints[a].rhs;
				m->constraints[a].rhs = v;
			}
		}
		for (r = 0; r < rows; r++) {
			matrix = m->constraints[r].coef;
			for (c = 0; c < cols / 2; c++) {
				double v = matrix[c];
				b = cols - 1 - c;
				matrix[c] = matrix[b];
				matrix[b] = v;
			}
			if (r % 2) {
				for (c = 0; c < cols; c++)
					matrix[c] = -matrix[c];
				m->constraints[r].rhs = -m->constraints[r].rhs;
			}
		}
	}
	return m;
}

static int check(int n, int permute, int refusal)
{
	struct lp_Model *m = conservation_model(n, permute);
	double *rhs;
	unsigned char *rows, *cols;
	int *degree;
	int i, j, k, removed, good = 1;
	size_t work = 100000000UL, bytes = 32UL * 1024UL * 1024UL;
	if (!m)
		return 0;
	if (refusal == 1)
		for (i = 0; i < m->m; i++)
			m->constraints[i].rhs += 1e-9;
	if (refusal == 2)
		for (i = 0; i < m->m; i++)
			for (j = 0; j < m->n; j++)
				if (m->constraints[i].coef[j] != 0.) {
					m->constraints[i].coef[j] *= 1. + 1e-12;
					break;
				}
	/* A RHS error equal to the modular prime must still be rejected by
	 * exact verification; residues alone cannot authorize deletion. */
	if (refusal == 6)
		m->constraints[0].rhs += 65521.;
	if (lp_model_build_sparse(m)) {
		lp_model_free(m);
		return 0;
	}
	rows = (unsigned char *)malloc((size_t)m->m);
	cols = (unsigned char *)malloc((size_t)m->n);
	rhs = (double *)malloc((size_t)m->m * sizeof(double));
	degree = (int *)malloc((size_t)m->n * sizeof(int));
	if (!rows || !cols || !rhs || !degree)
		good = 0;
	if (good) {
		for (i = 0; i < m->m; i++) {
			rows[i] = 1;
			rhs[i] = m->constraints[i].rhs;
		}
		for (j = 0; j < m->n; j++) {
			cols[j] = 1;
			degree[j] = m->column_start[j + 1] - m->column_start[j];
		}
		if (refusal == 3)
			bytes = 0;
		if (refusal == 4)
			work = 0;
		if (refusal == 5)
			work = 50;
		if (refusal == 7)
			work = 8000;
		removed = simplex_presolve_dependencies(
			m, rhs, rows, cols, degree, work, bytes);
		if (refusal == 6)
			good = removed < 3 * n * n - 3 * n + 2;
		else if (refusal == 7)
			good = removed >= 0 && removed <= 3 * n * n - 3 * n + 2;
		else
			good = removed == (refusal ? 0 : 3 * n * n - 3 * n + 2);
		good = good && affine_rank(m, rows) == affine_rank(m, NULL);
		if (!refusal) {
			int sums = 0;
			for (i = 0; i < m->m; i++) {
				int positive = 0, negative = 0;
				for (j = 0; j < m->n; j++) {
					positive +=
						m->constraints[i].coef[j] > 0.;
					negative +=
						m->constraints[i].coef[j] < 0.;
				}
				if (rows[i] &&
					((positive > 0) != (negative > 0)))
					sums++;
			}
			good = good && sums == 2 * n - 1;
		}
		for (j = 0; j < m->n; j++) {
			int count = 0;
			for (k = m->column_start[j]; k < m->column_start[j + 1];
				k++)
				count += rows[m->row_index[k]] != 0;
			good = good && degree[j] == count && cols[j] == 1;
		}
		if (!good)
			fprintf(stderr,
				"dependency n=%d permutation=%d refusal=%d removed=%d\n",
				n, permute, refusal, removed);
	}
	free(rows);
	free(cols);
	free(rhs);
	free(degree);
	lp_model_free(m);
	return good;
}

static int check_opt_in(void)
{
	struct lp_Model *m = conservation_model(4, 0);
	struct simplex_Presolve normal, candidate;
	int good;
	if (!m)
		return 0;
	if (lp_model_build_sparse(m)) {
		lp_model_free(m);
		return 0;
	}
	good = !simplex_presolve_run(&normal, m, 1e-7);
	if (good) {
		good = normal.stats.dependent_rows == 0;
		simplex_presolve_destroy(&normal);
	}
	if (!simplex_presolve_run_configured(&candidate, m, 1e-7, 1)) {
		good = good && candidate.stats.dependent_rows > 0;
		simplex_presolve_destroy(&candidate);
	} else
		good = 0;
	lp_model_free(m);
	return good;
}

int main(void)
{
	int n, permutation, refusal;
	if (!check_opt_in())
		return 1;
	for (n = 3; n <= 4; n++)
		for (permutation = 0; permutation < 2; permutation++)
			for (refusal = 0; refusal < 8; refusal++)
				if (!check(n, permutation, refusal))
					return 1;
	puts("Exact dependency and refusal checks passed.");
	return 0;
}
