/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#include <lp_simplex/lp_simplex.h>
#include "simplex_presolve_dominance.h"
#include "simplex_problem.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define CHECK(x) \
	do { \
		if (!(x)) { \
			printf("dominance failed at line %d: %s\n", __LINE__, \
				#x); \
			return 0; \
		} \
	} while (0)

/* Exercise the rule independently of preceding singleton/bound propagation. */
static int apply(struct lp_Model *m, int *removed, int *fixed, size_t budget)
{
	unsigned char *rows, *columns;
	int *rd, *cd, i, k;
	double *rhs;
	CHECK(lp_model_build_sparse(m) == 0);
	rows = (unsigned char *)malloc((size_t)m->m);
	columns = (unsigned char *)malloc((size_t)m->n);
	rd = (int *)calloc((size_t)m->m, sizeof(int));
	cd = (int *)calloc((size_t)m->n, sizeof(int));
	rhs = (double *)malloc((size_t)m->m * sizeof(double));
	CHECK(rows && columns && rd && cd && rhs);
	memset(rows, 1, (size_t)m->m);
	memset(columns, 1, (size_t)m->n);
	for (i = 0; i < m->m; i++) {
		rd[i] = m->row_start[i + 1] - m->row_start[i];
		rhs[i] = m->constraints[i].rhs;
		for (k = m->row_start[i]; k < m->row_start[i + 1]; k++)
			cd[m->column_index[k]]++;
	}
	simplex_presolve_dominance(m, m->bounds, rhs, rows, columns, rd, cd,
		&budget, removed, fixed);
	free(rows);
	free(columns);
	free(rd);
	free(cd);
	free(rhs);
	return 1;
}

static struct lp_Model *cover_model(void)
{
	static const double a[3][4] = {
		{1., 2., 0., 1.}, {1., 2., 1., 0.}, {0., 1., 1., 1.}};
	static const double c[4] = {3., 2., 4., 5.};
	struct lp_Model *m = lp_model_create(3, 4);
	int i, j;
	if (!m)
		return NULL;
	for (i = 0; i < 3; i++) {
		m->constraints[i].type = optm_CONS_T_GE;
		m->constraints[i].rhs = i == 2 ? 1. : 2.;
		for (j = 0; j < 4; j++)
			m->constraints[i].coef[j] = a[i][j];
	}
	for (j = 0; j < 4; j++) {
		m->objective[j] = c[j];
		m->bounds[j].b_type = optm_BOUND_T_BS;
		m->bounds[j].ub = 3.;
	}
	return m;
}

static int bounded_rules(void)
{
	struct lp_Model *m;
	int removed, fixed, mode, i, j;
	for (mode = 0; mode < 7; mode++) {
		m = cover_model();
		CHECK(m);
		if (mode == 1)
			m->bounds[1].ub =
				.25; /* cheaper cover has insufficient capacity */
		if (mode == 2) {
			m->objective[0] = -1.;
			m->objective[1] = -2.;
		}
		if (mode == 3)
			m->bounds[0].lb = -1.;
		if (mode == 4)
			m->constraints[0].type = optm_CONS_T_EQ;
		if (mode == 5) {
			for (i = 0; i < 3; i++) {
				m->constraints[i].type = optm_CONS_T_LE;
				m->constraints[i].rhs = -m->constraints[i].rhs;
				for (j = 0; j < 4; j++)
					m->constraints[i].coef[j] =
						-m->constraints[i].coef[j];
			}
		}
		CHECK(apply(m, &removed, &fixed, mode == 6 ? 0 : 100000));
		if (mode == 0 || mode == 5)
			CHECK(fixed > 0 && m->bounds[0].ub == 0.);
		else
			CHECK(m->bounds[0].ub > m->bounds[0].lb);
		lp_model_free(m);
	}
	/* Equality of a rounded product must not mask a capacity shortfall:
	 * 1.25 * nextafter(.8,0) is strictly below one. */
	m = cover_model();
	CHECK(m);
	for (i = 0; i < 3; i++) {
		m->constraints[i].rhs = 1.;
		m->constraints[i].coef[1] = 1.25;
	}
	m->bounds[1].ub = nextafter(.8, 0.);
	CHECK(apply(m, &removed, &fixed, 100000));
	CHECK(m->bounds[0].ub > 0.);
	lp_model_free(m);
	return 1;
}

static int general_rows(void)
{
	struct lp_Model *m;
	int mode, removed, fixed, i;
	for (mode = 0; mode < 4; mode++) {
		m = lp_model_create(2, 3);
		CHECK(m);
		m->constraints[0].type = optm_CONS_T_GE;
		m->constraints[0].rhs = 3.;
		m->constraints[0].coef[0] = 2.;
		m->constraints[0].coef[1] = 1.;
		m->constraints[1].type = optm_CONS_T_LE;
		m->constraints[1].rhs = -2.;
		m->constraints[1].coef[0] = -3.;
		m->constraints[1].coef[1] = -2.;
		m->constraints[1].coef[2] = -1.;
		for (i = 0; i < 3; i++) {
			m->objective[i] = 1.;
			m->bounds[i].b_type = optm_BOUND_T_BS;
			m->bounds[i].ub = 5.;
		}
		if (mode == 1)
			m->bounds[2].lb = -10.;
		if (mode == 2)
			m->constraints[1].coef[0] = -1.;
		if (mode == 3)
			m->constraints[1].rhs = -4.;
		CHECK(apply(m, &removed, &fixed, 10000));
		CHECK(removed == (mode == 0 ? 1 : 0));
		lp_model_free(m);
	}
	return 1;
}

static int unbounded_target(void)
{
	struct lp_Model *m;
	struct lp_simplex_Options options;
	struct lp_simplex_Result result;
	double x[3];
	int removed, fixed, mode;
	for (mode = 0; mode < 2; mode++) {
		m = lp_model_create(2, 3);
		CHECK(m);
		m->objective[0] = 3.;
		m->objective[1] = 1.;
		m->objective[2] = 4.;
		m->bounds[0].lb = -2.;
		m->constraints[0].type = optm_CONS_T_EQ;
		m->constraints[0].rhs = 3.;
		m->constraints[0].coef[0] = 1.;
		m->constraints[0].coef[1] = mode ? 2. : 1.;
		m->constraints[0].coef[2] = 1.;
		m->constraints[1].type = optm_CONS_T_GE;
		m->constraints[1].rhs = -3.;
		m->constraints[1].coef[0] = -2.;
		m->constraints[1].coef[1] = -1.;
		m->constraints[1].coef[2] = 1.;
		CHECK(apply(m, &removed, &fixed, 10000));
		if (!mode)
			CHECK(m->bounds[0].b_type == optm_BOUND_T_BS &&
				m->bounds[0].ub == -2.);
		else
			CHECK(m->bounds[0].b_type == optm_BOUND_T_LO);
		/* Restore the original bound and test the actual presolve/postsolve
		 * path, including a nonzero fixed value and objective offset. */
		m->bounds[0].b_type = optm_BOUND_T_LO;
		m->bounds[0].ub = HUGE_VAL;
		m->objective_offset = 3.25;
		lp_simplex_default_options(
			&options, lp_simplex_ALGORITHM_DUAL_REVISED);
		CHECK(lp_simplex_solve(m, &options, x, &result) == 0);
		CHECK(fabs(x[0] + 2.) < 1e-7);
		CHECK(fabs(result.objective - (mode ? -.25 : 2.25)) < 1e-7);
		CHECK(result.primal_infeasibility <= options.primal_tolerance);
		lp_model_free(m);
	}
	return 1;
}

/* Compare complete solves and original-model residuals, including postsolve.
 * Continuous random weighted covers are not tied to a named benchmark. */
static int fallback_and_storage_guards(void)
{
	struct lp_Model *small = cover_model(), *m;
	unsigned char *rows, *columns;
	int *rd, *cd, removed, fixed, i, k;
	double *rhs;
	size_t work = 100000;
	CHECK(small && lp_model_build_sparse(small) == 0);
	/* This shape exceeds the 32 MiB bit-index cap but stores only 9 entries.
	 * It exercises sparse posting-list fallback without allocating m*n. */
	m = lp_model_create_sparse(4097, 65537, 9);
	CHECK(m);
	m->nnz = 9;
	for (i = 0; i < m->n; i++) {
		m->column_start[i] = i < 4 ? small->column_start[i] : 9;
		if (i < 4) {
			m->bounds[i] = small->bounds[i];
			m->objective[i] = small->objective[i];
		}
	}
	m->column_start[m->n] = 9;
	for (i = 0; i < m->m; i++) {
		m->row_start[i] = i < 3 ? small->row_start[i] : 9;
		if (i < 3) {
			m->constraints[i].type = small->constraints[i].type;
			m->constraints[i].rhs = small->constraints[i].rhs;
		}
	}
	m->row_start[m->m] = 9;
	for (k = 0; k < 9; k++) {
		m->row_index[k] = small->row_index[k];
		m->value[k] = small->value[k];
		m->column_index[k] = small->column_index[k];
		m->row_value[k] = small->row_value[k];
	}
	rows = (unsigned char *)calloc((size_t)m->m, 1);
	columns = (unsigned char *)calloc((size_t)m->n, 1);
	rd = (int *)calloc((size_t)m->m, sizeof(int));
	cd = (int *)calloc((size_t)m->n, sizeof(int));
	rhs = (double *)calloc((size_t)m->m, sizeof(double));
	CHECK(rows && columns && rd && cd && rhs);
	for (i = 0; i < 3; i++) {
		rows[i] = 1;
		rd[i] = 3;
		rhs[i] = m->constraints[i].rhs;
	}
	for (i = 0; i < 4; i++) {
		columns[i] = 1;
		cd[i] = m->column_start[i + 1] - m->column_start[i];
	}
	simplex_presolve_dominance(m, m->bounds, rhs, rows, columns, rd, cd,
		&work, &removed, &fixed);
	CHECK(fixed > 0 && m->bounds[0].ub == 0.);
	/* An unsorted view is skipped, not misread as a valid merge ordering. */
	for (i = 0; i < 4; i++)
		m->bounds[i] = small->bounds[i];
	i = m->row_index[0];
	m->row_index[0] = m->row_index[1];
	m->row_index[1] = i;
	work = 100000;
	simplex_presolve_dominance(m, m->bounds, rhs, rows, columns, rd, cd,
		&work, &removed, &fixed);
	CHECK(!fixed && !removed && m->bounds[0].ub == 3.);
	/* Rows 34 and 104 collide in both 32-bit and 64-bit fingerprints.
	 * A hash hit must never replace the exact coefficient proof. */
	memset(rows, 0, (size_t)m->m);
	memset(columns, 0, (size_t)m->n);
	memset(rd, 0, (size_t)m->m * sizeof(int));
	memset(cd, 0, (size_t)m->n * sizeof(int));
	memset(rhs, 0, (size_t)m->m * sizeof(double));
	memset(m->row_start, 0, ((size_t)m->m + 1) * sizeof(int));
	m->nnz = 6;
	for (i = 0; i < m->n; i++)
		m->column_start[i] = i < 3 ? 2 * i : 6;
	m->column_start[m->n] = 6;
	m->row_index[0] = 0;
	m->row_index[1] = 34;
	m->row_index[2] = 0;
	m->row_index[3] = 104;
	m->row_index[4] = 34;
	m->row_index[5] = 104;
	for (k = 0; k < 6; k++) {
		m->value[k] = m->row_value[k] = 1.;
		m->row_start[m->row_index[k] + 1]++;
	}
	for (i = 0; i < m->m; i++)
		m->row_start[i + 1] += m->row_start[i];
	m->column_index[0] = 0;
	m->column_index[1] = 1;
	m->column_index[2] = 0;
	m->column_index[3] = 2;
	m->column_index[4] = 1;
	m->column_index[5] = 2;
	rows[0] = rows[34] = rows[104] = 1;
	rd[0] = rd[34] = rd[104] = 2;
	rhs[0] = rhs[34] = rhs[104] = 1.;
	m->constraints[0].type = m->constraints[34].type =
		m->constraints[104].type = optm_CONS_T_GE;
	for (i = 0; i < 3; i++) {
		columns[i] = 1;
		cd[i] = 2;
		m->bounds[i].lb = 0.;
		m->bounds[i].ub = 1.;
		m->bounds[i].b_type = optm_BOUND_T_BS;
	}
	m->objective[0] = 3.;
	m->objective[1] = 2.;
	m->objective[2] = 100.;
	work = 100000;
	simplex_presolve_dominance(m, m->bounds, rhs, rows, columns, rd, cd,
		&work, &removed, &fixed);
	CHECK(!fixed && !removed && m->bounds[0].ub == 1.);
	/* Signed unbounded dominance requires only restrictive rows in the
	 * fingerprint: the target may omit a negative GE coefficient. */
	m->constraints[0].type = optm_CONS_T_EQ;
	m->value[1] = m->row_value[2] = -2.;
	m->objective[1] = 1.;
	for (i = 0; i < 3; i++) {
		m->bounds[i].b_type = optm_BOUND_T_LO;
		m->bounds[i].ub = HUGE_VAL;
	}
	m->bounds[0].lb = -2.;
	work = 100000;
	simplex_presolve_dominance(m, m->bounds, rhs, rows, columns, rd, cd,
		&work, &removed, &fixed);
	CHECK(fixed && m->bounds[0].ub == -2.);
	free(rows);
	free(columns);
	free(rd);
	free(cd);
	free(rhs);
	lp_model_free(m);
	lp_model_free(small);
	m = lp_model_create(1, 2);
	CHECK(m);
	m->constraints[0].type = optm_CONS_T_GE;
	m->constraints[0].rhs = .75;
	for (i = 0; i < 2; i++) {
		m->constraints[0].coef[i] = 1.;
		m->objective[i] = i + 1.;
		m->bounds[i].b_type = optm_BOUND_T_BS;
		m->bounds[i].ub = .5;
	}
	CHECK(apply(m, &removed, &fixed, 10000));
	CHECK(!fixed); /* identical columns are both necessary because of capacity */
	lp_model_free(m);
	return 1;
}

static int equivalence(void)
{
	unsigned int rng = 17;
	int trial, i, j;
	for (trial = 0; trial < 60; trial++) {
		struct lp_Model *m = lp_model_create(4, 8);
		struct lp_simplex_Options o;
		struct lp_simplex_Result a, b;
		double x[8], y[8];
		CHECK(m);
		m->objective_offset = 3.25;
		for (j = 0; j < 8; j++) {
			rng = 1664525U * rng + 1013904223U;
			m->objective[j] = 1. + (rng % 7);
			m->bounds[j].b_type = optm_BOUND_T_BS;
			m->bounds[j].ub = 2.;
		}
		for (i = 0; i < 4; i++) {
			double sum = 0.;
			m->constraints[i].type =
				trial % 2 ? optm_CONS_T_LE : optm_CONS_T_GE;
			for (j = 0; j < 8; j++) {
				double value;
				rng = 1664525U * rng + 1013904223U;
				value = (double)((rng >> 8) % 4);
				sum += value;
				m->constraints[i].coef[j] =
					trial % 2 ? -value : value;
			}
			m->constraints[i].rhs =
				trial % 2 ? -.3 * sum : .3 * sum;
		}
		/* Force a reusable weighted dominance opportunity. */
		m->objective[0] = 7.;
		m->objective[1] = 1.;
		for (i = 0; i < 4; i++)
			m->constraints[i].coef[0] = m->constraints[i].coef[1];
		CHECK(lp_model_build_sparse(m) == 0);
		lp_simplex_default_options(
			&o, lp_simplex_ALGORITHM_DUAL_REVISED);
		o.presolve = 0;
		CHECK(lp_simplex_solve(m, &o, x, &a) == 0);
		o.presolve = 1;
		CHECK(lp_simplex_solve(m, &o, y, &b) == 0);
		CHECK(fabs(a.objective - b.objective) <=
			1e-6 * (1. + fabs(a.objective)));
		CHECK(b.primal_infeasibility <= o.primal_tolerance);
		lp_model_free(m);
	}
	return 1;
}

int main(void)
{
	if (!bounded_rules() || !general_rows() || !unbounded_target() ||
		!fallback_and_storage_guards() || !equivalence())
		return 1;
	puts("Weighted dominance, bound/sign/equality guards and 60 original-model equivalence checks passed");
	return 0;
}
