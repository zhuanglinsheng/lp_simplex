/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file test_infeasibility_certificate.c
 * @brief Regression tests for infeasibility certificate.
 */

#include "simplex_problem.h"
#include "simplex_pivot.h"
#include <lp_simplex/lp_simplex.h>
#include <math.h>
#include <stdio.h>

#define CHECK(condition) \
	do { \
		if (!(condition)) { \
			fprintf(stderr, \
				"certificate check failed at line %d\n", \
				__LINE__); \
			return 1; \
		} \
	} while (0)

static struct simplex_Problem *problem(void)
{
	struct simplex_Problem *p = simplex_problem_create_sparse(2, 1, 2);
	if (!p)
		return NULL;
	p->matrix.column_start[0] = 0;
	p->matrix.column_start[1] = 2;
	p->matrix.row_index[0] = 0;
	p->matrix.row_index[1] = 1;
	p->matrix.value[0] = 1.;
	p->matrix.value[1] = 1.;
	p->row_type[0] = optm_CONS_T_GE;
	p->row_type[1] = optm_CONS_T_LE;
	p->rhs[0] = 2.;
	p->rhs[1] = 1.;
	p->bounds[0].b_type = optm_BOUND_T_FR;
	p->objective[0] = -100.;
	return p;
}

/* A two-column contradiction with non-binary rational multipliers. Test
 * column and row permutations and power-of-two scaling, then make the same
 * structure feasible. No instance names or benchmark data enter the test. */
static int test_exact_transformations(void)
{
	int variant;
	for (variant = 0; variant < 16; variant++) {
		struct simplex_Problem *p =
			simplex_problem_create_sparse(3, 2, 5);
		int r0 = variant % 3, r1 = (r0 + 1) % 3, r2 = (r0 + 2) % 3;
		int first = variant % 2, second = 1 - first;
		int basis[3], j, k = 0;
		double scale[3];
		CHECK(p);
		scale[0] = ldexp(1., variant - 8);
		scale[1] = ldexp(1., 8 - variant);
		scale[2] = ldexp(1., variant % 4);
		p->row_type[r0] = optm_CONS_T_GE;
		p->row_type[r1] = optm_CONS_T_LE;
		p->row_type[r2] = optm_CONS_T_GE;
		p->rhs[r0] = 4. * scale[0];
		p->rhs[r1] = 7. * scale[1];
		p->rhs[r2] = 0.;
		for (j = 0; j < 2; j++) {
			p->bounds[j].b_type = optm_BOUND_T_FR;
			p->objective[j] = variant - j * 100.;
			p->matrix.column_start[j] = k;
			p->matrix.row_index[k] = r0;
			p->matrix.value[k++] =
				(j == first ? 1. : 3.) * scale[0];
			p->matrix.row_index[k] = r1;
			p->matrix.value[k++] =
				(j == first ? 2. : 7.1) * scale[1];
			if (j == second) {
				p->matrix.row_index[k] = r2;
				p->matrix.value[k++] = scale[2];
			}
		}
		p->matrix.column_start[2] = k;
		basis[0] = first;
		basis[1] = second;
		basis[2] = 2 + r2;
#ifdef TEST_HAVE_GMP
		CHECK(simplex_problem_verify_basis_infeasibility(p, basis, 2));
#endif
		p->rhs[r1] = 9. * scale[1];
		CHECK(!simplex_problem_verify_basis_infeasibility(p, basis, 2));
		basis[1] = basis[0];
		CHECK(!simplex_problem_verify_basis_infeasibility(p, basis, 2));
		simplex_problem_free(p);
	}
	return 0;
}

static int test_tableau_ray(void)
{
	struct lp_simplex_Options options;
	double table[6] = {1., 0., 0., 1e-20, 1., 1.};
	int basis[1] = {1}, iterations = 0;
	lp_simplex_default_options(&options, lp_simplex_ALGORITHM_TABLEAU);
	CHECK(simplex_run_pivots(
		      &iterations, table, 3, basis, 1, 2, 1, &options) == 9);
	table[3] = -1e-20;
	iterations = 0;
	CHECK(simplex_run_pivots(
		      &iterations, table, 3, basis, 1, 2, 1, &options) == 2);
	table[5] = -1.;
	iterations = 0;
	CHECK(simplex_run_pivots(
		      &iterations, table, 3, basis, 1, 2, 1, &options) == 9);
	return 0;
}

int main(void)
{
	struct simplex_Problem *p = problem();
	double y[2] = {1., -1.};
	int basis[2] = {0, 2}, i;
	CHECK(p);
	CHECK(test_tableau_ray() == 0);
	CHECK(test_exact_transformations() == 0);
	CHECK(simplex_problem_verify_infeasibility(p, y));
	/* The proof is independent of objective coefficients and orientation. */
	p->objective[0] = 1e30;
	y[0] = -1.;
	y[1] = 1.;
	CHECK(simplex_problem_verify_infeasibility(p, y));
#ifdef TEST_HAVE_GMP
	CHECK(simplex_problem_verify_basis_infeasibility(p, basis, 1));
#endif
	/* A tiny nonzero coefficient on a free variable is not zero: x can
	 * escape arbitrarily far. Do not turn a tolerance into a certificate. */
	y[0] = 1.;
	y[1] = -1.;
	p->matrix.value[1] = nextafter(1., 2.);
	CHECK(!simplex_problem_verify_infeasibility(p, y));
	p->matrix.value[1] = 1.;
	/* Feasible witnesses generated over changing scales and row directions
	 * must reject every proposed contradiction, including near boundaries. */
	for (i = 0; i < 200; i++) {
		double x = (i - 100) * .25;
		p->matrix.value[0] = (i % 7 + 1) * .125;
		p->matrix.value[1] = (i % 11 + 1) * 8.;
		p->rhs[0] = p->matrix.value[0] * x - (i % 3) * .125;
		p->rhs[1] = p->matrix.value[1] * x + (i % 5) * .25;
		y[0] = (i % 9 - 4) * .125;
		y[1] = (i % 13 - 6) * .25;
		CHECK(!simplex_problem_verify_infeasibility(p, y));
		CHECK(!simplex_problem_verify_basis_infeasibility(p, basis, 1));
	}
	/* Finite variable bounds can themselves contradict a single row. */
	p->bounds[0].b_type = optm_BOUND_T_BS;
	p->bounds[0].lb = 0.;
	p->bounds[0].ub = 1.;
	p->matrix.value[0] = 1.;
	p->rhs[0] = 2.;
	y[0] = 1.;
	y[1] = 0.;
	CHECK(simplex_problem_verify_infeasibility(p, y));
	y[0] = NAN;
	CHECK(!simplex_problem_verify_infeasibility(p, y));
	simplex_problem_free(p);
	return 0;
}
