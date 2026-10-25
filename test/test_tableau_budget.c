/**
 * @file test_tableau_budget.c
 * @brief Artificial-basis removal is numerical pivot work within the same
 * lease.
 */

#include "simplex_tableau.h"
#include <lp_simplex/lp_simplex.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#define CHECK(c) \
	do { \
		if (!(c)) { \
			fprintf(stderr, "tableau budget %d: %s\n", __LINE__, \
				#c); \
			return 1; \
		} \
	} while (0)

static void fixture(double *table, int *basis, int redundant)
{
	int i;
	memset(table, 0, 28 * sizeof(double));
	for (i = 0; i < 3; i++) {
		table[(i + 1) * 7 + i] = redundant ? 0. : -1.;
		table[(i + 1) * 7 + 3 + i] = 1.;
		basis[i] = 3 + i;
	}
}

int main(void)
{
	double table[28], saved[28], x[1];
	int basis[3], limit, work, code, i, pricing, budget;
	struct lp_Model *m = lp_model_create(1, 1);
	struct lp_simplex_Options o;
	struct lp_simplex_Result r;
	CHECK(m);
	for (limit = 0; limit <= 4; limit++) {
		fixture(table, basis, 0);
		memcpy(saved, table, sizeof(table));
		work = 0;
		code = simplex_tableau_remove_artificial_basis(
			table, 7, basis, 3, 3, 6, 1e-9, &work, limit);
		CHECK(code == (limit < 3 ? SIMPLEX_TABLEAU_CLEANUP_LIMIT
					 : SIMPLEX_TABLEAU_CLEANUP_DONE));
		CHECK(work == (limit < 3 ? limit : 3));
		for (i = 0; i < 3; i++)
			CHECK(basis[i] == (i < work ? i : 3 + i));
		if (!limit)
			CHECK(!memcmp(saved, table, sizeof(table)));
	}
	fixture(table, basis, 0);
	work = INT_MAX - 1;
	CHECK(simplex_tableau_remove_artificial_basis(table, 7, basis, 3, 3, 6,
		      1e-9, &work, INT_MAX) == SIMPLEX_TABLEAU_CLEANUP_LIMIT);
	CHECK(work == INT_MAX && basis[0] == 0 && basis[1] == 4);
	fixture(table, basis, 1);
	work = 2;
	CHECK(simplex_tableau_remove_artificial_basis(table, 7, basis, 3, 3, 6,
		      1e-9, &work, 2) == SIMPLEX_TABLEAU_CLEANUP_DONE);
	CHECK(work == 2 && basis[0] == 3 && basis[1] == 4);
	fixture(table, basis, 0);
	table[2 * 7 + 1] = NAN;
	work = 0;
	CHECK(simplex_tableau_remove_artificial_basis(table, 7, basis, 3, 3, 6,
		      1e-9, &work, 3) == SIMPLEX_TABLEAU_CLEANUP_NUMERICAL);
	CHECK(work == 1 && basis[0] == 3 && basis[1] == 4);
	m->constraints[0].type = optm_CONS_T_EQ;
	m->constraints[0].coef[0] = -1.;
	m->objective[0] = 1.;
	for (pricing = 0; pricing < 2; pricing++)
		for (budget = 1; budget <= 4; budget++) {
			lp_simplex_default_options(
				&o, lp_simplex_ALGORITHM_TABLEAU);
			o.presolve = 0;
			o.iteration_limit = budget;
			o.pricing = pricing ? lp_simplex_PRICING_DANTZIG
					    : lp_simplex_PRICING_BLAND;
			code = lp_simplex_solve(m, &o, x, &r);
			CHECK(r.iterations == (budget < 3 ? budget : 3));
			if (budget < 3)
				CHECK(code != 0 &&
					r.status ==
						lp_simplex_ExceedIterLimit &&
					!isfinite(x[0]));
			else
				CHECK(code == 0 &&
					r.status == lp_simplex_Success &&
					x[0] == 0. && r.objective == 0.);
		}
	lp_model_free(m);
	return 0;
}
