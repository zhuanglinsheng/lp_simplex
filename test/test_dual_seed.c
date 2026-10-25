/**
 * @file test_dual_seed.c
 * @brief Regression tests for dual seed.
 */

#include "simplex_dual.h"
#include "simplex_column_handoff.h"
#include "simplex_dual_internal.h"
#include <math.h>
#include <stdio.h>
#include <string.h>

#define CHECK(c) \
	do { \
		if (!(c)) { \
			printf("seed test line %d\n", __LINE__); \
			return 0; \
		} \
	} while (0)

static struct simplex_Problem *cover(int columns, int scaled)
{
	struct simplex_Problem *p =
		simplex_problem_create_sparse(2, columns, columns);
	int i, j, k = 0;
	if (!p)
		return NULL;
	for (i = 0; i < 2; i++) {
		p->rhs[i] = scaled ? (i == 0 ? 1e6 : 1e-6) : 1.;
		p->row_type[i] = optm_CONS_T_GE;
		p->matrix.row_start[i] = k;
		for (j = i; j < columns; j += 2) {
			p->matrix.column_index[k] = j;
			p->matrix.row_value[k++] = p->rhs[i];
		}
	}
	p->matrix.row_start[2] = k;
	for (j = 0; j < columns; j++) {
		p->objective[j] = j < 2 ? 2. : 1.;
		p->bounds[j].b_type = optm_BOUND_T_BS;
		p->bounds[j].v_type = optm_VAR_T_REAL;
		p->bounds[j].lb = 0.;
		p->bounds[j].ub = 1.;
		p->matrix.column_start[j] = j;
		p->matrix.row_index[j] = j % 2;
		p->matrix.value[j] = p->rhs[j % 2];
	}
	p->matrix.column_start[columns] = columns;
	return p;
}

static int appended_columns(int scaled, int budget, int mode)
{
	struct simplex_Problem *old = cover(2, scaled), *p = cover(4, scaled);
	struct simplex_DualBasisSeed seed;
	struct simplex_ColumnProjection projection;
	struct simplex_ColumnBasis parent_basis;
	int global_index[2], cursor[2], map[2] = {0, 1},
					reverse[4] = {0, 1, -1, -1},
					full[4] = {0, 1, 2, 3};
	unsigned char global_status[6], scratch[1];
	struct lp_simplex_Options o;
	struct lp_simplex_Result r;
	int index[2], code;
	unsigned char status[6];
	double x[4], y[2];

	CHECK(old && p);
	if (mode == 2)
		p->bounds[2].b_type = p->bounds[3].b_type = optm_BOUND_T_LO;
	simplex_dual_seed_init(&seed, old, index, status);
	seed.prefer_dual = mode != 0;
	simplex_column_basis_init(
		&parent_basis, p, global_index, global_status);
	CHECK(!simplex_column_projection_bind(
		&projection, p, old, map, reverse, cursor, 2));
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_DUAL_REVISED);
	CHECK(simplex_dual_solve_problem_seeded(old, &o, x, y, &r, 0, &seed) ==
		0);
	CHECK(seed.valid && !seed.accepted && fabs(r.objective - 4.) < 1e-10);
	CHECK(!simplex_column_basis_lift(
		&projection, &seed, &parent_basis, scratch, sizeof(scratch)));
	CHECK(!simplex_column_projection_bind(
		&projection, p, p, full, full, cursor, 2));
	CHECK(!simplex_column_basis_restrict(&projection, &parent_basis, &seed,
		mode != 0, scratch, sizeof(scratch)));
	if (budget)
		o.iteration_limit = 1;
	code = simplex_dual_solve_problem_seeded(p, &o, x, y, &r, 0, &seed);
	CHECK(seed.accepted);
	if (budget) {
		CHECK(code != 0 && r.status == lp_simplex_ExceedIterLimit);
		CHECK(!seed.valid && r.iterations <= 1);
	} else {
		CHECK(code == 0 && seed.valid &&
			fabs(r.objective - 2.) < 1e-10);
		CHECK(simplex_problem_primal_infeasibility(p, x) <=
			o.primal_tolerance);
		CHECK(simplex_problem_dual_infeasibility(p, x, y,
			      o.primal_tolerance, -1) <= o.dual_tolerance);
		if (mode == 1)
			CHECK(r.iterations == 0);
		if (!mode) {
			/* The same nonsingular basis is infeasible after changing the RHS.
			 * Reject it and rebuild every cached index from a cold basis. */
			simplex_problem_mark_changed(p);
			p->rhs[0] *= 1.5;
			CHECK(simplex_dual_solve_problem_seeded(
				      p, &o, x, y, &r, 0, &seed) == 0);
			CHECK(!seed.accepted && !seed.valid &&
				fabs(r.objective - 3.) < 1e-10);
			CHECK(!seed.import_attempted &&
				seed.rejection == SIMPLEX_BASIS_MODEL_REJECTED);
			simplex_dual_seed_init(&seed, p, index, status);
			CHECK(!simplex_dual_solve_problem_seeded(
				p, &o, x, y, &r, 0, &seed));
			CHECK(seed.valid);
			/* Distinct ids can still describe a singular basis. */
			index[0] = 0;
			index[1] = 2;
			memset(status, DUAL_STATUS_LOWER, sizeof(status));
			status[0] = status[2] = DUAL_STATUS_BASIC;
			CHECK(simplex_dual_solve_problem_seeded(
				      p, &o, x, y, &r, 0, &seed) == 0);
			CHECK(!seed.accepted && seed.valid &&
				fabs(r.objective - 3.) < 1e-10);
			CHECK(seed.import_attempted &&
				seed.rejection ==
					SIMPLEX_BASIS_FACTORIZATION_REJECTED);
			/* A duplicate basis is a rejected hint, not an infeasibility claim. */
			index[1] = index[0];
			CHECK(simplex_dual_solve_problem_seeded(
				      p, &o, x, y, &r, 0, &seed) == 0);
			CHECK(!seed.accepted && seed.valid &&
				fabs(r.objective - 3.) < 1e-10);
			CHECK(seed.import_attempted &&
				seed.rejection ==
					SIMPLEX_BASIS_LAYOUT_REJECTED);
			/* Wrong dimensions must neither consume pivots nor overwrite arrays. */
			seed.columns = 3;
			CHECK(simplex_dual_solve_problem_seeded(
				      p, &o, x, y, &r, 0, &seed) != 0);
			CHECK(!seed.valid && !seed.accepted &&
				r.iterations == 0);
		}
	}
	simplex_problem_free(old);
	simplex_problem_free(p);
	return 1;
}

int main(void)
{
	if (!appended_columns(0, 0, 0) || !appended_columns(1, 0, 0) ||
		!appended_columns(0, 1, 0) || !appended_columns(0, 0, 1) ||
		!appended_columns(1, 0, 1) || !appended_columns(0, 0, 2))
		return 1;
	puts("basis seeds preserve feasibility, remap logical ids and respect budgets");
	return 0;
}
