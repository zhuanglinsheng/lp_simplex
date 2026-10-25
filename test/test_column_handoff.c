/**
 * @file test_column_handoff.c
 * @brief Checked column restrictions distinguish model relation from LP
 * certificates.
 */

#include "simplex_solve_context.h"
#include "simplex_execution.h"
#include "simplex_column_handoff.h"
#include "simplex_dual_internal.h"
#include <lp_simplex/lp_simplex.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#define CHECK(c) \
	do { \
		if (!(c)) { \
			fprintf(stderr, "column handoff line %d: %s\n", \
				__LINE__, #c); \
			return 1; \
		} \
	} while (0)

static int model(struct simplex_Problem *p, const int *map, int n)
{
	struct lp_Model *m = lp_model_create(2, n);
	int i, j, code;
	if (!m)
		return -1;
	for (i = 0; i < 2; i++) {
		m->constraints[i].type = optm_CONS_T_GE;
		m->constraints[i].rhs = 1.;
	}
	for (j = 0; j < n; j++) {
		i = map ? map[j] : j;
		m->constraints[i % 2].coef[j] = 1.;
		m->objective[j] = i < 2 ? 2. : 1.;
		m->bounds[j].b_type = optm_BOUND_T_BS;
		m->bounds[j].lb = 0.;
		m->bounds[j].ub = 1.;
	}
	code = simplex_problem_from_model(p, m);
	lp_model_free(m);
	return code;
}

static int transfers(int tracing)
{
	struct simplex_Problem p, q, other, drop;
	struct simplex_SolveContext context;
	struct simplex_Execution journal;
	struct simplex_ColumnProjection v, full, dropped;
	struct simplex_ColumnBasis basis;
	struct simplex_DualBasisSeed seed;
	struct lp_simplex_Options o;
	struct lp_simplex_Result r;
	int map[2] = {1, 0}, reverse[4] = {1, 0, -1, -1}, all[4] = {0, 1, 2, 3},
	    dropmap[2] = {2, 3}, droprev[4] = {-1, -1, 0, 1};
	int cursor[2], index[2], global[2], saved[2], i;
	unsigned char status[6], global_status[6], saved_status[6], scratch[1];
	double x[4], y[2], staging[4];
	CHECK(!model(&p, NULL, 4) && !model(&q, map, 2));
	simplex_solve_context_init(
		&context, tracing ? &journal : NULL, &p.execution);
	q.execution = simplex_solve_context_derive(
		&p.execution, SIMPLEX_EXECUTION_COLUMN_MODEL);
	CHECK(!simplex_column_projection_bind(
		&v, &p, &q, map, reverse, cursor, 2));
	simplex_dual_seed_init(&seed, &q, index, status);
	simplex_column_basis_init(&basis, &p, global, global_status);
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_DUAL_REVISED);
	CHECK(!simplex_dual_solve_problem_seeded(&q, &o, x, y, &r, 0, &seed));
	CHECK(seed.valid && r.objective == 4.);
	CHECK(!simplex_column_point_lift(&v, x, staging, o.primal_tolerance));
	CHECK(staging[0] == 1. && staging[1] == 1. && staging[2] == 0. &&
		staging[3] == 0.);
	CHECK(!simplex_column_rows_lift(&v, y));
	/* Same row coordinates do not establish dual feasibility on omitted columns. */
	CHECK(simplex_problem_dual_infeasibility(&p, staging, y,
		      o.primal_tolerance, -1) > o.dual_tolerance);
	CHECK(!simplex_column_basis_lift(
		&v, &seed, &basis, scratch, sizeof(scratch)));
	memcpy(saved, index, sizeof(index));
	memcpy(saved_status, status, 4);
	/* Same rows, a new column dimension: a bare reused seed is not a remapping. */
	CHECK(!simplex_dual_solve_problem_seeded(&p, &o, x, y, &r, 0, &seed));
	CHECK(!seed.valid && !seed.accepted && !seed.import_attempted &&
		seed.rejection == SIMPLEX_BASIS_MODEL_REJECTED);
	CHECK(!memcmp(saved, index, sizeof(index)) &&
		!memcmp(saved_status, status, 4));
	CHECK(r.objective == 2.);
	CHECK(!model(&drop, dropmap, 2));
	drop.execution = simplex_solve_context_derive(
		&p.execution, SIMPLEX_EXECUTION_COLUMN_MODEL);
	CHECK(!simplex_column_projection_bind(
		&dropped, &p, &drop, dropmap, droprev, cursor, 2));
	CHECK(simplex_column_basis_restrict(&dropped, &basis, &seed, 1, scratch,
		      sizeof(scratch)) == SIMPLEX_COLUMN_LAYOUT);
	CHECK(!memcmp(saved, index, sizeof(index)) &&
		!memcmp(saved_status, status, 4));
	simplex_problem_destroy(&drop);
	/* The explicit parent-coordinate carrier survives destruction of the master. */
	v.valid = 0;
	simplex_problem_destroy(&q);
	CHECK(!simplex_column_projection_bind(
		&full, &p, &p, all, all, cursor, 2));
	CHECK(!simplex_column_basis_restrict(
		&full, &basis, &seed, 1, scratch, sizeof(scratch)));
	CHECK(simplex_problem_identity_matches(&seed.identity, &p));
	CHECK(!simplex_dual_solve_problem_seeded(&p, &o, x, y, &r, 0, &seed));
	CHECK(seed.accepted && seed.valid && r.objective == 2.);
	CHECK(!model(&other, NULL, 4));
	other.execution = p.execution;
	memcpy(saved, index, sizeof(index));
	memcpy(saved_status, status, sizeof(status));
	CHECK(!simplex_dual_solve_problem_seeded(
		&other, &o, x, y, &r, 0, &seed));
	CHECK(!seed.accepted && !seed.import_attempted &&
		seed.rejection == SIMPLEX_BASIS_MODEL_REJECTED);
	CHECK(!memcmp(saved, index, sizeof(index)) &&
		!memcmp(saved_status, status, sizeof(status)));
	if (tracing)
		for (i = 0; i < journal.count; i++) {
			const struct simplex_ExecutionEvent *e =
				journal.events + i;
			CHECK(e->complete && e->model > 0 && e->parent <= i);
			if (e->kind == SIMPLEX_EXECUTION_COLUMN_HANDOFF)
				CHECK(!e->result_available);
		}
	simplex_problem_destroy(&other);
	simplex_problem_destroy(&p);
	simplex_solve_context_destroy(&context);
	return 0;
}

static int refusals(void)
{
	struct simplex_Problem p, q;
	struct simplex_SolveContext context, other;
	struct simplex_Execution journal;
	struct simplex_ColumnProjection v, old;
	struct simplex_ColumnBasis basis;
	struct simplex_DualBasisSeed seed;
	struct lp_simplex_Options o;
	struct lp_simplex_Result r;
	int map[2] = {0, 1}, reverse[4] = {0, 1, -1, -1}, cursor[2], index[2],
	    global[2], saved[2], i;
	unsigned char status[4], global_status[6], saved_status[6], scratch[1];
	double x[4], y[2], stage[4], *saved_cost;
	CHECK(!model(&p, NULL, 4) && !model(&q, map, 2));
	simplex_solve_context_init(&context, &journal, &p.execution);
	q.execution = simplex_solve_context_derive(
		&p.execution, SIMPLEX_EXECUTION_COLUMN_MODEL);
	CHECK(!simplex_column_projection_bind(
		&v, &p, &q, map, reverse, cursor, 2));
	old = v;
	map[1] = 0;
	CHECK(simplex_column_projection_bind(&v, &p, &q, map, reverse, cursor,
		      2) == SIMPLEX_COLUMN_MAP);
	map[1] = 1;
	CHECK(!memcmp(&v, &old, sizeof(v)));
	reverse[3] = 0;
	CHECK(simplex_column_projection_bind(&v, &p, &q, map, reverse, cursor,
		      2) == SIMPLEX_COLUMN_MAP);
	reverse[3] = -1;
	for (i = 0; i < 6; i++) {
		if (i == 0)
			q.matrix.value[0] += 1.;
		if (i == 1)
			q.matrix.row_value[0] += 1.;
		if (i == 2)
			q.rhs[0] += 1.;
		if (i == 3)
			q.objective[0] += 1.;
		if (i == 4)
			q.bounds[0].ub += 1.;
		if (i == 5)
			p.bounds[3].lb = .1;
		CHECK(simplex_column_projection_bind(&v, &p, &q, map, reverse,
			      cursor, 2) == SIMPLEX_COLUMN_RELATION);
		if (i == 0)
			q.matrix.value[0] -= 1.;
		if (i == 1)
			q.matrix.row_value[0] -= 1.;
		if (i == 2)
			q.rhs[0] -= 1.;
		if (i == 3)
			q.objective[0] -= 1.;
		if (i == 4)
			q.bounds[0].ub -= 1.;
		if (i == 5)
			p.bounds[3].lb = 0.;
	}
	q.execution.context = &other;
	CHECK(simplex_column_projection_bind(&v, &p, &q, map, reverse, cursor,
		      2) == SIMPLEX_COLUMN_MODEL);
	q.execution.context = &context;
	q.revision = ULONG_MAX;
	CHECK(simplex_column_projection_bind(&v, &p, &q, map, reverse, cursor,
		      2) == SIMPLEX_COLUMN_MODEL);
	q.revision = 0;
	CHECK(simplex_column_projection_bind(&v, &p, &q, map, reverse, cursor,
		      1) == SIMPLEX_COLUMN_STORAGE);
	CHECK(!simplex_column_projection_bind(
		&v, &p, &q, map, reverse, cursor, 2));
	simplex_dual_seed_init(&seed, &q, index, status);
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_DUAL_REVISED);
	CHECK(!simplex_dual_solve_problem_seeded(&q, &o, x, y, &r, 0, &seed));
	simplex_column_basis_init(&basis, &p, global, global_status);
	memset(global, 0x5a, sizeof(global));
	memset(global_status, 0x5a, sizeof(global_status));
	memcpy(saved, global, sizeof(global));
	memcpy(saved_status, global_status, sizeof(global_status));
	CHECK(simplex_column_basis_lift(&v, &seed, &basis, scratch, 0) ==
		SIMPLEX_COLUMN_STORAGE);
	CHECK(!basis.valid && !memcmp(saved, global, sizeof(global)) &&
		!memcmp(saved_status, global_status, sizeof(global_status)));
	index[1] = index[0];
	CHECK(simplex_column_basis_lift(&v, &seed, &basis, scratch,
		      sizeof(scratch)) == SIMPLEX_COLUMN_LAYOUT);
	CHECK(!basis.valid && !memcmp(saved, global, sizeof(global)) &&
		!memcmp(saved_status, global_status, sizeof(global_status)));
	for (i = 0; i < 4; i++)
		stage[i] = 17.;
	x[0] = NAN;
	CHECK(simplex_column_point_lift(&v, x, stage, o.primal_tolerance) ==
		SIMPLEX_COLUMN_POINT);
	CHECK(stage[0] == 17.);
	y[0] = NAN;
	CHECK(simplex_column_rows_lift(&v, y) == SIMPLEX_COLUMN_MULTIPLIER);
	saved_cost = p.objective;
	p.objective = stage;
	CHECK(simplex_column_rows_lift(&v, y) == SIMPLEX_COLUMN_MODEL);
	p.objective = saved_cost;
	simplex_problem_mark_changed(&q);
	CHECK(simplex_column_point_lift(&v, x, stage, o.primal_tolerance) ==
		SIMPLEX_COLUMN_MODEL);
	CHECK(stage[0] == 17.);
	simplex_problem_destroy(&q);
	simplex_problem_destroy(&p);
	simplex_solve_context_destroy(&context);
	return 0;
}

static int integer_boundary(void)
{
	struct simplex_Problem p, q;
	struct simplex_ColumnProjection v;
	int starts[2] = {0, INT_MAX}, rows[1] = {1},
	    row_starts[3] = {0, INT_MAX, INT_MAX}, columns[1] = {0};
	int map[1] = {0}, cursor[2];
	double value[1] = {1.}, cost[1] = {0.}, rhs[2] = {0., 0.};
	unsigned char types[2] = {optm_CONS_T_GE, optm_CONS_T_GE};
	struct optm_VariableBound bounds[1];
	memset(&p, 0, sizeof(p));
	memset(bounds, 0, sizeof(bounds));
	bounds[0].b_type = optm_BOUND_T_LO;
	p.rows = p.matrix.rows = 2;
	p.columns = p.matrix.columns = 1;
	p.matrix.nonzeros = INT_MAX;
	p.matrix.column_start = starts;
	p.matrix.row_index = rows;
	p.matrix.value = p.matrix.row_value = value;
	p.matrix.row_start = row_starts;
	p.matrix.column_index = columns;
	p.objective = cost;
	p.rhs = rhs;
	p.row_type = types;
	p.bounds = bounds;
	q = p;
	/* Metadata is intentionally inconsistent. The first CSC entry targets a
	 * declared empty row whose cursor equals INT_MAX; reject before increment
	 * or access at that offset, without allocating huge metadata-sized arrays. */
	CHECK(simplex_column_projection_bind(&v, &p, &q, map, map, cursor, 2) ==
		SIMPLEX_COLUMN_RELATION);
	CHECK(cursor[1] == INT_MAX);
	return 0;
}

/* Numerical basis state outlives its diagnostic producer. Slot reuse must
 * not turn an unrelated event into the producer of that retained state. */
static int producer_lifetime(void)
{
	struct simplex_Problem p, q;
	struct simplex_SolveContext context;
	struct simplex_Execution h, other;
	struct simplex_ColumnProjection v;
	struct simplex_ColumnBasis basis;
	struct simplex_DualBasisSeed seed;
	struct lp_simplex_Options o;
	struct lp_simplex_Result r;
	int map[2] = {0, 1}, reverse[4] = {0, 1, -1, -1};
	int index[2], global[2], cursor[2], mode, lift, i, id, cause;
	unsigned char status[4], global_status[6], scratch[1];
	double x[4], y[2];
	for (mode = 0; mode < 6; mode++) {
		CHECK(!model(&p, NULL, 4) && !model(&q, map, 2));
		simplex_solve_context_init(&context, &h, &p.execution);
		simplex_solve_context_configure(&context, 0);
		q.execution = simplex_solve_context_derive(
			&p.execution, SIMPLEX_EXECUTION_COLUMN_MODEL);
		CHECK(!simplex_column_projection_bind(
			&v, &p, &q, map, reverse, cursor, 2));
		simplex_dual_seed_init(&seed, &q, index, status);
		simplex_column_basis_init(&basis, &p, global, global_status);
		lp_simplex_default_options(
			&o, lp_simplex_ALGORITHM_DUAL_REVISED);
		CHECK(!simplex_dual_solve_problem_seeded(
			&q, &o, x, y, &r, 0, &seed));
		if (mode == 3)
			context.execution_serial = ULONG_MAX;
		CHECK(!simplex_column_basis_lift(
			&v, &seed, &basis, scratch, sizeof(scratch)));
		lift = h.count;
		CHECK(h.events[lift - 1].kind ==
				SIMPLEX_EXECUTION_COLUMN_HANDOFF &&
			h.events[lift - 1].detail == SIMPLEX_COLUMN_BASIS_LIFT);
		if (mode == 1) {
			memset(&h, 0, sizeof(h));
			for (i = 0; i < lift; i++) {
				id = simplex_execution_begin(&p.execution,
					SIMPLEX_EXECUTION_COLUMN_HANDOFF,
					SIMPLEX_COLUMN_BASIS_LIFT, 0, 0);
				h.events[id - 1].complete = 1;
			}
		}
		if (mode == 2)
			h.events[lift - 1].complete = 0;
		if (mode == 4) {
			memset(&other, 0, sizeof(other));
			q.execution.history = &other;
		}
		if (mode == 5)
			while (h.count < SIMPLEX_EXECUTION_CAPACITY)
				(void)simplex_execution_begin(&p.execution,
					SIMPLEX_EXECUTION_SOLVE, 0, 0, 0);
		/* A surviving source is deliberately not the immediately prior action. */
		(void)simplex_execution_begin(
			&q.execution, SIMPLEX_EXECUTION_SOLVE, 0, 0, 0);
		CHECK(!simplex_column_basis_restrict(
			&v, &basis, &seed, 1, scratch, sizeof(scratch)));
		CHECK(basis.valid && seed.valid && seed.prefer_dual &&
			simplex_problem_identity_matches(&seed.identity, &q));
		if (mode != 5) {
			struct simplex_Execution *journal = q.execution.history;
			cause = journal->events[journal->count - 1].cause;
			if (cause != (mode == 0 ? lift : 0))
				fprintf(stderr,
					"column producer mode=%d stale-cause=%d lift=%d\n",
					mode, cause, lift);
			CHECK(cause == (mode == 0 ? lift : 0));
		}
		CHECK(!simplex_dual_solve_problem_seeded(
			&q, &o, x, y, &r, 0, &seed));
		CHECK(seed.accepted && r.status == lp_simplex_Success &&
			r.objective == 4.);
		simplex_problem_destroy(&q);
		simplex_problem_destroy(&p);
		simplex_solve_context_destroy(&context);
	}
	return 0;
}

int main(void)
{
	if (transfers(0) || transfers(1) || refusals() || integer_boundary() ||
		producer_lifetime())
		return 1;
	puts("column projections audit identities, exact relations, basis mappings, points and row candidates");
	return 0;
}
