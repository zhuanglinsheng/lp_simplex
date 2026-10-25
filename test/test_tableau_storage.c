/**
 * @file test_tableau_storage.c
 * @brief Compile the production Tableau call chain under a namespaced fault
 * allocator.
 *
 * No allocator hooks or alternate numerical path are added to the library.
 */

#include "simplex_stage.h"
#include "simplex_solve_context.h"
#include "simplex_tableau_solver.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>

static void *live[8];
static int requests, fail_at, active, bad_free;

void *storage_test_malloc(size_t bytes)
{
	void *p;
	int i;
	requests++;
	if (requests == fail_at)
		return NULL;
	p = malloc(bytes);
	if (!p)
		return NULL;
	for (i = 0; i < 8; i++)
		if (!live[i]) {
			live[i] = p;
			active++;
			return p;
		}
	free(p);
	return NULL;
}

void storage_test_free(void *p)
{
	int i;
	if (!p)
		return;
	for (i = 0; i < 8; i++)
		if (live[i] == p) {
			live[i] = NULL;
			active--;
			free(p);
			return;
		}
	bad_free++;
}

#define lp_simplex_malloc storage_test_malloc
#define lp_simplex_free storage_test_free
#define simplex_tableau_create storage_simplex_tableau_create
#define simplex_tableau_destroy storage_simplex_tableau_destroy
#define simplex_tableau_fill_constraint_types \
	storage_simplex_tableau_fill_constraint_types
#define simplex_tableau_fill_constraints \
	storage_simplex_tableau_fill_constraints
#define simplex_tableau_compute_size storage_simplex_tableau_compute_size
#define simplex_tableau_add_slack storage_simplex_tableau_add_slack
#define simplex_tableau_add_artificial storage_simplex_tableau_add_artificial
#define simplex_tableau_initialize_basis \
	storage_simplex_tableau_initialize_basis
#define simplex_tableau_initialize_phase_one_objective \
	storage_simplex_tableau_initialize_phase_one_objective
#define simplex_tableau_remove_artificial_basis \
	storage_simplex_tableau_remove_artificial_basis
#define simplex_tableau_remove_artificial_columns \
	storage_simplex_tableau_remove_artificial_columns
#define simplex_transform_alloc storage_simplex_transform_alloc
#define simplex_transform_free storage_simplex_transform_free
#define simplex_transform_size storage_simplex_transform_size
#define simplex_transform_problem storage_simplex_transform_problem
#define simplex_transform_recover storage_simplex_transform_recover
#define simplex_apply_pivot storage_simplex_apply_pivot
#define simplex_run_pivots storage_simplex_run_pivots
#define simplex_solve_standard_reported storage_simplex_solve_standard_reported
#define simplex_solve_standard storage_simplex_solve_standard
#define simplex_tableau_solve_model_reported \
	storage_simplex_tableau_solve_model_reported
#define simplex_tableau_solve_model storage_simplex_tableau_solve_model
#include "../src/tableau/simplex_tableau.c"
#include "../src/tableau/simplex_transform.c"
#include "../src/tableau/simplex_pivot.c"
#include "../src/tableau/simplex_phase.c"
#include "../src/tableau/simplex_tableau_solver.c"

#undef simplex_apply_pivot
#undef simplex_run_pivots
#undef simplex_solve_standard_reported
#undef simplex_solve_standard
#undef simplex_tableau_solve_model_reported
#undef simplex_tableau_solve_model
#undef lp_simplex_malloc
#undef lp_simplex_free
#undef simplex_tableau_create
#undef simplex_tableau_destroy
#undef simplex_tableau_fill_constraint_types
#undef simplex_tableau_fill_constraints
#undef simplex_tableau_compute_size
#undef simplex_tableau_add_slack
#undef simplex_tableau_add_artificial
#undef simplex_tableau_initialize_basis
#undef simplex_tableau_initialize_phase_one_objective
#undef simplex_tableau_remove_artificial_basis
#undef simplex_tableau_remove_artificial_columns
#undef simplex_transform_alloc
#undef simplex_transform_free
#undef simplex_transform_size
#undef simplex_transform_problem
#undef simplex_transform_recover
#include "simplex_tableau_solver.h"
#include "simplex_execution.h"

#define CHECK(c) \
	do { \
		if (!(c)) { \
			fprintf(stderr, "tableau storage %d: %s\n", __LINE__, \
				#c); \
			return 1; \
		} \
	} while (0)

static void reset(int failure)
{
	requests = active = bad_free = 0;
	fail_at = failure;
	memset(live, 0, sizeof(live));
}

static int allocations(void)
{
	double stale = 0., *table, *objective, *x, *coeff;
	int previous = 0, *basis, *types, failed, code;
	struct optm_LinearConstraint prior, *constraints;
	for (failed = 0; failed <= 3; failed++) {
		reset(failed);
		table = &stale;
		basis = types = &previous;
		code = storage_simplex_tableau_create(
			&table, &basis, &types, 2, 3, 6);
		if (failed)
			CHECK(code == lp_simplex_EXIT_FAILURE && !table &&
				!basis && !types && !active &&
				requests == failed);
		else
			CHECK(code == lp_simplex_EXIT_SUCCESS && table &&
				basis && types && active == 3 && requests == 3);
		storage_simplex_tableau_destroy(table, basis, types);
		CHECK(!active && !bad_free);
	}
	for (failed = 0; failed <= 4; failed++) {
		reset(failed);
		objective = x = coeff = &stale;
		constraints = &prior;
		code = storage_simplex_transform_alloc(
			2, 3, &objective, &x, &coeff, &constraints);
		if (failed)
			CHECK(code == lp_simplex_EXIT_FAILURE && !objective &&
				!x && !coeff && !constraints && !active &&
				requests == failed);
		else
			CHECK(code == lp_simplex_EXIT_SUCCESS && objective &&
				x && coeff && constraints && active == 4 &&
				requests == 4);
		storage_simplex_transform_free(
			objective, x, coeff, constraints);
		CHECK(!active && !bad_free);
	}
	return 0;
}

static int dimensions(void)
{
	double *table, *objective, *x, *coeff;
	int *basis, *types, m, n;
	size_t bytes;
	struct optm_LinearConstraint constraints[2], *allocated;
	struct optm_VariableBound bounds[2];
	memset(constraints, 0, sizeof(constraints));
	memset(bounds, 0, sizeof(bounds));
	reset(0);
	CHECK(simplex_tableau_storage_size(1, INT_MAX, 1, &bytes) &&
		bytes == (size_t)INT_MAX);
	CHECK(!simplex_tableau_storage_size(INT_MAX, 2, 1, &bytes) && !bytes);
	CHECK(!simplex_tableau_storage_size(2, 1, (size_t)-1, &bytes) &&
		!bytes);
	CHECK(!simplex_tableau_storage_size(-1, 1, 1, &bytes) && !bytes);
	CHECK(storage_simplex_tableau_create(&table, &basis, &types, 1, 2,
		      INT_MAX) == lp_simplex_EXIT_FAILURE &&
		!table && !basis && !types && !requests);
	CHECK(storage_simplex_tableau_create(&table, &basis, &types, 2, 2, 3) ==
			lp_simplex_EXIT_FAILURE &&
		!requests);
	CHECK(storage_simplex_transform_alloc(INT_MAX, 2, &objective, &x,
		      &coeff, &allocated) == lp_simplex_EXIT_FAILURE &&
		!objective && !x && !coeff && !allocated && !requests);
	CHECK(storage_simplex_transform_alloc(0, 2, &objective, &x, &coeff,
		      &allocated) == lp_simplex_EXIT_FAILURE &&
		!requests);
	CHECK(storage_simplex_tableau_compute_size(constraints, INT_MAX, 1, &m,
		      &n) == lp_simplex_EXIT_FAILURE &&
		!m && !n);
	CHECK(storage_simplex_tableau_compute_size(constraints, 1, INT_MAX, &m,
		      &n) == lp_simplex_EXIT_FAILURE &&
		!m && !n);
	constraints[0].type = optm_CONS_T_GE;
	CHECK(storage_simplex_tableau_compute_size(constraints, 1,
		      INT_MAX / 2 - 2, &m, &n) == lp_simplex_EXIT_FAILURE &&
		!m && !n);
	constraints[1].type = optm_CONS_T_LE;
	constraints[1].rhs = -1.;
	CHECK(storage_simplex_tableau_compute_size(constraints, 2, 3, &m, &n) ==
			lp_simplex_EXIT_SUCCESS &&
		m == 3 && n == 8);
	bounds[0].b_type = optm_BOUND_T_BS;
	CHECK(storage_simplex_transform_size(bounds, INT_MAX, 1, &m, &n) ==
			lp_simplex_EXIT_FAILURE &&
		!m && !n);
	bounds[1].b_type = optm_BOUND_T_BS;
	CHECK(storage_simplex_transform_size(bounds, INT_MAX / 2, 2, &m, &n) ==
			lp_simplex_EXIT_FAILURE &&
		!m && !n);
	bounds[1].b_type = optm_BOUND_T_FR;
	CHECK(storage_simplex_transform_size(bounds, 2, 2, &m, &n) ==
			lp_simplex_EXIT_SUCCESS &&
		m == 3 && n == 3);
	CHECK(!requests && !active && !bad_free);
	return 0;
}

static int reported_refusal(void)
{
	struct lp_Model model;
	struct optm_LinearConstraint constraint;
	struct optm_VariableBound bound;
	struct lp_simplex_Options options;
	struct simplex_StageReport report;
	struct simplex_SolveContext context;
	struct simplex_Execution history;
	struct simplex_ExecutionScope scope;
	double objective = 1., coefficient = 1., x, value;
	int kind, status, work, diagnostics;
	memset(&model, 0, sizeof(model));
	memset(&constraint, 0, sizeof(constraint));
	memset(&bound, 0, sizeof(bound));
	model.objective = &objective;
	model.constraints = &constraint;
	constraint.coef = &coefficient;
	constraint.type = optm_CONS_T_EQ;
	lp_simplex_default_options(&options, lp_simplex_ALGORITHM_TABLEAU);
	for (diagnostics = 0; diagnostics < 2; diagnostics++)
		for (kind = 0; kind < 2; kind++) {
			simplex_solve_context_init(&context,
				diagnostics ? &history : NULL, &scope);
			model.m = kind ? INT_MAX : 1;
			model.n = kind ? 1 : INT_MAX;
			bound.b_type = optm_BOUND_T_BS;
			model.bounds = kind ? &bound : NULL;
			x = value = NAN;
			work = 77;
			status = 0;
			CHECK(simplex_tableau_solve_model_reported(&model,
				      &options, &x, &value, &status, &work,
				      &report,
				      &scope) == lp_simplex_EXIT_FAILURE);
			CHECK(status == lp_simplex_MemoryAllocError &&
				work == 0 &&
				report.stop == SIMPLEX_STOP_RESOURCE);
			CHECK(!report.point_available &&
				!report.basis_available &&
				!report.factor_available &&
				!report.tableau.point_available &&
				!report.tableau.recovered_point);
			CHECK(isnan(x) && isnan(value) && context.models == 1);
			if (diagnostics)
				CHECK(history.events[0].kind ==
						SIMPLEX_EXECUTION_TABLEAU &&
					history.events[0].complete &&
					history.events[0].status ==
						lp_simplex_MemoryAllocError &&
					history.events[0].work == 0 &&
					history.events[0].stop ==
						SIMPLEX_STOP_RESOURCE &&
					!history.events[0].point_source);
			simplex_solve_context_destroy(&context);
		}
	return 0;
}

/* Exercise the full production call chain, including Phase-II scratch OOM
 * after Phase I has consumed work. No failed call can publish a partial point. */
static int solve_allocations(void)
{
	struct lp_Model *m = lp_model_create(1, 1);
	struct optm_VariableBound *bounds;
	struct simplex_SolveContext context;
	struct simplex_Execution h;
	struct simplex_ExecutionScope scope;
	struct simplex_StageReport report;
	struct lp_simplex_Options o;
	double x, value;
	int diagnostics, converted, failure, total, state, status, work, i;
	CHECK(m);
	bounds = m->bounds;
	m->constraints[0].type = optm_CONS_T_EQ;
	m->constraints[0].coef[0] = 1.;
	m->constraints[0].rhs = 2.;
	m->objective[0] = 1.;
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_TABLEAU);
	for (diagnostics = 0; diagnostics < 3; diagnostics++)
		for (converted = 0; converted < 2; converted++) {
			m->bounds = converted ? bounds : NULL;
			total = converted ? 9 : 5;
			for (failure = 0; failure <= total; failure++) {
				simplex_solve_context_init(&context,
					diagnostics ? &h : NULL, &scope);
				if (diagnostics == 2)
					while (h.count <
						SIMPLEX_EXECUTION_CAPACITY)
						(void)simplex_execution_begin(
							&scope,
							SIMPLEX_EXECUTION_SOLVE,
							0, 0, 0);
				reset(failure);
				x = value = -123.;
				state = storage_simplex_tableau_solve_model_reported(
					m, &o, &x, &value, &status, &work,
					&report, &scope);
				CHECK(!active && !bad_free);
				CHECK(requests == (failure ? failure : total));
				CHECK(report.tableau.phase_one_work +
						report.tableau.cleanup_work +
						report.tableau.phase_two_work ==
					work);
				CHECK(!report.basis_available &&
					!report.factor_available &&
					!report.dual.quality &&
					!report.lower_bound.quality);
				if (!failure) {
					CHECK(!state &&
						status == lp_simplex_Success &&
						report.point_available &&
						x == 2. && value == 2.);
					CHECK(report.tableau.recovered_point ==
						converted);
				} else {
					CHECK(state &&
						status ==
							lp_simplex_MemoryAllocError &&
						report.stop ==
							SIMPLEX_STOP_RESOURCE);
					CHECK(x == -123. && value == -123. &&
						!report.point_available &&
						!report.tableau
							.point_available &&
						!report.tableau
							.recovered_point);
					CHECK(work ==
						(failure == total ? 2 : 0));
					CHECK(report.tableau.phase ==
						(converted && failure <= 4
								? SIMPLEX_TABLEAU_INPUT
							: failure == total
								? SIMPLEX_TABLEAU_PHASE_TWO
								: SIMPLEX_TABLEAU_PHASE_ONE));
				}
				if (diagnostics == 1)
					for (i = 0; i < h.count; i++) {
						CHECK(h.events[i].complete);
						if (h.events[i].kind ==
								SIMPLEX_EXECUTION_TABLEAU ||
							h.events[i].kind ==
								SIMPLEX_EXECUTION_TABLEAU_STANDARD) {
							CHECK(h.events[i].status ==
									status &&
								h.events[i].work ==
									work &&
								h.events[i].stop ==
									report.stop);
							if (failure)
								CHECK(!h.events[i]
										.point_source);
						}
					}
				simplex_solve_context_destroy(&context);
			}
		}
	m->bounds = bounds;
	lp_model_free(m);
	return 0;
}

int main(void)
{
	return allocations() || dimensions() || reported_refusal() ||
	       solve_allocations();
}
