/**
 * @file test_checkpoint_fault.c
 * @brief Compile Dual boundaries with test-local startup and update failure seams.
 *
 * Production allocation and imported helper libraries remain unchanged.
 */

#include "simplex_stage.h"
#include "simplex_solve_context.h"
#include "simplex_execution.h"
#include <stdlib.h>

static int reject_allocation, reject_state, reject_scaling, reject_bridge;
static int reject_update, packed_updates, ordinary_updates;
static int reject_ratio, ratio_calls;
static int reported_update_resource, reported_factor_resource;
static int factor_calls, reject_factor_at, numerical_factor_at;
static int force_refactor, forced_refactors;

void *checkpoint_test_malloc(size_t n)
{
	return reject_allocation ? NULL : malloc(n);
}

#define lp_simplex_malloc checkpoint_test_malloc
#define simplex_dual_state_create checkpoint_test_state_create
#define simplex_scaling_create_problem checkpoint_test_scaling_create_problem
#define simplex_progress_bridge checkpoint_test_progress_bridge
#define simplex_dual_kernel_solve checkpoint_test_kernel_solve
#define simplex_dual_backend_solve checkpoint_test_backend_solve
#define simplex_basis_update checkpoint_test_basis_update
#define simplex_basis_update_packed checkpoint_test_basis_update_packed
#define simplex_basis_storage_failed checkpoint_test_storage_failed
#define simplex_basis_factorize checkpoint_test_factorize
#define simplex_dual_next_ratio_candidate checkpoint_test_next_ratio_candidate
#ifndef CHECKPOINT_TEST_SOURCE
#define CHECKPOINT_TEST_SOURCE "../src/dual/simplex_dual.c"
#include "../src/dual/simplex_dual_kernel.c"
#include "../src/dual/simplex_dual_backend.c"
#endif
#include CHECKPOINT_TEST_SOURCE

#undef lp_simplex_malloc
#undef simplex_dual_state_create
#undef simplex_scaling_create_problem
#undef simplex_progress_bridge
#undef simplex_dual_kernel_solve
#undef simplex_dual_backend_solve
#undef simplex_basis_update
#undef simplex_basis_update_packed
#undef simplex_basis_storage_failed
#undef simplex_basis_factorize
#undef simplex_dual_next_ratio_candidate
#include <stdio.h>
#include <string.h>

int simplex_basis_update(struct simplex_Basis *, int, const double *);
int simplex_basis_update_packed(
	struct simplex_Basis *, int, const double *, const int *, int);
int simplex_basis_storage_failed(const struct simplex_Basis *);
int simplex_basis_factorize(struct simplex_Basis *);

int checkpoint_test_factorize(struct simplex_Basis *basis)
{
	factor_calls++;
	if (factor_calls == reject_factor_at) {
		reported_factor_resource = 1;
		return lp_simplex_EXIT_FAILURE;
	}
	if (factor_calls == numerical_factor_at)
		return lp_simplex_EXIT_FAILURE;
	return simplex_basis_factorize(basis);
}

int simplex_dual_next_ratio_candidate(
	struct simplex_DualState *, double *, int, int, int);

int checkpoint_test_next_ratio_candidate(struct simplex_DualState *state,
	double *theta, int perturbation, int leaving, int kappa)
{
	ratio_calls++;
	if (reject_ratio) {
		state->ratio_failure_status = lp_simplex_MemoryAllocError;
		return SIMPLEX_DUAL_RATIO_ERROR;
	}
	return simplex_dual_next_ratio_candidate(
		state, theta, perturbation, leaving, kappa);
}

int checkpoint_test_basis_update(
	struct simplex_Basis *basis, int position, const double *direction)
{
	ordinary_updates++;
	if (reject_update)
		reported_update_resource = 1;
	if (force_refactor) {
		forced_refactors++;
		return 1;
	}
	return reject_update ? lp_simplex_EXIT_FAILURE
			     : simplex_basis_update(basis, position, direction);
}

int checkpoint_test_basis_update_packed(struct simplex_Basis *basis,
	int position, const double *direction, const int *index, int count)
{
	packed_updates++;
	if (reject_update)
		reported_update_resource = 1;
	if (force_refactor) {
		forced_refactors++;
		return 1;
	}
	return reject_update ? lp_simplex_EXIT_FAILURE
			     : simplex_basis_update_packed(basis, position,
				       direction, index, count);
}

int checkpoint_test_storage_failed(const struct simplex_Basis *basis)
{
	return reported_update_resource || reported_factor_resource ||
	       simplex_basis_storage_failed(basis);
}

int simplex_dual_state_create(struct simplex_DualState *,
	const struct simplex_Problem *, const struct lp_simplex_Options *, int);

int checkpoint_test_state_create(struct simplex_DualState *state,
	const struct simplex_Problem *p, const struct lp_simplex_Options *o,
	int backend)
{
	if (reject_state) {
		memset(state, 0, sizeof(*state));
		return lp_simplex_EXIT_FAILURE;
	}
	return simplex_dual_state_create(state, p, o, backend);
}

/* Return-category seams test the entry's resource exit, not the delegated
 * scaler/bridge allocators themselves. Normal calls use their real bodies. */
struct simplex_Problem *simplex_scaling_create_problem(
	const struct simplex_Problem *, struct simplex_Scaling *);
int simplex_progress_bridge(struct simplex_ProgressBridge *,
	const struct simplex_ProgressSink *, const struct simplex_Presolve *,
	const struct simplex_Scaling *, int, int,
	const struct simplex_Problem *, const struct simplex_ExecutionScope *,
	const struct simplex_Problem *);

struct simplex_Problem *checkpoint_test_scaling_create_problem(
	const struct simplex_Problem *p, struct simplex_Scaling *scaling)
{
	return reject_scaling ? NULL
			      : simplex_scaling_create_problem(p, scaling);
}

int checkpoint_test_progress_bridge(struct simplex_ProgressBridge *bridge,
	const struct simplex_ProgressSink *parent,
	const struct simplex_Presolve *presolve,
	const struct simplex_Scaling *scaling, int rows, int columns,
	const struct simplex_Problem *input,
	const struct simplex_ExecutionScope *target,
	const struct simplex_Problem *output)
{
	if (reject_bridge)
		return SIMPLEX_PROGRESS_BRIDGE_MEMORY;
	return simplex_progress_bridge(bridge, parent, presolve, scaling, rows,
		columns, input, target, output);
}

#define CHECK(c) \
	do { \
		if (!(c)) { \
			fprintf(stderr, "checkpoint fault line %d: %s\n", \
				__LINE__, #c); \
			return 1; \
		} \
	} while (0)

static int scenario(int tracing, int where)
{
	struct lp_Model *m = lp_model_create(1, 2);
	struct simplex_Problem p;
	struct simplex_DualCheckpoint cp;
	struct simplex_SolveContext context;
	struct simplex_Execution h;
	struct simplex_StageReport report;
	struct lp_simplex_Options o;
	struct lp_simplex_Result r;
	int index[1], saved_index[1], count;
	unsigned char status[3], saved_status[3];
	double x[2];
#ifndef OLD
	int producer;
#endif
	CHECK(m);
	m->constraints[0].type = optm_CONS_T_EQ;
	m->constraints[0].rhs = 1.;
	m->constraints[0].coef[0] = m->constraints[0].coef[1] = 1.;
	m->objective[1] = 1.;
	CHECK(!simplex_problem_from_model(&p, m));
	simplex_solve_context_init(&context, tracing ? &h : NULL, &p.execution);
	p.report = &report;
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_DUAL_REVISED);
	simplex_dual_checkpoint_init(&cp, &p, index, status);
	o.iteration_limit = 0;
	CHECK(simplex_dual_solve_problem_checkpoint(
		      &p, &o, x, NULL, &r, 0, &cp) != 0 &&
		cp.basis.valid);
	o.iteration_limit = 100;
	CHECK(!simplex_dual_solve_problem_checkpoint(
		&p, &o, x, NULL, &r, 0, &cp));
	CHECK(cp.basis.accepted && cp.basis.import_attempted);
#ifndef OLD
	CHECK(cp.basis.exported);
	producer = simplex_dual_seed_source_event(&cp.basis, &p.execution);
	CHECK((producer > 0) == tracing);
#endif
	memcpy(saved_index, index, sizeof(index));
	memcpy(saved_status, status, sizeof(status));
	count = tracing ? h.count : 0;
	reject_allocation = where == 0;
	reject_state = where == 1;
	reject_scaling = where == 2;
	reject_bridge = where == 3;
	factor_calls = 0;
	reject_factor_at = where == 4 ? 1 : 0;
	reported_factor_resource = 0;
	CHECK(simplex_dual_solve_problem_checkpoint(
		      &p, &o, x, NULL, &r, 0, &cp) != 0);
	reject_allocation = reject_state = reject_scaling = reject_bridge = 0;
	reject_factor_at = reported_factor_resource = 0;
	printf("{\"where\":%d,\"trace\":%d,\"status\":%d,\"work\":%d,\"valid\":%d,\"accepted\":%d,\"import_attempted\":%d}\n",
		where, tracing, r.status, r.iterations, cp.basis.valid,
		cp.basis.accepted, cp.basis.import_attempted);
	CHECK(r.status == lp_simplex_MemoryAllocError && r.iterations == 0 &&
		report.stop == SIMPLEX_STOP_RESOURCE);
	CHECK(!memcmp(index, saved_index, sizeof(index)) &&
		!memcmp(status, saved_status, sizeof(status)));
	/* State creation and checked-import failures enter a backend here. */
	if (tracing)
		CHECK(h.count == count + (where == 1 || where == 4));
#ifndef OLD
	CHECK(cp.basis.valid && !cp.basis.exported && !cp.basis.accepted &&
		cp.basis.import_attempted == (where == 4) &&
		!report.basis_available);
	if (where == 4)
		CHECK(cp.basis.rejection == SIMPLEX_BASIS_RESOURCE_REJECTED &&
			!simplex_dual_seed_retry_allowed(&cp.basis));
	CHECK(simplex_dual_seed_source_event(&cp.basis, &p.execution) ==
		producer);
#endif
	CHECK(!simplex_dual_solve_problem_checkpoint(
		&p, &o, x, NULL, &r, 0, &cp));
#ifndef OLD
	CHECK(cp.basis.exported && cp.basis.accepted && report.basis_available);
	if (tracing)
		CHECK(simplex_dual_seed_source_event(&cp.basis, &p.execution) >
			producer);
#endif
	simplex_solve_context_destroy(&context);
	simplex_problem_destroy(&p);
	lp_model_free(m);
	return 0;
}

static int kernel_storage_failure(int phase, int tracing)
{
	struct lp_Model *m = lp_model_create(1, 1);
	struct simplex_Problem p;
	struct simplex_SolveContext context;
	struct simplex_Execution history;
	struct simplex_StageReport report;
	struct simplex_DualBasisSeed seed;
	struct lp_simplex_Options options;
	struct lp_simplex_Result result;
	int index = 1;
	unsigned char status[] = {DUAL_STATUS_LOWER, DUAL_STATUS_BASIC};
	double x;
	CHECK(m);
	m->constraints[0].type =
		phase == 1 || phase == 3 ? optm_CONS_T_GE : optm_CONS_T_LE;
	m->constraints[0].rhs = 1.;
	m->constraints[0].coef[0] = 1.;
	m->objective[0] = phase == 1 || phase == 3 ? 1. : -1.;
	CHECK(!simplex_problem_from_model(&p, m));
	simplex_solve_context_init(
		&context, tracing ? &history : NULL, &p.execution);
	p.report = &report;
	lp_simplex_default_options(&options, lp_simplex_ALGORITHM_DUAL_REVISED);
	simplex_dual_seed_init(&seed, &p, &index, status);
	seed.valid = phase == 2;
	seed.prefer_dual = 0;
	ordinary_updates = packed_updates = ratio_calls = 0;
	reported_update_resource = reported_factor_resource = 0;
	reject_ratio = phase == 3;
	reject_update = phase != 3;
	CHECK(simplex_dual_solve_problem_seeded(
		&p, &options, &x, NULL, &result, 0, &seed));
	reject_update = reject_ratio = 0;
	printf("storage phase=%d trace=%d status=%d work=%d ordinary=%d packed=%d\n",
		phase, tracing, result.status, result.iterations,
		ordinary_updates, packed_updates);
	CHECK(result.status == lp_simplex_MemoryAllocError &&
		report.stop == SIMPLEX_STOP_RESOURCE);
	CHECK(!context.retries.count);
	if (tracing)
		CHECK(history.count == 1);
	CHECK(ordinary_updates == (phase == 0 || phase == 2) &&
		packed_updates == (phase == 1));
	if (phase == 3)
		CHECK(ratio_calls == 1);
	CHECK(result.iterations == (phase == 1));
	CHECK(!seed.exported && !report.basis_available);
	CHECK((seed.accepted != 0) == (phase == 2));
	reported_update_resource = reported_factor_resource = 0;
	CHECK(!simplex_dual_solve_problem_seeded(
		&p, &options, &x, NULL, &result, 0, &seed));
	CHECK(result.status == lp_simplex_Success && x == 1.);
	simplex_solve_context_destroy(&context);
	simplex_problem_destroy(&p);
	lp_model_free(m);
	return 0;
}

static int factor_failure(int mode, int tracing)
{
	struct lp_Model *m = lp_model_create(1, 1);
	struct simplex_Problem p;
	struct simplex_SolveContext context;
	struct simplex_Execution history;
	struct simplex_StageReport report;
	struct simplex_DualBasisSeed seed;
	struct lp_simplex_Options options;
	struct lp_simplex_Result result;
	int index = 1, saved_index = 1;
	unsigned char status[] = {DUAL_STATUS_LOWER, DUAL_STATUS_BASIC};
	unsigned char saved_status[2];
	double x;
	int cleanup = mode == 4 || mode == 5;
	CHECK(m);
	m->constraints[0].type = cleanup ? optm_CONS_T_LE : optm_CONS_T_GE;
	m->constraints[0].rhs = m->constraints[0].coef[0] = 1.;
	m->objective[0] = cleanup ? -1. : 1.;
	CHECK(!simplex_problem_from_model(&p, m));
	simplex_solve_context_init(
		&context, tracing ? &history : NULL, &p.execution);
	p.report = &report;
	lp_simplex_default_options(&options, lp_simplex_ALGORITHM_DUAL_REVISED);
	simplex_dual_seed_init(&seed, &p, &index, status);
	seed.valid = mode == 1 || mode == 5;
	seed.prefer_dual = mode == 1;
	memcpy(saved_status, status, sizeof(status));
	factor_calls = forced_refactors = 0;
	reported_factor_resource = reported_update_resource = 0;
	force_refactor = mode >= 3;
	reject_factor_at = mode == 3 || mode == 5 ? 3 : mode >= 2 ? 2 : 1;
	numerical_factor_at = mode == 2 ? 1 : 0;
	CHECK(simplex_dual_solve_problem_seeded(
		&p, &options, &x, NULL, &result, 0, &seed));
	printf("factor mode=%d trace=%d calls=%d status=%d work=%d rejection=%d\n",
		mode, tracing, factor_calls, result.status, result.iterations,
		seed.rejection);
	CHECK(factor_calls == reject_factor_at);
	CHECK(result.status == lp_simplex_MemoryAllocError &&
		report.stop == SIMPLEX_STOP_RESOURCE && !context.retries.count);
	CHECK(result.iterations == (mode == 3));
	CHECK(!report.basis_available && !seed.exported);
	CHECK(index == saved_index &&
		!memcmp(status, saved_status, sizeof(status)));
	CHECK(forced_refactors == (mode >= 3));
	if (mode == 1) {
		CHECK(seed.valid && seed.import_attempted && !seed.accepted &&
			seed.rejection == SIMPLEX_BASIS_RESOURCE_REJECTED &&
			!simplex_dual_seed_retry_allowed(&seed));
	}
	if (tracing) {
		CHECK(history.count == 1 + (mode == 2));
		if (mode == 2)
			CHECK(history.events[1].status ==
					lp_simplex_MemoryAllocError &&
				history.events[1].stop ==
					SIMPLEX_STOP_RESOURCE);
	}
	reject_factor_at = numerical_factor_at = force_refactor = 0;
	reported_factor_resource = reported_update_resource = 0;
	CHECK(!simplex_dual_solve_problem_seeded(
		&p, &options, &x, NULL, &result, 0, &seed));
	CHECK(result.status == lp_simplex_Success && x == 1.);
	simplex_solve_context_destroy(&context);
	simplex_problem_destroy(&p);
	lp_model_free(m);
	return 0;
}

int main(void)
{
	int tracing, where;
	for (tracing = 0; tracing < 2; tracing++)
		for (where = 0; where < 5; where++)
			if (scenario(tracing, where))
				return 1;
	for (tracing = 0; tracing < 2; tracing++)
		for (where = 0; where < 4; where++)
			if (kernel_storage_failure(where, tracing))
				return 1;
	for (tracing = 0; tracing < 2; tracing++)
		for (where = 0; where < 6; where++)
			if (factor_failure(where, tracing))
				return 1;
	return 0;
}
