/**
 * @file test_dual_storage.c
 * @brief Production Dual owned storage under a bounded fault allocator.
 */

#include "simplex_stage.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#define DUAL_STORAGE_SLOTS 512
static void *live[DUAL_STORAGE_SLOTS];
static size_t live_bytes[DUAL_STORAGE_SLOTS];
static int requests, fail_at, active, bad_free;
static size_t largest;
#ifdef LP_SIMPLEX_HAVE_KLU
static int klu_failure;
#define klu_analyze dual_storage_klu_analyze
#define klu_factor dual_storage_klu_factor
#endif

void *dual_storage_test_malloc(size_t bytes)
{
	void *p;
	int i;
	requests++;
	if (bytes > largest)
		largest = bytes;
	if (requests == fail_at || bytes > 1024 * 1024)
		return NULL;
	p = malloc(bytes ? bytes : 1);
	if (!p)
		return NULL;
	for (i = 0; i < DUAL_STORAGE_SLOTS; i++)
		if (!live[i]) {
			live[i] = p;
			live_bytes[i] = bytes;
			active++;
			return p;
		}
	free(p);
	return NULL;
}

void dual_storage_test_free(void *p)
{
	int i;
	if (!p)
		return;
	for (i = 0; i < DUAL_STORAGE_SLOTS; i++)
		if (live[i] == p) {
			live[i] = NULL;
			active--;
			free(p);
			return;
		}
	bad_free++;
}

void *dual_storage_test_realloc(void *p, size_t bytes)
{
	void *replacement;
	int i;
	if (!p)
		return dual_storage_test_malloc(bytes);
	for (i = 0; i < DUAL_STORAGE_SLOTS; i++)
		if (live[i] == p) {
			replacement = dual_storage_test_malloc(bytes);
			if (!replacement)
				return NULL;
			memcpy(replacement, p,
				bytes < live_bytes[i] ? bytes : live_bytes[i]);
			dual_storage_test_free(p);
			return replacement;
		}
	bad_free++;
	return NULL;
}

#define lp_simplex_malloc dual_storage_test_malloc
#define lp_simplex_realloc dual_storage_test_realloc
#define lp_simplex_free dual_storage_test_free
#define simplex_dual_solve_problem storage_simplex_dual_solve_problem
#define simplex_dual_solve_problem_seeded \
	storage_simplex_dual_solve_problem_seeded
#define simplex_dual_solve_problem_checkpoint \
	storage_simplex_dual_solve_problem_checkpoint
#define simplex_dual_find_feasible_reported \
	storage_simplex_dual_find_feasible_reported
#define simplex_dual_find_feasible storage_simplex_dual_find_feasible
#define simplex_dual_solve storage_simplex_dual_solve
#define simplex_dual_state_destroy storage_simplex_dual_state_destroy
#define simplex_dual_state_create storage_simplex_dual_state_create
#define simplex_dual_packed_work_prepare \
	storage_simplex_dual_packed_work_prepare
#define simplex_dual_pricing_create storage_simplex_dual_pricing_create
#define simplex_dual_pricing_destroy storage_simplex_dual_pricing_destroy
#define simplex_dual_nonbasic_add storage_simplex_dual_nonbasic_add
#define simplex_dual_nonbasic_remove storage_simplex_dual_nonbasic_remove
#define simplex_dual_nonbasic_initialize \
	storage_simplex_dual_nonbasic_initialize
#define simplex_dual_nonbasic_is_consistent \
	storage_simplex_dual_nonbasic_is_consistent
#define simplex_dual_ratio_remove storage_simplex_dual_ratio_remove
#define simplex_dual_prepare_ratio_test storage_simplex_dual_prepare_ratio_test
#define simplex_dual_next_ratio_candidate \
	storage_simplex_dual_next_ratio_candidate
#define simplex_dual_feasibility_create storage_simplex_dual_feasibility_create
#define simplex_dual_feasibility_destroy \
	storage_simplex_dual_feasibility_destroy
#define simplex_dual_feasibility_update storage_simplex_dual_feasibility_update
#define simplex_dual_feasibility_rebuild \
	storage_simplex_dual_feasibility_rebuild
#define simplex_dual_feasibility_update_packed \
	storage_simplex_dual_feasibility_update_packed
#define simplex_dual_feasibility_choose storage_simplex_dual_feasibility_choose
#define simplex_dual_feasibility_merit storage_simplex_dual_feasibility_merit
#define simplex_sparse_vector_create storage_simplex_sparse_vector_create
#define simplex_sparse_vector_destroy storage_simplex_sparse_vector_destroy
#define simplex_sparse_vector_clear storage_simplex_sparse_vector_clear
#define simplex_sparse_vector_set storage_simplex_sparse_vector_set
#define simplex_sparse_vector_add storage_simplex_sparse_vector_add
#define simplex_sparse_vector_pack storage_simplex_sparse_vector_pack
#define simplex_sparse_vector_pack_maximum \
	storage_simplex_sparse_vector_pack_maximum
#define simplex_sparse_vector_get storage_simplex_sparse_vector_get
#define simplex_sparse_vector_sort storage_simplex_sparse_vector_sort
#define simplex_sparse_vector_compact storage_simplex_sparse_vector_compact
#define simplex_basis_set_backend storage_simplex_basis_set_backend
#define simplex_basis_resolved_backend storage_simplex_basis_resolved_backend
#define simplex_basis_create storage_simplex_basis_create
#define simplex_basis_create_configured storage_simplex_basis_create_configured
#define simplex_basis_destroy storage_simplex_basis_destroy
#define simplex_basis_factorize storage_simplex_basis_factorize
#define simplex_basis_ftran_refined storage_simplex_basis_ftran_refined
#define simplex_basis_btran_refined storage_simplex_basis_btran_refined
#define simplex_basis_ftran_pair storage_simplex_basis_ftran_pair
#define simplex_basis_ftran storage_simplex_basis_ftran
#define simplex_basis_ftran_packed storage_simplex_basis_ftran_packed
#define simplex_basis_btran storage_simplex_basis_btran
#define simplex_basis_btran_packed storage_simplex_basis_btran_packed
#define simplex_basis_update storage_simplex_basis_update
#define simplex_basis_update_packed storage_simplex_basis_update_packed
#define simplex_basis_set_sparse_eta storage_simplex_basis_set_sparse_eta
#define simplex_basis_use_product_updates \
	storage_simplex_basis_use_product_updates
#define simplex_basis_update_count storage_simplex_basis_update_count
#define simplex_basis_storage_failed storage_simplex_basis_storage_failed
#define simplex_basis_compact_active storage_simplex_basis_compact_active
#define simplex_basis_compact_ever_active \
	storage_simplex_basis_compact_ever_active
#define simplex_basis_compact_requested storage_simplex_basis_compact_requested
#define simplex_basis_request_compact storage_simplex_basis_request_compact
#define simplex_basis_get_profile storage_simplex_basis_get_profile
#define simplex_sparse_lu_create storage_lu_create
#define simplex_sparse_lu_destroy storage_lu_destroy
#define simplex_sparse_lu_factorize storage_lu_factorize
#define simplex_sparse_lu_factorize_submatrix storage_lu_factorize_submatrix
#define simplex_sparse_lu_solve storage_lu_solve
#define simplex_sparse_lu_solve_pair storage_lu_solve_pair
#define simplex_sparse_lu_ft_update storage_lu_ft_update
#define simplex_sparse_lu_ft_active storage_lu_ft_active
#ifndef DUAL_STORAGE_SIMPLEX_DUAL_SOURCE
#define DUAL_STORAGE_SIMPLEX_DUAL_SOURCE "../src/dual/simplex_dual.c"
#endif
#include DUAL_STORAGE_SIMPLEX_DUAL_SOURCE
#ifndef DUAL_STORAGE_SIMPLEX_DUAL_STATE_SOURCE
#define DUAL_STORAGE_SIMPLEX_DUAL_STATE_SOURCE \
	"../src/dual/simplex_dual_state.c"
#endif
#include DUAL_STORAGE_SIMPLEX_DUAL_STATE_SOURCE
#ifndef DUAL_STORAGE_SIMPLEX_DUAL_PRICING_SOURCE
#define DUAL_STORAGE_SIMPLEX_DUAL_PRICING_SOURCE \
	"../src/dual/simplex_dual_pricing.c"
#endif
#include DUAL_STORAGE_SIMPLEX_DUAL_PRICING_SOURCE
#ifndef DUAL_STORAGE_SIMPLEX_DUAL_FEASIBILITY_SOURCE
#define DUAL_STORAGE_SIMPLEX_DUAL_FEASIBILITY_SOURCE \
	"../src/dual/simplex_dual_feasibility.c"
#endif
#include DUAL_STORAGE_SIMPLEX_DUAL_FEASIBILITY_SOURCE
#ifndef DUAL_STORAGE_SIMPLEX_SPARSE_VECTOR_SOURCE
#define DUAL_STORAGE_SIMPLEX_SPARSE_VECTOR_SOURCE \
	"../src/matrix/simplex_sparse_vector.c"
#endif
#include DUAL_STORAGE_SIMPLEX_SPARSE_VECTOR_SOURCE
#ifndef DUAL_STORAGE_SIMPLEX_BASIS_SOURCE
#define DUAL_STORAGE_SIMPLEX_BASIS_SOURCE "../src/basis/simplex_basis.c"
#endif
#include DUAL_STORAGE_SIMPLEX_BASIS_SOURCE
#include "../src/basis/simplex_sparse_lu.c"
#undef lp_simplex_malloc
#undef lp_simplex_realloc
#undef lp_simplex_free
#ifdef LP_SIMPLEX_HAVE_KLU
#undef klu_analyze
#undef klu_factor
klu_symbolic *klu_analyze(int32_t, int32_t *, int32_t *, klu_common *);
klu_numeric *klu_factor(
	int32_t *, int32_t *, double *, klu_symbolic *, klu_common *);

klu_symbolic *dual_storage_klu_analyze(
	int32_t n, int32_t *start, int32_t *index, klu_common *common)
{
	if (klu_failure >= 1 && klu_failure <= 3) {
		common->status = klu_failure == 1   ? KLU_OUT_OF_MEMORY
				 : klu_failure == 2 ? KLU_TOO_LARGE
						    : KLU_INVALID;
		return NULL;
	}
	return klu_analyze(n, start, index, common);
}

klu_numeric *dual_storage_klu_factor(int32_t *start, int32_t *index,
	double *value, klu_symbolic *symbolic, klu_common *common)
{
	if (klu_failure >= 4) {
		common->status = klu_failure == 4   ? KLU_OUT_OF_MEMORY
				 : klu_failure == 5 ? KLU_TOO_LARGE
						    : KLU_SINGULAR;
		return NULL;
	}
	return klu_factor(start, index, value, symbolic, common);
}
#endif
#define CHECK(c) \
	do { \
		if (!(c)) { \
			fprintf(stderr, "Dual storage line %d: %s\n", \
				__LINE__, #c); \
			return 1; \
		} \
	} while (0)

static void problem_shape(struct simplex_Problem *p, int rows, int columns)
{
	memset(p, 0, sizeof(*p));
	p->rows = p->matrix.rows = rows;
	p->columns = p->matrix.columns = columns;
}

static int extremes(const char *which)
{
	struct simplex_Problem p;
	struct simplex_DualState state;
	struct simplex_DualFeasibility heap;
	struct simplex_Basis basis;
	struct lp_simplex_Options options;
	struct lp_simplex_Result result;
	struct simplex_StageReport report;
	double x = 123.;
	int index = 0;
	requests = fail_at = 0;
	lp_simplex_default_options(&options, lp_simplex_ALGORITHM_DUAL_REVISED);
	if (!strcmp(which, "state")) {
		problem_shape(&p, 1, INT_MAX);
		CHECK(storage_simplex_dual_state_create(
			&state, &p, &options, 0));
		CHECK(!requests);
		storage_simplex_dual_state_destroy(&state);
	} else if (!strcmp(which, "frequency")) {
		problem_shape(&p, INT_MAX, 0);
		CHECK(storage_simplex_dual_state_create(
			&state, &p, &options, 0));
		storage_simplex_dual_state_destroy(&state);
	} else if (!strcmp(which, "heap")) {
		CHECK(storage_simplex_dual_feasibility_create(&heap, INT_MAX));
		storage_simplex_dual_feasibility_destroy(&heap);
	} else if (!strcmp(which, "heap-leaf")) {
		/* A leaf must return before computing 2*position+1 or reading rows. */
		feasibility_sift_down(NULL, NULL, NULL, INT_MAX, INT_MAX - 1);
	} else if (!strcmp(which, "basis")) {
		problem_shape(&p, INT_MAX, 0);
		CHECK(storage_simplex_basis_create_configured(
			&basis, &p.matrix, 0, &index, 0));
		storage_simplex_basis_destroy(&basis);
	} else if (!strcmp(which, "entry")) {
		problem_shape(&p, INT_MAX, 1);
		p.report = &report;
		CHECK(storage_simplex_dual_solve_problem(
			&p, &options, &x, NULL, &result, 0));
		CHECK(!requests && x == 123. && !result.iterations &&
			result.status == lp_simplex_MemoryAllocError);
		CHECK(report.stop == SIMPLEX_STOP_RESOURCE);
	} else if (!strcmp(which, "eta-growth") ||
		   !strcmp(which, "dense-growth") ||
		   !strcmp(which, "eta-count")) {
		struct simplex_BasisImpl impl;
		int start = INT_MAX;
		double direction[16] = {1.};
		memset(&impl, 0, sizeof(impl));
		basis.impl = &impl;
		impl.rows = 16;
		if (!strcmp(which, "eta-growth")) {
			impl.eta_capacity = INT_MAX - 1;
			CHECK(simplex_basis_reserve_eta(&basis, INT_MAX));
		} else if (!strcmp(which, "dense-growth")) {
			impl.dense_eta_capacity = INT_MAX - 1;
			CHECK(simplex_basis_reserve_dense_eta(&basis, INT_MAX));
		} else {
			impl.update_limit = 512;
			impl.eta_start = &start;
			impl.allow_sparse_eta = 1;
			CHECK(storage_simplex_basis_update_packed(
				&basis, 0, direction, &index, 1));
		}
	} else if (!strcmp(which, "ratio-leaf")) {
		memset(&state, 0, sizeof(state));
		state.candidate_count = INT_MAX;
		/* A leaf needs neither child arithmetic nor candidate access. This is
		 * a dimension-only probe, not a valid giant candidate array. */
		pricing_ratio_sift(&state, INT_MAX - 1);
	} else if (!strcmp(which, "lex-destroy")) {
		int position[] = {-1, 0}, candidate = 0;
		unsigned char status[] = {DUAL_STATUS_LOWER, DUAL_STATUS_BASIC};
		double breakpoint[] = {0., 0.};
		memset(&state, 0, sizeof(state));
		state.rows = 1;
		state.variables = 2;
		state.position = position;
		state.status = status;
		state.candidate_count = 1;
		state.candidate_index = &candidate;
		state.breakpoint = breakpoint;
		state.candidate_sign =
			dual_storage_test_malloc(2 * sizeof(int));
		CHECK(state.candidate_sign);
		state.candidate_sign[0] = 1;
		CHECK(pricing_lexicographic_candidate(&state, 0., 1, 1) == 0);
		storage_simplex_dual_pricing_destroy(&state);
		CHECK(!state.lex_storage && !state.lex_rho &&
			!state.lex_candidate && !state.lex_slot &&
			!state.lex_anchor_status && !state.lex_initialized);
		storage_simplex_dual_pricing_destroy(&state);
	} else
		return 1;
	CHECK(!active && !bad_free);
	return 0;
}

#ifndef DUAL_STORAGE_OLD
static int lex_layout(void)
{
	struct pricing_LexStorage layout;
	size_t numeric = 2 * sizeof(double), indices = 8 * sizeof(int);
	size_t total = numeric + indices + 4;
	CHECK(pricing_lex_storage(2, 4, total, &layout));
	CHECK(layout.numeric == numeric && layout.indices == indices &&
		layout.total == total);
	CHECK(!pricing_lex_storage(2, 4, total - 1, &layout) &&
		!layout.numeric && !layout.indices && !layout.total);
	CHECK(!pricing_lex_storage(2, 4, numeric + indices - 1, &layout));
	CHECK(!pricing_lex_storage(2, 4, numeric - 1, &layout));
	CHECK(!pricing_lex_storage(2, 4, indices - 1, &layout));
	CHECK(!pricing_lex_storage(-1, 4, (size_t)-1, &layout));
	CHECK(!pricing_lex_storage(2, 1, (size_t)-1, &layout));
	if (sizeof(double) == 8 && sizeof(int) == 4) {
		CHECK(!pricing_lex_storage(
			200000000, 400000000, (size_t)4294967295UL, &layout));
		CHECK(!layout.numeric && !layout.indices && !layout.total);
		if (sizeof(size_t) >= 8)
			CHECK(pricing_lex_storage(
				200000000, 400000000, (size_t)-1, &layout));
	}
	CHECK(!extremes("ratio-leaf") && !extremes("lex-destroy"));
	return 0;
}

static int lex_failure(int kappa, int ordered)
{
	struct simplex_Problem p;
	struct simplex_DualState state;
	struct lp_simplex_Options options;
	int starts[] = {0, 1, 3}, rows[] = {0, 0, 1};
	int row_starts[] = {0, 2, 3}, columns[] = {0, 1, 1};
	double values[] = {1., 1., 1.}, theta = 123.;
	unsigned char status[4];
	void *storage;
	int i, retained;
	problem_shape(&p, 2, 2);
	p.matrix.column_start = starts;
	p.matrix.row_index = rows;
	p.matrix.value = p.matrix.row_value = values;
	p.matrix.row_start = row_starts;
	p.matrix.column_index = columns;
	p.matrix.nonzeros = 3;
	lp_simplex_default_options(&options, lp_simplex_ALGORITHM_DUAL_REVISED);
	requests = fail_at = 0;
	CHECK(!storage_simplex_dual_state_create(&state, &p, &options, 0));
	for (i = 0; i < 4; i++) {
		state.position[i] = i < 2 ? -1 : i - 2;
		state.status[i] = i < 2 ? (kappa > 0 ? DUAL_STATUS_LOWER
						     : DUAL_STATUS_UPPER)
					: DUAL_STATUS_BASIC;
		state.reduced[i] = state.cost[i] = 0.;
	}
	state.basis[0] = 2;
	state.basis[1] = 3;
	CHECK(!storage_simplex_basis_factorize(&state.factor));
	state.rho[0] = -1.;
	state.rho[1] = 0.;
	storage_simplex_dual_nonbasic_initialize(&state);
	state.degeneracy.phase = SIMPLEX_PAN_LEXICOGRAPHIC;
	CHECK(storage_simplex_dual_prepare_ratio_test(&state, kappa) == 2);
	state.ratio_ordered = ordered;
	memcpy(status, state.status, sizeof(status));
	retained = active;
	requests = 0;
	fail_at = 1;
	CHECK(storage_simplex_dual_next_ratio_candidate(
		      &state, &theta, 1, 2, kappa) == SIMPLEX_DUAL_RATIO_ERROR);
	CHECK(requests == 1 && active == retained &&
		state.ratio_failure_status == lp_simplex_MemoryAllocError);
	CHECK(!state.lex_storage && !state.lex_rho && !state.lex_candidate &&
		!state.lex_slot && !state.lex_anchor_status &&
		!state.lex_initialized && theta == 123.);
	CHECK(!memcmp(status, state.status, sizeof(status)) &&
		state.basis[0] == 2 && state.basis[1] == 3);
	fail_at = 0;
	/* A stale failed candidate set cannot silently resume. A new preparation
	 * explicitly clears the failure before retrying the same symbolic tie. */
	CHECK(storage_simplex_dual_next_ratio_candidate(
		      &state, &theta, 1, 2, kappa) == SIMPLEX_DUAL_RATIO_ERROR);
	CHECK(requests == 1 && theta == 123.);
	CHECK(storage_simplex_dual_prepare_ratio_test(&state, kappa) == 2);
	state.ratio_ordered = ordered;
	CHECK(storage_simplex_dual_next_ratio_candidate(
		      &state, &theta, 1, 2, kappa) == 1 &&
		theta == 0.);
	CHECK(requests == 2 && active == retained + 1 &&
		state.lex_initialized &&
		!memcmp(status, state.lex_anchor_status, sizeof(status)));
	CHECK(state.lex_candidate ==
			(int *)((unsigned char *)state.lex_storage +
				2 * sizeof(double)) &&
		state.lex_slot == state.lex_candidate + 4 &&
		state.lex_anchor_status ==
			(unsigned char *)(state.lex_slot + 4));
	storage = state.lex_storage;
	CHECK(storage_simplex_dual_next_ratio_candidate(
		      &state, &theta, 1, 2, kappa) == 1);
	CHECK(requests == 2 && state.lex_storage == storage);
	for (i = 0; i < 4; i++)
		CHECK(state.cost[i] == 0.);
	storage_simplex_dual_pricing_destroy(&state);
	CHECK(!state.lex_storage && !state.lex_rho && !state.lex_candidate &&
		!state.lex_slot && !state.lex_anchor_status &&
		!state.lex_initialized);
	storage_simplex_dual_pricing_destroy(&state);
	storage_simplex_dual_state_destroy(&state);
	CHECK(!active && !bad_free);
	return 0;
}

static int layout(void)
{
	struct simplex_DualStorage d;
	struct basis_StartupStorage b;
	struct simplex_SparseVector vector;
	int variables, limit, previous = 0;
	size_t bytes, narrow = (size_t)4294967295UL;
	CHECK(!simplex_storage_bytes((size_t)-1, 2, (size_t)-1, &bytes) &&
		!bytes);
	CHECK(!simplex_dual_storage(1, INT_MAX, (size_t)-1, &d) &&
		!d.variables && !d.index && !d.status && !d.column && !d.row &&
		!d.multiplier && !d.frequency);
	CHECK(!simplex_dual_storage(-1, 2, (size_t)-1, &d));
	CHECK(!simplex_dual_storage(INT_MAX / 32, 1, narrow, &d));
	CHECK(simplex_dual_storage(2, 3, narrow, &d));
	CHECK(d.variables == 5 && d.index == 7 * sizeof(int) &&
		d.column == 35 * sizeof(double) &&
		d.row == 18 * sizeof(double));
	CHECK(!basis_startup_storage(INT_MAX, 0, narrow, &b));
	if (sizeof(size_t) >= 8) {
		CHECK(basis_startup_storage(INT_MAX, 0, (size_t)-1, &b) &&
			b.eta_capacity == INT_MAX);
		CHECK(simplex_dual_storage(0, INT_MAX, (size_t)-1, &d));
	}
	memset(&vector, 0, sizeof(vector));
	CHECK(storage_simplex_sparse_vector_create(&vector, -1));
	storage_simplex_sparse_vector_destroy(&vector);
	CHECK(!storage_simplex_sparse_vector_create(&vector, 0));
	storage_simplex_sparse_vector_destroy(&vector);
	CHECK(!active && !bad_free);
	/* Valid dimension metadata must not turn initialization into a negative
	 * or shrinking attempt budget. No giant matrix is allocated here. */
	CHECK(SIMPLEX_DUAL_CRASH_ATTEMPT_LIMIT(2, 5) == 22);
	for (variables = INT_MAX / 4 - 2; variables <= INT_MAX / 4 + 2;
		variables++) {
		limit = SIMPLEX_DUAL_CRASH_ATTEMPT_LIMIT(1, variables);
		CHECK(limit >= variables && limit >= previous);
		previous = limit;
	}
	CHECK(previous == INT_MAX &&
		SIMPLEX_DUAL_CRASH_ATTEMPT_LIMIT(0, INT_MAX) == INT_MAX &&
		SIMPLEX_DUAL_CRASH_ATTEMPT_LIMIT(INT_MAX, INT_MAX) == INT_MAX);
	return 0;
}

static int eta_solution(struct simplex_Basis *basis, int changed)
{
	double forward[16], transpose[16];
	int i;
	for (i = 0; i < 16; i++)
		forward[i] = transpose[i] = i + 1.;
	CHECK(!storage_simplex_basis_ftran(basis, forward));
	CHECK(!storage_simplex_basis_btran(basis, transpose));
	for (i = 0; i < 16; i++) {
		double f = -(i + 1.), t = f;
		/* B=-I before the update; afterwards column 0 is -e0-.5e3. */
		if (changed && i == 3)
			f += .5;
		if (changed && i == 0)
			t += 2.;
		CHECK(forward[i] == f && transpose[i] == t);
	}
	return 0;
}

static int eta_failure(int dense, int failure)
{
	struct simplex_CscMatrix matrix;
	struct simplex_Basis basis;
	struct simplex_BasisImpl *impl;
	int index[16], starts = 0, packed[] = {0, 3};
	int i, count, used, capacity;
	double direction[16] = {1.};
	double saved_value[64], saved_dense[64];
	int saved_index[64];
	memset(&matrix, 0, sizeof(matrix));
	matrix.rows = 16;
	matrix.column_start = &starts;
	for (i = 0; i < 16; i++)
		index[i] = i;
	requests = fail_at = 0;
	CHECK(!storage_simplex_basis_create(&basis, &matrix, 0, index));
	CHECK(!storage_simplex_basis_factorize(&basis));
	storage_simplex_basis_use_product_updates(&basis);
	storage_simplex_basis_set_sparse_eta(&basis, !dense);
	impl = basis.impl;
	impl->factor_work = 0;
	count = dense ? 4 : 64;
	for (i = 0; i < count; i++)
		CHECK(!storage_simplex_basis_update_packed(
			&basis, 0, direction, packed, 1));
	CHECK(!eta_solution(&basis, 0));
	used = impl->eta_start[count];
	capacity = dense ? impl->dense_eta_capacity : impl->eta_capacity;
	CHECK(capacity == (dense ? 4 : 64));
	if (dense)
		memcpy(saved_dense, impl->dense_eta, sizeof(saved_dense));
	else {
		memcpy(saved_value, impl->eta_value, sizeof(saved_value));
		memcpy(saved_index, impl->eta_index, sizeof(saved_index));
	}
	requests = 0;
	fail_at = failure;
	direction[3] = .5;
	CHECK(storage_simplex_basis_update_packed(
		&basis, 0, direction, packed, 2));
	CHECK(requests == failure &&
		storage_simplex_basis_storage_failed(&basis));
	CHECK(impl->update_count == count && impl->eta_start[count] == used &&
		impl->eta_capacity == 64 &&
		impl->dense_eta_count == (dense ? count : 0));
	if (dense)
		CHECK(impl->dense_eta_capacity == capacity &&
			!memcmp(saved_dense, impl->dense_eta,
				sizeof(saved_dense)));
	else
		CHECK(!memcmp(saved_value, impl->eta_value,
			      sizeof(saved_value)) &&
			!memcmp(saved_index, impl->eta_index,
				sizeof(saved_index)));
	CHECK(!eta_solution(&basis, 0));
	/* A numerical rejection must not inherit the preceding resource category. */
	direction[0] = 0.;
	CHECK(storage_simplex_basis_update(&basis, 0, direction));
	CHECK(!storage_simplex_basis_storage_failed(&basis));
	direction[0] = 1.;
	fail_at = 0;
	CHECK(!storage_simplex_basis_update_packed(
		&basis, 0, direction, packed, 2));
	CHECK(impl->update_count == count + 1 &&
		!storage_simplex_basis_storage_failed(&basis));
	CHECK(!eta_solution(&basis, 1));
	direction[3] = 0.;
	while (dense ? impl->dense_eta_count < impl->dense_eta_capacity
		     : impl->eta_start[impl->update_count] < impl->eta_capacity)
		CHECK(!storage_simplex_basis_update_packed(
			&basis, 0, direction, packed, 1));
	requests = 0;
	fail_at = 1;
	CHECK(storage_simplex_basis_update_packed(
		&basis, 0, direction, packed, 1));
	CHECK(storage_simplex_basis_storage_failed(&basis));
	CHECK(!eta_solution(&basis, 1));
	fail_at = 0;
	/* The caller owns the basis mapping; allocation rejection never rolls it
	 * back. Rebuilding this unchanged logical mapping discards the eta chain. */
	CHECK(!storage_simplex_basis_factorize(&basis));
	CHECK(!storage_simplex_basis_storage_failed(&basis) &&
		!impl->update_count);
	CHECK(!eta_solution(&basis, 0));
	storage_simplex_basis_destroy(&basis);
	CHECK(!active && !bad_free);
	return 0;
}

static int eta_storage(void)
{
	size_t bytes;
	CHECK(basis_eta_capacity(INT_MAX - 1, INT_MAX, 0) == INT_MAX);
	CHECK(basis_eta_capacity(INT_MAX - 1, INT_MAX, 1) == INT_MAX);
	CHECK(basis_eta_capacity(64, 65, 0) == 96);
	CHECK(basis_eta_capacity(4, 5, 1) == 7);
	CHECK(!basis_eta_bytes(INT_MAX, INT_MAX, (size_t)-1, &bytes) && !bytes);
	CHECK(!basis_eta_bytes(1024, 1024, 4294967, &bytes) && !bytes);
	CHECK(basis_eta_bytes(4, 16, 1024, &bytes) && bytes == 512);
	return eta_failure(0, 1) || eta_failure(0, 2) || eta_failure(1, 1) ||
	       extremes("eta-growth") || extremes("dense-growth") ||
	       extremes("eta-count");
}

static int ft_basis_failure(int failure)
{
	struct simplex_CscMatrix matrix;
	struct simplex_Basis basis;
	int starts[9], rows[8], index[8], i;
	double values[8], direction[8], forward[8], transpose[8];
	memset(&matrix, 0, sizeof(matrix));
	matrix.rows = matrix.columns = matrix.nonzeros = 8;
	matrix.column_start = starts;
	matrix.row_index = rows;
	matrix.value = values;
	for (i = 0; i < 8; i++) {
		starts[i] = rows[i] = index[i] = i;
		values[i] = 1.;
		direction[i] = i ? .5 : 1.;
		forward[i] = transpose[i] = i + 1.;
	}
	starts[8] = 8;
	requests = fail_at = 0;
	CHECK(!storage_simplex_basis_create(&basis, &matrix, 8, index));
	storage_simplex_basis_set_backend(
		&basis, SIMPLEX_BASIS_BACKEND_INTERNAL);
	CHECK(!storage_simplex_basis_factorize(&basis));
	CHECK(basis.impl->sparse.ft_active);
	requests = 0;
	fail_at = failure;
	CHECK(storage_simplex_basis_update(&basis, 0, direction) ==
		lp_simplex_EXIT_FAILURE);
	CHECK(storage_simplex_basis_storage_failed(&basis) &&
		!storage_simplex_basis_update_count(&basis));
	fail_at = 0;
	CHECK(!storage_simplex_basis_ftran(&basis, forward));
	CHECK(!storage_simplex_basis_btran(&basis, transpose));
	for (i = 0; i < 8; i++)
		CHECK(forward[i] == i + 1. && transpose[i] == i + 1.);
	CHECK(!storage_simplex_basis_update(&basis, 0, direction));
	CHECK(!storage_simplex_basis_storage_failed(&basis));
	storage_simplex_basis_destroy(&basis);
	CHECK(!active && !bad_free);
	return 0;
}

static int factor_storage(int backend)
{
	struct simplex_CscMatrix matrix;
	struct simplex_Basis basis;
	int starts[9], rows[8], index[8], i, failure, allocations;
	double values[8];
	memset(&matrix, 0, sizeof(matrix));
	matrix.rows = matrix.columns = matrix.nonzeros = 8;
	matrix.column_start = starts;
	matrix.row_index = rows;
	matrix.value = values;
	for (i = 0; i < 8; i++) {
		starts[i] = rows[i] = index[i] = i;
		values[i] = 1.;
	}
	starts[8] = 8;
	requests = fail_at = 0;
	CHECK(!storage_simplex_basis_create(&basis, &matrix, 8, index));
	storage_simplex_basis_set_backend(&basis, backend);
	requests = 0;
	CHECK(!storage_simplex_basis_factorize(&basis));
	allocations = requests;
	CHECK(allocations > 0 && !storage_simplex_basis_storage_failed(&basis));
	storage_simplex_basis_destroy(&basis);
	for (failure = 1; failure <= allocations; failure++) {
		fail_at = 0;
		CHECK(!storage_simplex_basis_create(&basis, &matrix, 8, index));
		storage_simplex_basis_set_backend(&basis, backend);
		requests = 0;
		fail_at = failure;
		CHECK(storage_simplex_basis_factorize(&basis));
		CHECK(storage_simplex_basis_storage_failed(&basis));
		fail_at = 0;
		/* Numerical singularity must clear the preceding resource category. */
		values[0] = 0.;
		CHECK(storage_simplex_basis_factorize(&basis));
		CHECK(!storage_simplex_basis_storage_failed(&basis));
		values[0] = 1.;
		CHECK(!storage_simplex_basis_factorize(&basis));
		CHECK(!storage_simplex_basis_storage_failed(&basis));
		storage_simplex_basis_destroy(&basis);
		CHECK(!active && !bad_free);
	}
#ifdef LP_SIMPLEX_HAVE_KLU
	if (backend == SIMPLEX_BASIS_BACKEND_DEFAULT)
		for (failure = 1; failure <= 6; failure++) {
			CHECK(!storage_simplex_basis_create(
				&basis, &matrix, 8, index));
			klu_failure = failure;
			CHECK(storage_simplex_basis_factorize(&basis));
			CHECK(storage_simplex_basis_storage_failed(&basis) ==
				(failure != 3 && failure != 6));
			klu_failure = 0;
			CHECK(!storage_simplex_basis_factorize(&basis));
			CHECK(!storage_simplex_basis_storage_failed(&basis));
			storage_simplex_basis_destroy(&basis);
			CHECK(!active && !bad_free);
		}
#endif
	printf("factor allocations=%d backend=%d\n", allocations, backend);
	return 0;
}
#endif

static int constructors(void)
{
	struct simplex_Problem p;
	struct simplex_DualState state, zero;
	struct lp_simplex_Options options;
	int starts[] = {0, 1, 2}, index[] = {0, 1};
	double values[] = {1., 3.};
	int baseline, failure, code, successes = 0;
	problem_shape(&p, 2, 2);
	p.matrix.column_start = starts;
	p.matrix.row_index = index;
	p.matrix.value = values;
	p.matrix.nonzeros = 2;
	lp_simplex_default_options(&options, lp_simplex_ALGORITHM_DUAL_REVISED);
	memset(&zero, 0, sizeof(zero));
	requests = fail_at = 0;
	CHECK(!storage_simplex_dual_state_create(&state, &p, &options, 0));
	baseline = requests;
	CHECK(state.regular_columns && state.modal_column_degree == 1);
	CHECK(state.variables == 4 && state.basis == state.workspace.index &&
		state.position == state.basis + 2);
	CHECK(state.breakpoint == state.lower + 6 * state.variables &&
		state.edge_weight == state.basic_value + 8 * state.rows);
	CHECK(state.nonbasic_slot ==
			state.candidate_sign + 3 * state.variables &&
		state.feasibility.structural_slot ==
			state.feasibility.heap + 3 * state.rows);
	CHECK(state.factor.impl->base_index ==
		state.factor.impl->core_basis + 5 * state.rows);
	storage_simplex_dual_state_destroy(&state);
	CHECK(!active && !bad_free);
	for (failure = 1; failure <= baseline; failure++) {
		requests = 0;
		fail_at = failure;
		code = storage_simplex_dual_state_create(
			&state, &p, &options, 0);
		if (!code) {
			/* The optional regular-column histogram can fail without rejecting a solve. */
			CHECK(failure == 1 && !state.regular_columns);
			successes++;
		} else
			CHECK(!memcmp(&state, &zero, sizeof(state)));
		storage_simplex_dual_state_destroy(&state);
		CHECK(!active && !bad_free);
	}
	CHECK(successes == 1);
	fail_at = 0;
	printf("startup allocations=%d optional=1 backend=%s\n", baseline,
#ifdef LP_SIMPLEX_HAVE_KLU
		"KLU"
#else
		"internal"
#endif
	);
	return 0;
}

#ifndef DUAL_STORAGE_OLD
static int compact_map_refusal(int backend, int failure)
{
	struct simplex_CscMatrix matrix;
	struct simplex_Basis basis;
	struct simplex_SparseVector packed;
	int starts[] = {0, 2}, rows[4097], columns[] = {0, 0};
	int indices[] = {0, 1}, index[4096], i;
	double values[] = {2., .5}, right[4096];
	memset(&matrix, 0, sizeof(matrix));
	rows[0] = 0;
	rows[1] = 1;
	for (i = 2; i <= 4096; i++)
		rows[i] = 2;
	for (i = 0; i < 4096; i++)
		index[i] = i + 1;
	index[0] = 0;
	matrix.rows = 4096;
	matrix.columns = 1;
	matrix.nonzeros = 2;
	matrix.column_start = starts;
	matrix.row_index = indices;
	matrix.value = values;
	matrix.row_start = rows;
	matrix.column_index = columns;
	matrix.row_value = values;
	requests = fail_at = 0;
	CHECK(!storage_simplex_basis_create_configured(
		&basis, &matrix, 1, index, 0));
	storage_simplex_basis_set_backend(&basis, backend);
	storage_simplex_basis_request_compact(&basis);
	/* The map is the first optional request after rebuilding the core. */
	fail_at = requests + failure;
	CHECK(!storage_simplex_basis_factorize(&basis));
	CHECK(storage_simplex_basis_compact_active(&basis));
	CHECK(!storage_simplex_basis_storage_failed(&basis));
	if (failure == 1)
		CHECK(!basis.impl->column_to_core);
#ifdef LP_SIMPLEX_HAVE_KLU
	else
		CHECK(basis.impl->column_to_core &&
			!basis.impl->sparse_lift.index);
#endif
	fail_at = 0;
	CHECK(!storage_simplex_sparse_vector_create(&packed, 4096));
	for (i = 0; i < 2; i++) {
		memset(right, 0, sizeof(right));
		right[0] = 3.;
		right[1] = 4.;
		CHECK(!storage_simplex_basis_btran(&basis, right));
		CHECK(fabs(right[0] - 2.5) < 1e-14 && right[1] == -4.);
		storage_simplex_sparse_vector_clear(&packed);
		CHECK(!storage_simplex_sparse_vector_set(&packed, 0, 3.));
		CHECK(!storage_simplex_sparse_vector_set(&packed, 1, 4.));
		CHECK(!storage_simplex_basis_btran_packed(&basis, &packed));
		CHECK(fabs(packed.dense[0] - 2.5) < 1e-14 &&
			packed.dense[1] == -4.);
		fail_at = 0;
		CHECK(!storage_simplex_basis_factorize(&basis));
		CHECK(basis.impl->column_to_core != NULL);
	}
	storage_simplex_basis_destroy(&basis);
	storage_simplex_sparse_vector_destroy(&packed);
	CHECK(!active && !bad_free);
	return 0;
}
#endif

#if defined(LP_SIMPLEX_HAVE_KLU) && !defined(DUAL_STORAGE_OLD)
static int packed_work_refusal(int failure)
{
	struct simplex_Problem p;
	struct simplex_DualState state;
	struct lp_simplex_Options options;
	int starts[] = {0, 1, 2}, rows[] = {0, 1}, before;
	double values[] = {1., 3.};
	problem_shape(&p, 2, 2);
	p.matrix.column_start = starts;
	p.matrix.row_index = rows;
	p.matrix.value = values;
	p.matrix.nonzeros = 2;
	lp_simplex_default_options(&options, lp_simplex_ALGORITHM_DUAL_REVISED);
	requests = fail_at = 0;
	CHECK(!storage_simplex_dual_state_create(&state, &p, &options, 0));
	requests = 0;
	fail_at = failure;
	CHECK(storage_simplex_dual_packed_work_prepare(&state) == !failure);
	CHECK(state.packed_work_attempted &&
		(!!state.packed_work.index == !failure));
	before = requests;
	fail_at = 0;
	CHECK(storage_simplex_dual_packed_work_prepare(&state) == !failure);
	CHECK(requests == before);
	storage_simplex_dual_state_destroy(&state);
	CHECK(!active && !bad_free);
	return 0;
}
#endif

int main(int argc, char **argv)
{
#ifdef LP_SIMPLEX_HAVE_KLU
	int failure;
#endif
	if (argc > 1)
		return extremes(argv[1]);
#ifndef DUAL_STORAGE_OLD
	if (compact_map_refusal(SIMPLEX_BASIS_BACKEND_INTERNAL, 1) ||
		compact_map_refusal(SIMPLEX_BASIS_BACKEND_DEFAULT, 1) ||
		layout() || eta_storage() || ft_basis_failure(1) ||
		ft_basis_failure(2) ||
		factor_storage(SIMPLEX_BASIS_BACKEND_INTERNAL) ||
		lex_layout() || lex_failure(1, 0) || lex_failure(-1, 0) ||
		lex_failure(1, 1) || lex_failure(-1, 1))
		return 1;
#ifdef LP_SIMPLEX_HAVE_KLU
	for (failure = 0; failure <= 3; failure++)
		if (packed_work_refusal(failure))
			return 1;
	for (failure = 4; failure <= 9; failure++)
		if (compact_map_refusal(SIMPLEX_BASIS_BACKEND_DEFAULT, failure))
			return 1;
	if (factor_storage(SIMPLEX_BASIS_BACKEND_DEFAULT))
		return 1;
#endif
#endif
	return constructors() || extremes("state") || extremes("frequency") ||
	       extremes("heap") || extremes("heap-leaf") || extremes("basis") ||
	       extremes("entry");
}
