/**
 * @file test_presolve_storage.c
 * @brief Presolve startup ownership and queue arithmetic under allocation faults.
 */

#include <lp_simplex/lp_simplex.h>
#include "utils.h"
#include <limits.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void *live[64];
static int requests, fail_at, active, bad_free;

void *presolve_storage_malloc(size_t bytes)
{
	void *p;
	int i;

	requests++;
	if (requests == fail_at || bytes > 1024 * 1024)
		return NULL;
	p = malloc(bytes ? bytes : 1);
	if (!p)
		return NULL;
	memset(p, 0xa5, bytes);
	for (i = 0; i < 64; i++)
		if (!live[i]) {
			live[i] = p;
			active++;
			return p;
		}
	free(p);
	return NULL;
}

void presolve_storage_free(void *p)
{
	int i;

	if (!p)
		return;
	for (i = 0; i < 64; i++)
		if (live[i] == p) {
			live[i] = NULL;
			active--;
			free(p);
			return;
		}
	bad_free++;
}

#define lp_simplex_malloc presolve_storage_malloc
#define lp_simplex_free presolve_storage_free
#define simplex_presolve_queue_init storage_queue_init
#define simplex_presolve_queue_clear storage_queue_clear
#define simplex_presolve_queue_push storage_queue_push
#define simplex_presolve_queue_pop storage_queue_pop
#define simplex_presolve_queue_destroy storage_queue_destroy
#define simplex_presolve_run storage_presolve_run
#define simplex_presolve_destroy storage_presolve_destroy
#define simplex_presolve_postsolve storage_presolve_postsolve
#ifndef PRESOLVE_QUEUE_SOURCE
#define PRESOLVE_QUEUE_SOURCE "../src/presolve/simplex_presolve_queue.c"
#endif
#include PRESOLVE_QUEUE_SOURCE
#include "../src/presolve/simplex_presolve.c"
#include "../src/presolve/simplex_presolve_postsolve.c"
#undef lp_simplex_malloc
#undef lp_simplex_free

static int refuse_reduced_point, reduced_point_requests;

static void *reduced_point_malloc(size_t bytes)
{
	reduced_point_requests++;
	return refuse_reduced_point ? NULL : lp_simplex_malloc(bytes);
}

#define lp_simplex_malloc reduced_point_malloc
#define lp_simplex_default_options storage_default_options
#define lp_simplex_solve storage_solve
#define lp_simplex_solve_experimental storage_solve_experimental
#define simplex_solve_reported storage_solve_reported
#define simplex_solve_reported_experimental storage_solve_reported_experimental
#ifndef PRESOLVE_CORE_SOURCE
#define PRESOLVE_CORE_SOURCE "../src/core/simplex.c"
#endif
#include PRESOLVE_CORE_SOURCE
#undef lp_simplex_malloc

#define CHECK(c) \
	do { \
		if (!(c)) { \
			fprintf(stderr, "presolve storage line %d: %s\n", \
				__LINE__, #c); \
			return 1; \
		} \
	} while (0)

static int sizes(void)
{
	struct simplex_PresolveStorage s;
	size_t narrow = (size_t)4294967295UL;

	CHECK(simplex_presolve_storage(3, 4, (size_t)-1, &s));
	CHECK(s.column_int == 4 * sizeof(int) && s.row_int == 3 * sizeof(int));
	CHECK(s.column_double == 4 * sizeof(double));
	CHECK(s.row_double == 3 * sizeof(double));
	CHECK(s.bounds == 4 * sizeof(struct optm_VariableBound));
	CHECK(s.column_flags == 4 && s.row_flags == 3);
	CHECK(!simplex_presolve_storage(-1, 4, (size_t)-1, &s));
	CHECK(!s.column_int && !s.bounds && !s.row_flags);
	CHECK(!simplex_presolve_storage(3, -1, (size_t)-1, &s));
	CHECK(!simplex_presolve_storage(INT_MAX, 1, narrow, &s));
	CHECK(!simplex_presolve_storage(1, INT_MAX, narrow, &s));
	CHECK(!simplex_presolve_storage(1, 1, sizeof(double) - 1, &s));
	CHECK(!simplex_presolve_storage(
		1, 1, sizeof(struct optm_VariableBound) - 1, &s));
	CHECK(simplex_presolve_search_work(3, 4, 12, (size_t)-1) == 1216);
	CHECK(simplex_presolve_search_work(3, 4, 12, 1216) == 1216);
	CHECK(!simplex_presolve_search_work(3, 4, 12, 1215));
	CHECK(!simplex_presolve_search_work(3, 4, 12, 6));
	CHECK(!simplex_presolve_search_work(3, 4, 12, 18));
	CHECK(!simplex_presolve_search_work(-1, 4, 12, narrow));
	CHECK(!simplex_presolve_search_work(3, -1, 12, narrow));
	CHECK(!simplex_presolve_search_work(3, 4, -1, narrow));
	CHECK(!simplex_presolve_search_work(INT_MAX, INT_MAX, INT_MAX, narrow));
	CHECK(!simplex_presolve_search_work(1, 1, INT_MAX, narrow));
	return 0;
}

/* Only the computed low tail is accessed. This probes representable ring
 * metadata without claiming to construct or allocate a giant real queue. */
static int large_tail(void)
{
	struct simplex_PresolveQueue q;
	int items[8] = {0};
	unsigned char present[8] = {0};

	memset(&q, 0, sizeof(q));
	q.item = items;
	q.present = present;
	q.capacity = INT_MAX;
	q.head = INT_MAX - 2;
	q.count = 4;
	CHECK(!storage_queue_push(&q, 3));
	CHECK(items[2] == 3 && present[3] && q.count == 5);
	return 0;
}

static int queues(void)
{
	struct simplex_PresolveQueue q;
	int i, f, value;

	for (f = 0; f <= 2; f++) {
		requests = 0;
		fail_at = f;
		memset(&q, 0xa5, sizeof(q));
		if (f) {
			CHECK(storage_queue_init(&q, 3));
			CHECK(!q.item && !q.present && !q.capacity);
		} else {
			CHECK(!storage_queue_init(&q, 3));
			for (i = 0; i < 1000; i++) {
				CHECK(!storage_queue_push(&q, 0));
				CHECK(!storage_queue_push(&q, 1));
				CHECK(!storage_queue_push(&q, 1));
				CHECK(!storage_queue_push(&q, 2));
				CHECK(q.count == 3);
				CHECK(storage_queue_push(&q, -1));
				CHECK(storage_queue_push(&q, 3));
				CHECK(storage_queue_pop(&q, &value) &&
					value == 0);
				CHECK(!storage_queue_push(&q, 0));
				CHECK(storage_queue_pop(&q, &value) &&
					value == 1);
				CHECK(storage_queue_pop(&q, &value) &&
					value == 2);
				CHECK(storage_queue_pop(&q, &value) &&
					value == 0);
				CHECK(!storage_queue_pop(&q, &value));
			}
			CHECK(!storage_queue_push(&q, 2));
			storage_queue_clear(&q);
			CHECK(!q.count && !q.head && !q.present[2]);
		}
		CHECK(requests == 2);
		storage_queue_destroy(&q);
		storage_queue_destroy(&q);
		CHECK(!active && !bad_free);
	}
	requests = fail_at = 0;
	CHECK(storage_queue_init(&q, -1));
	CHECK(!q.capacity && !q.item && !q.present && !requests);
	CHECK(!storage_queue_init(&q, 0));
	CHECK(!storage_queue_pop(&q, &value));
	CHECK(storage_queue_push(&q, 0));
	storage_queue_destroy(&q);
	return 0;
}

static int workspace(int cached)
{
	struct lp_Model *m = lp_model_create(3, 4);
	struct presolve_Workspace w;
	struct simplex_Presolve p;
	int f, i, j, item, trace;
	struct lp_simplex_Options options;
	struct lp_simplex_Result result;
	struct simplex_Execution history;
	double x[4];

	CHECK(m);
	lp_simplex_default_options(&options, lp_simplex_ALGORITHM_AUTO);
	for (j = 0; j < 4; j++) {
		m->bounds[j].b_type = optm_BOUND_T_BS;
		m->bounds[j].lb = -j - 1.;
		m->bounds[j].ub = j + 2.;
		for (i = 0; i < 3; i++)
			m->constraints[i].coef[j] = i + j + 1.;
	}
	for (i = 0; i < 3; i++)
		m->constraints[i].rhs = 10. + i;
	if (cached)
		CHECK(!lp_model_build_sparse(m));
	for (f = 0; f <= 16; f++) {
		requests = 0;
		fail_at = f;
		memset(&p, 0, sizeof(p));
		memset(&w, 0xa5, sizeof(w));
		if (f) {
			CHECK(storage_presolve_run(&p, m, 1e-8));
			CHECK(requests >= f);
			CHECK(!p.column_map && !p.row_map && !p.eliminated);
			CHECK(!p.eliminated_value && !p.reduced);
			CHECK(p.status == lp_simplex_MemoryAllocError);
			for (trace = 0; trace < 2; trace++) {
				requests = 0;
				CHECK(storage_solve_reported(m, &options, x,
					&result, trace ? &history : NULL));
				CHECK(result.status ==
					lp_simplex_MemoryAllocError);
				CHECK(result.iterations == 0 && requests >= f);
				CHECK(!active && !bad_free);
				for (i = 0; i < 4; i++)
					CHECK(isnan(x[i]));
				if (trace) {
					CHECK(history.count == 1);
					CHECK(history.events[0].kind ==
						SIMPLEX_EXECUTION_SOLVE);
					CHECK(history.events[0].complete);
					CHECK(history.events[0].status ==
						lp_simplex_MemoryAllocError);
				}
			}
		} else {
			CHECK(!presolve_workspace_create(&p, m, 1e-8, &w));
			CHECK(requests == 16);
			CHECK(w.dominance_work ==
				(size_t)(cached ? 1216 : 448));
			CHECK(w.propagation_tolerance == 1e-8);
			for (j = 0; j < 4; j++) {
				CHECK(w.column_kept[j] == -1 &&
					w.column_active[j]);
				CHECK(w.column_degree[j] == 3 &&
					!p.eliminated[j]);
				CHECK(p.eliminated_value[j] == 0.);
				CHECK(!memcmp(w.bounds + j, m->bounds + j,
					sizeof(*m->bounds)));
				CHECK(storage_queue_pop(
					      &w.column_queue, &item) &&
					item == j);
			}
			for (i = 0; i < 3; i++) {
				CHECK(w.row_kept[i] == -1 && w.row_active[i]);
				CHECK(w.row_degree[i] == 4 &&
					w.adjusted_rhs[i] == 10. + i);
				CHECK(storage_queue_pop(&w.row_queue, &item) &&
					item == i);
			}
			presolve_workspace_destroy(&w);
			presolve_workspace_destroy(&w);
		}
		storage_presolve_destroy(&p);
		storage_presolve_destroy(&p);
		CHECK(!active && !bad_free);
		for (j = 0; j < 4; j++) {
			CHECK(m->objective[j] == 0.);
			CHECK(m->bounds[j].b_type == optm_BOUND_T_BS);
			CHECK(m->bounds[j].lb == -j - 1. &&
				m->bounds[j].ub == j + 2.);
			if (cached) {
				CHECK(m->column_start[j] == 3 * j);
				for (i = 0; i < 3; i++) {
					CHECK(m->row_index[3 * j + i] == i);
					CHECK(m->value[3 * j + i] ==
						i + j + 1.);
				}
			}
		}
		if (cached) {
			CHECK(m->column_start[4] == 12 &&
				m->row_start[3] == 12);
			for (i = 0; i < 3; i++) {
				CHECK(m->row_start[i] == 4 * i);
				for (j = 0; j < 4; j++) {
					CHECK(m->column_index[4 * i + j] == j);
					CHECK(m->row_value[4 * i + j] ==
						i + j + 1.);
				}
			}
		}
		for (i = 0; i < 3; i++) {
			CHECK(m->constraints[i].rhs == 10. + i);
			for (j = 0; j < 4; j++)
				CHECK(m->constraints[i].coef[j] == i + j + 1.);
		}
	}
	/* Malformed private dimensions are rejected before any allocation or read. */
	requests = fail_at = 0;
	m->m = -1;
	CHECK(storage_presolve_run(&p, m, 1e-8));
	CHECK(!requests && !active && !bad_free);
	m->m = 3;
	m->n = -1;
	CHECK(storage_presolve_run(&p, m, 1e-8));
	CHECK(!requests && !active && !bad_free);
	m->n = 4;
	m->nnz = -1;
	CHECK(storage_presolve_run(&p, m, 1e-8));
	CHECK(!requests && !active && !bad_free);
	m->nnz = cached ? 12 : 0;
	lp_model_free(m);
	return 0;
}

static int exact_proportionality(void)
{
	double tiny = nextafter(0., 1.);
	CHECK(presolve_exact_scaled_coefficient(1., 1., 1.));
	CHECK(presolve_exact_scaled_coefficient(2., 6., 3.));
	CHECK(presolve_exact_scaled_coefficient(-3., -6., 2.));
	CHECK(presolve_exact_scaled_coefficient(1., -2., -2.));
	CHECK(presolve_exact_scaled_coefficient(DBL_MAX, DBL_MAX, 1.));
	CHECK(presolve_exact_scaled_coefficient(tiny, 2. * tiny, 2.));
	CHECK(presolve_exact_scaled_coefficient(2. * tiny, tiny, .5));
	CHECK(!presolve_exact_scaled_coefficient(tiny, tiny, .5));
	CHECK(!presolve_exact_scaled_coefficient(DBL_MAX, DBL_MAX, 2.));
	CHECK(!presolve_exact_scaled_coefficient(3., 1., 1. / 3.));
	CHECK(!presolve_exact_scaled_coefficient(1., 1. + ldexp(1., -40), 1.));
	CHECK(!presolve_exact_scaled_coefficient(1., 1., HUGE_VAL));
	CHECK(!presolve_exact_scaled_coefficient(1., 1., 0.));
	CHECK(!presolve_exact_scaled_coefficient(NAN, 1., 1.));
	return 0;
}

static int reduced_output_refusal(void)
{
	struct lp_Model *m;
	struct simplex_Presolve pre;
	struct lp_simplex_Options o;
	struct lp_simplex_Result r;
	struct simplex_Execution h;
	double x[4];
	int cached, tracing, algorithm, i, j, code;
	for (cached = 0; cached < 2; cached++) {
		m = lp_model_create(3, 4);
		CHECK(m);
		for (j = 0; j < 4; j++) {
			m->objective[j] = 1.;
			m->bounds[j].b_type = optm_BOUND_T_BS;
			m->bounds[j].ub = j == 3 ? 0. : 10.;
		}
		for (i = 0; i < 3; i++) {
			m->constraints[i].type = optm_CONS_T_GE;
			m->constraints[i].rhs = 1.;
			for (j = 0; j < 4; j++)
				m->constraints[i].coef[j] = i == j ? 2. : 1.;
		}
		if (cached)
			CHECK(!lp_model_build_sparse(m));
		fail_at = 0;
		CHECK(!storage_presolve_run(&pre, m, 1e-7));
		CHECK(pre.status == lp_simplex_CondUnsatisfied && pre.reduced &&
			pre.reduced->columns < m->n &&
			!simplex_pan_original_form_smaller(m, pre.reduced));
		storage_presolve_destroy(&pre);
		CHECK(!active && !bad_free);
		for (algorithm = lp_simplex_ALGORITHM_DUAL_REVISED;
			algorithm <= lp_simplex_ALGORITHM_AUTO; algorithm++)
			for (tracing = 0; tracing < 2; tracing++) {
				storage_default_options(&o, algorithm);
				refuse_reduced_point = 1;
				reduced_point_requests = 0;
				memset(&r, 0xff, sizeof(r));
				code = storage_solve_reported(
					m, &o, x, &r, tracing ? &h : NULL);
				refuse_reduced_point = 0;
				if (r.status != lp_simplex_MemoryAllocError)
					fprintf(stderr,
						"reduced output refusal: alg=%d status=%d requests=%d\n",
						algorithm, r.status,
						reduced_point_requests);
				CHECK(code &&
					r.status ==
						lp_simplex_MemoryAllocError &&
					!r.iterations &&
					reduced_point_requests == 1);
				CHECK(isnan(r.objective) &&
					isnan(r.dual_infeasibility) &&
					isinf(r.primal_infeasibility));
				for (i = 0; i < 4; i++)
					CHECK(isnan(x[i]));
				CHECK(!active && !bad_free);
				if (tracing)
					CHECK(h.count == 1 &&
						h.events[0].complete &&
						h.events[0].status ==
							lp_simplex_MemoryAllocError &&
						!h.events[0].work);
			}
		lp_model_free(m);
	}
	return 0;
}

int main(int argc, char **argv)
{
	if (argc > 1 && !strcmp(argv[1], "large-tail"))
		return large_tail();
	return reduced_output_refusal() || exact_proportionality() || sizes() ||
	       queues() || large_tail() || workspace(0) || workspace(1);
}
