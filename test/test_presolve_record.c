/**
 * @file test_presolve_record.c
 * @brief Presolve journal publication, moving realloc and compaction faults.
 */

#include <lp_simplex/lp_simplex.h>
#include "utils.h"
#include "simplex_model_audit.h"
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void *live[512];
static size_t lengths[512];
static int active, bad_free, requests, fail_at;

static void *record_allocate(size_t bytes)
{
	void *p;
	int i;
	if (bytes > 1024 * 1024)
		return NULL;
	p = malloc(bytes ? bytes : 1);
	if (!p)
		return NULL;
	memset(p, 0xa5, bytes);
	for (i = 0; i < 512; i++)
		if (!live[i]) {
			live[i] = p;
			lengths[i] = bytes;
			active++;
			return p;
		}
	free(p);
	return NULL;
}

void *record_malloc(size_t bytes)
{
	requests++;
	return requests == fail_at ? NULL : record_allocate(bytes);
}

void record_free(void *p)
{
	int i;
	if (!p)
		return;
	for (i = 0; i < 512; i++)
		if (live[i] == p) {
			memset(p, 0xdd, lengths[i]);
			free(p);
			live[i] = NULL;
			active--;
			return;
		}
	bad_free++;
}

void *record_realloc(void *p, size_t bytes)
{
	void *q;
	int i;

	requests++;
	if (requests == fail_at)
		return NULL;
	if (!p)
		return record_allocate(bytes);
	for (i = 0; i < 512; i++)
		if (live[i] == p) {
			q = record_allocate(bytes);
			if (!q)
				return NULL;
			memcpy(q, p, lengths[i] < bytes ? lengths[i] : bytes);
			record_free(p);
			return q;
		}
	bad_free++;
	return NULL;
}

#define lp_simplex_malloc record_malloc
#define lp_simplex_realloc record_realloc
#define lp_simplex_free record_free
#define simplex_presolve_queue_init record_queue_init
#define simplex_presolve_queue_clear record_queue_clear
#define simplex_presolve_queue_push record_queue_push
#define simplex_presolve_queue_pop record_queue_pop
#define simplex_presolve_queue_destroy record_queue_destroy
#define simplex_presolve_record_begin record_begin
#define simplex_presolve_record_add_term record_add_term
#define simplex_presolve_substitute_doubletons record_substitute
#define simplex_presolve_run record_run
#define simplex_presolve_postsolve record_postsolve
#define simplex_presolve_destroy record_destroy
#include "../src/presolve/simplex_presolve_queue.c"
#ifndef PRESOLVE_RECORD_SOURCE
#define PRESOLVE_RECORD_SOURCE "../src/presolve/simplex_presolve_record.c"
#endif
#include PRESOLVE_RECORD_SOURCE
#ifndef PRESOLVE_SUBSTITUTION_SOURCE
#define PRESOLVE_SUBSTITUTION_SOURCE \
	"../src/presolve/simplex_presolve_substitution.c"
#endif
#include PRESOLVE_SUBSTITUTION_SOURCE
#include "../src/presolve/simplex_presolve.c"
#include "../src/presolve/simplex_presolve_postsolve.c"
#undef lp_simplex_malloc
#undef lp_simplex_realloc
#undef lp_simplex_free
#define lp_simplex_default_options record_default_options
#define lp_simplex_solve record_solve
#define lp_simplex_solve_experimental record_solve_experimental
#define simplex_solve_reported record_solve_reported
#define simplex_solve_reported_experimental record_solve_reported_experimental
#include "../src/core/simplex.c"

#define CHECK(c) \
	do { \
		if (!(c)) { \
			fprintf(stderr, "presolve record line %d: %s\n", \
				__LINE__, #c); \
			return 1; \
		} \
	} while (0)

static int sizes(void)
{
	int cap;
	size_t bytes, narrow = (size_t)4294967295UL;
	CHECK(simplex_presolve_growth(
		4, 4, 1, 4, sizeof(double), (size_t)-1, &cap, &bytes));
	CHECK(cap == 8 && bytes == 8 * sizeof(double));
	CHECK(!simplex_presolve_growth(
		INT_MAX, INT_MAX, 1, 4, sizeof(int), (size_t)-1, &cap, &bytes));
	CHECK(!cap && !bytes);
	CHECK(!simplex_presolve_growth(
		4, -1, 1, 4, sizeof(int), narrow, &cap, &bytes));
	CHECK(!simplex_presolve_growth(
		4, 5, 1, 4, sizeof(int), narrow, &cap, &bytes));
	CHECK(!simplex_presolve_growth(
		4, 4, -1, 4, sizeof(int), narrow, &cap, &bytes));
	CHECK(!simplex_presolve_growth(INT_MAX / 2, INT_MAX / 2, 1, 4,
		sizeof(double), narrow, &cap, &bytes));
	CHECK(simplex_presolve_growth(
		4, 4, 1, 4, sizeof(double), 5 * sizeof(double), &cap, &bytes));
	CHECK(cap == 5 && bytes == 5 * sizeof(double));
	return 0;
}

static int record_overflow(void)
{
	struct simplex_Presolve p;
	memset(&p, 0, sizeof(p));
	p.journal.count = p.journal.capacity = INT_MAX;
	requests = fail_at = 0;
	CHECK(!record_begin(&p, 0, 0., 0));
	CHECK(!requests);
	return 0;
}

static int journals(void)
{
	struct simplex_Presolve p;
	int f, i, saved_count, saved_terms;
	struct simplex_PresolveSubstitution saved[32];
	int columns[64];
	double multipliers[64];

	for (f = 0; f <= 3; f++) {
		memset(&p, 0, sizeof(p));
		requests = 0;
		fail_at = f;
		if (f) {
			CHECK(!record_begin(&p, 0, 1., 2));
			CHECK(!p.journal.count && !p.journal.term_count);
			fail_at = 0;
		}
		CHECK(record_begin(&p, 0, 1., 2));
		CHECK(!p.journal.count && !p.journal.term_count);
		record_add_term(&p, 1, 2.);
		CHECK(!p.journal.count && p.journal.term_count == 1);
		record_add_term(&p, 2, -3.);
		CHECK(p.journal.count == 1 && p.journal.term_count == 2);
		CHECK(record_begin(&p, 3, 4., 0));
		CHECK(p.journal.count == 2 && p.journal.term_count == 2);
		record_destroy(&p);
		record_destroy(&p);
		CHECK(!active && !bad_free);
	}
	for (f = 1; f <= 3; f++) {
		memset(&p, 0, sizeof(p));
		fail_at = 0;
		for (i = 0; i < 32; i++) {
			CHECK(record_begin(&p, i, i + 1., 2));
			record_add_term(&p, i + 1, i + 2.);
			record_add_term(&p, i + 2, -i - 3.);
		}
		saved_count = p.journal.count;
		saved_terms = p.journal.term_count;
		memcpy(saved, p.journal.record, sizeof(saved));
		memcpy(columns, p.journal.term_column, sizeof(columns));
		memcpy(multipliers, p.journal.term_multiplier,
			sizeof(multipliers));
		requests = 0;
		fail_at = f;
		CHECK(!record_begin(&p, 99, 7., 2));
		CHECK(p.journal.count == saved_count &&
			p.journal.term_count == saved_terms);
		CHECK(!memcmp(saved, p.journal.record, sizeof(saved)));
		CHECK(!memcmp(columns, p.journal.term_column, sizeof(columns)));
		CHECK(!memcmp(multipliers, p.journal.term_multiplier,
			sizeof(multipliers)));
		requests = fail_at = 0;
		CHECK(!record_begin(&p, 99, 7., -1));
		CHECK(!requests && p.journal.count == saved_count);
		CHECK(record_begin(&p, 99, 7., 2));
		record_add_term(&p, 100, 8.);
		record_add_term(&p, 101, -9.);
		CHECK(p.journal.count == 33 && p.journal.term_count == 66);
		CHECK(!memcmp(saved, p.journal.record, sizeof(saved)));
		record_destroy(&p);
		CHECK(!active && !bad_free);
	}
	return 0;
}

static int replay_publication(void)
{
	struct simplex_Presolve p;
	struct lp_Model model;
	unsigned char eliminated[3] = {0, 0, 1};
	double eliminated_value[3] = {0., 0., 5.};
	double x[3];

	memset(&p, 0, sizeof(p));
	memset(&model, 0, sizeof(model));
	model.n = 3;
	p.original = &model;
	p.eliminated = eliminated;
	p.eliminated_value = eliminated_value;
	requests = fail_at = 0;
	CHECK(record_begin(&p, 0, 1., 2));
	record_add_term(&p, 1, 2.);
	record_postsolve(&p, NULL, x);
	CHECK(x[0] == 0. && x[1] == 0. && x[2] == 5.);
	record_add_term(&p, 2, -3.);
	CHECK(record_begin(&p, 1, 4., 1));
	record_postsolve(&p, NULL, x);
	CHECK(x[0] == -14. && x[1] == 0. && x[2] == 5.);
	record_add_term(&p, 2, 2.);
	record_postsolve(&p, NULL, x);
	CHECK(x[0] == 14. && x[1] == 14. && x[2] == 5.);
	p.eliminated = NULL;
	p.eliminated_value = NULL;
	record_destroy(&p);
	CHECK(!active && !bad_free);
	return 0;
}

static int mutable_growth(void)
{
	struct presolve_MutableRow row;
	struct presolve_Incidence incidence;
	int f, i, inserted, removed;
	for (f = 1; f <= 2; f++) {
		memset(&row, 0, sizeof(row));
		fail_at = 0;
		for (i = 0; i < 4; i++)
			CHECK(!substitution_add_to_row(
				&row, i, i + 1., &inserted, &removed));
		requests = 0;
		fail_at = f;
		CHECK(substitution_add_to_row(
			&row, 4, 5., &inserted, &removed));
		CHECK(!inserted && !removed && row.count == 4 &&
			row.capacity == 4);
		for (i = 0; i < 4; i++)
			CHECK(row.column[i] == i && row.value[i] == i + 1.);
		fail_at = 0;
		CHECK(!substitution_add_to_row(
			&row, 4, 5., &inserted, &removed));
		CHECK(inserted && !removed && row.count == 5);
		CHECK(!substitution_add_to_row(
			&row, 2, -3., &inserted, &removed));
		CHECK(removed && !inserted && row.count == 4);
		record_free(row.column);
		record_free(row.value);
		CHECK(!active && !bad_free);
	}
	memset(&incidence, 0, sizeof(incidence));
	for (i = 0; i < 4; i++)
		CHECK(!substitution_append_incidence(&incidence, i));
	requests = 0;
	fail_at = 1;
	CHECK(!substitution_append_incidence(&incidence, 2) && !requests);
	CHECK(substitution_append_incidence(&incidence, 4));
	CHECK(incidence.count == 4 && incidence.capacity == 4);
	for (i = 0; i < 4; i++)
		CHECK(incidence.row[i] == i);
	fail_at = 0;
	CHECK(!substitution_append_incidence(&incidence, 4));
	CHECK(incidence.count == 5);
	record_free(incidence.row);
	CHECK(!active && !bad_free);
	return 0;
}

static int compact(void)
{
	struct simplex_Problem *problem =
		simplex_problem_create_sparse(3, 4, 0);
	struct simplex_Problem *reduced;
	struct simplex_Presolve p;
	struct substitution_Workspace w;
	struct presolve_MutableRow rows[3];
	unsigned char columns[4] = {0, 1, 0, 1};
	int column_map[4], row_map[3], rcols[2] = {1, 3};
	double values[2] = {2., 3.};
	int f, i;

	CHECK(problem);
	memset(&p, 0, sizeof(p));
	memset(&w, 0, sizeof(w));
	memset(rows, 0, sizeof(rows));
	p.column_map = column_map;
	p.row_map = row_map;
	w.problem = problem;
	w.rows = rows;
	w.column_active = columns;
	w.active_rows = 2;
	w.active_columns = 2;
	for (i = 1; i < 3; i++) {
		rows[i].active = 1;
		rows[i].count = 2;
		rows[i].column = rcols;
		rows[i].value = values;
		problem->rhs[i] = i + 4.;
		problem->row_type[i] = optm_CONS_T_EQ;
	}
	for (f = 0; f <= 2; f++) {
		for (i = 0; i < 4; i++)
			column_map[i] = i + 10;
		for (i = 0; i < 3; i++)
			row_map[i] = i + 20;
		requests = 0;
		fail_at = f;
		reduced = substitution_compact_problem(&p, &w);
		if (f) {
			CHECK(!reduced);
			for (i = 0; i < 4; i++)
				CHECK(column_map[i] == i + 10);
			for (i = 0; i < 3; i++)
				CHECK(row_map[i] == i + 20);
		} else {
			CHECK(reduced && requests == 2);
			CHECK(column_map[0] == 11 && column_map[1] == 13);
			CHECK(row_map[0] == 21 && row_map[1] == 22);
			CHECK(reduced->matrix.nonzeros == 4);
			CHECK(reduced->rhs[0] == 5. && reduced->rhs[1] == 6.);
			for (i = 0; i < 4; i++) {
				CHECK(reduced->matrix.column_index[i] == i % 2);
				CHECK(reduced->matrix.row_value[i] ==
					values[i % 2]);
			}
		}
		simplex_problem_free(reduced);
		CHECK(!active && !bad_free);
	}
	simplex_problem_free(problem);
	return 0;
}

static struct lp_Model *fixture(void)
{
	struct lp_Model *m = lp_model_create(2, 3);
	int j;
	if (!m)
		return NULL;
	for (j = 0; j < 3; j++) {
		m->bounds[j].b_type = optm_BOUND_T_BS;
		m->bounds[j].ub = 3.;
	}
	m->objective[1] = 1.;
	m->constraints[0].rhs = 2.;
	m->constraints[0].coef[0] = m->constraints[0].coef[1] = 1.;
	m->constraints[1].rhs = 3.;
	m->constraints[1].coef[0] = m->constraints[1].coef[1] =
		m->constraints[1].coef[2] = 1.;
	if (lp_model_build_sparse(m)) {
		lp_model_free(m);
		return NULL;
	}
	return m;
}

static int fixture_unchanged(
	const struct lp_Model *m, const struct lp_Model *saved)
{
	static const int starts[4] = {0, 2, 4, 5};
	static const int row_starts[3] = {0, 2, 5};
	static const int rows[5] = {0, 1, 0, 1, 1};
	static const int columns[5] = {0, 1, 0, 1, 2};
	int i, j;

	if (memcmp(m, saved, sizeof(*m)) ||
		memcmp(m->column_start, starts, sizeof(starts)) ||
		memcmp(m->row_start, row_starts, sizeof(row_starts)) ||
		memcmp(m->row_index, rows, sizeof(rows)) ||
		memcmp(m->column_index, columns, sizeof(columns)))
		return 0;
	for (i = 0; i < 5; i++)
		if (m->value[i] != 1. || m->row_value[i] != 1.)
			return 0;
	for (j = 0; j < 3; j++)
		if (m->objective[j] != (j == 1 ? 1. : 0.) ||
			m->bounds[j].b_type != optm_BOUND_T_BS ||
			m->bounds[j].v_type != optm_VAR_T_REAL ||
			m->bounds[j].lb != 0. || m->bounds[j].ub != 3.)
			return 0;
	for (i = 0; i < 2; i++) {
		if (m->constraints[i].rhs != i + 2. ||
			m->constraints[i].type != optm_CONS_T_EQ)
			return 0;
		for (j = 0; j < 3; j++)
			if (m->constraints[i].coef[j] !=
				(i == 0 && j == 2 ? 0. : 1.))
				return 0;
	}
	return 1;
}

static int flow(void)
{
	struct lp_Model *m = fixture();
	struct lp_Model saved;
	struct simplex_Presolve p;
	struct lp_simplex_Options options;
	struct lp_simplex_Result result;
	struct simplex_Execution history;
	double x[3];
	int count, f, trace;

	CHECK(m);
	memcpy(&saved, m, sizeof(saved));
	requests = fail_at = 0;
	CHECK(!record_run(&p, m, 1e-8));
	count = requests;
	CHECK(p.status == lp_simplex_Success && p.journal.count > 0);
	record_postsolve(&p, NULL, x);
	CHECK(x[0] == 2. && x[1] == 0. && x[2] == 1.);
	CHECK(simplex_primal_infeasibility(m, x) == 0.);
	record_destroy(&p);
	CHECK(!active && !bad_free);
	record_default_options(&options, lp_simplex_ALGORITHM_AUTO);
	for (f = 1; f <= count; f++) {
		requests = 0;
		fail_at = f;
		CHECK(record_run(&p, m, 1e-8));
		CHECK(p.status == lp_simplex_MemoryAllocError);
		record_destroy(&p);
		CHECK(!active && !bad_free);
		for (trace = 0; trace < 2; trace++) {
			requests = 0;
			CHECK(record_solve_reported(m, &options, x, &result,
				trace ? &history : NULL));
			CHECK(result.status == lp_simplex_MemoryAllocError &&
				!result.iterations);
			CHECK(!active && !bad_free);
			if (trace) {
				CHECK(history.count == 1 &&
					history.events[0].complete);
				CHECK(history.events[0].status ==
					lp_simplex_MemoryAllocError);
			}
		}
		CHECK(fixture_unchanged(m, &saved));
	}
	printf("Presolve complete-path allocation points: %d\n", count);
	lp_model_free(m);
	return 0;
}

int main(int argc, char **argv)
{
	if (argc > 1 && !strcmp(argv[1], "overflow"))
		return record_overflow();
	if (argc > 1 && !strcmp(argv[1], "compact"))
		return compact();
	return sizes() || record_overflow() || journals() ||
	       replay_publication() || mutable_growth() || compact() || flow();
}
