/**
 * @file test_dominance_storage.c
 * @brief Optional dominance index limits and allocation-failure fallbacks.
 */

#include <lp_simplex/lp_simplex.h>
#include "utils.h"
#include "simplex_model_audit.h"
#include "policy/resources.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef TEST_DOMINANCE_BYTES
#undef DOMINANCE_BYTES
#define DOMINANCE_BYTES TEST_DOMINANCE_BYTES
#endif

static void *live[16];
static size_t length[16], owned, peak;
static unsigned long reject;
static int requests, bad_free;

void *dominance_test_malloc(size_t bytes)
{
	void *p;
	int i;
	requests++;
	if (requests < 32 && (reject & (1UL << requests)))
		return NULL;
	if (bytes > DOMINANCE_BYTES)
		return NULL;
	p = malloc(bytes ? bytes : 1);
	if (!p)
		return NULL;
	memset(p, 0xa5, bytes);
	for (i = 0; i < 16; i++)
		if (!live[i]) {
			live[i] = p;
			length[i] = bytes;
			owned += bytes;
			if (owned > peak)
				peak = owned;
			return p;
		}
	free(p);
	return NULL;
}

void dominance_test_free(void *p)
{
	int i;
	if (!p)
		return;
	for (i = 0; i < 16; i++)
		if (live[i] == p) {
			owned -= length[i];
			memset(p, 0xdd, length[i]);
			free(p);
			live[i] = NULL;
			return;
		}
	bad_free++;
}

#define lp_simplex_malloc dominance_test_malloc
#define lp_simplex_free dominance_test_free
#define simplex_presolve_dominance storage_dominance
#include "../src/presolve/simplex_presolve_dominance.c"
#undef lp_simplex_malloc
#undef lp_simplex_free

#define CHECK(c) \
	do { \
		if (!(c)) { \
			fprintf(stderr, "dominance storage line %d: %s\n", \
				__LINE__, #c); \
			return 1; \
		} \
	} while (0)

static int sizes(void)
{
	struct dominance_BitStorage b;
	struct dominance_PostingStorage p;
	size_t required, narrow = (size_t)4294967295UL;

	CHECK(!dominance_bit_storage(0, 4, 0, narrow, &b));
	CHECK(!dominance_bit_storage(3, -1, 0, narrow, &b));
	CHECK(!dominance_bit_storage(3, 4, 10, 9, &b));
	CHECK(dominance_bit_storage(3, 4, 3, (size_t)-1, &b));
	CHECK(b.words == 1 && b.bits == 3 * sizeof(unsigned long));
	CHECK(b.candidates == sizeof(unsigned long));
	required = 3 + b.bits + b.candidates;
	CHECK(dominance_bit_storage(3, 4, 3, required, &b));
	CHECK(!dominance_bit_storage(3, 4, 3, required - 1, &b));
	CHECK(!b.words && !b.bits && !b.candidates);
	CHECK(!dominance_bit_storage(INT_MAX, INT_MAX, 0, narrow, &b));
	CHECK(dominance_bit_storage(1, INT_MAX, 0, (size_t)-1, &b));
	CHECK(b.words == (size_t)INT_MAX / WORD_BITS +
				 ((size_t)INT_MAX % WORD_BITS != 0));
	CHECK(dominance_posting_storage(3, 4, 9, 67, (size_t)-1, &p));
	CHECK(p.order == 4 * sizeof(struct dominance_ColumnOrder));
	CHECK(p.starts == 4 * sizeof(int) && p.postings == 9 * sizeof(int));
	CHECK(p.cursor == 3 * sizeof(int));
	required = 67 + p.order + p.starts + p.postings + p.cursor;
	CHECK(dominance_posting_storage(3, 4, 9, 67, required, &p));
	CHECK(!dominance_posting_storage(3, 4, 9, 67, required - 1, &p));
	CHECK(!p.order && !p.starts && !p.postings && !p.cursor);
	CHECK(!dominance_posting_storage(INT_MAX, 1, 1, 0, narrow, &p));
	CHECK(!dominance_posting_storage(
		1, 1, (size_t)INT_MAX + 1, 0, (size_t)-1, &p));
	return 0;
}

static int empty(void)
{
	struct lp_Model m;
	int start[2] = {0, 0}, removed, fixed;
	size_t work = 123;

	memset(&m, 0, sizeof(m));
	m.n = 1;
	m.row_start = m.column_start = start;
	requests = 0;
	storage_dominance(&m, NULL, NULL, NULL, NULL, NULL, NULL, &work,
		&removed, &fixed);
	CHECK(!requests && !removed && !fixed && work == 123);
	m.m = 1;
	m.n = 0;
	storage_dominance(&m, NULL, NULL, NULL, NULL, NULL, NULL, &work,
		&removed, &fixed);
	CHECK(!requests && !removed && !fixed && work == 123);
	return 0;
}

static const double coefficients[3][4] = {
	{1., 2., 0., 1.}, {1., 2., 1., 0.}, {0., 1., 1., 1.}};
static const double costs[4] = {3., 2., 4., 5.};

static struct lp_Model *fixture(int columns)
{
	struct lp_Model *m = lp_model_create(3, columns);
	int i, j;
	if (!m)
		return NULL;
	for (i = 0; i < 3; i++) {
		m->constraints[i].type = optm_CONS_T_GE;
		m->constraints[i].rhs = i == 2 ? 1. : 2.;
		for (j = 0; j < columns; j++)
			m->constraints[i].coef[j] =
				j < 4 ? coefficients[i][j] : 0.;
	}
	for (j = 0; j < columns; j++) {
		m->objective[j] = j < 4 ? costs[j] : 3.;
		m->bounds[j].b_type = optm_BOUND_T_BS;
		m->bounds[j].ub = 3.;
	}
	if (lp_model_build_sparse(m)) {
		lp_model_free(m);
		return NULL;
	}
	return m;
}

static int unchanged(const struct lp_Model *m, const struct lp_Model *saved)
{
	static const int starts[5] = {0, 2, 5, 7, 9};
	static const int rstarts[4] = {0, 3, 6, 9};
	static const int rindex[9] = {0, 1, 0, 1, 2, 1, 2, 0, 2};
	static const int cindex[9] = {0, 1, 3, 0, 1, 2, 1, 2, 3};
	static const double value[9] = {1., 1., 2., 2., 1., 1., 1., 1., 1.};
	static const double rvalue[9] = {1., 2., 1., 1., 2., 1., 1., 1., 1.};
	int i, j;
	if (memcmp(m, saved, sizeof(*m)) ||
		memcmp(m->row_start, rstarts, sizeof(rstarts)) ||
		memcmp(m->row_index, rindex, sizeof(rindex)) ||
		memcmp(m->column_index, cindex, sizeof(cindex)) ||
		memcmp(m->value, value, sizeof(value)) ||
		memcmp(m->row_value, rvalue, sizeof(rvalue)))
		return 0;
	for (j = 0; j <= m->n; j++)
		if (m->column_start[j] != (j < 5 ? starts[j] : 9))
			return 0;
	for (j = 0; j < m->n; j++)
		if (m->objective[j] != (j < 4 ? costs[j] : 3.) ||
			m->bounds[j].b_type != optm_BOUND_T_BS ||
			m->bounds[j].lb != 0. || m->bounds[j].ub != 3.)
			return 0;
	for (i = 0; i < 3; i++) {
		if (m->constraints[i].rhs != (i == 2 ? 1. : 2.) ||
			m->constraints[i].type != optm_CONS_T_GE)
			return 0;
		for (j = 0; j < m->n; j++)
			if (m->constraints[i].coef[j] !=
				(j < 4 ? coefficients[i][j] : 0.))
				return 0;
	}
	return 1;
}

static int reductions(int wide)
{
	struct lp_Model *m = fixture(wide ? 512 : 4), *reduced;
	struct lp_Model saved;
	struct lp_simplex_Options options;
	struct lp_simplex_Result reference, result;
	struct optm_VariableBound *bounds;
	unsigned char rows[3], columns[512];
	int rd[3], cd[512], removed, fixed, scenario, i, j, k, changed;
	double rhs[3], x[512];
	size_t work;
	unsigned long masks[9] = {0, 1UL << 1, 1UL << 2, 1UL << 3,
		(1UL << 2) | (1UL << 4), (1UL << 2) | (1UL << 5),
		(1UL << 2) | (1UL << 6), (1UL << 2) | (1UL << 7),
		(1UL << 2) | (1UL << 8)};

	CHECK(m);
	memcpy(&saved, m, sizeof(saved));
	bounds = malloc((size_t)m->n * sizeof(*bounds));
	CHECK(bounds);
	lp_simplex_default_options(&options, lp_simplex_ALGORITHM_DUAL_REVISED);
	options.presolve = 0;
	CHECK(!lp_simplex_solve(m, &options, x, &reference));
	CHECK(reference.status == lp_simplex_Success);
	for (scenario = 0; scenario < (wide ? 1 : 11); scenario++) {
		memcpy(bounds, m->bounds, (size_t)m->n * sizeof(*bounds));
		memset(rows, 1, sizeof(rows));
		memset(columns, 1, sizeof(columns));
		memset(cd, 0, sizeof(cd));
		for (i = 0; i < 3; i++) {
			rd[i] = m->row_start[i + 1] - m->row_start[i];
			rhs[i] = m->constraints[i].rhs;
			for (k = m->row_start[i]; k < m->row_start[i + 1]; k++)
				cd[m->column_index[k]]++;
		}
		requests = bad_free = 0;
		owned = peak = 0;
		reject = scenario < 9 ? masks[scenario] : 0;
		work = scenario == 9 ? 1 : scenario == 10 ? 0 : 100000;
		storage_dominance(m, bounds, rhs, rows, columns, rd, cd, &work,
			&removed, &fixed);
		CHECK(work <= (size_t)(scenario == 9	? 1
				       : scenario == 10 ? 0
							: 100000));
		CHECK(!owned && !bad_free && peak <= DOMINANCE_BYTES);
		if (scenario == 1 || scenario == 10)
			CHECK(!removed && !fixed &&
				work == (size_t)(scenario == 10 ? 0 : 100000));
		if (scenario == 9)
			CHECK(!work && requests == 1);
		if (scenario == 10)
			CHECK(!requests);
		if (!wide && scenario == 0)
			CHECK(requests == 3);
		if (!wide && (scenario == 2 || scenario == 3))
			CHECK(requests == 8);
		if (!wide && scenario == 4)
			CHECK(requests == 4);
		if (!wide && scenario >= 5 && scenario <= 8)
			CHECK(requests == 8);
		if (wide && DOMINANCE_BYTES == 256)
			CHECK(requests == 1 && peak == 3);
		changed = 0;
		for (i = 0; i < 3; i++)
			changed += !rows[i];
		CHECK(removed == changed);
		changed = 0;
		for (j = 0; j < m->n; j++) {
			int degree = 0;
			changed += bounds[j].ub == bounds[j].lb;
			for (k = m->column_start[j]; k < m->column_start[j + 1];
				k++)
				degree += rows[m->row_index[k]] != 0;
			CHECK(cd[j] == degree);
		}
		CHECK(fixed == changed && unchanged(m, &saved));
		reduced = lp_model_create(3, m->n);
		CHECK(reduced);
		memcpy(reduced->objective, m->objective,
			(size_t)m->n * sizeof(double));
		memcpy(reduced->bounds, bounds, (size_t)m->n * sizeof(*bounds));
		for (i = 0; i < 3; i++)
			if (rows[i]) {
				reduced->constraints[i].rhs = rhs[i];
				reduced->constraints[i].type =
					m->constraints[i].type;
				memcpy(reduced->constraints[i].coef,
					m->constraints[i].coef,
					(size_t)m->n * sizeof(double));
			}
		CHECK(!lp_simplex_solve(reduced, &options, x, &result));
		CHECK(result.status == lp_simplex_Success);
		CHECK(fabs(result.objective - reference.objective) < 1e-7);
		CHECK(simplex_primal_infeasibility(m, x) <= 1e-7);
		lp_model_free(reduced);
	}
	free(bounds);
	lp_model_free(m);
	return 0;
}

int main(void)
{
	return sizes() || empty() || reductions(0) || reductions(1);
}
