/**
 * @file test_problem_storage.c
 * @brief Canonical/CSC/scaling production ownership under bounded allocation.
 */

#include "simplex_stage.h"
#include "simplex_solve_context.h"
#include "simplex_execution.h"
#include <lp_simplex/lp_simplex.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>

static void *live[128];
static int active, bad_free, requests, fail_at;
static size_t largest;

void *problem_storage_malloc(size_t bytes)
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
	memset(p, 0xa5, bytes);
	for (i = 0; i < 128; i++)
		if (!live[i]) {
			live[i] = p;
			active++;
			return p;
		}
	free(p);
	return NULL;
}

void problem_storage_free(void *p)
{
	int i;
	if (!p)
		return;
	for (i = 0; i < 128; i++)
		if (live[i] == p) {
			live[i] = NULL;
			active--;
			free(p);
			return;
		}
	bad_free++;
}

#define lp_simplex_malloc problem_storage_malloc
#define lp_simplex_free problem_storage_free
#define simplex_csc_from_model storage_csc_from_model
#define simplex_csc_destroy storage_csc_destroy
#define simplex_csc_column_dot storage_csc_column_dot
#define simplex_csc_column_to_dense storage_csc_column_to_dense
#define simplex_csc_column_axpy storage_csc_column_axpy
#define simplex_problem_mark_changed storage_problem_mark_changed
#define simplex_problem_from_model storage_problem_from_model
#define simplex_problem_create_sparse storage_problem_create_sparse
#define simplex_problem_create_with_pattern storage_problem_create_with_pattern
#define simplex_problem_destroy storage_problem_destroy
#define simplex_problem_free storage_problem_free
#define simplex_problem_primal_infeasibility \
	storage_problem_primal_infeasibility
#define simplex_problem_primal_summary storage_problem_primal_summary
#define simplex_problem_dual_infeasibility storage_problem_dual_infeasibility
#define simplex_scaling_create_problem storage_scaling_create_problem
#define simplex_scaling_recover_primal storage_scaling_recover_primal
#define simplex_scaling_recover_dual storage_scaling_recover_dual
#define simplex_scaling_destroy storage_scaling_destroy
#ifndef PROBLEM_STORAGE_CSC_SOURCE
#define PROBLEM_STORAGE_CSC_SOURCE "../src/matrix/simplex_csc.c"
#endif
#include PROBLEM_STORAGE_CSC_SOURCE
#ifndef PROBLEM_STORAGE_PROBLEM_SOURCE
#define PROBLEM_STORAGE_PROBLEM_SOURCE "../src/runtime/simplex_problem.c"
#endif
#include PROBLEM_STORAGE_PROBLEM_SOURCE
#ifndef PROBLEM_STORAGE_SCALING_SOURCE
#define PROBLEM_STORAGE_SCALING_SOURCE "../src/runtime/simplex_scaling.c"
#endif
#include PROBLEM_STORAGE_SCALING_SOURCE
#undef lp_simplex_malloc
#undef lp_simplex_free

#define CHECK(c) \
	do { \
		if (!(c)) { \
			fprintf(stderr, "problem storage line %d: %s\n", \
				__LINE__, #c); \
			return 1; \
		} \
	} while (0)

static struct lp_Model *fixture(int cached, int scaled)
{
	struct lp_Model *m = lp_model_create(2, 2);
	if (!m)
		return NULL;
	m->objective[0] = scaled == 2 ? nextafter(0., 1.) : 3.;
	m->objective[1] = -2.;
	m->constraints[0].coef[0] = scaled ? 1024. : 1.;
	m->constraints[1].coef[0] = 1.;
	m->constraints[1].coef[1] = scaled ? 1. / 1024. : 1.;
	m->constraints[0].type = optm_CONS_T_EQ;
	m->constraints[1].type = optm_CONS_T_GE;
	m->constraints[0].rhs = scaled ? 2048. : 2.;
	m->constraints[1].rhs = -2.;
	m->bounds[0].b_type = optm_BOUND_T_BS;
	m->bounds[0].lb = -4.;
	m->bounds[0].ub = 8.;
	m->bounds[1].b_type = optm_BOUND_T_LO;
	m->bounds[1].lb = .5;
	if (cached && lp_model_build_sparse(m)) {
		lp_model_free(m);
		return NULL;
	}
	return m;
}

static int fixture_unchanged(const struct lp_Model *m, int scaled)
{
	return m->objective[0] == (scaled == 2 ? nextafter(0., 1.) : 3.) &&
	       m->objective[1] == -2. &&
	       m->constraints[0].coef[0] == (scaled ? 1024. : 1.) &&
	       m->constraints[0].coef[1] == 0. &&
	       m->constraints[1].coef[0] == 1. &&
	       m->constraints[1].coef[1] == (scaled ? 1. / 1024. : 1.) &&
	       m->constraints[0].rhs == (scaled ? 2048. : 2.) &&
	       m->constraints[1].rhs == -2. &&
	       m->constraints[0].type == optm_CONS_T_EQ &&
	       m->constraints[1].type == optm_CONS_T_GE &&
	       m->bounds[0].b_type == optm_BOUND_T_BS &&
	       m->bounds[0].lb == -4. && m->bounds[0].ub == 8. &&
	       m->bounds[1].b_type == optm_BOUND_T_LO && m->bounds[1].lb == .5;
}

static int arrays_match(
	const struct simplex_CscMatrix *a, const struct simplex_CscMatrix *b)
{
	return a->rows == b->rows && a->columns == b->columns &&
	       a->nonzeros == b->nonzeros &&
	       !memcmp(a->column_start, b->column_start, 3 * sizeof(int)) &&
	       !memcmp(a->row_start, b->row_start, 3 * sizeof(int)) &&
	       !memcmp(a->row_index, b->row_index, 3 * sizeof(int)) &&
	       !memcmp(a->column_index, b->column_index, 3 * sizeof(int)) &&
	       !memcmp(a->value, b->value, 3 * sizeof(double)) &&
	       !memcmp(a->row_value, b->row_value, 3 * sizeof(double));
}

static int csc_and_model(int cached)
{
	struct lp_Model *m = fixture(cached, 0), saved;
	struct simplex_CscMatrix matrix, expected;
	struct simplex_Problem p;
	int allocations, failure, retained;
	CHECK(m);
	saved = *m;
	requests = fail_at = 0;
	CHECK(!storage_csc_from_model(m, &expected));
	CHECK(fixture_unchanged(m, 0));
	allocations = requests;
	retained = active;
	CHECK(allocations == (cached ? 0 : 7));
	if (cached)
		CHECK(!expected.owns_storage && expected.value == m->value);
	else
		CHECK(expected.owns_storage == SIMPLEX_CSC_OWNS_ALL);
	for (failure = 1; failure <= allocations; failure++) {
		requests = 0;
		fail_at = failure;
		CHECK(storage_csc_from_model(m, &matrix));
		storage_csc_destroy(&matrix);
		CHECK(active == retained && !bad_free &&
			!memcmp(&saved, m, sizeof(saved)) &&
			fixture_unchanged(m, 0));
	}
	fail_at = 0;
	CHECK(!storage_csc_from_model(m, &matrix));
	CHECK(arrays_match(&expected, &matrix));
	storage_csc_destroy(&matrix);
	requests = 0;
	CHECK(!storage_problem_from_model(&p, m));
	allocations = requests;
	CHECK(allocations == (cached ? 1 : 8) &&
		arrays_match(&expected, &p.matrix));
	CHECK(p.rhs == p.objective + 2 &&
		p.row_type == (unsigned char *)(p.bounds + 2));
	storage_problem_destroy(&p);
	for (failure = 1; failure <= allocations; failure++) {
		requests = 0;
		fail_at = failure;
		CHECK(storage_problem_from_model(&p, m));
		CHECK(!p.vector_storage && !p.matrix.column_start);
		storage_problem_destroy(&p);
		CHECK(active == retained && !bad_free &&
			!memcmp(&saved, m, sizeof(saved)) &&
			fixture_unchanged(m, 0));
	}
	fail_at = 0;
	storage_csc_destroy(&expected);
	lp_model_free(m);
	CHECK(!active && !bad_free);
	printf("canonical cached=%d allocations=%d\n", cached, allocations);
	return 0;
}

static int factories(void)
{
	struct simplex_Problem *p, *q;
	struct simplex_CscMatrix pattern, saved;
	int starts[] = {0, 1, 3}, row_starts[] = {0, 2, 3}, rows[] = {0, 1, 1},
	    columns[] = {0, 0, 1};
	double values[] = {1., 2., 3.};
	int allocations, failure;
	requests = fail_at = 0;
	p = storage_problem_create_sparse(2, 2, 3);
	CHECK(p);
	allocations = requests;
	CHECK(allocations == 8);
	storage_problem_free(p);
	for (failure = 1; failure <= allocations; failure++) {
		requests = 0;
		fail_at = failure;
		p = storage_problem_create_sparse(2, 2, 3);
		CHECK(!p && !active && !bad_free);
	}
	memset(&pattern, 0, sizeof(pattern));
	pattern.rows = pattern.columns = 2;
	pattern.nonzeros = 3;
	pattern.column_start = starts;
	pattern.row_start = row_starts;
	pattern.row_index = rows;
	pattern.column_index = columns;
	pattern.value = pattern.row_value = values;
	saved = pattern;
	requests = fail_at = 0;
	q = storage_problem_create_with_pattern(&pattern);
	CHECK(q);
	allocations = requests;
	CHECK(allocations == 3);
	CHECK(q->matrix.column_start == starts &&
		q->matrix.row_start == row_starts &&
		q->matrix.row_index == rows &&
		q->matrix.column_index == columns &&
		q->matrix.owns_storage == SIMPLEX_CSC_OWNS_VALUES &&
		q->matrix.value != values &&
		q->matrix.row_value == q->matrix.value + 3);
	storage_problem_free(q);
	for (failure = 1; failure <= allocations; failure++) {
		requests = 0;
		fail_at = failure;
		q = storage_problem_create_with_pattern(&pattern);
		CHECK(!q && !active && !bad_free &&
			!memcmp(&pattern, &saved, sizeof(saved)));
	}
	fail_at = 0;
	CHECK(values[0] == 1. && values[1] == 2. && values[2] == 3.);
	return 0;
}

static int scaling_result(const struct simplex_Problem *p,
	struct simplex_Problem *q, const struct simplex_Scaling *scaling,
	int expected_active)
{
	int i, j, k;
	double x[] = {2., 3.}, dual[] = {.25, 1.5}, expected[2];
	CHECK(scaling->active == expected_active && q->report == p->report);
	CHECK(q->matrix.column_start == p->matrix.column_start &&
		q->matrix.row_start == p->matrix.row_start &&
		q->matrix.row_index == p->matrix.row_index &&
		q->matrix.column_index == p->matrix.column_index);
	if (expected_active) {
		CHECK(q->owns_vectors &&
			q->matrix.owns_storage == SIMPLEX_CSC_OWNS_VALUES &&
			q->execution.context == p->execution.context &&
			q->execution.model != p->execution.model);
	} else {
		CHECK(!q->owns_vectors && !q->matrix.owns_storage &&
			!q->vector_storage && q->objective == p->objective &&
			q->bounds == p->bounds &&
			!memcmp(&q->execution, &p->execution,
				sizeof(q->execution)));
	}
	for (j = 0; j < 2; j++) {
		CHECK(q->objective[j] == p->objective[j] * scaling->column[j]);
		if (p->bounds[j].b_type == optm_BOUND_T_LO ||
			p->bounds[j].b_type == optm_BOUND_T_BS)
			CHECK(q->bounds[j].lb ==
				p->bounds[j].lb / scaling->column[j]);
		if (p->bounds[j].b_type == optm_BOUND_T_BS)
			CHECK(q->bounds[j].ub ==
				p->bounds[j].ub / scaling->column[j]);
		x[j] /= scaling->column[j];
	}
	for (i = 0; i < 2; i++) {
		CHECK(q->rhs[i] == p->rhs[i] * scaling->row[i] &&
			q->row_type[i] == p->row_type[i]);
		expected[i] = dual[i] * scaling->row[i];
	}
	for (j = 0; j < 2; j++)
		for (k = p->matrix.column_start[j];
			k < p->matrix.column_start[j + 1]; k++) {
			i = p->matrix.row_index[k];
			CHECK(q->matrix.value[k] == p->matrix.value[k] *
							    scaling->row[i] *
							    scaling->column[j]);
		}
	for (i = 0; i < 2; i++)
		for (k = p->matrix.row_start[i]; k < p->matrix.row_start[i + 1];
			k++) {
			j = p->matrix.column_index[k];
			CHECK(q->matrix.row_value[k] ==
				p->matrix.row_value[k] * scaling->row[i] *
					scaling->column[j]);
		}
	CHECK(storage_problem_primal_infeasibility(q, x) == 0.);
	storage_scaling_recover_primal(scaling, x);
	storage_scaling_recover_dual(scaling, dual);
	CHECK(x[0] == 2. && x[1] == 3. && dual[0] == expected[0] &&
		dual[1] == expected[1]);
	return 0;
}

static int scaling_failures(int enabled)
{
	struct lp_Model *m = fixture(1, enabled);
	struct simplex_Problem p, saved, *q;
	struct simplex_Scaling scaling, zero;
	struct simplex_SolveContext context;
	struct simplex_Execution history;
	struct simplex_StageReport report;
	struct simplex_CscMatrix frozen;
	int starts[6], indices[6], i;
	double values[6];
	unsigned char vectors[512];
	size_t bytes =
		4 * sizeof(double) + 2 * sizeof(struct optm_VariableBound) + 2;
	int allocations, failure, retained, models;
	CHECK(m && bytes <= sizeof(vectors));
	fail_at = 0;
	CHECK(!storage_problem_from_model(&p, m));
	simplex_solve_context_init(&context, &history, &p.execution);
	p.report = &report;
	p.revision = 7;
	saved = p;
	frozen = p.matrix;
	frozen.column_start = starts;
	frozen.row_start = starts + 3;
	frozen.row_index = indices;
	frozen.column_index = indices + 3;
	frozen.value = values;
	frozen.row_value = values + 3;
	for (i = 0; i < 3; i++) {
		starts[i] = p.matrix.column_start[i];
		starts[3 + i] = p.matrix.row_start[i];
		indices[i] = p.matrix.row_index[i];
		indices[3 + i] = p.matrix.column_index[i];
		values[i] = p.matrix.value[i];
		values[3 + i] = p.matrix.row_value[i];
	}
	memcpy(vectors, p.vector_storage, bytes);
	memset(&zero, 0, sizeof(zero));
	retained = active;
	requests = 0;
	q = storage_scaling_create_problem(&p, &scaling);
	CHECK(q);
	allocations = requests;
	CHECK(allocations == (enabled == 2 ? 6 : enabled ? 5 : 3));
	CHECK(!scaling_result(&p, q, &scaling, enabled == 1));
	CHECK(arrays_match(&p.matrix, &frozen) &&
		fixture_unchanged(m, enabled));
	storage_problem_free(q);
	storage_scaling_destroy(&scaling);
	for (failure = 1; failure <= allocations; failure++) {
		models = context.models;
		requests = 0;
		fail_at = failure;
		q = storage_scaling_create_problem(&p, &scaling);
		CHECK(!q && !memcmp(&scaling, &zero, sizeof(zero)) &&
			active == retained && !bad_free &&
			context.models == models);
		CHECK(!memcmp(&p, &saved, sizeof(p)) &&
			!memcmp(p.vector_storage, vectors, bytes));
		CHECK(arrays_match(&p.matrix, &frozen) &&
			fixture_unchanged(m, enabled));
		storage_scaling_destroy(&scaling);
		fail_at = 0;
		q = storage_scaling_create_problem(&p, &scaling);
		CHECK(q);
		CHECK(!scaling_result(&p, q, &scaling, enabled == 1));
		CHECK(arrays_match(&p.matrix, &frozen) &&
			fixture_unchanged(m, enabled));
		storage_problem_free(q);
		storage_scaling_destroy(&scaling);
	}
	storage_problem_destroy(&p);
	simplex_solve_context_destroy(&context);
	lp_model_free(m);
	CHECK(!active && !bad_free);
	printf("scaling active=%d allocations=%d\n", enabled, allocations);
	return 0;
}

#ifndef PROBLEM_STORAGE_OLD
static int sizes(void)
{
	struct problem_VectorStorage vector;
	struct simplex_CscStorage csc;
	struct simplex_Problem p, zero;
	struct simplex_CscMatrix pattern;
	struct lp_Model model;
	size_t bytes, numeric = 4 * sizeof(double),
		      bounds = 2 * sizeof(struct optm_VariableBound);
	CHECK(problem_vector_storage(2, 2, numeric + bounds + 2, &vector));
	CHECK(vector.numeric == numeric && vector.bounds == bounds &&
		vector.total == numeric + bounds + 2);
	CHECK(!problem_vector_storage(2, 2, numeric + bounds + 1, &vector) &&
		!vector.total && !vector.numeric && !vector.bounds);
	CHECK(!problem_vector_storage(2, 2, numeric + bounds - 1, &vector));
	CHECK(!problem_vector_storage(2, 2, 3 * sizeof(double), &vector));
	CHECK(!problem_vector_storage(-1, 2, (size_t)-1, &vector));
	CHECK(!problem_vector_storage(
		INT_MAX, INT_MAX, (size_t)4294967295UL, &vector));
	CHECK(!simplex_csc_storage(INT_MAX, 0, 0, (size_t)4294967295UL, &csc));
	CHECK(!simplex_csc_storage(0, INT_MAX, 0, (size_t)4294967295UL, &csc));
	CHECK(!simplex_csc_storage(1, 1, INT_MAX, (size_t)4294967295UL, &csc));
	CHECK(simplex_csc_storage(0, 0, 0, (size_t)-1, &csc) &&
		csc.rows == sizeof(int) && csc.columns == sizeof(int));
	CHECK(!scaling_storage(1, 2, 3 * sizeof(double) - 1, &bytes) && !bytes);
	CHECK(scaling_storage(1, 2, 3 * sizeof(double), &bytes) &&
		bytes == 3 * sizeof(double));
	CHECK(!scaling_storage(INT_MAX, 1, (size_t)4294967295UL, &bytes));
	if (sizeof(size_t) >= 8)
		CHECK(scaling_storage(INT_MAX, 1, (size_t)-1, &bytes));
	CHECK(!storage_problem_create_sparse(1, 1, -1));
	CHECK(!storage_problem_create_sparse(-1, 1, 0));
	memset(&pattern, 0, sizeof(pattern));
	pattern.rows = pattern.columns = 1;
	pattern.nonzeros = -1;
	CHECK(!storage_problem_create_with_pattern(&pattern));
	memset(&model, 0, sizeof(model));
	model.m = -1;
	memset(&zero, 0, sizeof(zero));
	CHECK(storage_problem_from_model(&p, &model));
	CHECK(!memcmp(&p, &zero, sizeof(p)));
	CHECK(!requests && !active && !bad_free);
	return 0;
}
#endif

static int extreme(const char *mode)
{
	struct simplex_Scaling scaling;
	struct simplex_CscMatrix matrix;
	struct lp_Model model;
	requests = fail_at = 0;
	if (!strcmp(mode, "scaling")) {
		CHECK(scaling_allocate(&scaling, INT_MAX, 1));
		storage_scaling_destroy(&scaling);
	} else if (!strcmp(mode, "csc-row")) {
		memset(&model, 0, sizeof(model));
		model.m = INT_MAX;
		CHECK(storage_csc_from_model(&model, &matrix));
		storage_csc_destroy(&matrix);
	} else
		return 1;
	CHECK(!active && !bad_free);
	return 0;
}

int main(int argc, char **argv)
{
	if (argc > 1)
		return extreme(argv[1]);
#ifndef PROBLEM_STORAGE_OLD
	CHECK(!sizes());
#endif
	return csc_and_model(0) || csc_and_model(1) || factories() ||
	       scaling_failures(0) || scaling_failures(1) ||
	       scaling_failures(2) || extreme("scaling") || extreme("csc-row");
}
