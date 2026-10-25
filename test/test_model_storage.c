/**
 * @file test_model_storage.c
 * @brief Production model bodies under a bounded allocator: no giant real
 * buffers.
 */

#include <lp_simplex/lp_simplex.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void *live[64];
static int requests, fail_at, active, bad_free;
static size_t cap = 1024 * 1024, largest;

void *model_storage_malloc(size_t bytes)
{
	void *p;
	int i;
	requests++;
	if (bytes > largest)
		largest = bytes;
	if (requests == fail_at || bytes > cap)
		return NULL;
	p = malloc(bytes);
	if (!p)
		return NULL;
	for (i = 0; i < 64; i++)
		if (!live[i]) {
			live[i] = p;
			active++;
			return p;
		}
	free(p);
	return NULL;
}

void model_storage_free(void *p)
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

void model_storage_destroy(struct lp_Model *model);
#define lp_simplex_malloc model_storage_malloc
#define lp_simplex_free model_storage_free
#define lp_model_create model_storage_create
#define lp_model_create_sparse model_storage_create_sparse
#define lp_model_build_sparse model_storage_build_sparse
#define lp_model_free model_storage_destroy
#ifndef MODEL_STORAGE_SOURCE
#define MODEL_STORAGE_SOURCE "../src/core/model.c"
#endif
#include MODEL_STORAGE_SOURCE

#undef lp_simplex_malloc
#undef lp_simplex_free
#undef lp_model_create
#undef lp_model_create_sparse
#undef lp_model_build_sparse
#undef lp_model_free
#include "simplex_model_storage.h"

#define CHECK(c) \
	do { \
		if (!(c)) { \
			fprintf(stderr, "model storage line %d: %s\n", \
				__LINE__, #c); \
			return 1; \
		} \
	} while (0)

static int sizes(void)
{
	struct simplex_ModelStorage storage;
	size_t bytes;
	CHECK(!simplex_model_storage(0, 1, 0, 1, (size_t)-1, &storage));
	CHECK(!simplex_model_storage(1, 1, -1, 0, (size_t)-1, &storage));
	CHECK(!simplex_model_storage(
		INT_MAX, INT_MAX, 0, 1, (size_t)-1, &storage));
	CHECK(!simplex_model_storage(
		65536, 65536, 0, 1, (size_t)4294967295UL, &storage));
	CHECK(!simplex_model_storage(
		1, INT_MAX, 0, 0, (size_t)4294967295UL, &storage));
	CHECK(!simplex_model_array_bytes(
		      (size_t)-1, sizeof(double), (size_t)-1, &bytes) &&
		!bytes);
	if (sizeof(size_t) >= 8) {
		CHECK(simplex_model_storage(
			65536, 65536, 0, 1, (size_t)-1, &storage));
		CHECK(storage.coefficients ==
			(size_t)65536 * 65536 * sizeof(double));
		CHECK(simplex_model_storage(
			1, INT_MAX, 0, 0, (size_t)-1, &storage));
		CHECK(storage.column_start ==
			((size_t)INT_MAX + 1) * sizeof(int));
	}
	requests = 0;
	CHECK(!model_storage_create(INT_MAX, INT_MAX) && !requests && !active);
	largest = 0;
	CHECK(!model_storage_create_sparse(1, INT_MAX, 0));
	CHECK(!active && !bad_free);
	CHECK(!model_storage_create_sparse(INT_MAX, 1, 0));
	CHECK(!active && !bad_free);
	return 0;
}

static int constructors(void)
{
	struct lp_Model *model;
	int kind, failure, calls, i;
	for (kind = 0; kind < 3; kind++) {
		calls = kind == 0 ? 5 : kind == 1 ? 6 : 10;
		for (failure = 1; failure <= calls + 1; failure++) {
			requests = 0;
			fail_at = failure <= calls ? failure : 0;
			model = kind == 0 ? model_storage_create(3, 4)
					  : model_storage_create_sparse(
						    3, 4, kind == 1 ? 0 : 5);
			if (fail_at)
				CHECK(!model && !active);
			else {
				CHECK(model && !model->nnz &&
					model->objective_sense == 1.);
				for (i = 0; i < 4; i++)
					CHECK(model->bounds[i].b_type ==
							optm_BOUND_T_LO &&
						model->bounds[i].lb == 0. &&
						model->objective[i] == 0.);
				for (i = 0; i < 3; i++)
					CHECK(model->constraints[i].coef ==
						(kind == 0 ? model->coefficients +
									(size_t)i *
										4
							   : NULL));
				model_storage_destroy(model);
				CHECK(!active);
			}
			CHECK(!bad_free);
		}
	}
	fail_at = 0;
	return 0;
}

static int sparse_matches(const struct lp_Model *model)
{
	int row, column, k, count = 0;
	for (row = 0; row < model->m; row++)
		for (column = 0; column < model->n; column++)
			count += model->constraints[row].coef[column] != 0.;
	if (model->nnz != count || model->column_start[model->n] != count ||
		model->row_start[model->m] != count)
		return 0;
	for (column = 0; column < model->n; column++)
		for (k = model->column_start[column];
			k < model->column_start[column + 1]; k++)
			if (model->value[k] !=
				model->constraints[model->row_index[k]]
					.coef[column])
				return 0;
	for (row = 0; row < model->m; row++)
		for (k = model->row_start[row]; k < model->row_start[row + 1];
			k++)
			if (model->row_value[k] !=
				model->constraints[row]
					.coef[model->column_index[k]])
				return 0;
	return 1;
}

static int rollback(void)
{
	struct lp_Model *model = model_storage_create(3, 4), before;
	int failure, existing;
	int starts[5], rows[4], columns[4], row_starts[4];
	double values[4], row_values[4];
	CHECK(model);
	model->constraints[0].coef[0] = 1.;
	model->constraints[1].coef[2] = 2.;
	model->constraints[2].coef[1] = 3.;
	model->constraints[2].coef[3] = 4.;
	CHECK(!model_storage_build_sparse(model) && sparse_matches(model));
	before = *model;
	memcpy(starts, model->column_start, sizeof(starts));
	memcpy(rows, model->row_index, sizeof(rows));
	memcpy(columns, model->column_index, sizeof(columns));
	memcpy(row_starts, model->row_start, sizeof(row_starts));
	memcpy(values, model->value, sizeof(values));
	memcpy(row_values, model->row_value, sizeof(row_values));
	model->constraints[0].coef[2] = 5.;
	existing = active;
	for (failure = 1; failure <= 7; failure++) {
		requests = 0;
		fail_at = failure;
		CHECK(model_storage_build_sparse(model) != 0);
		CHECK(active == existing && !bad_free &&
			!memcmp(&before, model, sizeof(before)));
		CHECK(!memcmp(starts, model->column_start, sizeof(starts)));
		CHECK(!memcmp(rows, model->row_index, sizeof(rows)));
		CHECK(!memcmp(columns, model->column_index, sizeof(columns)));
		CHECK(!memcmp(
			row_starts, model->row_start, sizeof(row_starts)));
		CHECK(!memcmp(values, model->value, sizeof(values)));
		CHECK(!memcmp(
			row_values, model->row_value, sizeof(row_values)));
	}
	fail_at = 0;
	CHECK(!model_storage_build_sparse(model) && model->nnz == 5 &&
		sparse_matches(model));
	memset(model->coefficients, 0, 12 * sizeof(double));
	before = *model;
	for (failure = 1; failure <= 3; failure++) {
		requests = 0;
		fail_at = failure;
		CHECK(model_storage_build_sparse(model) != 0);
		CHECK(active == existing &&
			!memcmp(&before, model, sizeof(before)) && !bad_free);
	}
	fail_at = 0;
	CHECK(!model_storage_build_sparse(model) && sparse_matches(model) &&
		!model->nnz);
	CHECK(!model->row_index && !model->value && !model->column_index &&
		!model->row_value);
	model_storage_destroy(model);
	CHECK(!active && !bad_free);
	return 0;
}

int main(void)
{
	return sizes() || constructors() || rollback();
}
