/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file model.c
 * @brief Dense and sparse linear-program models and their storage lifecycle.
 */

#include "utils.h"
#include "simplex_model_storage.h"
#include <limits.h>

#include <lp_simplex/model.h>

static void model_free_sparse(struct lp_Model *model)
{
	lp_simplex_free(model->column_start);
	lp_simplex_free(model->row_index);
	lp_simplex_free(model->value);
	lp_simplex_free(model->row_start);
	lp_simplex_free(model->column_index);
	lp_simplex_free(model->row_value);
	model->column_start = model->row_index = model->row_start =
		model->column_index = NULL;
	model->value = model->row_value = NULL;
	model->nnz = 0;
}

static struct lp_Model *model_allocate(int m, int n, int capacity, int dense)
{
	struct simplex_ModelStorage storage;
	struct lp_Model *model;
	int i;

	if (!simplex_model_storage(m, n, capacity, dense, (size_t)-1, &storage))
		return NULL;
	model = (struct lp_Model *)lp_simplex_malloc(sizeof(*model));
	if (!model)
		return NULL;
	lp_simplex_memset(model, 0, sizeof(*model));
	model->m = m;
	model->n = n;
	model->objective_sense = 1.;
	model->objective = (double *)lp_simplex_malloc(storage.objective);
	model->constraints = (struct optm_LinearConstraint *)lp_simplex_malloc(
		storage.constraints);
	model->bounds =
		(struct optm_VariableBound *)lp_simplex_malloc(storage.bounds);
	if (dense)
		model->coefficients =
			(double *)lp_simplex_malloc(storage.coefficients);
	else {
		model->column_start =
			(int *)lp_simplex_malloc(storage.column_start);
		model->row_start = (int *)lp_simplex_malloc(storage.row_start);
		if (capacity) {
			model->row_index =
				(int *)lp_simplex_malloc(storage.indices);
			model->column_index =
				(int *)lp_simplex_malloc(storage.indices);
			model->value =
				(double *)lp_simplex_malloc(storage.values);
			model->row_value =
				(double *)lp_simplex_malloc(storage.values);
		}
	}
	if (!model->objective || !model->constraints || !model->bounds ||
		(dense && !model->coefficients) ||
		(!dense && (!model->column_start || !model->row_start ||
				   (capacity && (!model->row_index ||
							!model->column_index ||
							!model->value ||
							!model->row_value))))) {
		lp_model_free(model);
		return NULL;
	}
	lp_simplex_memset(model->objective, 0, storage.objective);
	lp_simplex_memset(model->constraints, 0, storage.constraints);
	lp_simplex_memset(model->bounds, 0, storage.bounds);
	if (dense)
		lp_simplex_memset(model->coefficients, 0, storage.coefficients);
	else {
		lp_simplex_memset(model->column_start, 0, storage.column_start);
		lp_simplex_memset(model->row_start, 0, storage.row_start);
	}
	for (i = 0; i < m; i++)
		model->constraints[i].coef =
			dense ? model->coefficients + (size_t)i * (size_t)n
			      : NULL;
	for (i = 0; i < n; i++) {
		model->bounds[i].lb = 0.;
		model->bounds[i].ub = __lp_simplex_INF__;
		model->bounds[i].b_type = optm_BOUND_T_LO;
		model->bounds[i].v_type = optm_VAR_T_REAL;
	}
	return model;
}

struct lp_Model *lp_model_create(int m, int n)
{
	return model_allocate(m, n, 0, 1);
}

struct lp_Model *lp_model_create_sparse(int m, int n, int nonzero_capacity)
{
	return model_allocate(m, n, nonzero_capacity, 0);
}

int lp_model_build_sparse(struct lp_Model *model)
{
	struct simplex_ModelStorage storage;
	struct lp_Model staged;
	int i, j, k, nonzeros = 0;
	int *next = NULL;

	if (!model || !model->coefficients || !model->constraints ||
		!simplex_model_storage(
			model->m, model->n, 0, 1, (size_t)-1, &storage))
		return -1;
	for (i = 0; i < model->m; i++) {
		if (!model->constraints[i].coef)
			return -1;
		for (j = 0; j < model->n; j++)
			if (model->constraints[i].coef[j] != 0.) {
				if (nonzeros == INT_MAX)
					return -1;
				nonzeros++;
			}
	}
	if (!simplex_model_storage(
		    model->m, model->n, nonzeros, 1, (size_t)-1, &storage))
		return -1;
	lp_simplex_memset(&staged, 0, sizeof(staged));
	staged.column_start = (int *)lp_simplex_malloc(storage.column_start);
	staged.row_start = (int *)lp_simplex_malloc(storage.row_start);
	if (nonzeros) {
		staged.row_index = (int *)lp_simplex_malloc(storage.indices);
		staged.column_index = (int *)lp_simplex_malloc(storage.indices);
		staged.value = (double *)lp_simplex_malloc(storage.values);
		staged.row_value = (double *)lp_simplex_malloc(storage.values);
	}
	next = (int *)lp_simplex_malloc(storage.row_cursor);
	if (!staged.column_start || !staged.row_start || !next ||
		(nonzeros && (!staged.row_index || !staged.column_index ||
				     !staged.value || !staged.row_value))) {
		model_free_sparse(&staged);
		lp_simplex_free(next);
		return -1;
	}
	k = 0;
	for (j = 0; j < model->n; j++) {
		staged.column_start[j] = k;
		for (i = 0; i < model->m; i++) {
			double coefficient = model->constraints[i].coef[j];
			if (coefficient != 0.) {
				staged.row_index[k] = i;
				staged.value[k++] = coefficient;
			}
		}
	}
	staged.column_start[model->n] = k;
	lp_simplex_memset(staged.row_start, 0, storage.row_start);
	for (k = 0; k < nonzeros; k++)
		staged.row_start[staged.row_index[k] + 1]++;
	for (i = 0; i < model->m; i++) {
		staged.row_start[i + 1] += staged.row_start[i];
		next[i] = staged.row_start[i];
	}
	for (j = 0; j < model->n; j++)
		for (k = staged.column_start[j]; k < staged.column_start[j + 1];
			k++) {
			i = staged.row_index[k];
			staged.column_index[next[i]] = j;
			staged.row_value[next[i]++] = staged.value[k];
		}
	lp_simplex_free(next);
	/* Commit only a complete replacement. Allocation/count failures leave
	 * the previous sparse arrays and nnz available to the caller. */
	model_free_sparse(model);
	model->column_start = staged.column_start;
	model->row_start = staged.row_start;
	model->row_index = staged.row_index;
	model->column_index = staged.column_index;
	model->value = staged.value;
	model->row_value = staged.row_value;
	model->nnz = nonzeros;
	return 0;
}

void lp_model_free(struct lp_Model *model)
{
	if (model == NULL)
		return;
	lp_simplex_free(model->coefficients);
	lp_simplex_free(model->constraints);
	lp_simplex_free(model->bounds);
	lp_simplex_free(model->objective);
	model_free_sparse(model);
	lp_simplex_free(model);
}
