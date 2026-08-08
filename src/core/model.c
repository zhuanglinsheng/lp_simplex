/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#include "utils.h"

#include <lp_simplex/model.h>


struct lp_Model *lp_model_create(int m, int n)
{
	struct lp_Model *model;
	int i;

	if (m <= 0 || n <= 0)
		return NULL;
	model = (struct lp_Model *)lp_simplex_malloc(sizeof(*model));
	if (model == NULL)
		return NULL;
	model->m = m;
	model->n = n;
	model->nnz = 0;
	model->objective_offset = 0.;
	model->objective = NULL;
	model->coefficients = NULL;
	model->constraints = NULL;
	model->bounds = NULL;
	model->column_start = NULL;
	model->row_index = NULL;
	model->value = NULL;
	model->row_start = NULL;
	model->column_index = NULL;
	model->row_value = NULL;

	model->objective = (double *)lp_simplex_malloc((size_t)n * sizeof(double));
	model->coefficients = (double *)lp_simplex_malloc(
		(size_t)m * n * sizeof(double));
	model->constraints = (struct optm_LinearConstraint *)lp_simplex_malloc(
		(size_t)m * sizeof(*model->constraints));
	model->bounds = (struct optm_VariableBound *)lp_simplex_malloc(
		(size_t)n * sizeof(*model->bounds));
	if (model->objective == NULL || model->coefficients == NULL ||
	    model->constraints == NULL || model->bounds == NULL) {
		lp_model_free(model);
		return NULL;
	}

	lp_simplex_memset(model->objective, 0, (size_t)n * sizeof(double));
	lp_simplex_memset(model->coefficients, 0,
			  (size_t)m * n * sizeof(double));
	lp_simplex_memset(model->constraints, 0,
			  (size_t)m * sizeof(*model->constraints));
	lp_simplex_memset(model->bounds, 0,
			  (size_t)n * sizeof(*model->bounds));
	for (i = 0; i < m; i++)
		model->constraints[i].coef = model->coefficients + i * n;
	for (i = 0; i < n; i++) {
		model->bounds[i].lb = 0.;
		model->bounds[i].ub = __lp_simplex_INF__;
		model->bounds[i].b_type = optm_BOUND_T_LO;
		model->bounds[i].v_type = optm_VAR_T_REAL;
	}
	return model;
}


struct lp_Model *lp_model_create_sparse(
		int m, int n, int nonzero_capacity)
{
	struct lp_Model *model;
	int i;
	if (m <= 0 || n <= 0 || nonzero_capacity < 0)
		return NULL;
	model = (struct lp_Model *)lp_simplex_malloc(sizeof(*model));
	if (model == NULL)
		return NULL;
	model->m = m;
	model->n = n;
	model->nnz = 0;
	model->objective_offset = 0.;
	model->objective = (double *)lp_simplex_malloc((size_t)n * sizeof(double));
	model->coefficients = NULL;
	model->constraints = (struct optm_LinearConstraint *)lp_simplex_malloc(
		(size_t)m * sizeof(*model->constraints));
	model->bounds = (struct optm_VariableBound *)lp_simplex_malloc(
		(size_t)n * sizeof(*model->bounds));
	model->column_start = (int *)lp_simplex_malloc(
		(size_t)(n + 1) * sizeof(int));
	model->row_start = (int *)lp_simplex_malloc(
		(size_t)(m + 1) * sizeof(int));
	model->row_index = nonzero_capacity != 0 ? (int *)lp_simplex_malloc(
		(size_t)nonzero_capacity * sizeof(int)) : NULL;
	model->value = nonzero_capacity != 0 ? (double *)lp_simplex_malloc(
		(size_t)nonzero_capacity * sizeof(double)) : NULL;
	model->column_index = nonzero_capacity != 0 ? (int *)lp_simplex_malloc(
		(size_t)nonzero_capacity * sizeof(int)) : NULL;
	model->row_value = nonzero_capacity != 0 ? (double *)lp_simplex_malloc(
		(size_t)nonzero_capacity * sizeof(double)) : NULL;
	if (model->objective == NULL || model->constraints == NULL ||
	    model->bounds == NULL || model->column_start == NULL ||
	    model->row_start == NULL ||
	    (nonzero_capacity != 0 && (model->row_index == NULL ||
	     model->value == NULL || model->column_index == NULL ||
	     model->row_value == NULL))) {
		lp_model_free(model);
		return NULL;
	}
	lp_simplex_memset(model->objective, 0, (size_t)n * sizeof(double));
	lp_simplex_memset(model->constraints, 0,
		(size_t)m * sizeof(*model->constraints));
	lp_simplex_memset(model->bounds, 0, (size_t)n * sizeof(*model->bounds));
	lp_simplex_memset(model->column_start, 0, (size_t)(n + 1) * sizeof(int));
	lp_simplex_memset(model->row_start, 0, (size_t)(m + 1) * sizeof(int));
	for (i = 0; i < m; i++)
		model->constraints[i].coef = NULL;
	for (i = 0; i < n; i++) {
		model->bounds[i].lb = 0.;
		model->bounds[i].ub = __lp_simplex_INF__;
		model->bounds[i].b_type = optm_BOUND_T_LO;
		model->bounds[i].v_type = optm_VAR_T_REAL;
	}
	return model;
}


int lp_model_build_sparse(struct lp_Model *model)
{
	int i, j, k, nonzeros = 0;
	int *next = NULL;
	if (model == NULL || model->coefficients == NULL)
		return -1;
	for (i = 0; i < model->m; i++)
		for (j = 0; j < model->n; j++)
			if (model->constraints[i].coef[j] != 0.)
				nonzeros++;
	lp_simplex_free(model->column_start);
	lp_simplex_free(model->row_index);
	lp_simplex_free(model->value);
	lp_simplex_free(model->row_start);
	lp_simplex_free(model->column_index);
	lp_simplex_free(model->row_value);
	model->column_start = (int *)lp_simplex_malloc(
		(size_t)(model->n + 1) * sizeof(int));
	model->row_start = (int *)lp_simplex_malloc(
		(size_t)(model->m + 1) * sizeof(int));
	model->row_index = nonzeros > 0 ? (int *)lp_simplex_malloc(
		(size_t)nonzeros * sizeof(int)) : NULL;
	model->value = nonzeros > 0 ? (double *)lp_simplex_malloc(
		(size_t)nonzeros * sizeof(double)) : NULL;
	model->column_index = nonzeros > 0 ? (int *)lp_simplex_malloc(
		(size_t)nonzeros * sizeof(int)) : NULL;
	model->row_value = nonzeros > 0 ? (double *)lp_simplex_malloc(
		(size_t)nonzeros * sizeof(double)) : NULL;
	next = (int *)lp_simplex_malloc((size_t)model->m * sizeof(int));
	if (model->column_start == NULL || model->row_start == NULL || next == NULL ||
	    (nonzeros > 0 && (model->row_index == NULL || model->value == NULL ||
	     model->column_index == NULL || model->row_value == NULL))) {
		lp_simplex_free(next);
		return -1;
	}
	k = 0;
	for (j = 0; j < model->n; j++) {
		model->column_start[j] = k;
		for (i = 0; i < model->m; i++) {
			double coefficient = model->constraints[i].coef[j];
			if (coefficient != 0.) {
				model->row_index[k] = i;
				model->value[k++] = coefficient;
			}
		}
	}
	model->column_start[model->n] = k;
	lp_simplex_memset(model->row_start, 0,
		(size_t)(model->m + 1) * sizeof(int));
	for (k = 0; k < nonzeros; k++)
		model->row_start[model->row_index[k] + 1]++;
	for (i = 0; i < model->m; i++) {
		model->row_start[i + 1] += model->row_start[i];
		next[i] = model->row_start[i];
	}
	for (j = 0; j < model->n; j++)
		for (k = model->column_start[j]; k < model->column_start[j + 1]; k++) {
			i = model->row_index[k];
			model->column_index[next[i]] = j;
			model->row_value[next[i]++] = model->value[k];
		}
	lp_simplex_free(next);
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
	lp_simplex_free(model->column_start);
	lp_simplex_free(model->row_index);
	lp_simplex_free(model->value);
	lp_simplex_free(model->row_start);
	lp_simplex_free(model->column_index);
	lp_simplex_free(model->row_value);
	lp_simplex_free(model);
}
