/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#include <lp_simplex/model.h>
#include "utils.h"


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
	model->objective = NULL;
	model->coefficients = NULL;
	model->constraints = NULL;
	model->bounds = NULL;

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


void lp_model_free(struct lp_Model *model)
{
	if (model == NULL)
		return;
	lp_simplex_free(model->coefficients);
	lp_simplex_free(model->constraints);
	lp_simplex_free(model->bounds);
	lp_simplex_free(model->objective);
	lp_simplex_free(model);
}
