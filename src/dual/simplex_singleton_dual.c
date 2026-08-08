/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
/* Dualize equality models whose row residuals use paired singleton columns. */
#include "simplex_singleton_dual.h"
#include "simplex_dual.h"
#include "utils.h"

#include <lp_simplex/status.h>


static int singleton_dual_structure(
		const struct lp_Model *model, int *positive, int *negative,
		int *core_index, int *core_count)
{
	int i, j;
	for (i = 0; i < model->m; i++) {
		if (model->constraints[i].type != optm_CONS_T_EQ)
			return 0;
		positive[i] = -1;
		negative[i] = -1;
	}
	*core_count = 0;
	for (j = 0; j < model->n; j++) {
		int count = 0, row = -1;
		if (model->bounds[j].b_type != optm_BOUND_T_LO ||
		    model->bounds[j].lb != 0. ||
		    model->bounds[j].v_type != optm_VAR_T_REAL)
			return 0;
		if (model->column_start != NULL) {
			count = model->column_start[j + 1] - model->column_start[j];
			if (count == 1)
				row = model->row_index[model->column_start[j]];
		} else
			for (i = 0; i < model->m; i++)
				if (model->constraints[i].coef[j] != 0.) {
					count++;
					row = i;
				}
		if (count == 1 && model->objective[j] > 0.) {
			double coefficient = model->column_start != NULL
				? model->value[model->column_start[j]]
				: model->constraints[row].coef[j];
			if (coefficient > 0.) {
				if (positive[row] >= 0)
					return 0;
				positive[row] = j;
			} else {
				if (negative[row] >= 0)
					return 0;
				negative[row] = j;
			}
			core_index[j] = -1;
		} else {
			core_index[j] = (*core_count)++;
		}
	}
	if (*core_count <= 0)
		return 0;
	for (i = 0; i < model->m; i++)
		if (positive[i] < 0 || negative[i] < 0)
			return 0;
	return 1;
}


int simplex_singleton_dual_solve(
		const struct lp_Model *model,
		const struct lp_simplex_Options *options,
		double *x, struct lp_simplex_Result *result, int *applicable)
{
	int i, j, core_count = 0, state;
	int *positive = NULL, *negative = NULL, *core_index = NULL;
	double *dual_x = NULL, *row_dual = NULL;
	struct lp_Model *dual = NULL;
	*applicable = 0;
	positive = (int *)lp_simplex_malloc((size_t)model->m * sizeof(int));
	negative = (int *)lp_simplex_malloc((size_t)model->m * sizeof(int));
	core_index = (int *)lp_simplex_malloc((size_t)model->n * sizeof(int));
	if (positive == NULL || negative == NULL || core_index == NULL)
		goto MEMORY_FAILURE;
	if (!singleton_dual_structure(
		model, positive, negative, core_index, &core_count))
		goto NOT_APPLICABLE;
	*applicable = 1;
	dual = lp_model_create(core_count, model->m);
	dual_x = (double *)lp_simplex_malloc((size_t)model->m * sizeof(double));
	row_dual = (double *)lp_simplex_malloc((size_t)core_count * sizeof(double));
	if (dual == NULL || dual_x == NULL || row_dual == NULL)
		goto MEMORY_FAILURE;
	for (i = 0; i < model->m; i++) {
		int p = positive[i], n = negative[i];
		double ap = model->constraints[i].coef[p];
		double an = model->constraints[i].coef[n];
		dual->objective[i] = -model->constraints[i].rhs;
		dual->bounds[i].b_type = optm_BOUND_T_BS;
		dual->bounds[i].lb = model->objective[n] / an;
		dual->bounds[i].ub = model->objective[p] / ap;
		if (dual->bounds[i].lb > dual->bounds[i].ub)
			goto NOT_APPLICABLE;
	}
	for (j = 0; j < model->n; j++) {
		int row = core_index[j];
		if (row < 0)
			continue;
		dual->constraints[row].type = optm_CONS_T_LE;
		dual->constraints[row].rhs = model->objective[j];
		for (i = 0; i < model->m; i++)
			dual->constraints[row].coef[i] =
				model->constraints[i].coef[j];
	}
	state = simplex_dual_solve(dual, options, dual_x, row_dual, result, 0);
	if (state == lp_simplex_EXIT_SUCCESS &&
	    result->status == lp_simplex_Success) {
		double objective = 0., maximum = 0.;
		lp_simplex_memset(x, 0, (size_t)model->n * sizeof(double));
		for (j = 0; j < model->n; j++)
			if (core_index[j] >= 0) {
				x[j] = -row_dual[core_index[j]];
				if (x[j] < 0. && x[j] >= -options->primal_tolerance)
					x[j] = 0.;
			}
		for (i = 0; i < model->m; i++) {
			double residual = model->constraints[i].rhs;
			for (j = 0; j < model->n; j++)
				if (core_index[j] >= 0)
					residual -= model->constraints[i].coef[j] * x[j];
			if (residual >= 0.)
				x[positive[i]] = residual /
					model->constraints[i].coef[positive[i]];
			else
				x[negative[i]] = residual /
					model->constraints[i].coef[negative[i]];
		}
		for (j = 0; j < model->n; j++) {
			objective += model->objective[j] * x[j];
			if (x[j] < 0.)
				maximum = __lp_simplex_MAX__(maximum, -x[j]);
		}
		for (i = 0; i < model->m; i++) {
			double activity = 0.;
			for (j = 0; j < model->n; j++)
				activity += model->constraints[i].coef[j] * x[j];
			maximum = __lp_simplex_MAX__(maximum,
				__lp_simplex_ABS__(activity - model->constraints[i].rhs));
		}
		result->objective = objective;
		result->primal_infeasibility = maximum;
		if (maximum > 10. * options->primal_tolerance) {
			result->status = lp_simplex_PrecisionError;
			state = lp_simplex_EXIT_FAILURE;
		}
	}
	lp_model_free(dual);
	lp_simplex_free(dual_x);
	lp_simplex_free(row_dual);
	lp_simplex_free(positive);
	lp_simplex_free(negative);
	lp_simplex_free(core_index);
	return state;

MEMORY_FAILURE:
	*applicable = 1;
	result->status = lp_simplex_MemoryAllocError;
	lp_model_free(dual);
	lp_simplex_free(dual_x);
	lp_simplex_free(row_dual);
	lp_simplex_free(positive);
	lp_simplex_free(negative);
	lp_simplex_free(core_index);
	return lp_simplex_EXIT_FAILURE;

NOT_APPLICABLE:
	*applicable = 0;
	lp_model_free(dual);
	lp_simplex_free(dual_x);
	lp_simplex_free(row_dual);
	lp_simplex_free(positive);
	lp_simplex_free(negative);
	lp_simplex_free(core_index);
	return lp_simplex_EXIT_FAILURE;
}
