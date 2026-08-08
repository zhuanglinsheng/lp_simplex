/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#include "simplex_presolve.h"
#include "simplex_presolve_activity.h"
#include "simplex_presolve_queue.h"
#include "simplex_presolve_rules.h"
#include "simplex_presolve_substitution.h"
#include "utils.h"

#include <lp_simplex/status.h>

#include <math.h>
#include <stdlib.h>


/* A zero-cost column occurring in one inequality can be fixed at a finite
 * bound that relaxes that inequality.  This preserves at least every feasible
 * assignment of the remaining variables and cannot alter the objective. */
static int presolve_choose_singleton_inequality_column(
		const struct lp_Model *model,
		const struct optm_VariableBound *bounds,
		const unsigned char *row_active, const int column, double *value)
{
	int i, k, row = -1;
	double coefficient = 0.;
	const struct optm_VariableBound *bound = bounds + column;
	if (model->objective[column] != 0.)
		return 0;
	if (model->column_start != NULL) {
		for (k = model->column_start[column];
		     k < model->column_start[column + 1]; k++)
			if (row_active[model->row_index[k]]) {
				row = model->row_index[k];
				coefficient = model->value[k];
				break;
			}
	} else {
		for (i = 0; i < model->m; i++)
			if (row_active[i] &&
			    model->constraints[i].coef[column] != 0.) {
				row = i;
				coefficient = model->constraints[i].coef[column];
				break;
			}
	}
	if (row < 0 || coefficient == 0. ||
	    model->constraints[row].type == optm_CONS_T_EQ)
		return 0;
	if (model->constraints[row].type == optm_CONS_T_LE) {
		if (coefficient > 0. && simplex_presolve_has_lower(bound))
			*value = bound->lb;
		else if (coefficient < 0. && simplex_presolve_has_upper(bound))
			*value = bound->ub;
		else return 0;
	} else {
		if (coefficient > 0. && simplex_presolve_has_upper(bound))
			*value = bound->ub;
		else if (coefficient < 0. && simplex_presolve_has_lower(bound))
			*value = bound->lb;
		else return 0;
	}
	return 1;
}


static void presolve_deactivate_column(
		const struct lp_Model *model, const int column,
		const unsigned char *row_active, int *row_degree)
{
	int i, k;
	if (model->column_start != NULL) {
		for (k = model->column_start[column];
		     k < model->column_start[column + 1]; k++) {
			i = model->row_index[k];
			if (row_active[i])
				row_degree[i]--;
		}
		return;
	}
	for (i = 0; i < model->m; i++)
		if (row_active[i] && model->constraints[i].coef[column] != 0.)
			row_degree[i]--;
}


static void presolve_deactivate_row(
		const struct lp_Model *model, const int row,
		const unsigned char *column_active, int *column_degree)
{
	int j, k;
	if (model->row_start != NULL) {
		for (k = model->row_start[row]; k < model->row_start[row + 1]; k++) {
			j = model->column_index[k];
			if (column_active[j])
				column_degree[j]--;
		}
		return;
	}
	for (j = 0; j < model->n; j++)
		if (column_active[j] && model->constraints[row].coef[j] != 0.)
			column_degree[j]--;
}


static void presolve_remove_column_from_rhs(
		const struct lp_Model *model, const int column,
		const double value, const unsigned char *row_active, double *rhs)
{
	int i, k;
	if (model->column_start != NULL) {
		for (k = model->column_start[column];
		     k < model->column_start[column + 1]; k++)
			if (row_active[model->row_index[k]])
				rhs[model->row_index[k]] -= model->value[k] * value;
		return;
	}
	for (i = 0; i < model->m; i++)
		if (row_active[i])
			rhs[i] -= model->constraints[i].coef[column] * value;
}




struct presolve_RowSignature {
	int row;
	int degree;
	unsigned long hash;
};


static int presolve_compare_signature(const void *left, const void *right)
{
	const struct presolve_RowSignature *a =
		(const struct presolve_RowSignature *)left;
	const struct presolve_RowSignature *b =
		(const struct presolve_RowSignature *)right;
	if (a->degree != b->degree)
		return a->degree < b->degree ? -1 : 1;
	if (a->hash != b->hash)
		return a->hash < b->hash ? -1 : 1;
	return a->row < b->row ? -1 : a->row != b->row;
}


static unsigned long presolve_row_support_hash(
		const struct lp_Model *model, const int row,
		const unsigned char *column_active)
{
	unsigned long hash = 2166136261UL;
	int j, k;
	if (model->row_start != NULL) {
		for (k = model->row_start[row]; k < model->row_start[row + 1]; k++) {
			j = model->column_index[k];
			if (column_active[j])
				hash = (hash ^ (unsigned long)(j + 1)) * 16777619UL;
		}
	} else {
		for (j = 0; j < model->n; j++)
			if (column_active[j] && model->constraints[row].coef[j] != 0.)
				hash = (hash ^ (unsigned long)(j + 1)) * 16777619UL;
	}
	return hash;
}


static int presolve_parallel_rows(
		const struct lp_Model *model, const int first, const int second,
		const unsigned char *column_active, const double tolerance,
		double *scale)
{
	int j, ka, kb;
	double a, b;
	*scale = 0.;
	if (model->row_start == NULL) {
		for (j = 0; j < model->n; j++) {
			if (!column_active[j])
				continue;
			a = model->constraints[first].coef[j];
			b = model->constraints[second].coef[j];
			if ((a == 0.) != (b == 0.))
				return 0;
			if (a == 0.)
				continue;
			if (*scale == 0.)
				*scale = b / a;
			else if (__lp_simplex_ABS__(b - *scale * a) > tolerance *
				(1. + __lp_simplex_ABS__(b) +
				 __lp_simplex_ABS__(*scale * a)))
				return 0;
		}
		return *scale != 0.;
	}
	ka = model->row_start[first];
	kb = model->row_start[second];
	while (1) {
		while (ka < model->row_start[first + 1] &&
		       !column_active[model->column_index[ka]])
			ka++;
		while (kb < model->row_start[second + 1] &&
		       !column_active[model->column_index[kb]])
			kb++;
		if (ka == model->row_start[first + 1] ||
		    kb == model->row_start[second + 1])
			break;
		if (model->column_index[ka] != model->column_index[kb])
			return 0;
		a = model->row_value[ka++];
		b = model->row_value[kb++];
		if (*scale == 0.)
			*scale = b / a;
		else if (__lp_simplex_ABS__(b - *scale * a) > tolerance *
			(1. + __lp_simplex_ABS__(b) +
			 __lp_simplex_ABS__(*scale * a)))
			return 0;
	}
	while (ka < model->row_start[first + 1] &&
	       !column_active[model->column_index[ka]])
		ka++;
	while (kb < model->row_start[second + 1] &&
	       !column_active[model->column_index[kb]])
		kb++;
	return ka == model->row_start[first + 1] &&
		kb == model->row_start[second + 1] && *scale != 0.;
}


static int presolve_normalized_row_type(const int type, const double scale)
{
	if (scale >= 0. || type == optm_CONS_T_EQ)
		return type;
	return type == optm_CONS_T_GE ? optm_CONS_T_LE : optm_CONS_T_GE;
}


static int presolve_remove_parallel_rows(
		const struct lp_Model *model, unsigned char *row_active,
		const unsigned char *column_active, const int *row_degree,
		int *column_degree, double *rhs, const double tolerance,
		int *removed)
{
	struct presolve_RowSignature *signature;
	int count = 0, p, q;
	*removed = 0;
	signature = (struct presolve_RowSignature *)lp_simplex_malloc(
		(size_t)model->m * sizeof(*signature));
	if (signature == NULL)
		return lp_simplex_EXIT_FAILURE;
	for (p = 0; p < model->m; p++)
		if (row_active[p] && row_degree[p] >= 2) {
			signature[count].row = p;
			signature[count].degree = row_degree[p];
			signature[count].hash = presolve_row_support_hash(
				model, p, column_active);
			count++;
		}
	qsort(signature, (size_t)count, sizeof(*signature),
		presolve_compare_signature);
	for (p = 0; p < count; p++) {
		int first = signature[p].row;
		if (!row_active[first])
			continue;
		for (q = p + 1; q < count &&
		     signature[q].degree == signature[p].degree &&
		     signature[q].hash == signature[p].hash; q++) {
			int second = signature[q].row;
			int first_type, second_type;
			double scale, candidate, margin;
			if (!row_active[second] || !presolve_parallel_rows(model,
					first, second, column_active,
					__lp_simplex_MIN__(tolerance * 1e-4, 1e-12),
					&scale))
				continue;
			first_type = model->constraints[first].type;
			second_type = presolve_normalized_row_type(
				model->constraints[second].type, scale);
			if (first_type != second_type)
				continue;
			candidate = rhs[second] / scale;
			margin = tolerance * (1. + __lp_simplex_ABS__(rhs[first]) +
				__lp_simplex_ABS__(candidate));
			if (first_type == optm_CONS_T_EQ &&
			    __lp_simplex_ABS__(candidate - rhs[first]) > margin) {
				lp_simplex_free(signature);
				return lp_simplex_Infeasibility;
			}
			if (first_type == optm_CONS_T_GE && candidate > rhs[first])
				rhs[first] = candidate;
			else if (first_type == optm_CONS_T_LE && candidate < rhs[first])
				rhs[first] = candidate;
			row_active[second] = 0;
			presolve_deactivate_row(
				model, second, column_active, column_degree);
			(*removed)++;
		}
	}
	lp_simplex_free(signature);
	return lp_simplex_Success;
}


static struct simplex_Problem *presolve_build_sparse_model(
		const struct lp_Model *model, const int kept_rows,
		const int kept_columns, const int *column_map, const int *column_kept,
		const int *row_map, const int *row_kept, const double *adjusted_rhs,
		const struct optm_VariableBound *bounds)
{
	struct simplex_Problem *reduced;
	int i, j, k, nonzeros = 0, next = 0;
	for (i = 0; i < kept_rows; i++) {
		int original_row = row_map[i];
		if (model->row_start != NULL) {
			for (k = model->row_start[original_row];
			     k < model->row_start[original_row + 1]; k++)
				if (column_kept[model->column_index[k]] >= 0)
					nonzeros++;
		} else {
			for (j = 0; j < model->n; j++)
				if (column_kept[j] >= 0 &&
				    model->constraints[original_row].coef[j] != 0.)
					nonzeros++;
		}
	}
	reduced = simplex_problem_create_sparse(
		kept_rows, kept_columns, nonzeros);
	if (reduced == NULL)
		return NULL;
	for (j = 0; j < kept_columns; j++) {
		int original_column = column_map[j];
		reduced->objective[j] = model->objective[original_column];
		reduced->bounds[j] = bounds[original_column];
		reduced->matrix.column_start[j] = next;
		if (model->column_start != NULL) {
			for (k = model->column_start[original_column];
			     k < model->column_start[original_column + 1]; k++) {
				int row = row_kept[model->row_index[k]];
				if (row >= 0) {
					reduced->matrix.row_index[next] = row;
					reduced->matrix.value[next++] = model->value[k];
				}
			}
		} else {
			for (i = 0; i < kept_rows; i++) {
				double value = model->constraints[row_map[i]].coef[
					original_column];
				if (value != 0.) {
					reduced->matrix.row_index[next] = i;
					reduced->matrix.value[next++] = value;
				}
			}
		}
	}
	reduced->matrix.column_start[kept_columns] = next;
	next = 0;
	for (i = 0; i < kept_rows; i++) {
		int original_row = row_map[i];
		reduced->rhs[i] = adjusted_rhs[original_row];
		reduced->row_type[i] =
			(unsigned char)model->constraints[original_row].type;
		reduced->matrix.row_start[i] = next;
		if (model->row_start != NULL) {
			for (k = model->row_start[original_row];
			     k < model->row_start[original_row + 1]; k++) {
				int column = column_kept[model->column_index[k]];
				if (column >= 0) {
					reduced->matrix.column_index[next] = column;
					reduced->matrix.row_value[next++] = model->row_value[k];
				}
			}
		} else {
			for (j = 0; j < kept_columns; j++) {
				double value = model->constraints[original_row].coef[
					column_map[j]];
				if (value != 0.) {
					reduced->matrix.column_index[next] = j;
					reduced->matrix.row_value[next++] = value;
				}
			}
		}
	}
	reduced->matrix.row_start[kept_rows] = next;
	return reduced;
}


/* Temporary mutable state for the first presolve phase.  None of these arrays
 * survive into postsolve, so they do not belong in simplex_Presolve. */
struct presolve_Workspace {
	int *column_kept;
	int *row_kept;
	double *adjusted_rhs;
	struct optm_VariableBound *bounds;
	unsigned char *column_active;
	unsigned char *row_active;
	int *column_degree;
	int *row_degree;
	struct simplex_PresolveQueue column_queue;
	struct simplex_PresolveQueue row_queue;
	double propagation_tolerance;
};


static void presolve_workspace_destroy(struct presolve_Workspace *workspace)
{
	if (workspace == NULL)
		return;
	lp_simplex_free(workspace->column_kept);
	lp_simplex_free(workspace->row_kept);
	lp_simplex_free(workspace->adjusted_rhs);
	lp_simplex_free(workspace->bounds);
	lp_simplex_free(workspace->column_active);
	lp_simplex_free(workspace->row_active);
	lp_simplex_free(workspace->column_degree);
	lp_simplex_free(workspace->row_degree);
	simplex_presolve_queue_destroy(&workspace->column_queue);
	simplex_presolve_queue_destroy(&workspace->row_queue);
	lp_simplex_memset(workspace, 0, sizeof(*workspace));
}


static int presolve_workspace_create(
		struct simplex_Presolve *presolve, const struct lp_Model *model,
		const double tolerance, struct presolve_Workspace *workspace)
{
	int i, j, k;
	lp_simplex_memset(workspace, 0, sizeof(*workspace));
	presolve->stats.exact_bound_propagation = model->m <= 2048 &&
		model->n >= 5 * model->m;
	workspace->propagation_tolerance =
		presolve->stats.exact_bound_propagation ? 0. : tolerance;
	presolve->column_map = (int *)lp_simplex_malloc(
		(size_t)model->n * sizeof(int));
	presolve->row_map = (int *)lp_simplex_malloc(
		(size_t)model->m * sizeof(int));
	presolve->eliminated_value = (double *)lp_simplex_malloc(
		(size_t)model->n * sizeof(double));
	presolve->eliminated = (unsigned char *)lp_simplex_malloc(
		(size_t)model->n * sizeof(unsigned char));
	workspace->column_kept = (int *)lp_simplex_malloc(
		(size_t)model->n * sizeof(int));
	workspace->row_kept = (int *)lp_simplex_malloc(
		(size_t)model->m * sizeof(int));
	workspace->adjusted_rhs = (double *)lp_simplex_malloc(
		(size_t)model->m * sizeof(double));
	workspace->bounds = (struct optm_VariableBound *)lp_simplex_malloc(
		(size_t)model->n * sizeof(*workspace->bounds));
	workspace->column_active = (unsigned char *)lp_simplex_malloc(
		(size_t)model->n * sizeof(unsigned char));
	workspace->row_active = (unsigned char *)lp_simplex_malloc(
		(size_t)model->m * sizeof(unsigned char));
	workspace->column_degree = (int *)lp_simplex_malloc(
		(size_t)model->n * sizeof(int));
	workspace->row_degree = (int *)lp_simplex_malloc(
		(size_t)model->m * sizeof(int));
	if (presolve->column_map == NULL || presolve->row_map == NULL ||
	    presolve->eliminated_value == NULL || presolve->eliminated == NULL ||
	    workspace->column_kept == NULL || workspace->row_kept == NULL ||
	    workspace->adjusted_rhs == NULL || workspace->bounds == NULL ||
	    workspace->column_active == NULL || workspace->row_active == NULL ||
	    workspace->column_degree == NULL || workspace->row_degree == NULL)
		return lp_simplex_EXIT_FAILURE;
	if (simplex_presolve_queue_init(&workspace->column_queue, model->n) ==
	    lp_simplex_EXIT_FAILURE ||
	    simplex_presolve_queue_init(&workspace->row_queue, model->m) ==
	    lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	lp_simplex_memset(presolve->eliminated, 0,
		(size_t)model->n * sizeof(unsigned char));
	lp_simplex_memset(presolve->eliminated_value, 0,
		(size_t)model->n * sizeof(double));
	lp_simplex_memcpy(workspace->bounds, model->bounds,
		(size_t)model->n * sizeof(*workspace->bounds));
	lp_simplex_memset(workspace->column_active, 1,
		(size_t)model->n * sizeof(unsigned char));
	lp_simplex_memset(workspace->row_active, 1,
		(size_t)model->m * sizeof(unsigned char));
	lp_simplex_memset(workspace->column_degree, 0,
		(size_t)model->n * sizeof(int));
	lp_simplex_memset(workspace->row_degree, 0,
		(size_t)model->m * sizeof(int));
	if (model->row_start != NULL) {
		for (i = 0; i < model->m; i++) {
			workspace->row_degree[i] =
				model->row_start[i + 1] - model->row_start[i];
			for (k = model->row_start[i]; k < model->row_start[i + 1]; k++)
				workspace->column_degree[model->column_index[k]]++;
		}
	} else {
		for (i = 0; i < model->m; i++)
			for (j = 0; j < model->n; j++)
				if (model->constraints[i].coef[j] != 0.) {
					workspace->row_degree[i]++;
					workspace->column_degree[j]++;
				}
	}
	if (model->column_start != NULL)
		for (j = 0; j < model->n; j++)
			if (workspace->column_degree[j] == 1) {
				int row = model->row_index[model->column_start[j]];
				presolve->stats.singleton_column_candidates++;
				if (workspace->bounds[j].b_type == optm_BOUND_T_FR)
					presolve->stats.free_singleton_columns++;
				if (model->constraints[row].type == optm_CONS_T_EQ)
					presolve->stats.equality_singleton_columns++;
			}
	for (j = 0; j < model->n; j++)
		workspace->column_kept[j] = -1;
	for (i = 0; i < model->m; i++) {
		workspace->row_kept[i] = -1;
		workspace->adjusted_rhs[i] = model->constraints[i].rhs;
	}
	for (j = 0; j < model->n; j++)
		(void)simplex_presolve_queue_push(&workspace->column_queue, j);
	for (i = 0; i < model->m; i++)
		(void)simplex_presolve_queue_push(&workspace->row_queue, i);
	return lp_simplex_EXIT_SUCCESS;
}


static void presolve_enqueue_column_rows(
		const struct lp_Model *model, const int column,
		const unsigned char *row_active,
		struct simplex_PresolveQueue *row_queue)
{
	int i, k;
	if (model->column_start != NULL) {
		for (k = model->column_start[column];
		     k < model->column_start[column + 1]; k++) {
			i = model->row_index[k];
			if (row_active[i])
				(void)simplex_presolve_queue_push(row_queue, i);
		}
		return;
	}
	for (i = 0; i < model->m; i++)
		if (row_active[i] && model->constraints[i].coef[column] != 0.)
			(void)simplex_presolve_queue_push(row_queue, i);
}


static void presolve_enqueue_row_columns(
		const struct lp_Model *model, const int row,
		const unsigned char *column_active,
		struct simplex_PresolveQueue *column_queue)
{
	int j, k;
	if (model->row_start != NULL) {
		for (k = model->row_start[row]; k < model->row_start[row + 1]; k++) {
			j = model->column_index[k];
			if (column_active[j])
				(void)simplex_presolve_queue_push(column_queue, j);
		}
		return;
	}
	for (j = 0; j < model->n; j++)
		if (column_active[j] && model->constraints[row].coef[j] != 0.)
			(void)simplex_presolve_queue_push(column_queue, j);
}


static void presolve_enqueue_bound_neighbors(
		const struct lp_Model *model, const int row,
		const struct presolve_Workspace *workspace,
		struct simplex_PresolveQueue *row_queue,
		struct simplex_PresolveQueue *column_queue)
{
	int j, k;
	if (model->row_start != NULL) {
		for (k = model->row_start[row]; k < model->row_start[row + 1]; k++) {
			j = model->column_index[k];
			if (!workspace->column_active[j])
				continue;
			(void)simplex_presolve_queue_push(column_queue, j);
			presolve_enqueue_column_rows(
				model, j, workspace->row_active, row_queue);
		}
		return;
	}
	for (j = 0; j < model->n; j++)
		if (workspace->column_active[j] &&
		    model->constraints[row].coef[j] != 0.) {
			(void)simplex_presolve_queue_push(column_queue, j);
			presolve_enqueue_column_rows(
				model, j, workspace->row_active, row_queue);
		}
}


static int presolve_reduce_column(
		struct simplex_Presolve *presolve, const struct lp_Model *model,
		struct presolve_Workspace *workspace, const int column)
{
	int degree, fixed, singleton_inequality = 0, status;
	double value = 0.;
	if (!workspace->column_active[column])
		return lp_simplex_Success;
	degree = workspace->column_degree[column];
	fixed = workspace->bounds[column].b_type == optm_BOUND_T_BS &&
		workspace->bounds[column].lb == workspace->bounds[column].ub;
	if (fixed)
		value = workspace->bounds[column].lb;
	else if (degree == 1 && presolve_choose_singleton_inequality_column(
			model, workspace->bounds, workspace->row_active,
			column, &value)) {
		fixed = 1;
		singleton_inequality = 1;
	} else if (degree == 0) {
		status = simplex_presolve_choose_empty_column(
			model->objective[column], workspace->bounds + column, &value);
		if (status != lp_simplex_Success)
			return status;
		fixed = 1;
		presolve->stats.empty_columns++;
	}
	if (!fixed)
		return lp_simplex_Success;
	workspace->column_active[column] = 0;
	presolve->eliminated[column] = 1;
	presolve->eliminated_value[column] = value;
	presolve->stats.removed_columns++;
	if (degree != 0)
		presolve->stats.fixed_columns++;
	if (singleton_inequality)
		presolve->stats.singleton_inequality_columns++;
	presolve_remove_column_from_rhs(model, column, value,
		workspace->row_active, workspace->adjusted_rhs);
	presolve_deactivate_column(model, column, workspace->row_active,
		workspace->row_degree);
	presolve_enqueue_column_rows(
		model, column, workspace->row_active, &workspace->row_queue);
	return lp_simplex_Success;
}


static int presolve_find_singleton(
		const struct lp_Model *model, const int row,
		const unsigned char *column_active, int *column, double *coefficient)
{
	int j, k;
	if (model->row_start != NULL) {
		for (k = model->row_start[row]; k < model->row_start[row + 1]; k++) {
			j = model->column_index[k];
			if (column_active[j]) {
				*column = j;
				*coefficient = model->row_value[k];
				return 1;
			}
		}
	} else {
		for (j = 0; j < model->n; j++)
			if (column_active[j] && model->constraints[row].coef[j] != 0.) {
				*column = j;
				*coefficient = model->constraints[row].coef[j];
				return 1;
			}
	}
	return 0;
}


static int presolve_reduce_row(
		struct simplex_Presolve *presolve, const struct lp_Model *model,
		const double tolerance, struct presolve_Workspace *workspace,
		const int row)
{
	int entries, singleton = -1, classification;
	int tightened = 0, forcing = 0, status;
	double coefficient = 0.;
	if (!workspace->row_active[row])
		return lp_simplex_Success;
	entries = workspace->row_degree[row];
	if (entries == 1)
		presolve_find_singleton(model, row, workspace->column_active,
			&singleton, &coefficient);
	if (entries == 0) {
		if (!simplex_presolve_empty_row_feasible(
				model->constraints[row].type,
				workspace->adjusted_rhs[row], tolerance))
			return lp_simplex_Infeasibility;
		workspace->row_active[row] = 0;
		presolve_deactivate_row(model, row, workspace->column_active,
			workspace->column_degree);
		presolve_enqueue_row_columns(
			model, row, workspace->column_active, &workspace->column_queue);
		presolve->stats.removed_rows++;
		return lp_simplex_Success;
	}
	if (entries == 1 && coefficient != 0.) {
		status = simplex_presolve_tighten_singleton(
			workspace->bounds + singleton,
			model->constraints[row].type, coefficient,
			workspace->adjusted_rhs[row], tolerance, &tightened);
		if (status != lp_simplex_Success)
			return status;
		workspace->row_active[row] = 0;
		presolve_deactivate_row(model, row, workspace->column_active,
			workspace->column_degree);
		presolve_enqueue_row_columns(
			model, row, workspace->column_active, &workspace->column_queue);
		presolve->stats.removed_rows++;
		presolve->stats.singleton_rows++;
		presolve->stats.tightened_bounds += tightened;
		if (model->constraints[row].type == optm_CONS_T_EQ)
			presolve->stats.singleton_columns++;
		if (tightened != 0)
			presolve_enqueue_column_rows(model, singleton,
				workspace->row_active, &workspace->row_queue);
		return lp_simplex_Success;
	}
	classification = simplex_presolve_analyze_row(model, row,
		workspace->column_active, workspace->bounds,
		workspace->adjusted_rhs[row], tolerance,
		workspace->propagation_tolerance, 1, &tightened, &forcing);
	if (classification < 0)
		return lp_simplex_Infeasibility;
	if (classification > 0) {
		workspace->row_active[row] = 0;
		presolve_deactivate_row(model, row, workspace->column_active,
			workspace->column_degree);
		presolve_enqueue_row_columns(
			model, row, workspace->column_active, &workspace->column_queue);
		presolve->stats.removed_rows++;
		presolve->stats.redundant_rows++;
	} else if (tightened != 0) {
		presolve->stats.tightened_bounds += tightened;
		if (forcing) {
			presolve->stats.forcing_rows++;
			presolve->stats.forced_columns += tightened;
		}
		presolve_enqueue_bound_neighbors(model, row, workspace,
			&workspace->row_queue, &workspace->column_queue);
	}
	return lp_simplex_Success;
}


static int presolve_reduce_to_fixed_point(
		struct simplex_Presolve *presolve, const struct lp_Model *model,
		const double tolerance, struct presolve_Workspace *workspace)
{
	for (;;) {
		int column, row, duplicate_removed, status;
		presolve->stats.queue_epochs++;
		while (workspace->column_queue.count != 0 ||
		       workspace->row_queue.count != 0) {
			while (simplex_presolve_queue_pop(
					&workspace->column_queue, &column)) {
				status = presolve_reduce_column(
					presolve, model, workspace, column);
				if (status != lp_simplex_Success)
					return status;
			}
			if (simplex_presolve_queue_pop(&workspace->row_queue, &row)) {
				status = presolve_reduce_row(
					presolve, model, tolerance, workspace, row);
				if (status != lp_simplex_Success)
					return status;
			}
		}
		status = presolve_remove_parallel_rows(model,
			workspace->row_active, workspace->column_active,
			workspace->row_degree, workspace->column_degree,
			workspace->adjusted_rhs, tolerance, &duplicate_removed);
		if (status != lp_simplex_Success)
			return status;
		if (duplicate_removed == 0)
			return lp_simplex_Success;
		presolve->stats.removed_rows += duplicate_removed;
		presolve->stats.duplicate_rows += duplicate_removed;
		for (column = 0; column < model->n; column++)
			if (workspace->column_active[column])
				(void)simplex_presolve_queue_push(
					&workspace->column_queue, column);
		for (row = 0; row < model->m; row++)
			if (workspace->row_active[row])
				(void)simplex_presolve_queue_push(
					&workspace->row_queue, row);
	}
}


static void presolve_collect_active(
		struct simplex_Presolve *presolve, const struct lp_Model *model,
		struct presolve_Workspace *workspace,
		int *kept_rows, int *kept_columns)
{
	int i, j;
	*kept_rows = 0;
	*kept_columns = 0;
	for (j = 0; j < model->n; j++)
		if (workspace->column_active[j]) {
			workspace->column_kept[j] = *kept_columns;
			presolve->column_map[(*kept_columns)++] = j;
		}
	for (i = 0; i < model->m; i++)
		if (workspace->row_active[i]) {
			workspace->row_kept[i] = *kept_rows;
			presolve->row_map[(*kept_rows)++] = i;
		}
}


int simplex_presolve_run(
		struct simplex_Presolve *presolve, const struct lp_Model *model,
		const double tolerance)
{
	struct presolve_Workspace workspace;
	int kept_rows, kept_columns, status;
	if (presolve == NULL || model == NULL || model->bounds == NULL)
		return lp_simplex_EXIT_FAILURE;
	lp_simplex_memset(presolve, 0, sizeof(*presolve));
	presolve->original = model;
	presolve->status = lp_simplex_CondUnsatisfied;
	if (presolve_workspace_create(presolve, model, tolerance, &workspace) ==
	    lp_simplex_EXIT_FAILURE)
		goto failure;
	status = presolve_reduce_to_fixed_point(
		presolve, model, tolerance, &workspace);
	if (status == lp_simplex_EXIT_FAILURE)
		goto failure;
	if (status != lp_simplex_Success) {
		presolve->status = status;
		goto finish;
	}
	presolve_collect_active(presolve, model, &workspace,
		&kept_rows, &kept_columns);
	if (kept_rows == 0 || kept_columns == 0) {
		presolve->status = lp_simplex_Success;
		goto finish;
	}
	presolve->reduced = presolve_build_sparse_model(model,
		kept_rows, kept_columns, presolve->column_map,
		workspace.column_kept, presolve->row_map, workspace.row_kept,
		workspace.adjusted_rhs, workspace.bounds);
	if (presolve->reduced == NULL)
		goto failure;
	if (simplex_presolve_substitute_doubletons(presolve, tolerance) ==
	    lp_simplex_EXIT_FAILURE)
		goto failure;
finish:
	lp_simplex_free(presolve->row_map);
	presolve->row_map = NULL;
	if (presolve->status != lp_simplex_CondUnsatisfied) {
		lp_simplex_free(presolve->column_map);
		presolve->column_map = NULL;
	}
	presolve_workspace_destroy(&workspace);
	return lp_simplex_EXIT_SUCCESS;
failure:
	presolve_workspace_destroy(&workspace);
	simplex_presolve_destroy(presolve);
	return lp_simplex_EXIT_FAILURE;
}
