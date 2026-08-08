/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#include "simplex_presolve.h"
#include "simplex_presolve_substitution.h"
#include "utils.h"

#include <lp_simplex/status.h>

#include <math.h>
#include <stdlib.h>


static int presolve_finite_lower(const struct optm_VariableBound *bound)
{
	return bound->b_type == optm_BOUND_T_LO ||
		bound->b_type == optm_BOUND_T_BS;
}


static int presolve_finite_upper(const struct optm_VariableBound *bound)
{
	return bound->b_type == optm_BOUND_T_UP ||
		bound->b_type == optm_BOUND_T_BS;
}


static int presolve_empty_row_feasible(
		const struct optm_LinearConstraint *constraint,
		const double rhs, const double tolerance)
{
	if (constraint->type == optm_CONS_T_EQ)
		return __lp_simplex_ABS__(rhs) <= tolerance;
	if (constraint->type == optm_CONS_T_GE)
		return rhs <= tolerance;
	return rhs >= -tolerance;
}


static int presolve_choose_empty_column(
		const struct lp_Model *model,
		const struct optm_VariableBound *bounds,
		const int column, double *value)
{
	const struct optm_VariableBound *bound = bounds + column;
	double cost = model->objective[column];
	if (cost > 0.) {
		if (!presolve_finite_lower(bound))
			return lp_simplex_Unboundedness;
		*value = bound->lb;
	} else if (cost < 0.) {
		if (!presolve_finite_upper(bound))
			return lp_simplex_Unboundedness;
		*value = bound->ub;
	} else if (presolve_finite_lower(bound) && bound->lb > 0.) {
		*value = bound->lb;
	} else if (presolve_finite_upper(bound) && bound->ub < 0.) {
		*value = bound->ub;
	} else {
		*value = 0.;
	}
	return lp_simplex_Success;
}


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
		if (coefficient > 0. && presolve_finite_lower(bound))
			*value = bound->lb;
		else if (coefficient < 0. && presolve_finite_upper(bound))
			*value = bound->ub;
		else return 0;
	} else {
		if (coefficient > 0. && presolve_finite_upper(bound))
			*value = bound->ub;
		else if (coefficient < 0. && presolve_finite_lower(bound))
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


static void presolve_bound_set_lower(
		struct optm_VariableBound *bound, const double lower)
{
	bound->lb = lower;
	if (presolve_finite_upper(bound))
		bound->b_type = optm_BOUND_T_BS;
	else
		bound->b_type = optm_BOUND_T_LO;
}


static void presolve_bound_set_upper(
		struct optm_VariableBound *bound, const double upper)
{
	bound->ub = upper;
	if (presolve_finite_lower(bound))
		bound->b_type = optm_BOUND_T_BS;
	else
		bound->b_type = optm_BOUND_T_UP;
}


static int presolve_tighten_singleton(
		struct optm_VariableBound *bound, const int row_type,
		const double coefficient, const double rhs,
		const double tolerance, int *tightened)
{
	double implied = rhs / coefficient;
	double margin = tolerance * (1. + __lp_simplex_ABS__(implied));
	*tightened = 0;
	if (row_type == optm_CONS_T_EQ) {
		if ((presolve_finite_lower(bound) && implied < bound->lb - margin) ||
		    (presolve_finite_upper(bound) && implied > bound->ub + margin))
			return lp_simplex_Infeasibility;
		presolve_bound_set_lower(bound, implied);
		presolve_bound_set_upper(bound, implied);
		*tightened = 1;
		return lp_simplex_Success;
	}
	if ((row_type == optm_CONS_T_GE && coefficient > 0.) ||
	    (row_type == optm_CONS_T_LE && coefficient < 0.)) {
		if (presolve_finite_upper(bound) && implied > bound->ub + margin)
			return lp_simplex_Infeasibility;
		if (!presolve_finite_lower(bound) || implied > bound->lb) {
			presolve_bound_set_lower(bound, implied);
			*tightened = 1;
		}
	} else {
		if (presolve_finite_lower(bound) && implied < bound->lb - margin)
			return lp_simplex_Infeasibility;
		if (!presolve_finite_upper(bound) || implied < bound->ub) {
			presolve_bound_set_upper(bound, implied);
			*tightened = 1;
		}
	}
	return lp_simplex_Success;
}


static int presolve_apply_implied_lower(
		struct optm_VariableBound *bound, double candidate,
		const double tolerance)
{
	double margin;
	if (!isfinite(candidate))
		return 0;
	margin = tolerance * (1. + __lp_simplex_ABS__(candidate));
	candidate -= margin;
	if (presolve_finite_upper(bound) && candidate > bound->ub + margin)
		return -1;
	if (presolve_finite_upper(bound) && candidate > bound->ub)
		candidate = bound->ub;
	if (presolve_finite_lower(bound) &&
	    candidate <= bound->lb + margin)
		return 0;
	presolve_bound_set_lower(bound, candidate);
	return 1;
}


static int presolve_apply_implied_upper(
		struct optm_VariableBound *bound, double candidate,
		const double tolerance)
{
	double margin;
	if (!isfinite(candidate))
		return 0;
	margin = tolerance * (1. + __lp_simplex_ABS__(candidate));
	candidate += margin;
	if (presolve_finite_lower(bound) && candidate < bound->lb - margin)
		return -1;
	if (presolve_finite_lower(bound) && candidate < bound->lb)
		candidate = bound->lb;
	if (presolve_finite_upper(bound) &&
	    candidate >= bound->ub - margin)
		return 0;
	presolve_bound_set_upper(bound, candidate);
	return 1;
}


/* Analyze one sparse row and propagate bounds from its residual activity.
 * Return -1 for infeasible, 1 for redundant, and 0 when the row must stay. */
static int presolve_analyze_row(
		const struct lp_Model *model, const int row,
		const unsigned char *column_active,
		struct optm_VariableBound *bounds,
		const double rhs, const double tolerance,
		const double propagation_tolerance, const int propagate,
		int *tightened, int *forcing)
{
	long double minimum = 0., maximum = 0., magnitude = 0.;
	int minimum_infinite = 0, maximum_infinite = 0;
	int k, j, start, end;
	int type = model->constraints[row].type;
	*tightened = 0;
	*forcing = 0;
	if (model->row_start != NULL) {
		start = model->row_start[row];
		end = model->row_start[row + 1];
	} else {
		start = 0;
		end = model->n;
	}
	for (k = start; k < end; k++) {
		double coefficient = model->row_start != NULL
			? model->row_value[k] : model->constraints[row].coef[k];
		const struct optm_VariableBound *bound;
		double lower, upper;
		j = model->row_start != NULL ? model->column_index[k] : k;
		if (!column_active[j] || coefficient == 0.)
			continue;
		bound = bounds + j;
		if (coefficient > 0.) {
			if (presolve_finite_lower(bound)) {
				lower = coefficient * bound->lb;
				minimum += lower;
				magnitude += fabsl((long double)lower);
			} else minimum_infinite++;
			if (presolve_finite_upper(bound)) {
				upper = coefficient * bound->ub;
				maximum += upper;
				magnitude += fabsl((long double)upper);
			} else maximum_infinite++;
		} else {
			if (presolve_finite_upper(bound)) {
				lower = coefficient * bound->ub;
				minimum += lower;
				magnitude += fabsl((long double)lower);
			} else minimum_infinite++;
			if (presolve_finite_lower(bound)) {
				upper = coefficient * bound->lb;
				maximum += upper;
				magnitude += fabsl((long double)upper);
			} else maximum_infinite++;
		}
	}
	{
		long double margin = (long double)tolerance *
			(1. + fabsl((long double)rhs) + magnitude);
		int force_minimum, force_maximum;
		if ((type == optm_CONS_T_GE || type == optm_CONS_T_EQ) &&
		    maximum_infinite == 0 && maximum < (long double)rhs - margin)
			return -1;
		if ((type == optm_CONS_T_LE || type == optm_CONS_T_EQ) &&
		    minimum_infinite == 0 && minimum > (long double)rhs + margin)
			return -1;
		if (type == optm_CONS_T_GE && minimum_infinite == 0 &&
		    minimum >= (long double)rhs + margin)
			return 1;
		if (type == optm_CONS_T_LE && maximum_infinite == 0 &&
		    maximum <= (long double)rhs - margin)
			return 1;
		/* If a one-sided row is feasible only at an activity extreme, every
		 * participating variable is forced to the bound attaining that
		 * extreme.  The same certificate applies to either side of an
		 * equality.  Fixing the variables here lets the ordinary fixed-column
		 * loop perform RHS updates and postsolve bookkeeping. */
		force_minimum =
			(type == optm_CONS_T_LE || type == optm_CONS_T_EQ) &&
			minimum_infinite == 0 &&
			minimum == (long double)rhs;
		force_maximum =
			(type == optm_CONS_T_GE || type == optm_CONS_T_EQ) &&
			maximum_infinite == 0 &&
			maximum == (long double)rhs;
		if (force_minimum || force_maximum) {
			for (k = start; k < end; k++) {
				double coefficient = model->row_start != NULL
					? model->row_value[k]
					: model->constraints[row].coef[k];
				struct optm_VariableBound *bound;
				double value;
				j = model->row_start != NULL
					? model->column_index[k] : k;
				if (!column_active[j] || coefficient == 0.)
					continue;
				bound = bounds + j;
				if (bound->b_type == optm_BOUND_T_BS &&
				    bound->lb == bound->ub)
					continue;
				if (force_minimum)
					value = coefficient > 0. ? bound->lb : bound->ub;
				else
					value = coefficient > 0. ? bound->ub : bound->lb;
				presolve_bound_set_lower(bound, value);
				presolve_bound_set_upper(bound, value);
				(*tightened)++;
			}
			*forcing = *tightened != 0;
			return 0;
		}
	}
	if (!propagate)
		return 0;
	/* A finite residual activity gives a valid implied bound for one variable.
	 * If the row has exactly one infinite contribution, excluding that same
	 * variable also produces a finite residual. */
	for (k = start; k < end; k++) {
		double coefficient = model->row_start != NULL
			? model->row_value[k] : model->constraints[row].coef[k];
		struct optm_VariableBound *bound;
		int own_minimum_infinite, own_maximum_infinite;
		long double own_minimum = 0., own_maximum = 0.;
		long double residual;
		int state;
		j = model->row_start != NULL ? model->column_index[k] : k;
		if (!column_active[j] || coefficient == 0.)
			continue;
		bound = bounds + j;
		own_minimum_infinite = coefficient > 0.
			? !presolve_finite_lower(bound)
			: !presolve_finite_upper(bound);
		own_maximum_infinite = coefficient > 0.
			? !presolve_finite_upper(bound)
			: !presolve_finite_lower(bound);
		if (!own_minimum_infinite)
			own_minimum = (long double)coefficient *
				(coefficient > 0. ? bound->lb : bound->ub);
		if (!own_maximum_infinite)
			own_maximum = (long double)coefficient *
				(coefficient > 0. ? bound->ub : bound->lb);
		if ((type == optm_CONS_T_GE || type == optm_CONS_T_EQ) &&
		    (maximum_infinite == 0 ||
		     (maximum_infinite == 1 && own_maximum_infinite))) {
			residual = maximum - own_maximum;
			if (coefficient > 0.)
				state = presolve_apply_implied_lower(bound,
					(double)(((long double)rhs - residual) /
					coefficient), propagation_tolerance);
			else
				state = presolve_apply_implied_upper(bound,
					(double)(((long double)rhs - residual) /
					coefficient), propagation_tolerance);
			if (state < 0)
				return -1;
			*tightened += state;
		}
		if ((type == optm_CONS_T_LE || type == optm_CONS_T_EQ) &&
		    (minimum_infinite == 0 ||
		     (minimum_infinite == 1 && own_minimum_infinite))) {
			residual = minimum - own_minimum;
			if (coefficient > 0.)
				state = presolve_apply_implied_upper(bound,
					(double)(((long double)rhs - residual) /
					coefficient), propagation_tolerance);
			else
				state = presolve_apply_implied_lower(bound,
					(double)(((long double)rhs - residual) /
					coefficient), propagation_tolerance);
			if (state < 0)
				return -1;
			*tightened += state;
		}
	}
	return 0;
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


int simplex_presolve_run(
		struct simplex_Presolve *presolve, const struct lp_Model *model,
		const double tolerance)
{
	int i, j, k, kept_rows = 0, kept_columns = 0, changed;
	int duplicate_removed, duplicate_status;
	double propagation_tolerance;
	int *column_kept = NULL;
	int *row_kept = NULL;
	double *adjusted_rhs = NULL;
	struct optm_VariableBound *bounds = NULL;
	unsigned char *column_active = NULL;
	unsigned char *row_active = NULL;
	int *column_degree = NULL;
	int *row_degree = NULL;
	if (presolve == NULL || model == NULL || model->bounds == NULL)
		return lp_simplex_EXIT_FAILURE;
	lp_simplex_memset(presolve, 0, sizeof(*presolve));
	presolve->original = model;
	presolve->exact_bound_propagation = model->m <= 2048 &&
		model->n >= 5 * model->m;
	propagation_tolerance = presolve->exact_bound_propagation
		? 0. : tolerance;
	presolve->terminal_status = lp_simplex_CondUnsatisfied;
	presolve->column_map = (int *)lp_simplex_malloc(
		(size_t)model->n * sizeof(int));
	presolve->row_map = (int *)lp_simplex_malloc(
		(size_t)model->m * sizeof(int));
	presolve->eliminated_value = (double *)lp_simplex_malloc(
		(size_t)model->n * sizeof(double));
	presolve->eliminated = (unsigned char *)lp_simplex_malloc(
		(size_t)model->n * sizeof(unsigned char));
	column_kept = (int *)lp_simplex_malloc((size_t)model->n * sizeof(int));
	row_kept = (int *)lp_simplex_malloc((size_t)model->m * sizeof(int));
	adjusted_rhs = (double *)lp_simplex_malloc(
		(size_t)model->m * sizeof(double));
	bounds = (struct optm_VariableBound *)lp_simplex_malloc(
		(size_t)model->n * sizeof(*bounds));
	column_active = (unsigned char *)lp_simplex_malloc(
		(size_t)model->n * sizeof(unsigned char));
	row_active = (unsigned char *)lp_simplex_malloc(
		(size_t)model->m * sizeof(unsigned char));
	column_degree = (int *)lp_simplex_malloc(
		(size_t)model->n * sizeof(int));
	row_degree = (int *)lp_simplex_malloc(
		(size_t)model->m * sizeof(int));
	if (presolve->column_map == NULL || presolve->row_map == NULL ||
	    presolve->eliminated_value == NULL || presolve->eliminated == NULL ||
	    column_kept == NULL || row_kept == NULL || adjusted_rhs == NULL ||
	    bounds == NULL || column_active == NULL || row_active == NULL ||
	    column_degree == NULL || row_degree == NULL)
		goto failure;
	lp_simplex_memset(presolve->eliminated, 0,
		(size_t)model->n * sizeof(unsigned char));
	lp_simplex_memset(presolve->eliminated_value, 0,
		(size_t)model->n * sizeof(double));
	lp_simplex_memcpy(bounds, model->bounds,
		(size_t)model->n * sizeof(*bounds));
	lp_simplex_memset(column_active, 1,
		(size_t)model->n * sizeof(unsigned char));
	lp_simplex_memset(row_active, 1,
		(size_t)model->m * sizeof(unsigned char));
	lp_simplex_memset(column_degree, 0, (size_t)model->n * sizeof(int));
	lp_simplex_memset(row_degree, 0, (size_t)model->m * sizeof(int));
	if (model->row_start != NULL) {
		for (i = 0; i < model->m; i++) {
			row_degree[i] = model->row_start[i + 1] - model->row_start[i];
			for (k = model->row_start[i]; k < model->row_start[i + 1]; k++)
				column_degree[model->column_index[k]]++;
		}
	} else {
		for (i = 0; i < model->m; i++)
			for (j = 0; j < model->n; j++)
				if (model->constraints[i].coef[j] != 0.) {
					row_degree[i]++;
					column_degree[j]++;
				}
	}
	if (model->column_start != NULL)
		for (j = 0; j < model->n; j++)
			if (column_degree[j] == 1) {
				int incident_row = model->row_index[
					model->column_start[j]];
				presolve->singleton_column_candidates++;
				if (bounds[j].b_type == optm_BOUND_T_FR)
					presolve->free_singleton_columns++;
				if (model->constraints[incident_row].type == optm_CONS_T_EQ)
					presolve->equality_singleton_columns++;
			}
	for (j = 0; j < model->n; j++)
		column_kept[j] = -1;
	for (i = 0; i < model->m; i++)
		row_kept[i] = -1;
	for (i = 0; i < model->m; i++)
		adjusted_rhs[i] = model->constraints[i].rhs;

	/* Reductions are iterated to a fixed point.  Removing a row can expose an
	 * empty column; fixing that column can in turn expose singleton rows. */
reduce_again:
	do {
		changed = 0;
		presolve->passes++;
		for (j = 0; j < model->n; j++) {
			int degree, fixed, singleton_inequality = 0, status;
			double value = 0.;
			if (!column_active[j])
				continue;
			degree = column_degree[j];
			fixed = bounds[j].b_type == optm_BOUND_T_BS &&
				bounds[j].lb == bounds[j].ub;
			if (fixed)
				value = bounds[j].lb;
			else if (degree == 1 &&
			    presolve_choose_singleton_inequality_column(model, bounds,
				row_active, j, &value)) {
				fixed = 1;
				singleton_inequality = 1;
			}
			else if (degree == 0) {
				status = presolve_choose_empty_column(
					model, bounds, j, &value);
				if (status != lp_simplex_Success) {
					presolve->terminal = 1;
					presolve->terminal_status = status;
					goto finish;
				}
				fixed = 1;
				presolve->empty_columns++;
			}
			if (!fixed)
				continue;
			column_active[j] = 0;
			presolve->eliminated[j] = 1;
			presolve->eliminated_value[j] = value;
			presolve->removed_columns++;
			if (degree != 0)
				presolve->fixed_columns++;
			if (singleton_inequality)
				presolve->singleton_inequality_columns++;
			presolve_remove_column_from_rhs(
				model, j, value, row_active, adjusted_rhs);
			presolve_deactivate_column(
				model, j, row_active, row_degree);
			changed = 1;
		}
		for (i = 0; i < model->m; i++) {
			int entries, singleton = -1;
			double coefficient = 0.;
			int classification, tightened = 0, forcing = 0, status;
			if (!row_active[i])
				continue;
			entries = row_degree[i];
			if (entries == 1 && model->row_start != NULL) {
				for (k = model->row_start[i];
				     k < model->row_start[i + 1]; k++) {
					j = model->column_index[k];
					if (!column_active[j])
						continue;
					singleton = j;
					coefficient = model->row_value[k];
					break;
				}
			} else if (entries == 1) {
				for (j = 0; j < model->n; j++) {
					if (!column_active[j] ||
					    model->constraints[i].coef[j] == 0.)
						continue;
					singleton = j;
					coefficient = model->constraints[i].coef[j];
					break;
				}
			}
			if (entries == 0) {
				if (!presolve_empty_row_feasible(
						model->constraints + i,
						adjusted_rhs[i], tolerance)) {
					presolve->terminal = 1;
					presolve->terminal_status = lp_simplex_Infeasibility;
					goto finish;
				}
				row_active[i] = 0;
				presolve_deactivate_row(
					model, i, column_active, column_degree);
				presolve->removed_rows++;
				changed = 1;
				continue;
			}
			if (entries == 1 && coefficient != 0.) {
				status = presolve_tighten_singleton(
					bounds + singleton, model->constraints[i].type,
					coefficient, adjusted_rhs[i], tolerance,
					&tightened);
				if (status != lp_simplex_Success) {
					presolve->terminal = 1;
					presolve->terminal_status = status;
					goto finish;
				}
				row_active[i] = 0;
				presolve_deactivate_row(
					model, i, column_active, column_degree);
				presolve->removed_rows++;
				presolve->singleton_rows++;
				presolve->tightened_bounds += tightened;
				if (model->constraints[i].type == optm_CONS_T_EQ)
					presolve->singleton_columns++;
				changed = 1;
				continue;
			}
			classification = presolve_analyze_row(model, i,
				column_active, bounds, adjusted_rhs[i], tolerance,
				propagation_tolerance, 1, &tightened, &forcing);
			if (classification < 0) {
				presolve->terminal = 1;
				presolve->terminal_status = lp_simplex_Infeasibility;
				goto finish;
			}
			if (classification > 0) {
				row_active[i] = 0;
				presolve_deactivate_row(
					model, i, column_active, column_degree);
				presolve->removed_rows++;
				presolve->redundant_rows++;
				changed = 1;
			} else if (tightened != 0) {
				presolve->tightened_bounds += tightened;
				if (forcing) {
					presolve->forcing_rows++;
					presolve->forced_columns += tightened;
				}
				changed = 1;
			}
		}
	} while (changed);
	duplicate_status = presolve_remove_parallel_rows(model, row_active,
		column_active, row_degree, column_degree, adjusted_rhs,
		tolerance, &duplicate_removed);
	if (duplicate_status == lp_simplex_EXIT_FAILURE)
		goto failure;
	if (duplicate_status == lp_simplex_Infeasibility) {
		presolve->terminal = 1;
		presolve->terminal_status = duplicate_status;
		goto finish;
	}
	if (duplicate_removed != 0) {
		presolve->removed_rows += duplicate_removed;
		presolve->duplicate_rows += duplicate_removed;
		goto reduce_again;
	}
	for (j = 0; j < model->n; j++)
		if (column_active[j]) {
			column_kept[j] = kept_columns;
			presolve->column_map[kept_columns++] = j;
		}
	for (i = 0; i < model->m; i++) {
		if (row_active[i]) {
			row_kept[i] = kept_rows;
			presolve->row_map[kept_rows++] = i;
		} else row_kept[i] = -1;
	}
	if (kept_columns == 0) {
		presolve->terminal = 1;
		presolve->terminal_status = lp_simplex_Success;
		goto finish;
	}
	if (kept_rows == 0) {
		presolve->terminal = 1;
		presolve->terminal_status = lp_simplex_Success;
		goto finish;
	}
	presolve->reduced = presolve_build_sparse_model(model, kept_rows,
		kept_columns, presolve->column_map, column_kept,
		presolve->row_map, row_kept, adjusted_rhs, bounds);
	if (presolve->reduced == NULL)
		goto failure;
	if (simplex_presolve_substitute_doubletons(presolve, tolerance) ==
	    lp_simplex_EXIT_FAILURE)
		goto failure;
finish:
	lp_simplex_free(column_kept);
	lp_simplex_free(row_kept);
	lp_simplex_free(adjusted_rhs);
	lp_simplex_free(bounds);
	lp_simplex_free(column_active);
	lp_simplex_free(row_active);
	lp_simplex_free(column_degree);
	lp_simplex_free(row_degree);
	return lp_simplex_EXIT_SUCCESS;
failure:
	lp_simplex_free(column_kept);
	lp_simplex_free(row_kept);
	lp_simplex_free(adjusted_rhs);
	lp_simplex_free(bounds);
	lp_simplex_free(column_active);
	lp_simplex_free(row_active);
	lp_simplex_free(column_degree);
	lp_simplex_free(row_degree);
	simplex_presolve_destroy(presolve);
	return lp_simplex_EXIT_FAILURE;
}


void simplex_presolve_postsolve(
		const struct simplex_Presolve *presolve,
		const double *reduced_x, double *original_x)
{
	int j;
	for (j = 0; j < presolve->original->n; j++)
		original_x[j] = presolve->eliminated[j]
			? presolve->eliminated_value[j] : 0.;
	if (reduced_x != NULL && presolve->reduced != NULL)
		for (j = 0; j < presolve->reduced->columns; j++)
			original_x[presolve->column_map[j]] = reduced_x[j];
	for (j = presolve->substitution_count - 1; j >= 0; j--) {
		const struct simplex_PresolveSubstitution *record =
			presolve->substitution + j;
		double value = record->constant;
		int k;
		for (k = 0; k < record->term_count; k++) {
			int position = record->term_start + k;
			value += presolve->substitution_term_multiplier[position] *
				original_x[presolve->substitution_term_column[position]];
		}
		original_x[record->column] = value;
	}
}


void simplex_presolve_destroy(struct simplex_Presolve *presolve)
{
	if (presolve == NULL)
		return;
	simplex_problem_free(presolve->reduced);
	lp_simplex_free(presolve->column_map);
	lp_simplex_free(presolve->row_map);
	lp_simplex_free(presolve->eliminated_value);
	lp_simplex_free(presolve->eliminated);
	lp_simplex_free(presolve->substitution);
	lp_simplex_free(presolve->substitution_term_column);
	lp_simplex_free(presolve->substitution_term_multiplier);
	presolve->reduced = NULL;
	presolve->column_map = NULL;
	presolve->row_map = NULL;
	presolve->eliminated_value = NULL;
	presolve->eliminated = NULL;
	presolve->substitution = NULL;
	presolve->substitution_term_column = NULL;
	presolve->substitution_term_multiplier = NULL;
}
