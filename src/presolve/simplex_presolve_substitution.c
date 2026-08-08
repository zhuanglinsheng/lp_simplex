/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#include "simplex_presolve_substitution.h"
#include "simplex_presolve.h"
#include "simplex_presolve_queue.h"
#include "simplex_presolve_record.h"
#include "simplex_presolve_rules.h"
#include "utils.h"

#include <lp_simplex/status.h>

#include <math.h>


struct presolve_MutableRow {
	int *column;
	double *value;
	int count;
	int capacity;
	unsigned char active;
};


struct presolve_Incidence {
	int *row;
	int count;
	int capacity;
};


static int substitution_append_incidence(
		struct presolve_Incidence *incidence, const int row)
{
	int capacity, k;
	int *grown;
	/* Removed coefficients leave a harmless tombstone in the incidence list.
	 * Reusing that row/column pair must not append another copy: duplicates
	 * inflate scans and can repeatedly schedule the same structural rule. */
	for (k = 0; k < incidence->count; k++)
		if (incidence->row[k] == row)
			return lp_simplex_EXIT_SUCCESS;
	if (incidence->count < incidence->capacity) {
		incidence->row[incidence->count++] = row;
		return lp_simplex_EXIT_SUCCESS;
	}
	capacity = incidence->capacity == 0 ? 4 : incidence->capacity * 2;
	grown = (int *)lp_simplex_realloc(incidence->row,
		(size_t)capacity * sizeof(int));
	if (grown == NULL)
		return lp_simplex_EXIT_FAILURE;
	incidence->row = grown;
	incidence->capacity = capacity;
	incidence->row[incidence->count++] = row;
	return lp_simplex_EXIT_SUCCESS;
}


static int substitution_find_column(
		const struct presolve_MutableRow *row, const int column)
{
	int k;
	for (k = 0; k < row->count; k++)
		if (row->column[k] == column)
			return k;
	return -1;
}


static int substitution_add_to_row(
		struct presolve_MutableRow *row, const int column,
		const double delta, int *inserted, int *removed)
{
	int position = substitution_find_column(row, column);
	*inserted = 0;
	*removed = 0;
	if (position >= 0) {
		row->value[position] += delta;
		if (row->value[position] == 0.) {
			row->count--;
			row->column[position] = row->column[row->count];
			row->value[position] = row->value[row->count];
			*removed = 1;
		}
		return lp_simplex_EXIT_SUCCESS;
	}
	if (delta == 0.)
		return lp_simplex_EXIT_SUCCESS;
	if (row->count == row->capacity) {
		int capacity = row->capacity == 0 ? 4 : row->capacity * 2;
		int *columns = (int *)lp_simplex_realloc(row->column,
			(size_t)capacity * sizeof(int));
		double *values;
		if (columns == NULL)
			return lp_simplex_EXIT_FAILURE;
		row->column = columns;
		values = (double *)lp_simplex_realloc(row->value,
			(size_t)capacity * sizeof(double));
		if (values == NULL)
			return lp_simplex_EXIT_FAILURE;
		row->value = values;
		row->capacity = capacity;
	}
	row->column[row->count] = column;
	row->value[row->count++] = delta;
	*inserted = 1;
	return lp_simplex_EXIT_SUCCESS;
}


static void substitution_remove_from_row(
		struct presolve_MutableRow *row, const int position)
{
	row->count--;
	row->column[position] = row->column[row->count];
	row->value[position] = row->value[row->count];
}


static int substitution_transfer_bounds(
		const struct optm_VariableBound *eliminated,
		struct optm_VariableBound *kept, const double constant,
		const double multiplier, const double tolerance)
{
	double lower = __lp_simplex_NINF__, upper = __lp_simplex_INF__;
	double margin;
	if (simplex_presolve_has_lower(eliminated)) {
		double value = (eliminated->lb - constant) / multiplier;
		if (multiplier > 0.) lower = value;
		else upper = value;
	}
	if (simplex_presolve_has_upper(eliminated)) {
		double value = (eliminated->ub - constant) / multiplier;
		if (multiplier > 0.) upper = value;
		else lower = value;
	}
	if (isfinite(lower) &&
	    (!simplex_presolve_has_lower(kept) || lower > kept->lb))
		simplex_presolve_set_lower(kept, lower);
	if (isfinite(upper) &&
	    (!simplex_presolve_has_upper(kept) || upper < kept->ub))
		simplex_presolve_set_upper(kept, upper);
	margin = tolerance * (1. +
		(simplex_presolve_has_lower(kept) ? __lp_simplex_ABS__(kept->lb) : 0.) +
		(simplex_presolve_has_upper(kept) ? __lp_simplex_ABS__(kept->ub) : 0.));
	return simplex_presolve_has_lower(kept) &&
		simplex_presolve_has_upper(kept) && kept->lb > kept->ub + margin
		? lp_simplex_Infeasibility : lp_simplex_Success;
}


static void substitution_destroy_rows(
		struct presolve_MutableRow *rows, const int count,
		struct presolve_Incidence *incidence, const int columns)
{
	int i;
	if (rows != NULL)
		for (i = 0; i < count; i++) {
			lp_simplex_free(rows[i].column);
			lp_simplex_free(rows[i].value);
		}
	if (incidence != NULL)
		for (i = 0; i < columns; i++)
			lp_simplex_free(incidence[i].row);
	lp_simplex_free(rows);
	lp_simplex_free(incidence);
}


static void substitution_expression_bound_flags(
		const struct simplex_Problem *problem,
		const struct presolve_MutableRow *row,
		const int eliminated_position, const double constant,
		const double eliminated_coefficient, const double tolerance,
		int *lower_implied, int *upper_implied)
{
	const struct optm_VariableBound *eliminated =
		problem->bounds + row->column[eliminated_position];
	long double minimum = constant, maximum = constant, magnitude = 0.;
	int minimum_finite = 1, maximum_finite = 1;
	int k;
	*lower_implied = !simplex_presolve_has_lower(eliminated);
	*upper_implied = !simplex_presolve_has_upper(eliminated);
	for (k = 0; k < row->count; k++) {
		const struct optm_VariableBound *bound;
		double multiplier;
		if (k == eliminated_position)
			continue;
		bound = problem->bounds + row->column[k];
		multiplier = -row->value[k] / eliminated_coefficient;
		if (multiplier > 0.) {
			if (simplex_presolve_has_lower(bound)) {
				long double contribution = multiplier * bound->lb;
				minimum += contribution;
				magnitude += fabsl(contribution);
			} else minimum_finite = 0;
			if (simplex_presolve_has_upper(bound)) {
				long double contribution = multiplier * bound->ub;
				maximum += contribution;
				magnitude += fabsl(contribution);
			} else maximum_finite = 0;
		} else {
			if (simplex_presolve_has_upper(bound)) {
				long double contribution = multiplier * bound->ub;
				minimum += contribution;
				magnitude += fabsl(contribution);
			} else minimum_finite = 0;
			if (simplex_presolve_has_lower(bound)) {
				long double contribution = multiplier * bound->lb;
				maximum += contribution;
				magnitude += fabsl(contribution);
			} else maximum_finite = 0;
		}
	}
	{
		long double margin = (long double)tolerance *
			(1. + fabsl((long double)constant) + magnitude);
		if (simplex_presolve_has_lower(eliminated) && minimum_finite &&
		    minimum >= eliminated->lb - margin)
			*lower_implied = 1;
		if (simplex_presolve_has_upper(eliminated) && maximum_finite &&
		    maximum <= eliminated->ub + margin)
			*upper_implied = 1;
	}
}


static int substitution_remove_singleton_columns(
		struct simplex_Presolve *presolve,
		struct simplex_Problem *problem,
		struct presolve_MutableRow *rows,
	struct presolve_Incidence *incidence,
	unsigned char *column_active, int *column_degree,
	struct simplex_PresolveQueue *column_queue,
	const double tolerance,
		int *active_rows, int *active_columns, int *substitutions,
		int *projections)
{
	int j;
	while (*active_rows > 1 && *active_columns > 1 &&
	       simplex_presolve_queue_pop(column_queue, &j)) {
			int degree, position, row_index = -1, t;
			int lower_implied, upper_implied, project_lower;
			struct presolve_MutableRow *row;
			double coefficient, constant, maximum = 0.;
			if (!column_active[j] ||
			    column_degree[j] < 1 || column_degree[j] > 2)
				continue;
			degree = column_degree[j];
			/* A degree-two column can be eliminated through either incident
			 * equality without increasing the total matrix nonzero count.  Pick
			 * the shortest equality whose expression implies both column bounds.
			 * Degree-one columns also allow one-sided projection. */
			for (t = 0; t < incidence[j].count; t++) {
				int candidate = incidence[j].row[t];
				struct presolve_MutableRow *candidate_row;
				int candidate_position, candidate_lower, candidate_upper, u;
				double candidate_coefficient, candidate_constant;
				double candidate_maximum = 0.;
				if (!rows[candidate].active ||
				    problem->row_type[candidate] != optm_CONS_T_EQ)
					continue;
				candidate_row = rows + candidate;
				candidate_position = substitution_find_column(
					candidate_row, j);
				if (candidate_position < 0)
					continue;
				candidate_coefficient =
					candidate_row->value[candidate_position];
				for (u = 0; u < candidate_row->count; u++)
					candidate_maximum = __lp_simplex_MAX__(
						candidate_maximum,
						__lp_simplex_ABS__(candidate_row->value[u]));
				if (__lp_simplex_ABS__(candidate_coefficient) <
				    1e-10 * __lp_simplex_MAX__(1., candidate_maximum))
					continue;
				candidate_constant = problem->rhs[candidate] /
					candidate_coefficient;
				if (!isfinite(candidate_constant))
					continue;
				substitution_expression_bound_flags(problem, candidate_row,
					candidate_position, candidate_constant,
					candidate_coefficient, tolerance,
					&candidate_lower, &candidate_upper);
				if (degree == 2 &&
				    !(candidate_lower && candidate_upper))
					continue;
				if (!candidate_lower && !candidate_upper)
					continue;
				if (row_index < 0 ||
				    candidate_row->count < rows[row_index].count)
					row_index = candidate;
			}
			if (row_index < 0)
				continue;
			row = rows + row_index;
			position = substitution_find_column(row, j);
			coefficient = row->value[position];
			for (t = 0; t < row->count; t++)
				maximum = __lp_simplex_MAX__(maximum,
					__lp_simplex_ABS__(row->value[t]));
			if (__lp_simplex_ABS__(coefficient) <
			    1e-10 * __lp_simplex_MAX__(1., maximum))
				continue;
			constant = problem->rhs[row_index] / coefficient;
			if (!isfinite(constant))
				continue;
			substitution_expression_bound_flags(problem, row, position,
				constant, coefficient, tolerance,
				&lower_implied, &upper_implied);
			if (!lower_implied && !upper_implied)
				continue;
			/* A one-sided projection replaces the singleton basis column by a
			 * zero-cost logical slack.  For a costed singleton this can destroy
			 * the natural dual-feasible crash even though the primal projection
			 * is algebraically valid.  Keep such columns unless both bounds are
			 * implied and the row can be removed by full substitution. */
			if (!(lower_implied && upper_implied) &&
			    problem->objective[j] != 0.)
				continue;
			if (simplex_presolve_record_begin(presolve,
					presolve->column_map[j], constant,
					row->count - 1) == NULL)
				return lp_simplex_EXIT_FAILURE;
			for (t = 0; t < row->count; t++)
				if (t != position) {
					int kept = row->column[t];
					double multiplier = -row->value[t] / coefficient;
					simplex_presolve_record_add_term(presolve,
						presolve->column_map[kept], multiplier);
					problem->objective[kept] +=
						problem->objective[j] * multiplier;
				}
			project_lower = !lower_implied;
			if (lower_implied && upper_implied) {
				/* Substitute the pivot-row expression into every other row
				 * containing the low-degree column. */
				for (t = 0; t < incidence[j].count; t++) {
					int affected = incidence[j].row[t];
					struct presolve_MutableRow *target = rows + affected;
					int target_position, u;
					double target_coefficient;
					if (!target->active || affected == row_index)
						continue;
					target_position = substitution_find_column(target, j);
					if (target_position < 0)
						continue;
					target_coefficient = target->value[target_position];
					substitution_remove_from_row(target, target_position);
					column_degree[j]--;
					problem->rhs[affected] -= target_coefficient * constant;
					for (u = 0; u < row->count; u++)
						if (u != position) {
							int kept = row->column[u];
							int inserted, removed;
							double multiplier =
								-row->value[u] / coefficient;
							if (substitution_add_to_row(target, kept,
									target_coefficient * multiplier,
									&inserted, &removed) ==
							    lp_simplex_EXIT_FAILURE)
								return lp_simplex_EXIT_FAILURE;
							column_degree[kept] += inserted - removed;
							(void)simplex_presolve_queue_push(
								column_queue, kept);
							if (inserted && substitution_append_incidence(
									incidence + kept, affected) ==
							    lp_simplex_EXIT_FAILURE)
								return lp_simplex_EXIT_FAILURE;
						}
				}
				for (t = 0; t < row->count; t++)
					if (column_active[row->column[t]] &&
					    column_degree[row->column[t]] > 0) {
						column_degree[row->column[t]]--;
						(void)simplex_presolve_queue_push(
							column_queue, row->column[t]);
					}
				row->active = 0;
				(*active_rows)--;
				(*substitutions)++;
				if (degree == 1)
					presolve->stats.singleton_column_rows++;
				else
					presolve->stats.implied_free_columns++;
			} else {
				substitution_remove_from_row(row, position);
				for (t = 0; t < row->count; t++)
					row->value[t] = -row->value[t] / coefficient;
				problem->rhs[row_index] = project_lower
					? problem->bounds[j].lb - constant
					: problem->bounds[j].ub - constant;
				problem->row_type[row_index] = project_lower
					? optm_CONS_T_GE : optm_CONS_T_LE;
				(*projections)++;
				presolve->stats.singleton_projection_columns++;
			}
			column_active[j] = 0;
			column_degree[j] = 0;
			presolve->eliminated[presolve->column_map[j]] = 1;
			(*active_columns)--;
	}
	return lp_simplex_EXIT_SUCCESS;
}


static void substitution_eliminate_fixed_column(
		struct simplex_Presolve *presolve,
		struct simplex_Problem *problem,
		struct presolve_MutableRow *rows,
		struct presolve_Incidence *incidence,
		unsigned char *column_active, int *column_degree,
		struct simplex_PresolveQueue *queue, const int column,
		const double value, int *active_columns, int *removed_columns)
{
	int k;
	for (k = 0; k < incidence[column].count; k++) {
		int row_index = incidence[column].row[k];
		struct presolve_MutableRow *row = rows + row_index;
		int position;
		if (!row->active)
			continue;
		position = substitution_find_column(row, column);
		if (position < 0)
			continue;
		problem->rhs[row_index] -= row->value[position] * value;
		substitution_remove_from_row(row, position);
		if (column_degree[column] > 0)
			column_degree[column]--;
		if (row->count <= 1 && queue != NULL)
			(void)simplex_presolve_queue_push(queue, row_index);
	}
	column_active[column] = 0;
	column_degree[column] = 0;
	presolve->eliminated[presolve->column_map[column]] = 1;
	presolve->eliminated_value[presolve->column_map[column]] = value;
	(*active_columns)--;
	(*removed_columns)++;
}


/* Close the cheap presolve rules exposed by substitution without rebuilding
 * an intermediate model.  Row/column degrees are updated in place and every
 * newly singleton row is queued at most once when it crosses degree one. */
static int substitution_close_fixed_point(
		struct simplex_Presolve *presolve,
		struct simplex_Problem *problem,
		struct presolve_MutableRow *rows,
		struct presolve_Incidence *incidence,
		unsigned char *column_active, int *column_degree,
		const double tolerance,
		int *active_rows, int *active_columns,
		int *removed_rows, int *removed_columns)
{
	struct simplex_PresolveQueue queue;
	int i, j, row_index, activity_changed;
	int result = lp_simplex_Success;
	if (simplex_presolve_queue_init(&queue, problem->rows) ==
	    lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
restart:
	simplex_presolve_queue_clear(&queue);
	for (j = 0; j < problem->columns; j++) {
		double value;
		int fixed, status;
		if (!column_active[j])
			continue;
		fixed = problem->bounds[j].b_type == optm_BOUND_T_BS &&
			problem->bounds[j].lb == problem->bounds[j].ub;
		if (fixed)
			value = problem->bounds[j].lb;
		else if (column_degree[j] == 0) {
			status = simplex_presolve_choose_empty_column(
				problem->objective[j], problem->bounds + j, &value);
			if (status != lp_simplex_Success) {
				result = status;
				goto finish;
			}
			presolve->stats.empty_columns++;
			fixed = 1;
		}
		if (fixed) {
			presolve->stats.fixed_columns += column_degree[j] != 0;
			substitution_eliminate_fixed_column(presolve, problem, rows,
				incidence, column_active, column_degree, NULL,
				j, value, active_columns, removed_columns);
		}
	}
	for (i = 0; i < problem->rows; i++)
		if (rows[i].active && rows[i].count <= 1)
			(void)simplex_presolve_queue_push(&queue, i);
	while (simplex_presolve_queue_pop(&queue, &row_index)) {
		struct presolve_MutableRow *row;
		double coefficient;
		struct optm_VariableBound *bound;
		int status, tightened;
		if (!rows[row_index].active || rows[row_index].count > 1)
			continue;
		row = rows + row_index;
		if (row->count == 0) {
			if (!simplex_presolve_empty_row_feasible(
					problem->row_type[row_index],
					problem->rhs[row_index], tolerance)) {
				result = lp_simplex_Infeasibility;
				goto finish;
			}
			row->active = 0;
			(*active_rows)--;
			(*removed_rows)++;
			continue;
		}
		j = row->column[0];
		coefficient = row->value[0];
		if (!column_active[j] || coefficient == 0.)
			continue;
		bound = problem->bounds + j;
		status = simplex_presolve_tighten_singleton(bound,
			problem->row_type[row_index], coefficient,
			problem->rhs[row_index], tolerance, &tightened);
		if (status != lp_simplex_Success) {
			result = status;
			goto finish;
		}
		if (problem->row_type[row_index] == optm_CONS_T_EQ)
			presolve->stats.singleton_columns++;
		row->active = 0;
		if (column_degree[j] > 0)
			column_degree[j]--;
		(*active_rows)--;
		(*removed_rows)++;
		presolve->stats.singleton_rows++;
		presolve->stats.tightened_bounds += tightened;
		if (bound->b_type == optm_BOUND_T_BS && bound->lb == bound->ub) {
			presolve->stats.fixed_columns++;
			substitution_eliminate_fixed_column(presolve, problem, rows,
				incidence, column_active, column_degree, &queue,
				j, bound->lb, active_columns,
				removed_columns);
		}
	}
	activity_changed = 0;
	for (i = 0; i < problem->rows; i++) {
		struct presolve_MutableRow *row = rows + i;
		long double minimum = 0., maximum = 0., magnitude = 0.;
		int minimum_infinite = 0, maximum_infinite = 0;
		int force_minimum, force_maximum, redundant = 0, fixed = 0;
		long double margin;
		if (!row->active || row->count < 2)
			continue;
		for (j = 0; j < row->count; j++) {
			const struct optm_VariableBound *bound =
				problem->bounds + row->column[j];
			double coefficient = row->value[j];
			double contribution;
			if (coefficient > 0.) {
				if (simplex_presolve_has_lower(bound)) {
					contribution = coefficient * bound->lb;
					minimum += contribution;
					magnitude += fabsl((long double)contribution);
				} else minimum_infinite++;
				if (simplex_presolve_has_upper(bound)) {
					contribution = coefficient * bound->ub;
					maximum += contribution;
					magnitude += fabsl((long double)contribution);
				} else maximum_infinite++;
			} else {
				if (simplex_presolve_has_upper(bound)) {
					contribution = coefficient * bound->ub;
					minimum += contribution;
					magnitude += fabsl((long double)contribution);
				} else minimum_infinite++;
				if (simplex_presolve_has_lower(bound)) {
					contribution = coefficient * bound->lb;
					maximum += contribution;
					magnitude += fabsl((long double)contribution);
				} else maximum_infinite++;
			}
		}
		margin = (long double)tolerance *
			(1. + fabsl((long double)problem->rhs[i]) + magnitude);
		if ((problem->row_type[i] == optm_CONS_T_GE ||
		     problem->row_type[i] == optm_CONS_T_EQ) &&
		    maximum_infinite == 0 &&
		    maximum < (long double)problem->rhs[i] - margin) {
			result = lp_simplex_Infeasibility;
			goto finish;
		}
		if ((problem->row_type[i] == optm_CONS_T_LE ||
		     problem->row_type[i] == optm_CONS_T_EQ) &&
		    minimum_infinite == 0 &&
		    minimum > (long double)problem->rhs[i] + margin) {
			result = lp_simplex_Infeasibility;
			goto finish;
		}
		if (problem->row_type[i] == optm_CONS_T_GE &&
		    minimum_infinite == 0 &&
		    minimum >= (long double)problem->rhs[i] + margin)
			redundant = 1;
		if (problem->row_type[i] == optm_CONS_T_LE &&
		    maximum_infinite == 0 &&
		    maximum <= (long double)problem->rhs[i] - margin)
			redundant = 1;
		if (redundant) {
			for (j = 0; j < row->count; j++)
				if (column_degree[row->column[j]] > 0)
					column_degree[row->column[j]]--;
			row->active = 0;
			(*active_rows)--;
			(*removed_rows)++;
			presolve->stats.redundant_rows++;
			activity_changed = 1;
			continue;
		}
		force_minimum =
			(problem->row_type[i] == optm_CONS_T_LE ||
			 problem->row_type[i] == optm_CONS_T_EQ) &&
			minimum_infinite == 0 &&
			minimum == (long double)problem->rhs[i];
		force_maximum =
			(problem->row_type[i] == optm_CONS_T_GE ||
			 problem->row_type[i] == optm_CONS_T_EQ) &&
			maximum_infinite == 0 &&
			maximum == (long double)problem->rhs[i];
		if (!force_minimum && !force_maximum)
			continue;
		for (j = 0; j < row->count; j++) {
			int column = row->column[j];
			struct optm_VariableBound *bound = problem->bounds + column;
			double value;
			if (bound->b_type == optm_BOUND_T_BS && bound->lb == bound->ub)
				continue;
			if (force_minimum)
				value = row->value[j] > 0. ? bound->lb : bound->ub;
			else
				value = row->value[j] > 0. ? bound->ub : bound->lb;
			simplex_presolve_set_lower(bound, value);
			simplex_presolve_set_upper(bound, value);
			fixed++;
		}
		if (fixed != 0) {
			presolve->stats.forcing_rows++;
			presolve->stats.forced_columns += fixed;
			presolve->stats.tightened_bounds += fixed;
			activity_changed = 1;
		}
	}
	if (activity_changed)
		goto restart;
	for (j = 0; j < problem->columns; j++)
		if (column_active[j] && column_degree[j] == 0) {
			double value;
			int status = simplex_presolve_choose_empty_column(
				problem->objective[j], problem->bounds + j, &value);
			if (status != lp_simplex_Success) {
				result = status;
				goto finish;
			}
			presolve->stats.empty_columns++;
			substitution_eliminate_fixed_column(presolve, problem, rows,
				incidence, column_active, column_degree, NULL,
				j, value, active_columns, removed_columns);
		}
finish:
	simplex_presolve_queue_destroy(&queue);
	return result;
}


struct substitution_Workspace {
	struct simplex_Problem *problem;
	struct presolve_MutableRow *rows;
	struct presolve_Incidence *incidence;
	unsigned char *column_active;
	int *column_degree;
	struct simplex_PresolveQueue doubleton_queue;
	struct simplex_PresolveQueue singleton_queue;
	int active_rows;
	int active_columns;
	int substitutions;
	int doubletons;
	int projections;
	int closed_rows;
	int closed_columns;
};


static void substitution_workspace_destroy(
		struct substitution_Workspace *workspace)
{
	if (workspace == NULL)
		return;
	if (workspace->problem != NULL)
		substitution_destroy_rows(workspace->rows,
			workspace->problem->rows, workspace->incidence,
			workspace->problem->columns);
	else
		substitution_destroy_rows(workspace->rows, 0,
			workspace->incidence, 0);
	lp_simplex_free(workspace->column_active);
	lp_simplex_free(workspace->column_degree);
	simplex_presolve_queue_destroy(&workspace->doubleton_queue);
	simplex_presolve_queue_destroy(&workspace->singleton_queue);
	lp_simplex_memset(workspace, 0, sizeof(*workspace));
}


static int substitution_workspace_create(
		struct simplex_Problem *problem,
		struct substitution_Workspace *workspace)
{
	int i, j, k;
	lp_simplex_memset(workspace, 0, sizeof(*workspace));
	workspace->problem = problem;
	workspace->rows = (struct presolve_MutableRow *)lp_simplex_malloc(
		(size_t)problem->rows * sizeof(*workspace->rows));
	workspace->incidence = (struct presolve_Incidence *)lp_simplex_malloc(
		(size_t)problem->columns * sizeof(*workspace->incidence));
	workspace->column_active = (unsigned char *)lp_simplex_malloc(
		(size_t)problem->columns * sizeof(unsigned char));
	workspace->column_degree = (int *)lp_simplex_malloc(
		(size_t)problem->columns * sizeof(int));
	if (workspace->rows != NULL)
		lp_simplex_memset(workspace->rows, 0,
			(size_t)problem->rows * sizeof(*workspace->rows));
	if (workspace->incidence != NULL)
		lp_simplex_memset(workspace->incidence, 0,
			(size_t)problem->columns * sizeof(*workspace->incidence));
	if (workspace->rows == NULL || workspace->incidence == NULL ||
	    workspace->column_active == NULL || workspace->column_degree == NULL)
		return lp_simplex_EXIT_FAILURE;
	if (simplex_presolve_queue_init(&workspace->doubleton_queue,
		problem->rows) == lp_simplex_EXIT_FAILURE ||
	    simplex_presolve_queue_init(&workspace->singleton_queue,
		problem->columns) == lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	lp_simplex_memset(workspace->column_active, 1,
		(size_t)problem->columns * sizeof(unsigned char));
	lp_simplex_memset(workspace->column_degree, 0,
		(size_t)problem->columns * sizeof(int));
	/* The immutable input already provides exact column degrees.  Reserve the
	 * initial incidence storage once instead of growing every column through
	 * several reallocations while copying rows. */
	for (j = 0; j < problem->columns; j++) {
		int count = problem->matrix.column_start[j + 1] -
			problem->matrix.column_start[j];
		workspace->incidence[j].capacity = count;
		workspace->incidence[j].row = (int *)lp_simplex_malloc(
			(size_t)count * sizeof(int));
		if (count != 0 && workspace->incidence[j].row == NULL)
			return lp_simplex_EXIT_FAILURE;
	}
	for (i = 0; i < problem->rows; i++) {
		int count = problem->matrix.row_start[i + 1] -
			problem->matrix.row_start[i];
		workspace->rows[i].active = 1;
		workspace->rows[i].count = workspace->rows[i].capacity = count;
		workspace->rows[i].column = (int *)lp_simplex_malloc(
			(size_t)count * sizeof(int));
		workspace->rows[i].value = (double *)lp_simplex_malloc(
			(size_t)count * sizeof(double));
		if (count != 0 && (workspace->rows[i].column == NULL ||
		    workspace->rows[i].value == NULL))
			return lp_simplex_EXIT_FAILURE;
		for (k = 0; k < count; k++) {
			j = problem->matrix.column_index[
				problem->matrix.row_start[i] + k];
			workspace->rows[i].column[k] = j;
			workspace->rows[i].value[k] = problem->matrix.row_value[
				problem->matrix.row_start[i] + k];
			workspace->column_degree[j]++;
			if (substitution_append_incidence(
					workspace->incidence + j, i) ==
			    lp_simplex_EXIT_FAILURE)
				return lp_simplex_EXIT_FAILURE;
		}
	}
	workspace->active_rows = problem->rows;
	workspace->active_columns = problem->columns;
	for (j = 0; j < problem->columns; j++)
		(void)simplex_presolve_queue_push(&workspace->singleton_queue, j);
	return lp_simplex_EXIT_SUCCESS;
}


static void substitution_rebuild_doubleton_queue(
		struct substitution_Workspace *workspace)
{
	int i;
	simplex_presolve_queue_clear(&workspace->doubleton_queue);
	for (i = 0; i < workspace->problem->rows; i++)
		if (workspace->rows[i].active && workspace->rows[i].count == 2 &&
		    workspace->problem->row_type[i] == optm_CONS_T_EQ)
			(void)simplex_presolve_queue_push(
				&workspace->doubleton_queue, i);
}


static int substitution_eliminate_doubletons(
		struct simplex_Presolve *presolve,
		struct substitution_Workspace *workspace, const double tolerance)
{
	struct simplex_Problem *problem = workspace->problem;
	int row_index;
	while (workspace->active_rows > 1 && workspace->active_columns > 1 &&
	       simplex_presolve_queue_pop(
		       &workspace->doubleton_queue, &row_index)) {
		struct presolve_MutableRow *row = workspace->rows + row_index;
		int first, second, eliminate, keep, eliminate_position;
		int first_degree, second_degree, k;
		double eliminate_coefficient, keep_coefficient;
		double constant, multiplier;
		if (!row->active || row->count != 2 ||
		    problem->row_type[row_index] != optm_CONS_T_EQ)
			continue;
		first = row->column[0];
		second = row->column[1];
		if (!workspace->column_active[first] ||
		    !workspace->column_active[second])
			continue;
		first_degree = workspace->column_degree[first];
		second_degree = workspace->column_degree[second];
		eliminate = first_degree <= second_degree ? first : second;
		keep = eliminate == first ? second : first;
		eliminate_position = eliminate == first ? 0 : 1;
		eliminate_coefficient = row->value[eliminate_position];
		keep_coefficient = row->value[1 - eliminate_position];
		if (__lp_simplex_ABS__(eliminate_coefficient) < 1e-10 *
		    __lp_simplex_MAX__(1., __lp_simplex_ABS__(keep_coefficient)))
			continue;
		constant = problem->rhs[row_index] / eliminate_coefficient;
		multiplier = -keep_coefficient / eliminate_coefficient;
		if (!isfinite(constant) || !isfinite(multiplier) || multiplier == 0.)
			continue;
		if (substitution_transfer_bounds(problem->bounds + eliminate,
			problem->bounds + keep, constant, multiplier, tolerance) !=
		    lp_simplex_Success)
			return lp_simplex_Infeasibility;
		if (simplex_presolve_record_begin(presolve,
			presolve->column_map[eliminate], constant, 1) == NULL)
			return lp_simplex_EXIT_FAILURE;
		simplex_presolve_record_add_term(presolve,
			presolve->column_map[keep], multiplier);
		problem->objective[keep] += problem->objective[eliminate] * multiplier;
		for (k = 0; k < workspace->incidence[eliminate].count; k++) {
			int affected = workspace->incidence[eliminate].row[k];
			struct presolve_MutableRow *target = workspace->rows + affected;
			int position, inserted, removed;
			double coefficient;
			if (!target->active || affected == row_index)
				continue;
			position = substitution_find_column(target, eliminate);
			if (position < 0)
				continue;
			coefficient = target->value[position];
			substitution_remove_from_row(target, position);
			workspace->column_degree[eliminate]--;
			problem->rhs[affected] -= coefficient * constant;
			if (substitution_add_to_row(target, keep,
				coefficient * multiplier, &inserted, &removed) ==
			    lp_simplex_EXIT_FAILURE)
				return lp_simplex_EXIT_FAILURE;
			workspace->column_degree[keep] += inserted - removed;
			(void)simplex_presolve_queue_push(
				&workspace->singleton_queue, keep);
			if (inserted && substitution_append_incidence(
					workspace->incidence + keep, affected) ==
			    lp_simplex_EXIT_FAILURE)
				return lp_simplex_EXIT_FAILURE;
			if (problem->row_type[affected] == optm_CONS_T_EQ &&
			    target->count == 2)
				(void)simplex_presolve_queue_push(
					&workspace->doubleton_queue, affected);
		}
		if (workspace->column_degree[first] > 0) {
			workspace->column_degree[first]--;
			(void)simplex_presolve_queue_push(
				&workspace->singleton_queue, first);
		}
		if (workspace->column_degree[second] > 0) {
			workspace->column_degree[second]--;
			(void)simplex_presolve_queue_push(
				&workspace->singleton_queue, second);
		}
		row->active = 0;
		workspace->column_active[eliminate] = 0;
		workspace->column_degree[eliminate] = 0;
		presolve->eliminated[presolve->column_map[eliminate]] = 1;
		workspace->active_rows--;
		workspace->active_columns--;
		workspace->substitutions++;
		workspace->doubletons++;
	}
	return lp_simplex_Success;
}


static struct simplex_Problem *substitution_compact_problem(
		struct simplex_Presolve *presolve,
		const struct substitution_Workspace *workspace)
{
	struct simplex_Problem *problem = workspace->problem;
	struct simplex_Problem *reduced = NULL;
	int *column_position = NULL, *cursor = NULL;
	int i, j, k, nonzeros = 0, next = 0, new_row = 0, new_column = 0;
	column_position = (int *)lp_simplex_malloc(
		(size_t)problem->columns * sizeof(int));
	if (column_position == NULL)
		return NULL;
	for (j = 0; j < problem->columns; j++) {
		column_position[j] = -1;
		if (workspace->column_active[j]) {
			column_position[j] = new_column++;
			presolve->column_map[new_column - 1] = presolve->column_map[j];
		}
	}
	for (i = 0; i < problem->rows; i++)
		if (workspace->rows[i].active) {
			presolve->row_map[new_row++] = presolve->row_map[i];
			nonzeros += workspace->rows[i].count;
		}
	reduced = simplex_problem_create_sparse(workspace->active_rows,
		workspace->active_columns, nonzeros);
	if (reduced == NULL)
		goto failure;
	for (j = 0; j < problem->columns; j++)
		if (workspace->column_active[j]) {
			int destination = column_position[j];
			reduced->objective[destination] = problem->objective[j];
			reduced->bounds[destination] = problem->bounds[j];
		}
	new_row = 0;
	for (i = 0; i < problem->rows; i++)
		if (workspace->rows[i].active) {
			reduced->rhs[new_row] = problem->rhs[i];
			reduced->row_type[new_row] = problem->row_type[i];
			reduced->matrix.row_start[new_row] = next;
			for (k = 0; k < workspace->rows[i].count; k++) {
				reduced->matrix.column_index[next] =
					column_position[workspace->rows[i].column[k]];
				reduced->matrix.row_value[next++] = workspace->rows[i].value[k];
			}
			new_row++;
		}
	reduced->matrix.row_start[workspace->active_rows] = next;
	lp_simplex_memset(reduced->matrix.column_start, 0,
		(size_t)(workspace->active_columns + 1) * sizeof(int));
	for (k = 0; k < nonzeros; k++)
		reduced->matrix.column_start[
			reduced->matrix.column_index[k] + 1]++;
	for (j = 0; j < workspace->active_columns; j++)
		reduced->matrix.column_start[j + 1] +=
			reduced->matrix.column_start[j];
	cursor = (int *)lp_simplex_malloc(
		(size_t)workspace->active_columns * sizeof(int));
	if (cursor == NULL)
		goto failure;
	lp_simplex_memcpy(cursor, reduced->matrix.column_start,
		(size_t)workspace->active_columns * sizeof(int));
	for (i = 0; i < workspace->active_rows; i++)
		for (k = reduced->matrix.row_start[i];
		     k < reduced->matrix.row_start[i + 1]; k++) {
			int column = reduced->matrix.column_index[k];
			int destination = cursor[column]++;
			reduced->matrix.row_index[destination] = i;
			reduced->matrix.value[destination] =
				reduced->matrix.row_value[k];
		}
	lp_simplex_free(cursor);
	lp_simplex_free(column_position);
	return reduced;
failure:
	lp_simplex_free(cursor);
	lp_simplex_free(column_position);
	simplex_problem_free(reduced);
	return NULL;
}


static int substitution_reduce_to_fixed_point(
		struct simplex_Presolve *presolve,
		struct substitution_Workspace *workspace, const double tolerance)
{
	struct simplex_Problem *problem = workspace->problem;
	for (;;) {
		int before = workspace->substitutions + workspace->projections +
			workspace->closed_rows + workspace->closed_columns;
		int j, status;
		if (substitution_remove_singleton_columns(presolve, problem,
			workspace->rows, workspace->incidence,
			workspace->column_active, workspace->column_degree,
			&workspace->singleton_queue, tolerance,
			&workspace->active_rows, &workspace->active_columns,
			&workspace->substitutions, &workspace->projections) ==
		    lp_simplex_EXIT_FAILURE)
			return lp_simplex_EXIT_FAILURE;
		substitution_rebuild_doubleton_queue(workspace);
		status = substitution_eliminate_doubletons(
			presolve, workspace, tolerance);
		if (status != lp_simplex_Success)
			return status;
		if (substitution_remove_singleton_columns(presolve, problem,
			workspace->rows, workspace->incidence,
			workspace->column_active, workspace->column_degree,
			&workspace->singleton_queue, tolerance,
			&workspace->active_rows, &workspace->active_columns,
			&workspace->substitutions, &workspace->projections) ==
		    lp_simplex_EXIT_FAILURE)
			return lp_simplex_EXIT_FAILURE;
		status = substitution_close_fixed_point(presolve, problem,
			workspace->rows, workspace->incidence,
			workspace->column_active, workspace->column_degree,
			tolerance, &workspace->active_rows,
			&workspace->active_columns, &workspace->closed_rows,
			&workspace->closed_columns);
		if (status != lp_simplex_Success)
			return status;
		if (before == workspace->substitutions + workspace->projections +
		    workspace->closed_rows + workspace->closed_columns ||
		    workspace->active_rows <= 1 || workspace->active_columns <= 1)
			return lp_simplex_Success;
		/* Cheap closure can expose fresh low-degree columns.  Begin another
		 * structural epoch without rebuilding the mutable matrix. */
		simplex_presolve_queue_clear(&workspace->singleton_queue);
		for (j = 0; j < problem->columns; j++)
			if (workspace->column_active[j])
				(void)simplex_presolve_queue_push(
					&workspace->singleton_queue, j);
	}
}


int simplex_presolve_substitute_doubletons(
		struct simplex_Presolve *presolve, const double tolerance)
{
	struct substitution_Workspace workspace;
	struct simplex_Problem *problem = presolve->reduced;
	struct simplex_Problem *reduced;
	int status;
	if (problem == NULL)
		return lp_simplex_EXIT_SUCCESS;
	if (substitution_workspace_create(problem, &workspace) ==
	    lp_simplex_EXIT_FAILURE)
		goto failure;
	status = substitution_reduce_to_fixed_point(
		presolve, &workspace, tolerance);
	if (status == lp_simplex_EXIT_FAILURE)
		goto failure;
	if (status != lp_simplex_Success) {
		presolve->status = status;
		goto success;
	}
	if (workspace.active_rows == 0 || workspace.active_columns == 0) {
		presolve->stats.removed_rows +=
			workspace.substitutions + workspace.closed_rows;
		presolve->stats.removed_columns += workspace.substitutions +
			workspace.projections + workspace.closed_columns;
		presolve->status = lp_simplex_Success;
		goto success;
	}
	if (workspace.substitutions == 0 && workspace.projections == 0 &&
	    workspace.closed_rows == 0 && workspace.closed_columns == 0)
		goto success;
	reduced = substitution_compact_problem(presolve, &workspace);
	if (reduced == NULL)
		goto failure;
	presolve->stats.removed_columns += workspace.substitutions +
		workspace.projections + workspace.closed_columns;
	presolve->stats.removed_rows +=
		workspace.substitutions + workspace.closed_rows;
	presolve->stats.doubleton_rows += workspace.doubletons;
	substitution_workspace_destroy(&workspace);
	simplex_problem_free(presolve->reduced);
	presolve->reduced = reduced;
	return lp_simplex_EXIT_SUCCESS;
success:
	substitution_workspace_destroy(&workspace);
	return lp_simplex_EXIT_SUCCESS;
failure:
	substitution_workspace_destroy(&workspace);
	return lp_simplex_EXIT_FAILURE;
}
