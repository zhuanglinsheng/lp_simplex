/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#include "simplex_presolve_substitution.h"
#include "simplex_presolve.h"
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


static int substitution_finite_lower(
		const struct optm_VariableBound *bound)
{
	return bound->b_type == optm_BOUND_T_LO ||
		bound->b_type == optm_BOUND_T_BS;
}


static int substitution_finite_upper(
		const struct optm_VariableBound *bound)
{
	return bound->b_type == optm_BOUND_T_UP ||
		bound->b_type == optm_BOUND_T_BS;
}


static void substitution_set_lower(
		struct optm_VariableBound *bound, const double value)
{
	bound->lb = value;
	bound->b_type = substitution_finite_upper(bound)
		? optm_BOUND_T_BS : optm_BOUND_T_LO;
}


static void substitution_set_upper(
		struct optm_VariableBound *bound, const double value)
{
	bound->ub = value;
	bound->b_type = substitution_finite_lower(bound)
		? optm_BOUND_T_BS : optm_BOUND_T_UP;
}


static int substitution_append_incidence(
		struct presolve_Incidence *incidence, const int row)
{
	int capacity;
	int *grown;
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


static int substitution_reserve_terms(
		struct simplex_Presolve *presolve, const int additional)
{
	int capacity;
	int *columns;
	double *multipliers;
	if (presolve->substitution_term_count + additional <=
	    presolve->substitution_term_capacity)
		return lp_simplex_EXIT_SUCCESS;
	capacity = presolve->substitution_term_capacity == 0
		? 64 : presolve->substitution_term_capacity;
	while (capacity < presolve->substitution_term_count + additional)
		capacity *= 2;
	columns = (int *)lp_simplex_realloc(presolve->substitution_term_column,
		(size_t)capacity * sizeof(int));
	if (columns == NULL)
		return lp_simplex_EXIT_FAILURE;
	presolve->substitution_term_column = columns;
	multipliers = (double *)lp_simplex_realloc(
		presolve->substitution_term_multiplier,
		(size_t)capacity * sizeof(double));
	if (multipliers == NULL)
		return lp_simplex_EXIT_FAILURE;
	presolve->substitution_term_multiplier = multipliers;
	presolve->substitution_term_capacity = capacity;
	return lp_simplex_EXIT_SUCCESS;
}


static struct simplex_PresolveSubstitution *substitution_new_record(
		struct simplex_Presolve *presolve, const int column,
		const double constant, const int terms)
{
	struct simplex_PresolveSubstitution *grown;
	int capacity;
	if (presolve->substitution_count == presolve->substitution_capacity) {
		capacity = presolve->substitution_capacity == 0
			? 32 : presolve->substitution_capacity * 2;
		grown = (struct simplex_PresolveSubstitution *)lp_simplex_realloc(
			presolve->substitution,
			(size_t)capacity * sizeof(*grown));
		if (grown == NULL)
			return NULL;
		presolve->substitution = grown;
		presolve->substitution_capacity = capacity;
	}
	if (substitution_reserve_terms(presolve, terms) ==
	    lp_simplex_EXIT_FAILURE)
		return NULL;
	grown = presolve->substitution + presolve->substitution_count++;
	grown->column = column;
	grown->constant = constant;
	grown->term_start = presolve->substitution_term_count;
	grown->term_count = terms;
	return grown;
}


static void substitution_add_record_term(
		struct simplex_Presolve *presolve, const int column,
		const double multiplier)
{
	int position = presolve->substitution_term_count++;
	presolve->substitution_term_column[position] = column;
	presolve->substitution_term_multiplier[position] = multiplier;
}


static int substitution_transfer_bounds(
		const struct optm_VariableBound *eliminated,
		struct optm_VariableBound *kept, const double constant,
		const double multiplier, const double tolerance)
{
	double lower = __lp_simplex_NINF__, upper = __lp_simplex_INF__;
	double margin;
	if (substitution_finite_lower(eliminated)) {
		double value = (eliminated->lb - constant) / multiplier;
		if (multiplier > 0.) lower = value;
		else upper = value;
	}
	if (substitution_finite_upper(eliminated)) {
		double value = (eliminated->ub - constant) / multiplier;
		if (multiplier > 0.) upper = value;
		else lower = value;
	}
	if (isfinite(lower) &&
	    (!substitution_finite_lower(kept) || lower > kept->lb))
		substitution_set_lower(kept, lower);
	if (isfinite(upper) &&
	    (!substitution_finite_upper(kept) || upper < kept->ub))
		substitution_set_upper(kept, upper);
	margin = tolerance * (1. +
		(substitution_finite_lower(kept) ? __lp_simplex_ABS__(kept->lb) : 0.) +
		(substitution_finite_upper(kept) ? __lp_simplex_ABS__(kept->ub) : 0.));
	return substitution_finite_lower(kept) &&
		substitution_finite_upper(kept) && kept->lb > kept->ub + margin
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
	*lower_implied = !substitution_finite_lower(eliminated);
	*upper_implied = !substitution_finite_upper(eliminated);
	for (k = 0; k < row->count; k++) {
		const struct optm_VariableBound *bound;
		double multiplier;
		if (k == eliminated_position)
			continue;
		bound = problem->bounds + row->column[k];
		multiplier = -row->value[k] / eliminated_coefficient;
		if (multiplier > 0.) {
			if (substitution_finite_lower(bound)) {
				long double contribution = multiplier * bound->lb;
				minimum += contribution;
				magnitude += fabsl(contribution);
			} else minimum_finite = 0;
			if (substitution_finite_upper(bound)) {
				long double contribution = multiplier * bound->ub;
				maximum += contribution;
				magnitude += fabsl(contribution);
			} else maximum_finite = 0;
		} else {
			if (substitution_finite_upper(bound)) {
				long double contribution = multiplier * bound->ub;
				minimum += contribution;
				magnitude += fabsl(contribution);
			} else minimum_finite = 0;
			if (substitution_finite_lower(bound)) {
				long double contribution = multiplier * bound->lb;
				maximum += contribution;
				magnitude += fabsl(contribution);
			} else maximum_finite = 0;
		}
	}
	{
		long double margin = (long double)tolerance *
			(1. + fabsl((long double)constant) + magnitude);
		if (substitution_finite_lower(eliminated) && minimum_finite &&
		    minimum >= eliminated->lb - margin)
			*lower_implied = 1;
		if (substitution_finite_upper(eliminated) && maximum_finite &&
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
		const double tolerance,
		int *active_rows, int *active_columns, int *substitutions,
		int *projections)
{
	int j, progress;
	do {
		progress = 0;
		for (j = 0; j < problem->columns; j++) {
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
			if (substitution_new_record(presolve,
					presolve->column_map[j], constant,
					row->count - 1) == NULL)
				return lp_simplex_EXIT_FAILURE;
			for (t = 0; t < row->count; t++)
				if (t != position) {
					int kept = row->column[t];
					double multiplier = -row->value[t] / coefficient;
					substitution_add_record_term(presolve,
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
							if (inserted && substitution_append_incidence(
									incidence + kept, affected) ==
							    lp_simplex_EXIT_FAILURE)
								return lp_simplex_EXIT_FAILURE;
						}
				}
				for (t = 0; t < row->count; t++)
					if (column_active[row->column[t]] &&
					    column_degree[row->column[t]] > 0)
						column_degree[row->column[t]]--;
				row->active = 0;
				(*active_rows)--;
				(*substitutions)++;
				if (degree == 1)
					presolve->singleton_column_rows++;
				else
					presolve->implied_free_columns++;
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
				presolve->singleton_projection_columns++;
			}
			column_active[j] = 0;
			column_degree[j] = 0;
			presolve->eliminated[presolve->column_map[j]] = 1;
			(*active_columns)--;
			progress = 1;
		}
	} while (progress && *active_rows > 1 && *active_columns > 1);
	return lp_simplex_EXIT_SUCCESS;
}


static int substitution_empty_row_feasible(
		const int type, const double rhs, const double tolerance)
{
	if (type == optm_CONS_T_EQ)
		return __lp_simplex_ABS__(rhs) <= tolerance;
	if (type == optm_CONS_T_GE)
		return rhs <= tolerance;
	return rhs >= -tolerance;
}


static int substitution_choose_empty_column(
		const struct simplex_Problem *problem, const int column,
		double *value)
{
	const struct optm_VariableBound *bound = problem->bounds + column;
	double cost = problem->objective[column];
	if (cost > 0.) {
		if (!substitution_finite_lower(bound))
			return lp_simplex_Unboundedness;
		*value = bound->lb;
	} else if (cost < 0.) {
		if (!substitution_finite_upper(bound))
			return lp_simplex_Unboundedness;
		*value = bound->ub;
	} else if (substitution_finite_lower(bound) && bound->lb > 0.) {
		*value = bound->lb;
	} else if (substitution_finite_upper(bound) && bound->ub < 0.) {
		*value = bound->ub;
	} else *value = 0.;
	return lp_simplex_Success;
}


static void substitution_eliminate_fixed_column(
		struct simplex_Presolve *presolve,
		struct simplex_Problem *problem,
		struct presolve_MutableRow *rows,
		struct presolve_Incidence *incidence,
		unsigned char *column_active, int *column_degree,
		int *queue, int *queue_count, const int column,
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
		if (row->count <= 1 && *queue_count < problem->rows)
			queue[(*queue_count)++] = row_index;
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
		int *queue, const double tolerance,
		int *active_rows, int *active_columns,
		int *removed_rows, int *removed_columns)
{
	int i, j, queue_head, queue_count, activity_changed;
restart:
	queue_head = 0;
	queue_count = 0;
	for (i = 0; i < problem->rows; i++)
		if (rows[i].active && rows[i].count <= 1)
			queue[queue_count++] = i;
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
			status = substitution_choose_empty_column(problem, j, &value);
			if (status != lp_simplex_Success)
				return status;
			presolve->empty_columns++;
			fixed = 1;
		}
		if (fixed) {
			presolve->fixed_columns += column_degree[j] != 0;
			substitution_eliminate_fixed_column(presolve, problem, rows,
				incidence, column_active, column_degree, queue,
				&queue_count, j, value, active_columns, removed_columns);
		}
	}
	while (queue_head < queue_count) {
		struct presolve_MutableRow *row;
		double coefficient, implied, margin;
		int row_index = queue[queue_head++];
		struct optm_VariableBound *bound;
		int lower, status = lp_simplex_Success;
		if (!rows[row_index].active || rows[row_index].count > 1)
			continue;
		row = rows + row_index;
		if (row->count == 0) {
			if (!substitution_empty_row_feasible(
					problem->row_type[row_index],
					problem->rhs[row_index], tolerance))
				return lp_simplex_Infeasibility;
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
		implied = problem->rhs[row_index] / coefficient;
		margin = tolerance * (1. + __lp_simplex_ABS__(implied));
		if (problem->row_type[row_index] == optm_CONS_T_EQ) {
			if ((substitution_finite_lower(bound) &&
			     implied < bound->lb - margin) ||
			    (substitution_finite_upper(bound) &&
			     implied > bound->ub + margin))
				return lp_simplex_Infeasibility;
			substitution_set_lower(bound, implied);
			substitution_set_upper(bound, implied);
			presolve->singleton_columns++;
		} else {
			lower = (problem->row_type[row_index] == optm_CONS_T_GE &&
				 coefficient > 0.) ||
				(problem->row_type[row_index] == optm_CONS_T_LE &&
				 coefficient < 0.);
			if (lower) {
				if (substitution_finite_upper(bound) &&
				    implied > bound->ub + margin)
					status = lp_simplex_Infeasibility;
				else if (!substitution_finite_lower(bound) ||
					 implied > bound->lb)
					substitution_set_lower(bound, implied);
			} else {
				if (substitution_finite_lower(bound) &&
				    implied < bound->lb - margin)
					status = lp_simplex_Infeasibility;
				else if (!substitution_finite_upper(bound) ||
					 implied < bound->ub)
					substitution_set_upper(bound, implied);
			}
			if (status != lp_simplex_Success)
				return status;
		}
		row->active = 0;
		if (column_degree[j] > 0)
			column_degree[j]--;
		(*active_rows)--;
		(*removed_rows)++;
		presolve->singleton_rows++;
		presolve->tightened_bounds++;
		if (bound->b_type == optm_BOUND_T_BS && bound->lb == bound->ub) {
			presolve->fixed_columns++;
			substitution_eliminate_fixed_column(presolve, problem, rows,
				incidence, column_active, column_degree, queue,
				&queue_count, j, bound->lb, active_columns,
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
				if (substitution_finite_lower(bound)) {
					contribution = coefficient * bound->lb;
					minimum += contribution;
					magnitude += fabsl((long double)contribution);
				} else minimum_infinite++;
				if (substitution_finite_upper(bound)) {
					contribution = coefficient * bound->ub;
					maximum += contribution;
					magnitude += fabsl((long double)contribution);
				} else maximum_infinite++;
			} else {
				if (substitution_finite_upper(bound)) {
					contribution = coefficient * bound->ub;
					minimum += contribution;
					magnitude += fabsl((long double)contribution);
				} else minimum_infinite++;
				if (substitution_finite_lower(bound)) {
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
		    maximum < (long double)problem->rhs[i] - margin)
			return lp_simplex_Infeasibility;
		if ((problem->row_type[i] == optm_CONS_T_LE ||
		     problem->row_type[i] == optm_CONS_T_EQ) &&
		    minimum_infinite == 0 &&
		    minimum > (long double)problem->rhs[i] + margin)
			return lp_simplex_Infeasibility;
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
			presolve->redundant_rows++;
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
			substitution_set_lower(bound, value);
			substitution_set_upper(bound, value);
			fixed++;
		}
		if (fixed != 0) {
			presolve->forcing_rows++;
			presolve->forced_columns += fixed;
			presolve->tightened_bounds += fixed;
			activity_changed = 1;
		}
	}
	if (activity_changed)
		goto restart;
	for (j = 0; j < problem->columns; j++)
		if (column_active[j] && column_degree[j] == 0) {
			double value;
			int status = substitution_choose_empty_column(problem, j, &value);
			if (status != lp_simplex_Success)
				return status;
			presolve->empty_columns++;
			substitution_eliminate_fixed_column(presolve, problem, rows,
				incidence, column_active, column_degree, queue,
				&queue_count, j, value, active_columns, removed_columns);
		}
	return lp_simplex_Success;
}


int simplex_presolve_substitute_doubletons(
		struct simplex_Presolve *presolve, const double tolerance)
{
	struct simplex_Problem *problem = presolve->reduced;
	struct presolve_MutableRow *rows = NULL;
	struct presolve_Incidence *incidence = NULL;
	unsigned char *column_active = NULL;
	int *column_degree = NULL;
	int *queue = NULL, *new_column_map = NULL, *new_row_map = NULL;
	int queue_head = 0, queue_count = 0;
	int i, j, k, active_rows, active_columns, substitutions = 0;
	int doubletons = 0, projections = 0, closed_rows = 0, closed_columns = 0;
	int closure_status;
	struct simplex_Problem *reduced = NULL;
	if (problem == NULL)
		return lp_simplex_EXIT_SUCCESS;
	rows = (struct presolve_MutableRow *)lp_simplex_malloc(
		(size_t)problem->rows * sizeof(*rows));
	incidence = (struct presolve_Incidence *)lp_simplex_malloc(
		(size_t)problem->columns * sizeof(*incidence));
	column_active = (unsigned char *)lp_simplex_malloc(
		(size_t)problem->columns * sizeof(unsigned char));
	column_degree = (int *)lp_simplex_malloc(
		(size_t)problem->columns * sizeof(int));
	queue = (int *)lp_simplex_malloc((size_t)problem->rows * sizeof(int));
	if (rows == NULL || incidence == NULL || column_active == NULL ||
	    column_degree == NULL || queue == NULL)
		goto failure;
	lp_simplex_memset(rows, 0, (size_t)problem->rows * sizeof(*rows));
	lp_simplex_memset(incidence, 0,
		(size_t)problem->columns * sizeof(*incidence));
	lp_simplex_memset(column_active, 1,
		(size_t)problem->columns * sizeof(unsigned char));
	lp_simplex_memset(column_degree, 0,
		(size_t)problem->columns * sizeof(int));
	for (i = 0; i < problem->rows; i++) {
		int count = problem->matrix.row_start[i + 1] -
			problem->matrix.row_start[i];
		rows[i].active = 1;
		rows[i].count = rows[i].capacity = count;
		rows[i].column = (int *)lp_simplex_malloc(
			(size_t)count * sizeof(int));
		rows[i].value = (double *)lp_simplex_malloc(
			(size_t)count * sizeof(double));
		if (count != 0 && (rows[i].column == NULL || rows[i].value == NULL))
			goto failure;
		for (k = 0; k < count; k++) {
			j = problem->matrix.column_index[
				problem->matrix.row_start[i] + k];
			rows[i].column[k] = j;
			rows[i].value[k] = problem->matrix.row_value[
				problem->matrix.row_start[i] + k];
			column_degree[j]++;
			if (substitution_append_incidence(incidence + j, i) ==
			    lp_simplex_EXIT_FAILURE)
				goto failure;
		}
		if (problem->row_type[i] == optm_CONS_T_EQ && count == 2)
			queue[queue_count++] = i;
	}
	active_rows = problem->rows;
	active_columns = problem->columns;
	/* A singleton column in an equality can be substituted without matrix
	 * fill.  Its bounds must already be implied by the remaining row activity;
	 * otherwise eliminating it would create a new ranged constraint. */
	if (substitution_remove_singleton_columns(presolve, problem, rows,
			incidence, column_active, column_degree, tolerance, &active_rows,
			&active_columns, &substitutions, &projections) ==
		    lp_simplex_EXIT_FAILURE)
		goto failure;
	/* Low-degree column aggregation can expose new doubleton equalities.
	 * Rebuild the work queue from the mutated rows instead of keeping only
	 * the doubletons present in the input model. */
	queue_head = 0;
	queue_count = 0;
	for (i = 0; i < problem->rows; i++)
		if (rows[i].active && rows[i].count == 2 &&
		    problem->row_type[i] == optm_CONS_T_EQ)
			queue[queue_count++] = i;
	while (queue_head < queue_count && active_rows > 1 && active_columns > 1) {
		int row_index = queue[queue_head++];
		struct presolve_MutableRow *row = rows + row_index;
		int first, second, eliminate, keep, eliminate_position;
		int first_degree, second_degree;
		double eliminate_coefficient, keep_coefficient;
		double constant, multiplier;
		if (!row->active || row->count != 2 ||
		    problem->row_type[row_index] != optm_CONS_T_EQ)
			continue;
		first = row->column[0];
		second = row->column[1];
		if (!column_active[first] || !column_active[second])
			continue;
		first_degree = column_degree[first];
		second_degree = column_degree[second];
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
				problem->bounds + keep, constant, multiplier,
				tolerance) != lp_simplex_Success)
			goto infeasible;
		if (substitution_new_record(presolve,
				presolve->column_map[eliminate], constant, 1) == NULL)
			goto failure;
		substitution_add_record_term(presolve,
			presolve->column_map[keep], multiplier);
		problem->objective[keep] += problem->objective[eliminate] * multiplier;
		for (k = 0; k < incidence[eliminate].count; k++) {
			int affected = incidence[eliminate].row[k];
			struct presolve_MutableRow *target = rows + affected;
			int position, inserted, removed;
			double coefficient;
			if (!target->active || affected == row_index)
				continue;
			position = substitution_find_column(target, eliminate);
			if (position < 0)
				continue;
			coefficient = target->value[position];
			substitution_remove_from_row(target, position);
			column_degree[eliminate]--;
			problem->rhs[affected] -= coefficient * constant;
			if (substitution_add_to_row(target, keep,
					coefficient * multiplier, &inserted, &removed) ==
				    lp_simplex_EXIT_FAILURE)
				goto failure;
			column_degree[keep] += inserted - removed;
			if (inserted && substitution_append_incidence(
					incidence + keep, affected) == lp_simplex_EXIT_FAILURE)
				goto failure;
			if (problem->row_type[affected] == optm_CONS_T_EQ &&
			    target->count == 2 && queue_count < problem->rows)
				queue[queue_count++] = affected;
		}
		if (column_degree[first] > 0)
			column_degree[first]--;
		if (column_degree[second] > 0)
			column_degree[second]--;
		row->active = 0;
		column_active[eliminate] = 0;
		column_degree[eliminate] = 0;
		presolve->eliminated[presolve->column_map[eliminate]] = 1;
		active_rows--;
		active_columns--;
		substitutions++;
		doubletons++;
	}
	/* Removing doubleton rows lowers adjacent column degrees.  Run the
	 * singleton-column fixed point again to capture the resulting cascade. */
	if (substitution_remove_singleton_columns(presolve, problem, rows,
			incidence, column_active, column_degree, tolerance, &active_rows,
			&active_columns, &substitutions, &projections) ==
		    lp_simplex_EXIT_FAILURE)
		goto failure;
	closure_status = substitution_close_fixed_point(presolve, problem, rows,
		incidence, column_active, column_degree, queue, tolerance,
		&active_rows, &active_columns, &closed_rows, &closed_columns);
	if (closure_status == lp_simplex_EXIT_FAILURE)
		goto failure;
	if (closure_status != lp_simplex_Success) {
		presolve->terminal = 1;
		presolve->terminal_status = closure_status;
		substitution_destroy_rows(rows, problem->rows,
			incidence, problem->columns);
		lp_simplex_free(column_active);
		lp_simplex_free(column_degree);
		lp_simplex_free(queue);
		return lp_simplex_EXIT_SUCCESS;
	}
	if (active_rows == 0 || active_columns == 0) {
		presolve->removed_rows += substitutions + closed_rows;
		presolve->removed_columns += substitutions + projections + closed_columns;
		presolve->terminal = 1;
		presolve->terminal_status = lp_simplex_Success;
		substitution_destroy_rows(rows, problem->rows,
			incidence, problem->columns);
		lp_simplex_free(column_active);
		lp_simplex_free(column_degree);
		lp_simplex_free(queue);
		return lp_simplex_EXIT_SUCCESS;
	}
	if (substitutions == 0 && projections == 0 &&
	    closed_rows == 0 && closed_columns == 0) {
		substitution_destroy_rows(rows, problem->rows,
			incidence, problem->columns);
		lp_simplex_free(column_active);
		lp_simplex_free(column_degree);
		lp_simplex_free(queue);
		return lp_simplex_EXIT_SUCCESS;
	}
	{
		int nonzeros = 0, next = 0, new_row = 0, new_column = 0;
		int *column_position = (int *)lp_simplex_malloc(
			(size_t)problem->columns * sizeof(int));
		new_column_map = (int *)lp_simplex_malloc(
			(size_t)active_columns * sizeof(int));
		new_row_map = (int *)lp_simplex_malloc(
			(size_t)active_rows * sizeof(int));
		if (column_position == NULL || new_column_map == NULL ||
		    new_row_map == NULL) {
			lp_simplex_free(column_position);
			goto failure;
		}
		for (j = 0; j < problem->columns; j++) {
			column_position[j] = -1;
			if (column_active[j]) {
				column_position[j] = new_column++;
				new_column_map[new_column - 1] = presolve->column_map[j];
			}
		}
		for (i = 0; i < problem->rows; i++)
			if (rows[i].active) {
				new_row_map[new_row++] = presolve->row_map[i];
				nonzeros += rows[i].count;
			}
		reduced = simplex_problem_create_sparse(
			active_rows, active_columns, nonzeros);
		if (reduced == NULL) {
			lp_simplex_free(column_position);
			goto failure;
		}
		for (j = 0; j < problem->columns; j++)
			if (column_active[j]) {
				int destination = column_position[j];
				reduced->objective[destination] = problem->objective[j];
				reduced->bounds[destination] = problem->bounds[j];
			}
		new_row = 0;
		for (i = 0; i < problem->rows; i++)
			if (rows[i].active) {
				reduced->rhs[new_row] = problem->rhs[i];
				reduced->row_type[new_row] = problem->row_type[i];
				reduced->matrix.row_start[new_row] = next;
				for (k = 0; k < rows[i].count; k++) {
					reduced->matrix.column_index[next] =
						column_position[rows[i].column[k]];
					reduced->matrix.row_value[next++] = rows[i].value[k];
				}
				new_row++;
			}
		reduced->matrix.row_start[active_rows] = next;
		lp_simplex_memset(reduced->matrix.column_start, 0,
			(size_t)(active_columns + 1) * sizeof(int));
		for (k = 0; k < nonzeros; k++)
			reduced->matrix.column_start[
				reduced->matrix.column_index[k] + 1]++;
		for (j = 0; j < active_columns; j++)
			reduced->matrix.column_start[j + 1] +=
				reduced->matrix.column_start[j];
		{
			int *cursor = (int *)lp_simplex_malloc(
				(size_t)active_columns * sizeof(int));
			if (cursor == NULL) {
				lp_simplex_free(column_position);
				goto failure;
			}
			lp_simplex_memcpy(cursor, reduced->matrix.column_start,
				(size_t)active_columns * sizeof(int));
			for (i = 0; i < active_rows; i++)
				for (k = reduced->matrix.row_start[i];
				     k < reduced->matrix.row_start[i + 1]; k++) {
					int column = reduced->matrix.column_index[k];
					int destination = cursor[column]++;
					reduced->matrix.row_index[destination] = i;
					reduced->matrix.value[destination] =
						reduced->matrix.row_value[k];
				}
			lp_simplex_free(cursor);
		}
		lp_simplex_free(column_position);
	}
	for (j = 0; j < active_columns; j++)
		presolve->column_map[j] = new_column_map[j];
	for (i = 0; i < active_rows; i++)
		presolve->row_map[i] = new_row_map[i];
	presolve->removed_columns += substitutions + projections + closed_columns;
	presolve->removed_rows += substitutions + closed_rows;
	presolve->doubleton_rows += doubletons;
	substitution_destroy_rows(rows, problem->rows, incidence, problem->columns);
	simplex_problem_free(presolve->reduced);
	presolve->reduced = reduced;
	lp_simplex_free(column_active);
	lp_simplex_free(column_degree);
	lp_simplex_free(queue);
	lp_simplex_free(new_column_map);
	lp_simplex_free(new_row_map);
	return lp_simplex_EXIT_SUCCESS;
infeasible:
	presolve->terminal = 1;
	presolve->terminal_status = lp_simplex_Infeasibility;
	substitution_destroy_rows(rows, problem->rows, incidence, problem->columns);
	lp_simplex_free(column_active);
	lp_simplex_free(column_degree);
	lp_simplex_free(queue);
	return lp_simplex_EXIT_SUCCESS;
failure:
	simplex_problem_free(reduced);
	substitution_destroy_rows(rows,
		problem != NULL ? problem->rows : 0, incidence,
		problem != NULL ? problem->columns : 0);
	lp_simplex_free(column_active);
	lp_simplex_free(column_degree);
	lp_simplex_free(queue);
	lp_simplex_free(new_column_map);
	lp_simplex_free(new_row_map);
	return lp_simplex_EXIT_FAILURE;
}
