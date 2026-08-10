/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 *
 * Sparse conversion to equality form for Pan's generalized simplex method.
 */
#include "simplex_pan_standard.h"
#include "utils.h"

#include <lp_simplex/status.h>

#include <math.h>


static int pan_bound_has_width(const struct optm_VariableBound *bound)
{
	return bound->b_type == optm_BOUND_T_BS && bound->ub > bound->lb;
}


static int pan_transformed_column_count(
		const struct simplex_Problem *problem)
{
	int columns = 0;
	int j;
	for (j = 0; j < problem->columns; j++) {
		const struct optm_VariableBound *bound = problem->bounds + j;
		if (bound->b_type == optm_BOUND_T_BS && bound->ub == bound->lb)
			continue;
		columns++;
		if (bound->b_type == optm_BOUND_T_FR)
			columns++;
	}
	return columns;
}


static int pan_slack_column_count(const struct simplex_Problem *problem)
{
	int columns = 0;
	int i, j;
	for (i = 0; i < problem->rows; i++)
		if (problem->row_type[i] != optm_CONS_T_EQ)
			columns++;
	for (j = 0; j < problem->columns; j++)
		if (pan_bound_has_width(problem->bounds + j))
			columns++;
	return columns;
}


static int pan_extra_bound_rows(const struct simplex_Problem *problem)
{
	int rows = 0;
	int j;
	for (j = 0; j < problem->columns; j++)
		if (pan_bound_has_width(problem->bounds + j))
			rows++;
	return rows;
}


static void pan_variable_transform(
		const struct optm_VariableBound *bound,
		double *shift, double *scale)
{
	*shift = 0.;
	*scale = 1.;
	if (bound->b_type == optm_BOUND_T_UP) {
		*shift = bound->ub;
		*scale = -1.;
	} else if (bound->b_type == optm_BOUND_T_LO ||
		   bound->b_type == optm_BOUND_T_BS) {
		*shift = bound->lb;
	}
}


static int pan_standard_nonzeros(
		const struct simplex_Problem *problem)
{
	int nonzeros = 0;
	int j;
	for (j = 0; j < problem->columns; j++) {
		const struct optm_VariableBound *bound = problem->bounds + j;
		int count = problem->matrix.column_start[j + 1] -
			problem->matrix.column_start[j];
		if (bound->b_type == optm_BOUND_T_BS && bound->ub == bound->lb)
			continue;
		nonzeros += count;
		if (pan_bound_has_width(bound))
			nonzeros++;
		if (bound->b_type == optm_BOUND_T_FR)
			nonzeros += count;
	}
	nonzeros += pan_slack_column_count(problem);
	return nonzeros;
}


static int pan_standard_allocate(
		struct simplex_PanStandard *standard,
		const int rows, const int columns, const int original_columns,
		const int nonzeros)
{
	standard->rows = rows;
	standard->columns = columns;
	standard->original_columns = original_columns;
	standard->nonzeros = nonzeros;
	standard->column_start = (int *)lp_simplex_malloc(
		(size_t)(columns + 1) * sizeof(int));
	standard->row_index = nonzeros > 0 ? (int *)lp_simplex_malloc(
		(size_t)nonzeros * sizeof(int)) : NULL;
	standard->value = nonzeros > 0 ? (double *)lp_simplex_malloc(
		(size_t)nonzeros * sizeof(double)) : NULL;
	standard->normalized_value = nonzeros > 0 ? (double *)lp_simplex_malloc(
		(size_t)nonzeros * sizeof(double)) : NULL;
	standard->column_norm = (double *)lp_simplex_malloc(
		(size_t)(columns > 0 ? columns : 1) * sizeof(double));
	standard->rhs = (double *)lp_simplex_malloc((size_t)rows * sizeof(double));
	standard->objective = (double *)lp_simplex_malloc(
		(size_t)(columns > 0 ? columns : 1) * sizeof(double));
	standard->original_column = (int *)lp_simplex_malloc(
		(size_t)(columns > 0 ? columns : 1) * sizeof(int));
	standard->original_scale = (double *)lp_simplex_malloc(
		(size_t)(columns > 0 ? columns : 1) * sizeof(double));
	standard->original_shift = (double *)lp_simplex_malloc(
		(size_t)original_columns * sizeof(double));
	if (standard->column_start == NULL || standard->rhs == NULL ||
	    standard->column_norm == NULL || standard->objective == NULL ||
	    standard->original_column == NULL ||
	    standard->original_scale == NULL || standard->original_shift == NULL ||
	    (nonzeros > 0 && (standard->row_index == NULL ||
	     standard->value == NULL || standard->normalized_value == NULL)))
		return lp_simplex_EXIT_FAILURE;
	return lp_simplex_EXIT_SUCCESS;
}


static void pan_append_structural_column(
		struct simplex_PanStandard *standard,
		const struct simplex_Problem *problem,
		const int original, const double scale, const int upper_row,
		int *column, int *position)
{
	int k;
	standard->column_start[*column] = *position;
	standard->objective[*column] = scale * problem->objective[original];
	standard->original_column[*column] = original;
	standard->original_scale[*column] = scale;
	for (k = problem->matrix.column_start[original];
	     k < problem->matrix.column_start[original + 1]; k++) {
		standard->row_index[*position] = problem->matrix.row_index[k];
		standard->value[*position] = scale * problem->matrix.value[k];
		(*position)++;
	}
	if (upper_row >= 0) {
		standard->row_index[*position] = upper_row;
		standard->value[*position] = 1.;
		(*position)++;
	}
	(*column)++;
}


static void pan_append_slack_column(
		struct simplex_PanStandard *standard,
		const int row, const double coefficient,
		int *column, int *position)
{
	standard->column_start[*column] = *position;
	standard->row_index[*position] = row;
	standard->value[*position] = coefficient;
	standard->objective[*column] = 0.;
	standard->original_column[*column] = -1;
	standard->original_scale[*column] = 0.;
	(*position)++;
	(*column)++;
}


int simplex_pan_standard_create(
		struct simplex_PanStandard *standard,
		const struct simplex_Problem *problem)
{
	int bound_row, column = 0, position = 0;
	int transformed, slacks, rows, nonzeros;
	int i, j, k;
	lp_simplex_memset(standard, 0, sizeof(*standard));
	transformed = pan_transformed_column_count(problem);
	slacks = pan_slack_column_count(problem);
	rows = problem->rows + pan_extra_bound_rows(problem);
	nonzeros = pan_standard_nonzeros(problem);
	if (pan_standard_allocate(standard, rows, transformed + slacks,
		problem->columns, nonzeros) == lp_simplex_EXIT_FAILURE) {
		simplex_pan_standard_destroy(standard);
		return lp_simplex_EXIT_FAILURE;
	}
	lp_simplex_memcpy(standard->rhs, problem->rhs,
		(size_t)problem->rows * sizeof(double));
	lp_simplex_memset(standard->original_shift, 0,
		(size_t)problem->columns * sizeof(double));
	standard->objective_offset = 0.;
	bound_row = problem->rows;
	for (j = 0; j < problem->columns; j++) {
		const struct optm_VariableBound *bound = problem->bounds + j;
		double shift, scale;
		int upper_row = -1;
		if (bound->b_type == optm_BOUND_T_BS && bound->ub == bound->lb) {
			shift = bound->lb;
			standard->original_shift[j] = shift;
			standard->objective_offset += problem->objective[j] * shift;
			for (k = problem->matrix.column_start[j];
			     k < problem->matrix.column_start[j + 1]; k++)
				standard->rhs[problem->matrix.row_index[k]] -=
					problem->matrix.value[k] * shift;
			continue;
		}
		pan_variable_transform(bound, &shift, &scale);
		standard->original_shift[j] = shift;
		standard->objective_offset += problem->objective[j] * shift;
		for (k = problem->matrix.column_start[j];
		     k < problem->matrix.column_start[j + 1]; k++)
			standard->rhs[problem->matrix.row_index[k]] -=
				problem->matrix.value[k] * shift;
		if (pan_bound_has_width(bound)) {
			upper_row = bound_row++;
			standard->rhs[upper_row] = bound->ub - bound->lb;
		}
		pan_append_structural_column(standard, problem, j, scale,
			upper_row, &column, &position);
		if (bound->b_type == optm_BOUND_T_FR)
			pan_append_structural_column(standard, problem, j, -1.,
				-1, &column, &position);
	}
	for (i = 0; i < problem->rows; i++) {
		if (problem->row_type[i] == optm_CONS_T_LE)
			pan_append_slack_column(standard, i, 1., &column, &position);
		else if (problem->row_type[i] == optm_CONS_T_GE)
			pan_append_slack_column(standard, i, -1., &column, &position);
	}
	for (i = problem->rows; i < rows; i++)
		pan_append_slack_column(standard, i, 1., &column, &position);
	standard->column_start[column] = position;
	standard->columns = column;
	standard->nonzeros = position;
	if (column != transformed + slacks || position != nonzeros) {
		simplex_pan_standard_destroy(standard);
		return lp_simplex_EXIT_FAILURE;
	}
	for (j = 0; j < standard->columns; j++) {
		double square = 0.;
		for (k = standard->column_start[j];
		     k < standard->column_start[j + 1]; k++)
			square += standard->value[k] * standard->value[k];
		standard->column_norm[j] = square > 0. ? sqrt(square) : 0.;
		for (k = standard->column_start[j];
		     k < standard->column_start[j + 1]; k++)
			standard->normalized_value[k] = square > 0.
				? standard->value[k] / standard->column_norm[j] : 0.;
	}
	return lp_simplex_EXIT_SUCCESS;
}


void simplex_pan_standard_destroy(struct simplex_PanStandard *standard)
{
	if (standard == NULL)
		return;
	lp_simplex_free(standard->column_start);
	lp_simplex_free(standard->row_index);
	lp_simplex_free(standard->value);
	lp_simplex_free(standard->normalized_value);
	lp_simplex_free(standard->column_norm);
	lp_simplex_free(standard->rhs);
	lp_simplex_free(standard->objective);
	lp_simplex_free(standard->original_column);
	lp_simplex_free(standard->original_scale);
	lp_simplex_free(standard->original_shift);
	lp_simplex_memset(standard, 0, sizeof(*standard));
}


void simplex_pan_standard_recover(
		const struct simplex_PanStandard *standard,
		const double *z, double *x)
{
	int j;
	lp_simplex_memcpy(x, standard->original_shift,
		(size_t)standard->original_columns * sizeof(double));
	for (j = 0; j < standard->columns; j++)
		if (standard->original_column[j] >= 0)
			x[standard->original_column[j]] +=
				standard->original_scale[j] * z[j];
}
