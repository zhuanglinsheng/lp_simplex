/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_pan_standard.c
 * @brief Sparse conversion to equality form for Pan's generalized simplex
 * method.
 */

#include "policy/pan.h"
#include "policy/controls.h"
#include "simplex_pan_standard.h"
#include "simplex_pan_storage.h"
#include "utils.h"
#include "linalg.h"

#include <lp_simplex/status.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int pan_bound_has_width(const struct optm_VariableBound *bound)
{
	return bound->b_type == optm_BOUND_T_BS && bound->ub > bound->lb;
}

/* Count the actual materialized representation before allocating or copying
 * matrix entries. In particular, free columns duplicate their support. */
static int pan_standard_dimensions(
	const struct simplex_Problem *p, int *rows, int *columns, int *nonzeros)
{
	int i, j;
	if (!p || p->rows < 0 || p->columns < 0 || p->matrix.nonzeros < 0 ||
		!p->matrix.column_start || (p->columns && !p->bounds) ||
		(p->rows && !p->row_type))
		return 0;
	*rows = p->rows;
	*columns = *nonzeros = 0;
	for (i = 0; i < p->rows; i++)
		if (p->row_type[i] != optm_CONS_T_EQ &&
			(!pan_storage_add(columns, 1) ||
				!pan_storage_add(nonzeros, 1)))
			return 0;
	if (p->matrix.column_start[0] != 0)
		return 0;
	for (j = 0; j < p->columns; j++) {
		const struct optm_VariableBound *bound = p->bounds + j;
		int start = p->matrix.column_start[j];
		int end = p->matrix.column_start[j + 1];
		int count;
		if (start < 0 || end < start || end > p->matrix.nonzeros)
			return 0;
		count = end - start;
		if (bound->b_type == optm_BOUND_T_BS && bound->ub == bound->lb)
			continue;
		if (!pan_storage_add(columns, 1) ||
			!pan_storage_add(nonzeros, count))
			return 0;
		if (bound->b_type == optm_BOUND_T_FR &&
			(!pan_storage_add(columns, 1) ||
				!pan_storage_add(nonzeros, count)))
			return 0;
		if (pan_bound_has_width(bound) &&
			(!pan_storage_add(rows, 1) ||
				!pan_storage_add(columns, 1) ||
				!pan_storage_add(nonzeros, 2)))
			return 0;
	}
	return p->matrix.column_start[p->columns] == p->matrix.nonzeros;
}

static void pan_variable_transform(
	const struct optm_VariableBound *bound, double *shift, double *scale)
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

static int pan_standard_allocate(struct simplex_PanStandard *standard,
	const int rows, const int columns, const int original_columns,
	const int nonzeros)
{
	if (!pan_storage_fits(
		    rows, columns, original_columns, nonzeros, 0, (size_t)-1))
		return lp_simplex_EXIT_FAILURE;
	standard->rows = rows;

	standard->columns = columns;
	standard->original_columns = original_columns;
	standard->nonzeros = nonzeros;
	standard->column_start =
		(int *)lp_simplex_malloc(((size_t)columns + 1) * sizeof(int));
	standard->row_index = nonzeros > 0
				      ? (int *)lp_simplex_malloc(
						(size_t)nonzeros * sizeof(int))
				      : NULL;
	standard->value = nonzeros > 0
				  ? (double *)lp_simplex_malloc(
					    (size_t)nonzeros * sizeof(double))
				  : NULL;
	standard->normalized_value =
		nonzeros > 0 ? (double *)lp_simplex_malloc(
				       (size_t)nonzeros * sizeof(double))
			     : NULL;
	standard->column_norm = (double *)lp_simplex_malloc(
		(size_t)(columns > 0 ? columns : 1) * sizeof(double));
	standard->rhs =
		(double *)lp_simplex_malloc((size_t)rows * sizeof(double));
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
		standard->original_scale == NULL ||
		standard->original_shift == NULL ||
		(nonzeros > 0 && (standard->row_index == NULL ||
					 standard->value == NULL ||
					 standard->normalized_value == NULL)))
		return lp_simplex_EXIT_FAILURE;
	return lp_simplex_EXIT_SUCCESS;
}

static void pan_append_structural_column(struct simplex_PanStandard *standard,
	const struct simplex_Problem *problem, const int original,
	const double scale, const int upper_row, int *column, int *position)
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

static void pan_append_slack_column(struct simplex_PanStandard *standard,
	const int row, const double coefficient, int *column, int *position)
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

int simplex_pan_standard_create(struct simplex_PanStandard *standard,
	const struct simplex_Problem *problem)
{
	int bound_row, column = 0, position = 0;
	int columns, rows, nonzeros;
	int i, j, k;

	lp_simplex_memset(standard, 0, sizeof(*standard));
	standard->report = problem->report;
	if (!pan_standard_dimensions(problem, &rows, &columns, &nonzeros) ||
		pan_standard_allocate(standard, rows, columns, problem->columns,
			nonzeros) == lp_simplex_EXIT_FAILURE) {
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

		if (bound->b_type == optm_BOUND_T_BS &&
			bound->ub == bound->lb) {
			shift = bound->lb;
			standard->original_shift[j] = shift;
			standard->objective_offset +=
				problem->objective[j] * shift;
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
			pan_append_slack_column(
				standard, i, 1., &column, &position);
		else if (problem->row_type[i] == optm_CONS_T_GE)
			pan_append_slack_column(
				standard, i, -1., &column, &position);
	}
	for (i = problem->rows; i < rows; i++)
		pan_append_slack_column(standard, i, 1., &column, &position);
	standard->column_start[column] = position;
	standard->columns = column;
	standard->nonzeros = position;
	if (column != columns || position != nonzeros) {
		simplex_pan_standard_destroy(standard);
		return lp_simplex_EXIT_FAILURE;
	}
	for (j = 0; j < standard->columns; j++) {
		double square = 0.;

		for (k = standard->column_start[j];
			k < standard->column_start[j + 1]; k++)
			square = hypot(square, standard->value[k]);
		standard->column_norm[j] = square;
		for (k = standard->column_start[j];
			k < standard->column_start[j + 1]; k++)
			standard->normalized_value[k] =
				square > 0. ? standard->value[k] /
						      standard->column_norm[j]
					    : 0.;
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
	const struct simplex_PanStandard *standard, const double *z, double *x)
{
	int j;

	lp_simplex_memcpy(x, standard->original_shift,
		(size_t)standard->original_columns * sizeof(double));
	for (j = 0; j < standard->columns; j++)
		if (standard->original_column[j] >= 0)
			x[standard->original_column[j]] +=
				standard->original_scale[j] * z[j];
}

int simplex_pan_standard_encode(const struct simplex_PanStandard *standard,
	const double *x, double *z, double *residual, const double tolerance)
{
	int i, j, k;
	double *recovered;

	lp_simplex_memcpy(residual, standard->rhs,
		(size_t)standard->rows * sizeof(double));
	for (j = 0; j < standard->columns; j++) {
		int original = standard->original_column[j];

		z[j] = 0.;
		if (original < 0)
			continue;
		if (!isfinite(x[original]))
			return lp_simplex_EXIT_FAILURE;
		z[j] = (x[original] - standard->original_shift[original]) *
		       standard->original_scale[j];
		/* For a free variable, exactly one of its opposite columns is used.
		 * Other negative values are caught by the residual/nonnegative audit. */
		if (z[j] < 0.)
			z[j] = 0.;
		if (z[j] <= tolerance * SIMPLEX_PAN_INNER_TOLERANCE_FRACTION)
			z[j] = 0.;
		for (k = standard->column_start[j];
			k < standard->column_start[j + 1]; k++)
			residual[standard->row_index[k]] -=
				standard->value[k] * z[j];
	}
	for (j = 0; j < standard->columns; j++) {
		if (standard->original_column[j] >= 0)
			continue;
		k = standard->column_start[j];
		if (standard->column_start[j + 1] - k != 1 ||
			standard->value[k] == 0.)
			return lp_simplex_EXIT_FAILURE;
		i = standard->row_index[k];
		z[j] = residual[i] / standard->value[k];
		if (!isfinite(z[j]) || z[j] < -tolerance)
			return lp_simplex_EXIT_FAILURE;
		if (z[j] < 0.)
			z[j] = 0.;
		if (z[j] <= tolerance * SIMPLEX_PAN_INNER_TOLERANCE_FRACTION)
			z[j] = 0.;
		residual[i] -= standard->value[k] * z[j];
	}
	for (i = 0; i < standard->rows; i++)
		if (!isfinite(residual[i]) || fabs(residual[i]) > tolerance)
			return lp_simplex_EXIT_FAILURE;
	recovered = (double *)lp_simplex_malloc(
		(size_t)(standard->original_columns > 0
				 ? standard->original_columns
				 : 1) *
		sizeof(double));
	if (recovered == NULL)
		return lp_simplex_EXIT_FAILURE;
	simplex_pan_standard_recover(standard, z, recovered);
	for (j = 0; j < standard->original_columns; j++)
		if (!isfinite(x[j]) || fabs(recovered[j] - x[j]) > tolerance) {
			lp_simplex_free(recovered);
			return lp_simplex_EXIT_FAILURE;
		}
	lp_simplex_free(recovered);
	return lp_simplex_EXIT_SUCCESS;
}

double simplex_pan_column_dot(const struct simplex_PanStandard *standard,
	const int column, const double *vector)
{
	return lp_simplex_linalg_sparse_dot(standard->column_start[column],
		standard->column_start[column + 1], standard->row_index,
		standard->value, vector);
}

double simplex_pan_residual(const struct simplex_PanStandard *standard,
	const double *value, double *residual, double *correction)
{
	double maximum = 0.;
	int i, j, k;

	lp_simplex_memcpy(residual, standard->rhs,
		(size_t)standard->rows * sizeof(double));
	lp_simplex_memset(
		correction, 0, (size_t)standard->rows * sizeof(double));
	for (j = 0; j < standard->columns; j++) {
		if (value[j] == 0.)
			continue;
		for (k = standard->column_start[j];
			k < standard->column_start[j + 1]; k++) {
			double a = -standard->value[k], product = a * value[j],
			       next;

			i = standard->row_index[k];
			next = residual[i] + product;
			correction[i] +=
				(fabs(residual[i]) >= fabs(product)
						? (residual[i] - next) + product
						: (product - next) +
							  residual[i]) +
				fma(a, value[j], -product);
			residual[i] = next;
		}
	}
	for (i = 0; i < standard->rows; i++) {
		residual[i] += correction[i];
		if (!isfinite(residual[i]))
			return HUGE_VAL;
		maximum = fmax(maximum, fabs(residual[i]));
	}
	return maximum;
}

/* Presolve can remove structural rows yet add many implied finite bounds.
 * Pan currently expands every width into an equality and a slack. Compare
 * the actual equality-form dimensions before choosing its input model. */
int simplex_pan_original_form_smaller(
	const struct lp_Model *model, const struct simplex_Problem *reduced)
{
	size_t original_rows = (size_t)model->m,
	       original_columns = (size_t)model->n;
	size_t reduced_rows = (size_t)reduced->rows,
	       reduced_columns = (size_t)reduced->columns;
	int i, j;

	for (i = 0; i < model->m; i++)
		if (model->constraints[i].type != optm_CONS_T_EQ)
			original_columns++;
	for (i = 0; i < reduced->rows; i++)
		if (reduced->row_type[i] != optm_CONS_T_EQ)
			reduced_columns++;
	for (j = 0; j < model->n; j++) {
		const struct optm_VariableBound *bound = model->bounds + j;

		if (bound->b_type == optm_BOUND_T_FR)
			original_columns++;
		if (bound->b_type == optm_BOUND_T_BS) {
			if (bound->ub == bound->lb)
				original_columns--;
			else if (bound->ub > bound->lb) {
				original_rows++;
				original_columns++;
			}
		}
	}
	for (j = 0; j < reduced->columns; j++) {
		const struct optm_VariableBound *bound = reduced->bounds + j;

		if (bound->b_type == optm_BOUND_T_FR)
			reduced_columns++;
		if (bound->b_type == optm_BOUND_T_BS) {
			if (bound->ub == bound->lb)
				reduced_columns--;
			else if (bound->ub > bound->lb) {
				reduced_rows++;
				reduced_columns++;
			}
		}
	}
	if (reduced_rows > original_rows &&
		reduced_columns >= original_columns) {
		if (getenv(SIMPLEX_ENV_PROFILE))
			fprintf(stderr,
				"presolve: Pan equality form original=%lu/%lu reduced=%lu/%lu selecting original\n",
				(unsigned long)original_rows,
				(unsigned long)original_columns,
				(unsigned long)reduced_rows,
				(unsigned long)reduced_columns);
		return 1;
	}
	return 0;
}
