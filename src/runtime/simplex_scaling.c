/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_scaling.c
 * @brief Global row/column equilibration for the revised-simplex kernel.
 */

#include "simplex_solve_context.h"
#include "simplex_execution.h"
#include "policy/resources.h"
#include "policy/controls.h"
#include "simplex_scaling.h"
#include "utils.h"
#include "simplex_storage.h"

#include <lp_simplex/status.h>

#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>

static double scaling_balancing_factor(
	const double maximum, const double current)
{
	int current_exponent, maximum_exponent;
	int delta, exponent;

	if (maximum <= 0.)
		return 1.;
	(void)frexp(maximum, &maximum_exponent);
	(void)frexp(current, &current_exponent);
	current_exponent--;
	/* The coefficient is already scaled by current.  Move half of its binary
	 * exponent into this row/column, as in iterative geometric equilibration. */
	delta = -maximum_exponent / 2;
	exponent = current_exponent + delta;
	exponent = __lp_simplex_MAX__(SIMPLEX_SCALING_MIN_EXPONENT,
		__lp_simplex_MIN__(SIMPLEX_SCALING_MAX_EXPONENT, exponent));
	return ldexp(1., exponent - current_exponent);
}

static int scaling_storage(int rows, int columns, size_t limit, size_t *bytes)
{
	size_t row_bytes, column_bytes;
	*bytes = 0;
	if (rows < 0 || columns < 0 ||
		!simplex_storage_bytes(
			(size_t)rows, sizeof(double), limit, &row_bytes) ||
		!simplex_storage_bytes((size_t)columns, sizeof(double), limit,
			&column_bytes) ||
		row_bytes > limit - column_bytes)
		return 0;
	*bytes = row_bytes + column_bytes;
	return 1;
}

static int scaling_allocate(
	struct simplex_Scaling *scaling, const int rows, const int columns)
{
	size_t i, bytes;
	lp_simplex_memset(scaling, 0, sizeof(*scaling));
	if (!scaling_storage(rows, columns, (size_t)-1, &bytes))
		return lp_simplex_EXIT_FAILURE;
	if (bytes > 0) {
		scaling->storage = (double *)lp_simplex_malloc(bytes);
		if (scaling->storage == NULL)
			return lp_simplex_EXIT_FAILURE;
	}
	scaling->rows = rows;
	scaling->columns = columns;
	scaling->row = scaling->storage;
	scaling->column =
		scaling->storage != NULL ? scaling->storage + rows : NULL;
	for (i = 0; i < bytes / sizeof(double); i++)
		scaling->storage[i] = 1.;
	return lp_simplex_EXIT_SUCCESS;
}

static int scaling_compute(
	const struct simplex_Problem *problem, struct simplex_Scaling *scaling)
{
	double *maximum;
	int pass, i, j, k;
	int capacity = __lp_simplex_MAX__(problem->rows, problem->columns);
	size_t bytes;

	if (capacity == 0)
		return lp_simplex_EXIT_SUCCESS;
	if (capacity < 0 || !simplex_storage_bytes((size_t)capacity,
				    sizeof(double), (size_t)-1, &bytes))
		return lp_simplex_EXIT_FAILURE;
	maximum = (double *)lp_simplex_malloc(bytes);
	if (maximum == NULL)
		return lp_simplex_EXIT_FAILURE;
	for (pass = 0; pass < SIMPLEX_SCALING_PASSES; pass++) {
		lp_simplex_memset(
			maximum, 0, (size_t)problem->rows * sizeof(double));
		for (j = 0; j < problem->columns; j++)
			for (k = problem->matrix.column_start[j];
				k < problem->matrix.column_start[j + 1]; k++) {
				i = problem->matrix.row_index[k];
				maximum[i] = __lp_simplex_MAX__(maximum[i],
					__lp_simplex_ABS__(
						problem->matrix.value[k] *
						scaling->row[i] *
						scaling->column[j]));
			}
		for (i = 0; i < problem->rows; i++)
			scaling->row[i] *= scaling_balancing_factor(
				maximum[i], scaling->row[i]);

		lp_simplex_memset(
			maximum, 0, (size_t)problem->columns * sizeof(double));
		for (j = 0; j < problem->columns; j++)
			for (k = problem->matrix.column_start[j];
				k < problem->matrix.column_start[j + 1]; k++) {
				i = problem->matrix.row_index[k];
				maximum[j] = __lp_simplex_MAX__(maximum[j],
					__lp_simplex_ABS__(
						problem->matrix.value[k] *
						scaling->row[i] *
						scaling->column[j]));
			}
		for (j = 0; j < problem->columns; j++)
			scaling->column[j] *= scaling_balancing_factor(
				maximum[j], scaling->column[j]);
	}
	lp_simplex_free(maximum);
	for (i = 0; i < problem->rows; i++)
		if (scaling->row[i] != 1.)
			scaling->active = 1;
	for (j = 0; j < problem->columns; j++)
		if (scaling->column[j] != 1.)
			scaling->active = 1;
	return lp_simplex_EXIT_SUCCESS;
}

/* Powers of two are exact in the normal range. A subnormal result must
 * round-trip: silent underflow would change the LP used for certification. */
static int scaling_product(double value, double factor, double *scaled)
{
	*scaled = value * factor;
	return isfinite(*scaled) && (value == 0. || fabs(*scaled) >= DBL_MIN ||
					    *scaled / factor == value);
}

static int scaling_copy_vectors(const struct simplex_Problem *problem,
	struct simplex_Problem *scaled, const struct simplex_Scaling *scaling)
{
	int i, j;

	for (j = 0; j < problem->columns; j++) {
		double column = scaling->column[j];

		if (!scaling_product(problem->objective[j], column,
			    &scaled->objective[j]))
			return 0;
		scaled->bounds[j] = problem->bounds[j];
		if (scaled->bounds[j].b_type == optm_BOUND_T_LO ||
			scaled->bounds[j].b_type == optm_BOUND_T_BS)
			if (!scaling_product(problem->bounds[j].lb, 1. / column,
				    &scaled->bounds[j].lb))
				return 0;
		if (scaled->bounds[j].b_type == optm_BOUND_T_UP ||
			scaled->bounds[j].b_type == optm_BOUND_T_BS)
			if (!scaling_product(problem->bounds[j].ub, 1. / column,
				    &scaled->bounds[j].ub))
				return 0;
	}
	for (i = 0; i < problem->rows; i++) {
		if (!scaling_product(
			    problem->rhs[i], scaling->row[i], &scaled->rhs[i]))
			return 0;
		scaled->row_type[i] = problem->row_type[i];
	}
	return 1;
}

static int scaling_copy_matrix(const struct simplex_Problem *problem,
	struct simplex_Problem *scaled, const struct simplex_Scaling *scaling)
{
	int i, j, k;

	for (j = 0; j < problem->columns; j++)
		for (k = problem->matrix.column_start[j];
			k < problem->matrix.column_start[j + 1]; k++) {
			i = problem->matrix.row_index[k];
			if (!scaling_product(problem->matrix.value[k],
				    scaling->row[i] * scaling->column[j],
				    &scaled->matrix.value[k]))
				return 0;
		}
	for (i = 0; i < problem->rows; i++)
		for (k = problem->matrix.row_start[i];
			k < problem->matrix.row_start[i + 1]; k++) {
			j = problem->matrix.column_index[k];
			if (!scaling_product(problem->matrix.row_value[k],
				    scaling->row[i] * scaling->column[j],
				    &scaled->matrix.row_value[k]))
				return 0;
		}
	return 1;
}

struct simplex_Problem *simplex_scaling_create_problem(
	const struct simplex_Problem *problem, struct simplex_Scaling *scaling)
{
	struct simplex_Problem *scaled = NULL;
	int i;

	if (scaling_allocate(scaling, problem->rows, problem->columns) ==
			lp_simplex_EXIT_FAILURE ||
		scaling_compute(problem, scaling) == lp_simplex_EXIT_FAILURE) {
		simplex_scaling_destroy(scaling);
		return NULL;
	}
	if (scaling->active) {
		scaled = simplex_problem_create_with_pattern(&problem->matrix);
		if (!scaled) {
			simplex_scaling_destroy(scaling);
			return NULL;
		}
		if (!scaling_copy_vectors(problem, scaled, scaling) ||
			!scaling_copy_matrix(problem, scaled, scaling)) {
			/* Equilibration is optional. Keep the exact input rather than
			 * give the kernel an overflowed or rounded-away coefficient. */
			simplex_problem_free(scaled);
			for (i = 0; i < scaling->rows; i++)
				scaling->row[i] = 1.;
			for (i = 0; i < scaling->columns; i++)
				scaling->column[i] = 1.;
			scaling->active = 0;
			if (getenv(SIMPLEX_ENV_PROFILE))
				fprintf(stderr,
					"scaling: non-reversible numeric conversion; using original coordinates\n");
		}
	}
	/* Scaling changes only values.  Identity scaling can borrow the entire
	 * immutable problem; otherwise only numeric arrays need new storage. */
	if (!scaling->active) {
		scaled = (struct simplex_Problem *)lp_simplex_malloc(
			sizeof(*scaled));
		if (scaled != NULL) {
			*scaled = *problem;
			scaled->owns_vectors = 0;
			scaled->vector_storage = NULL;
			scaled->matrix.owns_storage = 0;
		}
	}
	if (scaled == NULL) {
		simplex_scaling_destroy(scaling);
		return NULL;
	}
	scaled->report = problem->report;
	scaled->execution = scaling->active ? simplex_solve_context_derive(
						      &problem->execution,
						      SIMPLEX_EXECUTION_SCALING)
					    : problem->execution;
	if (getenv(SIMPLEX_ENV_PROFILE) != NULL) {
		double row_min = 1., row_max = 1.;
		double column_min = 1., column_max = 1.;
		int i, j;

		for (i = 0; i < problem->rows; i++) {
			row_min = __lp_simplex_MIN__(row_min, scaling->row[i]);
			row_max = __lp_simplex_MAX__(row_max, scaling->row[i]);
		}
		for (j = 0; j < problem->columns; j++) {
			column_min = __lp_simplex_MIN__(
				column_min, scaling->column[j]);
			column_max = __lp_simplex_MAX__(
				column_max, scaling->column[j]);
		}
		fprintf(stderr,
			"scaling: passes=%d active=%d row=[%.3g,%.3g] column=[%.3g,%.3g]\n",
			SIMPLEX_SCALING_PASSES, scaling->active, row_min,
			row_max, column_min, column_max);
	}
	return scaled;
}

void simplex_scaling_recover_primal(
	const struct simplex_Scaling *scaling, double *x)
{
	int j;

	for (j = 0; j < scaling->columns; j++)
		x[j] *= scaling->column[j];
}

void simplex_scaling_recover_dual(
	const struct simplex_Scaling *scaling, double *dual)
{
	int i;

	for (i = 0; i < scaling->rows; i++)
		dual[i] *= scaling->row[i];
}

void simplex_scaling_destroy(struct simplex_Scaling *scaling)
{
	if (scaling == NULL)
		return;
	lp_simplex_free(scaling->storage);
	lp_simplex_memset(scaling, 0, sizeof(*scaling));
}
