/* Global row/column equilibration for the revised-simplex kernel. */
#include "simplex_scaling.h"
#include "utils.h"

#include <lp_simplex/status.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>


#define SIMPLEX_SCALING_PASSES 4
#define SIMPLEX_SCALING_MIN_EXPONENT (-48)
#define SIMPLEX_SCALING_MAX_EXPONENT 48


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


static int scaling_allocate(
		struct simplex_Scaling *scaling, const int rows, const int columns)
{
	int i;
	lp_simplex_memset(scaling, 0, sizeof(*scaling));
	if (rows + columns > 0) {
		scaling->storage = (double *)lp_simplex_malloc(
			(size_t)(rows + columns) * sizeof(double));
		if (scaling->storage == NULL)
			return lp_simplex_EXIT_FAILURE;
	}
	scaling->rows = rows;
	scaling->columns = columns;
	scaling->row = scaling->storage;
	scaling->column = scaling->storage != NULL
		? scaling->storage + rows : NULL;
	for (i = 0; i < rows + columns; i++)
		scaling->storage[i] = 1.;
	return lp_simplex_EXIT_SUCCESS;
}


static int scaling_compute(
		const struct simplex_Problem *problem,
		struct simplex_Scaling *scaling)
{
	double *maximum;
	int pass, i, j, k;
	int capacity = __lp_simplex_MAX__(problem->rows, problem->columns);
	if (capacity == 0)
		return lp_simplex_EXIT_SUCCESS;
	maximum = (double *)lp_simplex_malloc((size_t)capacity * sizeof(double));
	if (maximum == NULL)
		return lp_simplex_EXIT_FAILURE;
	for (pass = 0; pass < SIMPLEX_SCALING_PASSES; pass++) {
		lp_simplex_memset(maximum, 0,
			(size_t)problem->rows * sizeof(double));
		for (j = 0; j < problem->columns; j++)
			for (k = problem->matrix.column_start[j];
			     k < problem->matrix.column_start[j + 1]; k++) {
				i = problem->matrix.row_index[k];
				maximum[i] = __lp_simplex_MAX__(maximum[i],
					__lp_simplex_ABS__(problem->matrix.value[k] *
						scaling->row[i] * scaling->column[j]));
			}
		for (i = 0; i < problem->rows; i++)
			scaling->row[i] *= scaling_balancing_factor(
				maximum[i], scaling->row[i]);

		lp_simplex_memset(maximum, 0,
			(size_t)problem->columns * sizeof(double));
		for (j = 0; j < problem->columns; j++)
			for (k = problem->matrix.column_start[j];
			     k < problem->matrix.column_start[j + 1]; k++) {
				i = problem->matrix.row_index[k];
				maximum[j] = __lp_simplex_MAX__(maximum[j],
					__lp_simplex_ABS__(problem->matrix.value[k] *
						scaling->row[i] * scaling->column[j]));
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


static void scaling_copy_vectors(
		const struct simplex_Problem *problem,
		struct simplex_Problem *scaled,
		const struct simplex_Scaling *scaling)
{
	int i, j;
	for (j = 0; j < problem->columns; j++) {
		double column = scaling->column[j];
		scaled->objective[j] = problem->objective[j] * column;
		scaled->bounds[j] = problem->bounds[j];
		if (scaled->bounds[j].b_type == optm_BOUND_T_LO ||
		    scaled->bounds[j].b_type == optm_BOUND_T_BS)
			scaled->bounds[j].lb /= column;
		if (scaled->bounds[j].b_type == optm_BOUND_T_UP ||
		    scaled->bounds[j].b_type == optm_BOUND_T_BS)
			scaled->bounds[j].ub /= column;
	}
	for (i = 0; i < problem->rows; i++) {
		scaled->rhs[i] = problem->rhs[i] * scaling->row[i];
		scaled->row_type[i] = problem->row_type[i];
	}
}


static void scaling_copy_matrix(
		const struct simplex_Problem *problem,
		struct simplex_Problem *scaled,
		const struct simplex_Scaling *scaling)
{
	int i, j, k;
	lp_simplex_memcpy(scaled->matrix.column_start,
		problem->matrix.column_start,
		(size_t)(problem->columns + 1) * sizeof(int));
	lp_simplex_memcpy(scaled->matrix.row_index, problem->matrix.row_index,
		(size_t)problem->matrix.nonzeros * sizeof(int));
	for (j = 0; j < problem->columns; j++)
		for (k = problem->matrix.column_start[j];
		     k < problem->matrix.column_start[j + 1]; k++) {
			i = problem->matrix.row_index[k];
			scaled->matrix.value[k] = problem->matrix.value[k] *
				scaling->row[i] * scaling->column[j];
		}
	lp_simplex_memcpy(scaled->matrix.row_start, problem->matrix.row_start,
		(size_t)(problem->rows + 1) * sizeof(int));
	lp_simplex_memcpy(scaled->matrix.column_index,
		problem->matrix.column_index,
		(size_t)problem->matrix.nonzeros * sizeof(int));
	for (i = 0; i < problem->rows; i++)
		for (k = problem->matrix.row_start[i];
		     k < problem->matrix.row_start[i + 1]; k++) {
			j = problem->matrix.column_index[k];
			scaled->matrix.row_value[k] = problem->matrix.row_value[k] *
				scaling->row[i] * scaling->column[j];
		}
}


struct simplex_Problem *simplex_scaling_create_problem(
		const struct simplex_Problem *problem,
		struct simplex_Scaling *scaling)
{
	struct simplex_Problem *scaled;
	if (scaling_allocate(scaling, problem->rows, problem->columns) ==
		lp_simplex_EXIT_FAILURE ||
	    scaling_compute(problem, scaling) == lp_simplex_EXIT_FAILURE) {
		simplex_scaling_destroy(scaling);
		return NULL;
	}
	scaled = simplex_problem_create_sparse(problem->rows, problem->columns,
		problem->matrix.nonzeros);
	if (scaled == NULL) {
		simplex_scaling_destroy(scaling);
		return NULL;
	}
	scaling_copy_vectors(problem, scaled, scaling);
	scaling_copy_matrix(problem, scaled, scaling);
	if (getenv("LP_SIMPLEX_PROFILE") != NULL) {
		double row_min = 1., row_max = 1.;
		double column_min = 1., column_max = 1.;
		int i, j;
		for (i = 0; i < problem->rows; i++) {
			row_min = __lp_simplex_MIN__(row_min, scaling->row[i]);
			row_max = __lp_simplex_MAX__(row_max, scaling->row[i]);
		}
		for (j = 0; j < problem->columns; j++) {
			column_min = __lp_simplex_MIN__(column_min,
				scaling->column[j]);
			column_max = __lp_simplex_MAX__(column_max,
				scaling->column[j]);
		}
		fprintf(stderr,
			"scaling: passes=%d active=%d row=[%.3g,%.3g] column=[%.3g,%.3g]\n",
			SIMPLEX_SCALING_PASSES, scaling->active,
			row_min, row_max, column_min, column_max);
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
