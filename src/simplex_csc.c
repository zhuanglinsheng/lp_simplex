/* Sparse immutable column storage used by revised simplex. */
#include "simplex_csc.h"
#include "utils.h"
#include <lp_simplex/status.h>


int simplex_csc_from_model(
		const struct lp_Model *model, struct simplex_CscMatrix *matrix)
{
	int i, j, k, nonzeros = 0;
	matrix->rows = model->m;
	matrix->columns = model->n;
	matrix->nonzeros = 0;
	matrix->column_start = NULL;
	matrix->row_index = NULL;
	matrix->value = NULL;
	for (j = 0; j < model->n; j++) {
		for (i = 0; i < model->m; i++) {
			if (model->constraints[i].coef[j] != 0.)
				nonzeros++;
		}
	}
	matrix->column_start = (int *)lp_simplex_malloc(
		(size_t)(model->n + 1) * sizeof(int));
	if (nonzeros > 0) {
		matrix->row_index = (int *)lp_simplex_malloc(
			(size_t)nonzeros * sizeof(int));
		matrix->value = (double *)lp_simplex_malloc(
			(size_t)nonzeros * sizeof(double));
	}
	if (matrix->column_start == NULL ||
	    (nonzeros > 0 && (matrix->row_index == NULL || matrix->value == NULL))) {
		simplex_csc_destroy(matrix);
		return lp_simplex_EXIT_FAILURE;
	}
	k = 0;
	for (j = 0; j < model->n; j++) {
		matrix->column_start[j] = k;
		for (i = 0; i < model->m; i++) {
			double coefficient = model->constraints[i].coef[j];
			if (coefficient == 0.)
				continue;
			matrix->row_index[k] = i;
			matrix->value[k] = coefficient;
			k++;
		}
	}
	matrix->column_start[model->n] = k;
	matrix->nonzeros = k;
	return lp_simplex_EXIT_SUCCESS;
}


void simplex_csc_destroy(struct simplex_CscMatrix *matrix)
{
	if (matrix == NULL)
		return;
	lp_simplex_free(matrix->column_start);
	lp_simplex_free(matrix->row_index);
	lp_simplex_free(matrix->value);
	matrix->column_start = NULL;
	matrix->row_index = NULL;
	matrix->value = NULL;
	matrix->nonzeros = 0;
}


double simplex_csc_column_dot(
		const struct simplex_CscMatrix *matrix,
		const int column, const double *vector)
{
	int k;
	double result = 0.;
	for (k = matrix->column_start[column];
	     k < matrix->column_start[column + 1]; k++)
		result += vector[matrix->row_index[k]] * matrix->value[k];
	return result;
}


void simplex_csc_column_to_dense(
		const struct simplex_CscMatrix *matrix,
		const int structural_columns, const int variable,
		double *dense)
{
	int k;
	lp_simplex_memset(dense, 0, (size_t)matrix->rows * sizeof(double));
	if (variable >= structural_columns) {
		dense[variable - structural_columns] = -1.;
		return;
	}
	for (k = matrix->column_start[variable];
	     k < matrix->column_start[variable + 1]; k++)
		dense[matrix->row_index[k]] = matrix->value[k];
}


void simplex_csc_column_axpy(
		const struct simplex_CscMatrix *matrix,
		const int structural_columns, const int variable,
		const double scale, double *dense)
{
	int k;
	if (variable >= structural_columns) {
		dense[variable - structural_columns] -= scale;
		return;
	}
	for (k = matrix->column_start[variable];
	     k < matrix->column_start[variable + 1]; k++)
		dense[matrix->row_index[k]] += scale * matrix->value[k];
}
