#ifndef LP_SIMPLEX_CSC_INTERNAL_H
#define LP_SIMPLEX_CSC_INTERNAL_H

#include <lp_simplex/model.h>

struct simplex_CscMatrix {
	int rows;
	int columns;
	int nonzeros;
	int *column_start;
	int *row_index;
	double *value;
	int *row_start;
	int *column_index;
	double *row_value;
	int owns_storage;
};

int simplex_csc_from_model(const struct lp_Model *model,
		struct simplex_CscMatrix *matrix);
void simplex_csc_destroy(struct simplex_CscMatrix *matrix);
double simplex_csc_column_dot(const struct simplex_CscMatrix *matrix,
		int column, const double *vector);
void simplex_csc_column_to_dense(const struct simplex_CscMatrix *matrix,
		int structural_columns, int variable, double *dense);
void simplex_csc_column_axpy(const struct simplex_CscMatrix *matrix,
		int structural_columns, int variable, double scale,
		double *dense);

#endif
