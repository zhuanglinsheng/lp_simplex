#ifndef LP_SIMPLEX_BASIS_INTERNAL_H
#define LP_SIMPLEX_BASIS_INTERNAL_H

#include "simplex_csc.h"
#include "simplex_sparse_lu.h"
#ifdef LP_SIMPLEX_HAVE_KLU
#include <klu.h>
#endif

struct simplex_Basis {
	int rows;
	int structural_columns;
	const struct simplex_CscMatrix *matrix;
	int *index;
#ifdef LP_SIMPLEX_HAVE_KLU
	int *base_column_start;
	int *base_row_index;
	double *base_value;
	klu_symbolic *symbolic;
	klu_numeric *numeric;
	klu_common common;
#else
	struct simplex_SparseLu sparse;
#endif
	double *eta;
	int *eta_pivot;
	int update_count;
	int update_limit;
	int profile_enabled;
	double profile_factor_seconds;
	double profile_ftran_seconds;
	double profile_btran_seconds;
	long profile_factor_calls;
	long profile_ftran_calls;
	long profile_btran_calls;
};

int simplex_basis_create(struct simplex_Basis *basis,
		const struct simplex_CscMatrix *matrix,
		int structural_columns, int *index);
void simplex_basis_destroy(struct simplex_Basis *basis);
int simplex_basis_factorize(struct simplex_Basis *basis);
int simplex_basis_ftran(const struct simplex_Basis *basis, double *vector);
int simplex_basis_btran(const struct simplex_Basis *basis, double *vector);
int simplex_basis_edge_weights(const struct simplex_Basis *basis,
		double *weights);
int simplex_basis_update(struct simplex_Basis *basis,
		int leaving_position, const double *direction);

#endif
