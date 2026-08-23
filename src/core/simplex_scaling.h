/* Reversible power-of-two equilibration for sparse LP kernels. */
#ifndef LP_SIMPLEX_SCALING_INTERNAL_H
#define LP_SIMPLEX_SCALING_INTERNAL_H

#include "simplex_problem.h"


struct simplex_Scaling {
	int rows;
	int columns;
	double *storage;
	double *row;
	double *column;
	int active;
};


struct simplex_Problem *simplex_scaling_create_problem(
		const struct simplex_Problem *problem,
		struct simplex_Scaling *scaling);

void simplex_scaling_recover_primal(
		const struct simplex_Scaling *scaling, double *x);

void simplex_scaling_recover_dual(
		const struct simplex_Scaling *scaling, double *dual);

void simplex_scaling_destroy(struct simplex_Scaling *scaling);

#endif
