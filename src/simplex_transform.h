#ifndef LP_SIMPLEX_TRANSFORM_H
#define LP_SIMPLEX_TRANSFORM_H

#include <lp_simplex/model.h>

int simplex_transform_alloc(
		int m, int n, double **objective, double **solution,
		double **coefficients,
		struct optm_LinearConstraint **constraints);

void simplex_transform_free(double *objective, double *solution,
		double *coefficients,
		struct optm_LinearConstraint *constraints);

void simplex_transform_size(
		const struct optm_VariableBound *bounds,
		int m, int n, int *new_m, int *new_n);

void simplex_transform_problem(
		const double *objective,
		const struct optm_LinearConstraint *constraints,
		const struct optm_VariableBound *bounds,
		int m, int n, int new_n, double *new_objective,
		double *objective_offset, double *new_coefficients,
		struct optm_LinearConstraint *new_constraints);

void simplex_transform_recover(
		const struct optm_VariableBound *bounds, int n,
		const double *transformed_x,
		double transformed_objective,
		double objective_offset, double *x, double *objective);

#endif
