/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#include "simplex_phase.h"
#include "simplex_pivot.h"
#include "simplex_tableau.h"
#include "linalg.h"
#include "utils.h"
#include <lp_simplex/status.h>

#define SIMPLEX_FEASIBILITY_TOLERANCE 1e-5


static int simplex_run_phase_one(
		double **table, int *ld, int **basis, int **constraint_types,
		int *variable_count, int *iteration, int *status,
		const struct optm_LinearConstraint *constraints,
		int m, int n, const char *criteria, int iteration_limit)
{
	int rows, columns;
	int slack_count, artificial_count;

	simplex_tableau_compute_size(constraints, m, n, &rows, &columns);
	*ld = columns;
	if (simplex_tableau_create(table, basis, constraint_types, m, rows, *ld) ==
	    lp_simplex_EXIT_FAILURE) {
		*status = lp_simplex_MemoryAllocError;
		return lp_simplex_EXIT_FAILURE;
	}
	simplex_tableau_fill_constraint_types(constraints, *constraint_types, m);
	simplex_tableau_fill_constraints(*table, *ld, constraints, rows, columns, m, n);
	slack_count = simplex_tableau_add_slack(*table, *ld, *constraint_types, m, n);
	artificial_count = simplex_tableau_add_artificial(
		*table, *ld, *constraint_types, m, n, slack_count);
	*variable_count = n + slack_count + artificial_count;
	if (m > *variable_count) {
		*status = lp_simplex_OverDetermination;
		goto FAILURE;
	}
	simplex_tableau_initialize_basis(*basis, *constraint_types, m, n, slack_count);
	simplex_tableau_initialize_phase_one_objective(
		*table, *ld, *constraint_types, m, columns);

	switch (simplex_run_pivots(iteration, *table, *ld, *basis, m,
				   *variable_count, n + slack_count,
				   criteria, iteration_limit)) {
	case 0:
		*status = lp_simplex_ExceedIterLimit;
		goto FAILURE;
	case 1:
		if ((*table)[columns - 1] > SIMPLEX_FEASIBILITY_TOLERANCE) {
			*status = lp_simplex_Infeasibility;
			goto FAILURE;
		}
		simplex_tableau_remove_artificial_basis(
			*table, *ld, *basis, m, n + slack_count, *variable_count);
		simplex_tableau_remove_artificial_columns(
			*table, *ld, m, n + slack_count, artificial_count);
		*variable_count = n + slack_count;
		return lp_simplex_EXIT_SUCCESS;
	case 2:
		*status = lp_simplex_Unboundedness;
		goto FAILURE;
	case 4:
		*status = lp_simplex_MemoryAllocError;
		goto FAILURE;
	case 9:
		*status = lp_simplex_PrecisionError;
		goto FAILURE;
	}

FAILURE:
	simplex_tableau_destroy(*table, *basis, *constraint_types);
	return lp_simplex_EXIT_FAILURE;
}


static int simplex_run_phase_two(
		double *table, int ld, int *basis, int *constraint_types,
		int *iteration, int *status, int m, int variable_count,
		const char *criteria, int iteration_limit)
{
	switch (simplex_run_pivots(iteration, table, ld, basis, m, variable_count,
				   variable_count, criteria, iteration_limit)) {
	case 0:
		*status = lp_simplex_ExceedIterLimit;
		break;
	case 1:
		*status = lp_simplex_Success;
		return lp_simplex_EXIT_SUCCESS;
	case 2:
		*status = lp_simplex_Unboundedness;
		break;
	case 4:
		*status = lp_simplex_MemoryAllocError;
		break;
	case 9:
		*status = lp_simplex_PrecisionError;
		break;
	}
	simplex_tableau_destroy(table, basis, constraint_types);
	return lp_simplex_EXIT_FAILURE;
}


int simplex_solve_standard(
		const double *objective,
		const struct optm_LinearConstraint *constraints,
		int m, int n, const char *criteria, int iteration_limit,
		double *x, double *value, int *status)
{
	int i, j, k;
	int ld, variable_count;
	int iteration = 0;
	int *basis = NULL;
	double *table = NULL;
	int *constraint_types = NULL;

	if (status == NULL)
		return lp_simplex_EXIT_FAILURE;
	if (objective == NULL || constraints == NULL || x == NULL || value == NULL ||
	    m <= 0 || n <= 0 || iteration_limit <= 0) {
		*status = lp_simplex_CondUnsatisfied;
		return lp_simplex_EXIT_FAILURE;
	}
	if (criteria == NULL)
		criteria = "";
	for (k = 0; k < m; k++) {
		if (constraints[k].coef == NULL || constraints[k].type < optm_CONS_T_EQ ||
		    constraints[k].type > optm_CONS_T_LE) {
			*status = lp_simplex_CondUnsatisfied;
			return lp_simplex_EXIT_FAILURE;
		}
	}
	if (simplex_run_phase_one(&table, &ld, &basis, &constraint_types,
				  &variable_count, &iteration, status,
				  constraints, m, n, criteria, iteration_limit) ==
	    lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;

	for (j = 0; j < n; j++)
		table[j] = -objective[j];
	for (i = 0; i < m; i++) {
		int row = (i + 1) * ld;
		double ratio = -table[basis[i]];
		lp_simplex_linalg_daxpy(variable_count + 1, ratio,
					table + row, 1, table, 1);
	}
	if (simplex_run_phase_two(table, ld, basis, constraint_types, &iteration,
				  status, m, variable_count, criteria,
				  iteration_limit) == lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;

	*value = table[variable_count];
	lp_simplex_memset(x, 0, (size_t)n * sizeof(double));
	for (i = 0; i < m; i++) {
		if (basis[i] < n)
			x[basis[i]] = table[variable_count + (i + 1) * ld];
	}
	simplex_tableau_destroy(table, basis, constraint_types);
	return lp_simplex_EXIT_SUCCESS;
}
