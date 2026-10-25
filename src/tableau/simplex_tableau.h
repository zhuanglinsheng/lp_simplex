/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_tableau.h
 * @brief Tableau allocation, initialization and basis bookkeeping.
 */

#ifndef LP_SIMPLEX_TABLEAU_H
#define LP_SIMPLEX_TABLEAU_H

#include <lp_simplex/model.h>

int simplex_tableau_create(double **table, int **basis, int **constraint_types,
	int m, int rows, int columns);

void simplex_tableau_destroy(double *table, int *basis, int *constraint_types);

void simplex_tableau_fill_constraint_types(
	const struct optm_LinearConstraint *constraints, int *types, int m);

void simplex_tableau_fill_constraints(double *table, int ld,
	const struct optm_LinearConstraint *constraints, int rows, int columns,
	int m, int n);

/* Zero output dimensions on failure; integer offsets and byte counts must fit. */
int simplex_tableau_compute_size(
	const struct optm_LinearConstraint *constraints, int m, int n,
	int *rows, int *columns);

int simplex_tableau_add_slack(
	double *table, int ld, const int *types, int m, int n);

int simplex_tableau_add_artificial(
	double *table, int ld, const int *types, int m, int n, int slack_count);

void simplex_tableau_initialize_basis(
	int *basis, const int *types, int m, int n, int slack_count);

void simplex_tableau_initialize_phase_one_objective(
	double *table, int ld, const int *types, int m, int columns);

enum simplex_TableauCleanup {
	SIMPLEX_TABLEAU_CLEANUP_DONE = 0,
	SIMPLEX_TABLEAU_CLEANUP_NUMERICAL = -1,
	SIMPLEX_TABLEAU_CLEANUP_LIMIT = 1
};

/* Only an exactly zero real row and RHS can be discarded. Each actual cleanup
 * pivot consumes one unit from the same cumulative Phase-I/II lease. */
int simplex_tableau_remove_artificial_basis(double *table, int ld, int *basis,
	int m, int real_columns, int variable_count, double pivot_tolerance,
	int *iterations, int limit);

void simplex_tableau_remove_artificial_columns(
	double *table, int ld, int m, int real_columns, int artificial_count);

#endif
