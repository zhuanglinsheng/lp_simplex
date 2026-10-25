/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_tableau.c
 * @brief Tableau allocation, initialization and basis bookkeeping.
 */

#include "simplex_tableau.h"
#include "simplex_pivot.h"
#include "linalg.h"
#include "utils.h"
#include "simplex_tableau_storage.h"

#include <lp_simplex/status.h>
#include <float.h>
#include <math.h>

/* To create in heap (need to be released) simplex table, index set of basis
 * and constraint type recorder
 */
int simplex_tableau_create(double **table, int **basis, int **constypes,
	const int m, const int nrow, const int ncol)
{
	size_t table_bytes, index_bytes;
	*table = NULL;
	*basis = NULL;
	*constypes = NULL;
	if (m <= 0 || m == INT_MAX || nrow != m + 1 ||
		!simplex_tableau_storage_size(
			nrow, ncol, sizeof(double), &table_bytes) ||
		!simplex_tableau_storage_size(1, m, sizeof(int), &index_bytes))
		return lp_simplex_EXIT_FAILURE;
	*table = (double *)lp_simplex_malloc(table_bytes);
	if (!*table)
		goto FAILURE;
	*basis = (int *)lp_simplex_malloc(index_bytes);
	if (!*basis)
		goto FAILURE;
	*constypes = (int *)lp_simplex_malloc(index_bytes);
	if (!*constypes)
		goto FAILURE;
	return lp_simplex_EXIT_SUCCESS;
FAILURE:
	simplex_tableau_destroy(*table, *basis, *constypes);
	*table = NULL;
	*basis = *constypes = NULL;
	return lp_simplex_EXIT_FAILURE;
}

void simplex_tableau_destroy(double *table, int *basis, int *constypes)
{
	if (table)
		lp_simplex_free(table);
	if (basis)
		lp_simplex_free(basis);
	if (constypes)
		lp_simplex_free(constypes);
}

/* Fill in constraint type array from "constraints"
 *
 * Constraints rhs are transformed to be nonnegative,
 * "LE" and "GE" types are transformed respectively
 */
void simplex_tableau_fill_constraint_types(
	const struct optm_LinearConstraint *constraints, int *constypes,
	const int m)
{
	int i;

	for (i = 0; i < m; i++) {
		const struct optm_LinearConstraint *cons = constraints + i;

		if (cons->rhs >= 0)
			constypes[i] = cons->type;
		else {
			switch (cons->type) {
			case optm_CONS_T_EQ:
				constypes[i] = optm_CONS_T_EQ;
				break;
			case optm_CONS_T_GE:
				constypes[i] = optm_CONS_T_LE;
				break;
			case optm_CONS_T_LE:
				constypes[i] = optm_CONS_T_GE;
				break;
			}
		}
	}
}

/* Fill in coef and rhs of constraints
 *
 * Constraints rhs are transformed to be nonnegative
 */
void simplex_tableau_fill_constraints(double *table, const int ldtable,
	const struct optm_LinearConstraint *constraints, const int nrow,
	const int ncol, const int m, const int n)
{
	int i, j;

	lp_simplex_memset(
		table, 0., (size_t)nrow * (size_t)ldtable * sizeof(double));

	for (i = 0; i < m; i++) {
		const struct optm_LinearConstraint *cons = constraints + i;
		int row = (i + 1) * ldtable;

		if (cons->rhs >= 0) {
			table[ncol - 1 + row] = cons->rhs;
			lp_simplex_memcpy(
				table + row, cons->coef, n * sizeof(double));
		} else {
			table[ncol - 1 + row] = -cons->rhs;
			for (j = 0; j < n; j++)
				table[j + row] = -cons->coef[j];
		}
	}
}

/* Determine the size of (basic) simplex table
 *
 * "GE" constraint has a slack var and an artificial var, hence will generate
 * an additional variable than usual
 */
int simplex_tableau_compute_size(
	const struct optm_LinearConstraint *constraints, const int m,
	const int n, int *nrow, int *ncol)
{
	int i, rows, columns;
	size_t bytes;
	*nrow = *ncol = 0;
	if (!constraints || m <= 0 || n <= 0 || m == INT_MAX ||
		n > INT_MAX - m - 1)
		return lp_simplex_EXIT_FAILURE;
	rows = m + 1;
	columns = m + n + 1;
	/* Even the minimum layout must fit before reading row metadata. */
	if (!simplex_tableau_storage_size(
		    rows, columns, sizeof(double), &bytes))
		return lp_simplex_EXIT_FAILURE;
	for (i = 0; i < m; i++) {
		const struct optm_LinearConstraint *cons = constraints + i;
		if ((optm_CONS_T_GE == cons->type && cons->rhs >= 0) ||
			(optm_CONS_T_LE == cons->type && cons->rhs < 0)) {
			if (columns == INT_MAX)
				return lp_simplex_EXIT_FAILURE;
			columns++;
		}
	}
	if (!simplex_tableau_storage_size(
		    rows, columns, sizeof(double), &bytes))
		return lp_simplex_EXIT_FAILURE;
	*nrow = rows;
	*ncol = columns;
	return lp_simplex_EXIT_SUCCESS;
}

/* Add slack variables (GE, LE) to simplex table
 * Return the number of slack variables
 */
int simplex_tableau_add_slack(double *table, const int ldtable,
	const int *constypes, const int m, const int n)
{
	int i, nslack = 0;

	for (i = 0; i < m; i++) {
		if (optm_CONS_T_GE == constypes[i]) {
			table[n + nslack + (i + 1) * ldtable] = -1.;
			nslack++;
		}
		if (optm_CONS_T_LE == constypes[i]) {
			table[n + nslack + (i + 1) * ldtable] = 1.;
			nslack++;
		}
	}
	return nslack;
}

/* Add artificial variables (GE, EQ) to simplex table
 * Return the number of artificial variables
 */
int simplex_tableau_add_artificial(double *table, const int ldtable,
	const int *constypes, const int m, const int n, const int nslack)
{
	int i, nartif = 0;

	for (i = 0; i < m; i++) {
		if (optm_CONS_T_LE != constypes[i]) {
			table[n + nslack + nartif + (i + 1) * ldtable] = 1.;
			table[n + nslack + nartif] = -1.;
			nartif++;
		}
	}
	return nartif;
}

/* Fill in the basis index set of artificial LP
 */
void simplex_tableau_initialize_basis(int *basis, const int *constypes,
	const int m, const int n, const int nslack)
{
	int i, tmp_nbasis = 0, tmp_nslack = 0, tmp_nartif = 0;

	for (i = 0; i < m; i++) {
		switch (constypes[i]) {
		case optm_CONS_T_EQ:
			*(basis + tmp_nbasis) = n + nslack + tmp_nartif;
			tmp_nartif++;
			break;
		case optm_CONS_T_GE:
			*(basis + tmp_nbasis) = n + nslack + tmp_nartif;
			tmp_nartif++;
			tmp_nslack++;
			break;
		case optm_CONS_T_LE:
			*(basis + tmp_nbasis) = n + tmp_nslack;
			tmp_nslack++;
			break;
		}
		tmp_nbasis++;
	}
}

void simplex_tableau_initialize_phase_one_objective(double *table,
	const int ldtable, const int *constypes, const int m, const int ncol)
{
	int i, rowi;

	for (i = 0; i < m; i++) {
		if (optm_CONS_T_LE == constypes[i])
			continue;
		rowi = (i + 1) * ldtable;
		lp_simplex_linalg_daxpy(ncol, 1, table + rowi, 1, table, 1);
	}
}

int simplex_tableau_remove_artificial_basis(double *table, int ldtable,
	int *basis, const int m, const int nreal, int nvar,
	double pivot_tolerance, int *iterations, int limit)
{
	int i, j, r, q;
	double ele, maxv, row_scale;

	if (!iterations || *iterations < 0 || limit < *iterations)
		return SIMPLEX_TABLEAU_CLEANUP_NUMERICAL;
	for (i = 0; i < m; i++) {
		if (basis[i] < nreal)
			continue;
		maxv = 0.;
		row_scale = 0.;
		q = nvar;
		for (j = 0; j < nreal; j++) {
			ele = __lp_simplex_ABS__(table[j + (i + 1) * ldtable]);
			if (!isfinite(ele))
				return SIMPLEX_TABLEAU_CLEANUP_NUMERICAL;
			row_scale = __lp_simplex_MAX__(row_scale, ele);
			if (is_in_arri(j, basis, m))
				continue;
			if (ele > maxv) {
				maxv = ele;
				q = j;
			}
		}
		/* Small is not redundant. Compare the candidate with this row's
		 * real coefficients. An unresolved basic residual must never
		 * silently erase an equation. Keep removed artificial indices
		 * only on exactly redundant rows, skipped during phase II. */
		if (row_scale == 0. && table[nvar + (i + 1) * ldtable] == 0.)
			continue;
		if (maxv < DBL_MIN || maxv / row_scale <= pivot_tolerance ||
			!isfinite(table[nvar + (i + 1) * ldtable]))
			return SIMPLEX_TABLEAU_CLEANUP_NUMERICAL;
		if (*iterations >= limit)
			return SIMPLEX_TABLEAU_CLEANUP_LIMIT;
		(*iterations)++;
		simplex_apply_pivot(table, ldtable, m, nvar, i, q, 1, 1, 0);
		for (r = 0; r <= m; r++)
			for (j = 0; j <= nvar; j++)
				if (!isfinite(table[r * ldtable + j]))
					return SIMPLEX_TABLEAU_CLEANUP_NUMERICAL;
		basis[i] = q;
	}
	return SIMPLEX_TABLEAU_CLEANUP_DONE;
}

void simplex_tableau_remove_artificial_columns(double *table, const int ldtable,
	const int m, const int nreal, const int nartif)
{
	int i, rowi;

	if (nartif <= 0) /* Delete artificial columns */
		return;
	for (i = 0; i < m + 1; i++) {
		rowi = i * ldtable;
		table[nreal + rowi] = table[nreal + nartif + rowi];
	}
}
