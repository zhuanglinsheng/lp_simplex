/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#include "simplex_tableau.h"
#include "simplex_pivot.h"
#include "linalg.h"
#include "utils.h"
#include <lp_simplex/status.h>

/* To create in heap (need to be released) simplex table, index set of basis
 * and constraint type recorder
 */
int simplex_tableau_create(
		double **table, int **basis, int **constypes,
		const int m, const int nrow, const int ncol)
{
	*table = NULL;
	*basis = NULL;
	*constypes = NULL;

	*table = (double *)lp_simplex_malloc(nrow * ncol * sizeof(double));
	if (*table == NULL)
		return lp_simplex_EXIT_FAILURE;
	*basis = (int *)lp_simplex_malloc(m * sizeof(int));
	if (*basis == NULL) {
		lp_simplex_free(*table);
		return lp_simplex_EXIT_FAILURE;
	}
	*constypes = (int *)lp_simplex_malloc(m * sizeof(int));
	if (*constypes == NULL) {
		lp_simplex_free(*table);
		lp_simplex_free(*basis);
		return lp_simplex_EXIT_FAILURE;
	}
	return lp_simplex_EXIT_SUCCESS;
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
		const struct optm_LinearConstraint *constraints,
		int *constypes, const int m)
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
void simplex_tableau_fill_constraints(
		double *table, const int ldtable,
		const struct optm_LinearConstraint *constraints,
		const int nrow, const int ncol, const int m, const int n)
{
	int i, j;

	lp_simplex_memset(table, 0., nrow * ldtable * sizeof(double));

	for (i = 0; i < m; i++) {
		const struct optm_LinearConstraint *cons = constraints + i;
		int row = (i + 1) * ldtable;

		if (cons->rhs >= 0) {
			table[ncol - 1 + row] = cons->rhs;
			lp_simplex_memcpy(table + row, cons->coef, n * sizeof(double));
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
void simplex_tableau_compute_size(
		const struct optm_LinearConstraint *constraints,
		const int m, const int n, int *nrow, int *ncol)
{
	int i;
	*nrow = m + 1;
	*ncol = m + n + 1;

	for (i = 0; i < m; i++) {
		const struct optm_LinearConstraint *cons = constraints + i;

		if (optm_CONS_T_GE == cons->type && cons->rhs >= 0)
			(*ncol)++;
		if (optm_CONS_T_LE == cons->type && cons->rhs < 0)
			(*ncol)++;
	}
}


/* Add slack variables (GE, LE) to simplex table
 * Return the number of slack variables
 */
int simplex_tableau_add_slack(
		double *table, const int ldtable, const int *constypes,
		const int m, const int n)
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
int simplex_tableau_add_artificial(
		double *table, const int ldtable, const int *constypes,
		const int m, const int n, const int nslack)
{
	int i, nartif = 0;

	for (i = 0; i < m; i++) {
		if (optm_CONS_T_LE != constypes[i]) {
			table[n + nslack + nartif + (i + 1) * ldtable] =  1.;
			table[n + nslack + nartif] = -1.;
			nartif++;
		}
	}
	return nartif;
}


/* Fill in the basis index set of artificial LP
 */
void simplex_tableau_initialize_basis(
		int *basis, const int *constypes,
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


void simplex_tableau_initialize_phase_one_objective(
		double *table, const int ldtable, const int *constypes,
		const int m, const int ncol)
{
	int i, rowi;

	for (i = 0; i < m; i++) {
		if (optm_CONS_T_LE == constypes[i])
			continue;
		rowi = (i + 1) * ldtable;
		lp_simplex_linalg_daxpy(ncol, 1, table + rowi, 1, table, 1);
	}
}


void simplex_tableau_remove_artificial_basis(
		double *table, int ldtable, int *basis,
		const int m, const int nreal, int nvar)
{
	int i, j, q = nvar;
	double ele, maxv;

	/* basis index should not exceeds number of variables */
	if (maxabs_arri(basis, m, 1) < nreal)
		return;
	for (i = 0; i < m; i++) {
		if (basis[i] <nreal)
			continue;
		maxv = 0.;
		q = nvar;
		for (j = 0; j < nreal; j++) {
			if (is_in_arri(j, basis, m))
				continue;
			ele = __lp_simplex_ABS__(table[j + (i + 1) * ldtable]);
			if (ele > maxv) {
				maxv = ele;
				q = j;
			}
		}
		if (maxv < 1e-9) {
			table[nvar + (i + 1) * ldtable] = 0.;
			lp_simplex_memset(table + (i + 1) * ldtable, 0, nreal * sizeof(double));
			continue;
		}
		simplex_apply_pivot(table, ldtable, m, nvar, i, q, 1, 1, 0);
		basis[i] = q;
	}
}


void simplex_tableau_remove_artificial_columns(
		double *table, const int ldtable,
		const int m, const int nreal, const int nartif)
{
	int i, rowi;

	if (nartif <= 0)  /* Delete artificial columns */
		return;
	for (i = 0; i < m + 1; i++) {
		rowi = i * ldtable;
		table[nreal + rowi] = table[nreal + nartif + rowi];
	}
}
