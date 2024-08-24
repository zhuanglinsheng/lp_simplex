/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#include <lp_simplex/lp_simplex_utils.h>
#include <lp_simplex/lp_simplex.h>

/* Checker of the checking "LP is feasible" */
#define __lp_simplex_FEASIBLE__			1e-5

/* To create in heap (need to be released) simplex table, index set of basis
 * and constraint type recorder
 */
static int create_buffer(double **table, int **basis, int **constypes,
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

static void free_buffer(double *table, int *basis, int *constypes)
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
static void fill_constypes(const struct optm_LinearConstraint *constraints, int *constypes, const int m)
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
static void fill_conscoefs(double *table, const int ldtable, const struct optm_LinearConstraint *constraints,
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
static void table_size_usul(const struct optm_LinearConstraint *constraints,
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
static int add_slack(double *table, const int ldtable, const int *constypes, const int m, const int n)
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
static int add_artif(double *table, const int ldtable, const int *constypes,
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
static void fill_artiflp_basis(int *basis, const int *constypes,
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

static void fill_artiflp_nrcost(double *table, const int ldtable, const int *constypes,
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

static void transf_artif_basis(double *table, int ldtable, int *basis, const int m, const int nreal, int nvar)
{
	int i, j, q = nvar;
	double ele, maxv = __lp_simplex_NINF__;

	/* basis index should not exceeds number of variables */
	if (maxabs_arri(basis, m, 1) <= nreal)
		return;
	for (i = 0; i < m; i++) {
		if (basis[i] <nreal)
			continue;
		for (j = 0; j < nreal; j++) {
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
		lp_simplex_pivot_core(table, ldtable, m, nvar, i, q, 1, 1, 0);
		basis[i] = q;
		q = nvar;
	}
}

static void delete_artif_cols(double *table, const int ldtable, const int m, const int nreal, const int nartif)
{
	int i, rowi;

	if (nartif <= 0)  /* Delete artificial columns */
		return;
	for (i = 0; i < m + 1; i++) {
		rowi = i * ldtable;
		table[nreal + rowi] = table[nreal + nartif + rowi];
	}
}

/* Phase 1: get a BFS for the original problem using the usual way - artificial LP
 *
 * Work:
 * 	1. allocate memory for table, basis, constypes
 * 	2. form a basic feasible solution (BSF)
 * 	3. assign ldtable and nvar, the number of vars in BSF
 */
static int simplex_phase_1_usul(double **table, int *ldtable, int **basis, int **constypes,
				int *nvar, int *epoch, int *code,
				const struct optm_LinearConstraint *constraints,
				const int m, const int n, const char *criteria, const int niter)
{
	int nrow, ncol;
	int nslack, nartif;

	table_size_usul(constraints, m, n, &nrow, &ncol);
	*ldtable = ncol;  /* leading dimension of table in memory */
	if (create_buffer(table, basis, constypes, m, nrow, *ldtable) == lp_simplex_EXIT_FAILURE) {
		*code = lp_simplex_MemoryAllocError;
		return lp_simplex_EXIT_FAILURE;
	}

	fill_constypes(constraints, *constypes, m);
	fill_conscoefs(*table, *ldtable, constraints, nrow, ncol, m, n);
	nslack = add_slack(*table, *ldtable, *constypes, m, n);
	nartif = add_artif(*table, *ldtable, *constypes, m, n, nslack);
	*nvar = n + nslack + nartif;  /* will be recovered to `n + nslack` upon success */
	if (m > (*nvar)) {
		*code = lp_simplex_OverDetermination;
		goto END;
	}
	fill_artiflp_basis(*basis, *constypes, m, n, nslack);
	fill_artiflp_nrcost(*table, *ldtable, *constypes, m, ncol);

	switch (lp_simplex_bsc(epoch, *table, *ldtable, *basis, m, *nvar, n + nslack, criteria, niter)) {
	case 0:
		*code = lp_simplex_ExceedIterLimit;
		goto END;
	case 1:
		if ((*table)[ncol - 1] > __lp_simplex_FEASIBLE__) {
			*code = lp_simplex_Infeasibility;
			goto END;
		}
		transf_artif_basis(*table, *ldtable, *basis, m, n + nslack, *nvar);
		delete_artif_cols(*table, *ldtable, m, n + nslack, nartif);
		*nvar = n + nslack;
		return lp_simplex_EXIT_SUCCESS;
	case 2:
		*code = lp_simplex_Unboundedness;
		goto END;
	case 3:
		*code = lp_simplex_Degeneracy;
		goto END;
	case 9:
		*code = lp_simplex_PrecisionError;
		goto END;
	}
END:
	free_buffer(*table, *basis, *constypes);
	return lp_simplex_EXIT_FAILURE;
}

/* Phase 2: Solve the original problem
 */
static int simplex_phase_2_usul(double *table, int ldtable, int *basis, int *constypes,
				int *epoch, int *code, const int m, const int n,
				const int nvar, const char *criteria, const int niter)
{
	switch (lp_simplex_bsc(epoch, table, ldtable, basis, m, nvar, nvar, criteria, niter)) {
	case 0:
		*code = lp_simplex_ExceedIterLimit;
		goto END;
	case 1:
		*code = lp_simplex_Success;
		return lp_simplex_EXIT_SUCCESS;
	case 2:
		*code = lp_simplex_Unboundedness;
		goto END;
	case 3:
		*code = lp_simplex_Degeneracy;
		goto END;
	case 9:
		*code = lp_simplex_PrecisionError;
		goto END;
	}
END:
	free_buffer(table, basis, constypes);
	return lp_simplex_EXIT_FAILURE;  /* error code already updated */
}

int lp_simplex_std(const double *objective, const struct optm_LinearConstraint *constraints,
			const int m, const int n, const char *criteria, const int niter,
			double *x, double *value, int *code)
{
	int i, j;
	int ldtable;
	int nvar;
	int epoch = 0;
	int *basis = NULL;
	double *table = NULL;
	int *constypes = NULL;

	assert(objective != NULL);
	assert(constraints != NULL);
	assert(x != NULL);
	assert(value != NULL);
	assert(code != NULL);

	if (simplex_phase_1_usul(&table, &ldtable, &basis, &constypes, &nvar, &epoch, code,
				 constraints, m, n, criteria, niter) == lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;

	for (j = 0; j < n; j++)  /* Fill in original objective coefficients */
		table[j] = -objective[j];
	for (i = 0; i < m; i++) {  /* row_0 = row_0 - ratio * row_{i+1} */
		int rowi = (i + 1) * ldtable;
		double ratio = -table[basis[i]];

		lp_simplex_linalg_daxpy(nvar + 1, ratio, table + rowi, 1, table, 1);
	}
	if (simplex_phase_2_usul(table, ldtable, basis, constypes, &epoch, code,
				 m, n, nvar, criteria, niter) == lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;

	*value = table[nvar];
	lp_simplex_memset(x, 0., n * sizeof(double));
	for (i = 0; i < m; i++) {
		if (basis[i] < n)
			x[basis[i]] = table[nvar + (i + 1) * ldtable];
	}
	free_buffer(table, basis, constypes);
	return lp_simplex_EXIT_SUCCESS;
}
