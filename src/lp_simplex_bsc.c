/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

#include <lp_simplex/lp_simplex_magics.h>
#include <lp_simplex/lp_simplex_utils.h>
#include <lp_simplex/lp_simplex.h>
#include <assert.h>
#include <stddef.h>

/* Check simplex optimality: all zero row coefficients are non-positive
 * Return:
 * 	1 if is optimal
 * 	0 if is not optimal
 */
static int is_simplex_optimal(const double *table, const int n)
{
	int j;

	for (j = 0; j < n; j++) {
		if (table[j] > __lp_simplex_CTR_SPLX_OPTIMAL__)
			return 0;
	}
	return 1;
}

/* In simplex iteration, check whether value is improved.
 *
 * Return:
 *     0    not circled
 *     1    circled but optimal already
 *     2    circled
 */
static int check_simplex_degenerated(const double *table, const int n, const double old_value)
{
	if (old_value <= table[n] + __lp_simplex_CHC_SPLX_DEGENERATED__) {
		if (is_simplex_optimal(table, n))
			return 1;
		return 2;
	} else
		return 0;
}

/* Choose the variable to leave basis
 * Return the index of the variable and check weather LP is "bounded"
 */
static int simplex_pivot_leave_rule(const double *table, const int ldtable,
				    const int m, const int n, const int q, int *bounded)
{
	int i, p = n;
	double y_i_0, y_i_q, x_iq, min_x_iq = __lp_simplex_INF__;
	*bounded = 0;

	for (i = 0; i < m; i++) {
		y_i_0 = table[n + (i + 1) * ldtable];
		y_i_q = table[q + (i + 1) * ldtable];

		if (y_i_q <= __lp_simplex_CTR_SPLX_PIV_LEV__)
			continue;
		else {
			x_iq = y_i_0 / y_i_q;

			if (x_iq < min_x_iq) {
				min_x_iq = x_iq;
				p = i;
			}
			*bounded = 1;
		}
	}
	return p;
}

/* Fast pivot rule: choosing the variable to enter basis
 * Return the index of the variable (< n)
 *
 * Note:
 * 	On failure, the algorithm returns `n`. Logically, this NEVER happens,
 *	but numerically, there are many criteriors reporting optimality,
 *	leading to unpredicted results
 */
static int simplex_pivot_enter_rule_datzig(const double *table, const int *basis, const int m, const int n)
{
	int j, q = n;
	double beta_j = 0.;
	double beta_q = 0.;

	for (j = 0; j < n; j++) {
		if (is_in_arri(j, basis, m))
			continue;
		beta_j = table[j];

		if (beta_j > beta_q) {
			q = j;
			beta_q = beta_j;
		}
	}
	return q; /* return n if no p is found (optimal already) */
}

/* Bland's rule: choosing the variable to enter basis
 * Return the index of the variable (< n)
 *
 * Note: on failure, the algorithm returns n
 */
static int simplex_pivot_enter_rule_bland(const double *table, const int *basis, const int m, const int n)
{
	int j;
	double epsilon = __lp_simplex_CTR_SPLX_BLAND_EPS__;
BLAND_BEGIN:
	for (j = 0; j < n; j++) {
		if (!is_in_arri(j, basis, m)) {
			if (table[j] > epsilon)
				return j;
		}
	}
	if (epsilon >= __lp_simplex_CTR_SPLX_BLAND_EPS_MIN__) {
		epsilon /= 10.;
		goto BLAND_BEGIN;
	} else
		return n;
}

/* Key subroutine of pivoting
 *
 * Parameter:
 *	p	idx of variable to leave basis
 *	q	idx of variable to enter basis
 *
 * Work:
 *	rule 1. row_p normalized by dividing y_p_q
 *	rule 2. row_i -= row_p * y_i_q
 *	rule 3. row_0 -= row_p * beta_q
 */
void simplex_pivot_core(double *table, const int ldtable,
			const int m, const int n, const int p, const int q,
			const int rule1, const int rule2, const int rule3)
{
	int i, ncol = n + 1, rowp = (p + 1) * ldtable;
	double y_p_q = table[q + rowp];

	if (rule1)
		lp_simplex_linalg_dscal(ncol, 1 / y_p_q, table + rowp, 1);
	if (rule2) {
		for (i = 0; i < m; i++) {
			int rowi = (i + 1) * ldtable;
			double rto =  -table[q + rowi];

			if (i == p)
				continue;
			lp_simplex_linalg_daxpy(ncol, rto, table + rowp, 1, table + rowi, 1);
		}
	}
	if (rule3)
		lp_simplex_linalg_daxpy(ncol, -table[q], table + rowp, 1, table, 1);
}


#define __PAN_97_INVALID_BASIS_CRIT 1e-12

void simplex_pan97_trsf(const double *table, const int ldtable, const int *basis,
			const int m, const int n, const int p, const int q)
{
	int *non_basis_1 = (int *)lp_simplex_malloc(sizeof(int) * m);

	/* H = I - tau * u * u^T */
	double *vec_u = (double *)lp_simplex_malloc(sizeof(double) * m);
	double tau = 0, beta = 0;

	int i, j, k = 0;
	int m1 = m;

	/* calculate m1 and general basis, copy column "q" */
	for (i = 0; i < m; i++) {
		if (table[(i + 1) * (n + 1) + n] < __PAN_97_INVALID_BASIS_CRIT) {
			m1--;
			vec_u[k] = table[(i + 1) * (n + 1) + q];
			non_basis_1[k] = i;
			k++;
		}
	}
	if (m1 == m)
		goto END;
/*
	lp_simplex_linalg_dlarfg(m - m1, vec_u, vec_u + 1, 1, &tau);
	beta = vec_u[0];
	vec_u[0] = 1.0;
*/
END:
	lp_simplex_free(non_basis_1);
	lp_simplex_free(vec_u);
}


/* Pivot starting from a basic representation for one round
 *
 * Return:
 *	0: current BFS is NOT optimal (stop before converged)
 *	1: current BSF is optimal
 *	2: LP is unbounded
 *	9: numerical precision error
 */
static int simplex_pivot_on(double *table, const int ldtable, int *basis,
			    const int m, const int n, const char *criteria)
{
	int bounded;
	int q, p;

	if (is_simplex_optimal(table, n))
		return 1;
	if (7 == lp_simplex_strlen(criteria) && 0 == lp_simplex_memcmp("dantzig", criteria, 7)) {
		q = simplex_pivot_enter_rule_datzig(table, basis, m, n);
		p = simplex_pivot_leave_rule(table, ldtable, m, n, q, &bounded);
	}
	else if (5 == lp_simplex_strlen(criteria) && 0 == lp_simplex_memcmp("bland", criteria, 5)) {
		q = simplex_pivot_enter_rule_bland(table, basis, m, n);
		p = simplex_pivot_leave_rule(table, ldtable, m, n, q, &bounded);
	}
	else {  /* default method: "pan97" */
	/*
		q = simplex_pivot_enter_rule_datzig(table, basis, m, n);
		simplex_pan97_trsf(table, ldtable, basis, m, n, p, q);
		p = simplex_pivot_leave_rule(table, ldtable, m, n, q, &bounded);
	*/
	}
	if (n <= q) {
		return 9;
	}
	if (bounded == 0)
		return 2;
	basis[p] = q;
	simplex_pivot_core(table, ldtable, m, n, p, q, 1, 1, 1);
	return 0;
}

int simplex_pivot_bsc(int *epoch, double *table, const int ldtable, int *basis,
			const int m, const int n, const int nreal,
			const char *criteria, const int niter)
{
	double old_value = __lp_simplex_INF__;
	int degen_iter = 0;

	assert(table != NULL);
	assert(basis != NULL);
	assert(epoch != NULL);

	while (*epoch < niter) {
		(*epoch)++;
		switch (simplex_pivot_on(table, ldtable, basis, m, n, criteria)) {
		case 0:
			break;
		case 1:
			return 1;
		case 2:
			return 2;
		case 9:
			return 9;
		}
		if (check_simplex_degenerated(table, n, old_value) == 2) {
			degen_iter++;
			if (degen_iter > 5)
				return 3;
		} else
			degen_iter = 0;
		old_value = table[n];
	}
	return 0;
}
