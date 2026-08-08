/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#include "simplex_pivot.h"
#include "linalg.h"
#include "utils.h"


/* Identifier of non-zero beta */
#define __lp_simplex_ZEROS_BETA__           1e-9

/* Checker of the general checking "LP is optimal" */
#define __lp_simplex_CTR_SPLX_OPTIMAL__     1e-9

/* Numerically stable staged threshold used by Bland's entering rule. */
#define __lp_simplex_BLAND_EPS__            1e-3

/* Controller for pivot leaving rule */
#define __lp_simplex_PIV_LEV__              1e-9


/* Choose the variable to leave basis
 * Return the index of the variable and check weather LP is "bounded"
 */
static int simplex_pivot_leave_rule(
		const double *table, const int ldtable,
		const int *basis, const int m, const int n,
		const int q, int *bounded)
{
	int i, p = n;
	double y_i_0, y_i_q, x_iq, min_x_iq = __lp_simplex_INF__;
	*bounded = 0;

	for (i = 0; i < m; i++) {
		y_i_0 = table[n + (i + 1) * ldtable];
		y_i_q = table[q + (i + 1) * ldtable];

		if (y_i_q <= __lp_simplex_PIV_LEV__)
			continue;
		else {
			x_iq = y_i_0 / y_i_q;

			if (x_iq < min_x_iq - __lp_simplex_PIV_LEV__ ||
			    (__lp_simplex_ABS__(x_iq - min_x_iq) <= __lp_simplex_PIV_LEV__ &&
			     (p == n || basis[i] < basis[p]))) {
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
static int simplex_pivot_enter_rule_datzig(
		const double *table, const unsigned char *is_basic, const int n)
{
	int j, q = n;
	double beta_j = 0.;
	double beta_q = 0.;

	for (j = 0; j < n; j++) {
		if (is_basic[j])
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
static int simplex_pivot_enter_rule_bland(
		const double *table, const unsigned char *is_basic, const int n)
{
	int j;
	double epsilon = __lp_simplex_BLAND_EPS__;

	for (;;) {
		for (j = 0; j < n; j++) {
			if (!is_basic[j] && table[j] > epsilon)
				return j;
		}
		if (epsilon < __lp_simplex_CTR_SPLX_OPTIMAL__)
			return n;
		epsilon /= 10.;
	}
}


/* Key subroutine of pivoting given p and q
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
void simplex_apply_pivot(
		double *table, const int ldtable,
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
			if (rto == 0.)
				continue;
			lp_simplex_linalg_daxpy(ncol, rto, table + rowp, 1, table + rowi, 1);
		}
	}
	if (rule3)
		lp_simplex_linalg_daxpy(ncol, -table[q], table + rowp, 1, table, 1);
}


/* Pivot starting from a basic representation for one round
 *
 * Return:
 *	0: current BFS is NOT optimal (stop before converged)
 *	1: current BSF is optimal
 *	2: LP is unbounded
 *	9: numerical precision error
 */
static int lp_simplex_pivot_on(
		double *table, const int ldtable, int *basis,
		unsigned char *is_basic, const int m, const int n,
		const int use_dantzig)
{
	int bounded = 0;
	int q = 0, p = 0;
	int i, residual = 0;

	/* Rebuild from the authoritative basis on every iteration.  This is O(m+n),
	 * but makes each entering-column scan O(n) instead of O(m*n), while also
	 * tolerating phase-I basis rewrites and duplicate entries exactly. */
	lp_simplex_memset(is_basic, 0, (size_t)n * sizeof(unsigned char));
	for (i = 0; i < m; i++) {
		if (basis[i] < 0)
			return 9;
		/* A redundant phase-I row may deliberately retain an artificial
		 * column index after those columns have been removed from phase II. */
		if (basis[i] < n)
			is_basic[basis[i]] = 1;
	}

	if (use_dantzig) {
		q = simplex_pivot_enter_rule_datzig(table, is_basic, n);
	}
	else {
		q = simplex_pivot_enter_rule_bland(table, is_basic, n);
	}
	if (n <= q) {
		/* No eligible nonbasic column remains.  A sizeable positive basic
		 * residual signals a damaged tableau rather than a valid optimum. */
		for (i = 0; i < n; i++) {
			if (table[i] > __lp_simplex_CTR_SPLX_OPTIMAL__) {
				residual = 1;
				break;
			}
		}
		return residual ? 9 : 1;
	}
	p = simplex_pivot_leave_rule(table, ldtable, basis, m, n, q, &bounded);
	if (bounded == 0)
		return 2;
	basis[p] = q;
	simplex_apply_pivot(table, ldtable, m, n, p, q, 1, 1, 1);
	return 0;
}


int simplex_run_pivots(
		int *epoch, double *table, const int ldtable, int *basis,
		const int m, const int n, const int nreal,
		const char *criteria, const int niter)
{
	unsigned char *is_basic;
	int result = 0;
	int use_dantzig;

	(void)nreal;
	assert(table != NULL);
	assert(basis != NULL);
	assert(epoch != NULL);
	if (criteria == NULL)
		criteria = "";
	use_dantzig = 7 == lp_simplex_strlen(criteria) &&
		0 == lp_simplex_memcmp("dantzig", criteria, 7);
	is_basic = (unsigned char *)lp_simplex_malloc((size_t)n * sizeof(unsigned char));
	if (is_basic == NULL)
		return 4;

	while (*epoch < niter) {
		(*epoch)++;
		switch (lp_simplex_pivot_on(table, ldtable, basis, is_basic,
					   m, n, use_dantzig)) {
		case 0:
			break;
		case 1:
			result = 1;
			goto END;
		case 2:
			result = 2;
			goto END;
		case 9:
			result = 9;
			goto END;
		}
	}
END:
	lp_simplex_free(is_basic);
	return result;
}
