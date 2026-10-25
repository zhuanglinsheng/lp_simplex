/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_pivot.c
 * @brief Tableau entering and leaving variable selection and pivot operations.
 */

#include "simplex_pivot.h"
#include "linalg.h"
#include "utils.h"

#include <math.h>

/* Harris ratio test: bound the step by every row's allowed primal error,
 * then prefer a strong pivot inside that window. A positive coefficient
 * below the pivot floor still constrains the window; it is never dropped
 * as evidence of unboundedness. Return m when no safe pivot is available. */
static int simplex_pivot_leave_rule(const double *table, const int ldtable,
	const int *basis, const int m, const int n, const int q,
	const double pivot_tolerance, const double primal_tolerance,
	int *bounded)
{
	int i, p = m;
	double rhs, direction, step, lower = -primal_tolerance,
				     upper = HUGE_VAL;
	double strongest = 0.;
	*bounded = 0;
	for (i = 0; i < m; i++) {
		rhs = table[n + (i + 1) * ldtable];
		direction = table[q + (i + 1) * ldtable];
		if (!isfinite(rhs) || !isfinite(direction) ||
			rhs < -primal_tolerance)
			return m;
		if (direction > 0.) {
			step = (rhs + primal_tolerance) / direction;
			if (step < upper)
				upper = step;
		} else if (direction < 0.) {
			/* A tiny negative basic value can give a backward step. Its
			 * actual ratio must also respect the rows with negative direction. */
			step = (rhs + primal_tolerance) / direction;
			if (step > lower)
				lower = step;
		}
	}
	for (i = 0; i < m; i++) {
		rhs = table[n + (i + 1) * ldtable];
		direction = table[q + (i + 1) * ldtable];
		if (direction <= pivot_tolerance)
			continue;
		step = rhs / direction;
		if (!isfinite(step) || step < lower || step > upper)
			continue;
		if (direction > strongest ||
			(direction == strongest &&
				(p == m || basis[i] < basis[p]))) {
			strongest = direction;
			p = i;
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
static int simplex_pivot_enter_rule_datzig(const double *table,
	const unsigned char *is_basic, const int n,
	const double optimality_tolerance)
{
	int j, q = n;
	double beta_j = 0.;
	double beta_q = optimality_tolerance;

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
static int simplex_pivot_enter_rule_bland(const double *table,
	const unsigned char *is_basic, const int n,
	const double optimality_tolerance)
{
	int j;
	/* Bland selects the smallest eligible index at the requested tolerance.
	 * A staged absolute threshold changes that rule and depends on units. */
	for (j = 0; j < n; j++)
		if (!is_basic[j] && table[j] > optimality_tolerance)
			return j;
	return n;
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
void simplex_apply_pivot(double *table, const int ldtable, const int m,
	const int n, const int p, const int q, const int rule1, const int rule2,
	const int rule3)
{
	int i, ncol = n + 1, rowp = (p + 1) * ldtable;
	double y_p_q = table[q + rowp];

	if (rule1)
		lp_simplex_linalg_dscal(ncol, 1 / y_p_q, table + rowp, 1);
	if (rule2) {
		for (i = 0; i < m; i++) {
			int rowi = (i + 1) * ldtable;
			double rto = -table[q + rowi];

			if (i == p)
				continue;
			if (rto == 0.)
				continue;
			lp_simplex_linalg_daxpy(
				ncol, rto, table + rowp, 1, table + rowi, 1);
		}
	}
	if (rule3)
		lp_simplex_linalg_daxpy(
			ncol, -table[q], table + rowp, 1, table, 1);
	/* Keep the new unit column exact. Leaving cancellation residue in a
	 * basic column lets later degenerate pivots amplify that residue. */
	if (rule1 && isfinite(table[q + rowp]))
		table[q + rowp] = 1.;
	if (rule2 && table[q + rowp] == 1.)
		for (i = 0; i < m; i++)
			if (i != p && isfinite(table[q + (i + 1) * ldtable]))
				table[q + (i + 1) * ldtable] = 0.;
	if (rule3 && table[q + rowp] == 1. && isfinite(table[q]))
		table[q] = 0.;
}

/* Pivot starting from a basic representation for one round
 *
 * Return:
 *	0: current BFS is NOT optimal (stop before converged)
 *	1: current BSF is optimal
 *	2: LP is unbounded
 *	9: numerical precision error
 */
static int lp_simplex_pivot_on(double *table, const int ldtable, int *basis,
	unsigned char *is_basic, const int m, const int n,
	const int use_dantzig, const double optimality_tolerance,
	const double pivot_tolerance, const double primal_tolerance)
{
	int bounded = 0;
	int q = 0, p = 0;
	int i, residual = 0;

	/* The artificial objective can overflow while every constraint RHS and
	 * reduced cost remains finite. Its value cannot certify either phase. */
	if (!isfinite(table[n]))
		return 9;

	/* Rebuild from the authoritative basis on every iteration.  This is O(m+n),
	 * but makes each entering-column scan O(n) instead of O(m*n), while also
	 * tolerating phase-I basis rewrites and duplicate entries exactly. */
	lp_simplex_memset(is_basic, 0, (size_t)n * sizeof(unsigned char));
	for (i = 0; i < m; i++) {
		if (basis[i] < 0 || !isfinite(table[n + (i + 1) * ldtable]) ||
			table[n + (i + 1) * ldtable] < -primal_tolerance)
			return 9;
		/* A redundant phase-I row may deliberately retain an artificial
		 * column index after those columns have been removed from phase II. */
		if (basis[i] < n)
			is_basic[basis[i]] = 1;
	}

	for (i = 0; i < n; i++)
		if (!isfinite(table[i]))
			return 9;
	if (use_dantzig) {
		q = simplex_pivot_enter_rule_datzig(
			table, is_basic, n, optimality_tolerance);
	} else {
		q = simplex_pivot_enter_rule_bland(
			table, is_basic, n, optimality_tolerance);
	}
	if (n <= q) {
		/* No eligible nonbasic column remains.  A sizeable positive basic
		 * residual signals a damaged tableau rather than a valid optimum. */
		for (i = 0; i < n; i++) {
			if (table[i] > optimality_tolerance) {
				residual = 1;
				break;
			}
		}
		return residual ? 9 : 1;
	}
	p = simplex_pivot_leave_rule(table, ldtable, basis, m, n, q,
		pivot_tolerance, primal_tolerance, &bounded);
	if (bounded == 0) {
		/* A positive coefficient below the pivot threshold still bounds
		 * the ray. Refusing that pivot is numerical failure, not unboundedness.
		 * Phase II also requires a feasible finite starting basis. */
		for (i = 0; i < m; i++) {
			double direction = table[q + (i + 1) * ldtable];
			double value = table[n + (i + 1) * ldtable];
			if (!isfinite(direction) || !isfinite(value) ||
				direction > 0. || value < -primal_tolerance)
				return 9;
		}
		return 2;
	}
	basis[p] = q;
	simplex_apply_pivot(table, ldtable, m, n, p, q, 1, 1, 1);
	return 0;
}

int simplex_run_pivots(int *epoch, double *table, const int ldtable, int *basis,
	const int m, const int n, const int nreal,
	const struct lp_simplex_Options *options)
{
	unsigned char *is_basic;
	int result = 0;
	int use_dantzig;

	(void)nreal;
	assert(table != NULL);
	assert(basis != NULL);
	assert(epoch != NULL);
	use_dantzig = options->pricing == lp_simplex_PRICING_DANTZIG;
	is_basic = (unsigned char *)lp_simplex_malloc(
		(size_t)n * sizeof(unsigned char));
	if (is_basic == NULL)
		return 4;

	while (*epoch < options->iteration_limit) {
		(*epoch)++;
		switch (lp_simplex_pivot_on(table, ldtable, basis, is_basic, m,
			n, use_dantzig, options->dual_tolerance,
			options->pivot_tolerance, options->primal_tolerance)) {
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
