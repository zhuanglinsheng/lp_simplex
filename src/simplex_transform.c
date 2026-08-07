/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#include "simplex_transform.h"
#include "utils.h"
#include <lp_simplex/status.h>


int simplex_transform_alloc(
		const int M, const int N,
		double **obj2, double **x2, double **coef2,
		struct optm_LinearConstraint **constraints2)
{
	*obj2 = NULL;
	*x2 = NULL;
	*coef2 = NULL;
	*constraints2 = NULL;

	*obj2 = lp_simplex_malloc(N * sizeof(double));
	if (*obj2 == NULL)
		return lp_simplex_EXIT_FAILURE;
	*x2 = lp_simplex_malloc(N * sizeof(double));
	if (*x2 == NULL) {
		lp_simplex_free(*obj2);
		return lp_simplex_EXIT_FAILURE;
	}
	*coef2 = lp_simplex_malloc(M * N * sizeof(double));
	if (*coef2 == NULL) {
		lp_simplex_free(*obj2);
		lp_simplex_free(*x2);
		return lp_simplex_EXIT_FAILURE;
	}
	*constraints2 = lp_simplex_malloc(M * sizeof(struct optm_LinearConstraint));
	if (*constraints2 == NULL) {
		lp_simplex_free(*obj2);
		lp_simplex_free(*x2);
		lp_simplex_free(*coef2);
		return lp_simplex_EXIT_FAILURE;
	}
	return lp_simplex_EXIT_SUCCESS;
}


void simplex_transform_free(
		double *obj2, double *x2, double *coef2,
		struct optm_LinearConstraint *constraints2)
{
	if (obj2)
		lp_simplex_free(obj2);
	if (x2)
		lp_simplex_free(x2);
	if (coef2)
		lp_simplex_free(coef2);
	if (constraints2)
		lp_simplex_free(constraints2);
}


/* Get the size of standard form LP.
 *
 * A free variable needs two nonnegative variables.  Only a variable bounded
 * on both sides needs an extra constraint: after x = lb + y, y <= ub - lb.
 * An upper-only variable is represented as x = ub - y and needs no extra
 * constraint.
 */
void simplex_transform_size(
		const struct optm_VariableBound *bounds,
		const int m, const int n, int *_M, int *_N)
{
	int j;
	*_M = m;
	*_N = n;

	for (j = 0; j < n; j++) {
		const struct optm_VariableBound *bd = bounds + j;

		if (optm_BOUND_T_FR == bd->b_type)
			(*_N)++;
		if (optm_BOUND_T_BS == bd->b_type)
			(*_M)++;
	}
}


/* Transform original LP into standard form
 *
 * Note:
 *	Original LP: allow for more variable bounds
 *	Standard LP: x >= 0
 */
void simplex_transform_problem(
		const double *objective,
		const struct optm_LinearConstraint *constraints,
		const struct optm_VariableBound *bounds,
		const int m, const int n, const int _N,
		double *obj2, double *obj_diff, double *coef2,
		struct optm_LinearConstraint *constraints2)
{
	int i, j;
	int ctr_var = 0;
	int ctr_ubcons = 0;

	for (i = 0; i < m; i++) {
		constraints2[i].coef = coef2 + i * _N;
		constraints2[i].rhs = (constraints + i)->rhs;
		constraints2[i].type = (constraints + i)->type;
	}
	for (j = 0; j < n; j++) {
		const struct optm_VariableBound *bd = bounds + j;
		double shift = 0.;
		double scale = 1.;

		if (optm_BOUND_T_UP == bd->b_type) {
			shift = bd->ub;
			scale = -1.;
		} else if (optm_BOUND_T_LO == bd->b_type || optm_BOUND_T_BS == bd->b_type) {
			shift = bd->lb;
		}

		obj2[ctr_var] = objective[j] * scale;
		*obj_diff += objective[j] * shift;
		for (i = 0; i < m; i++) {
			double aij = constraints[i].coef[j];
			coef2[ctr_var + i * _N] = aij * scale;
			constraints2[i].rhs -= aij * shift;
		}
		ctr_var++;

		if (optm_BOUND_T_FR == bd->b_type) {
			obj2[ctr_var] = -objective[j];
			for (i = 0; i < m; i++)
				coef2[ctr_var + i * _N] = -constraints[i].coef[j];
			ctr_var++;
		} else if (optm_BOUND_T_BS == bd->b_type) {
			int idx = m + ctr_ubcons;
			constraints2[idx].coef = coef2 + idx * _N;
			constraints2[idx].rhs = bd->ub - bd->lb;
			constraints2[idx].type = optm_CONS_T_LE;
			constraints2[idx].coef[ctr_var - 1] = 1.;
			ctr_ubcons++;
		}
	}
}


/* Recover original LP solution and value from the standard form
 */
void simplex_transform_recover(
		const struct optm_VariableBound *bounds, const int n,
		const double *x2, const double value2, const double obj_diff,
		double *x, double *value)
{
	int ctr_var;
	int j;

	ctr_var = 0;

	for (j = 0; j < n; j++) {
		const struct optm_VariableBound *bd = bounds + j;

		if (optm_BOUND_T_FR == bd->b_type) {
			x[j] = x2[ctr_var] - x2[ctr_var + 1];
			ctr_var++;
		} else if (optm_BOUND_T_LO == bd->b_type || optm_BOUND_T_BS == bd->b_type) {
			x[j] = x2[ctr_var] + bd->lb;
		} else if (optm_BOUND_T_UP == bd->b_type) {
			x[j] = bd->ub - x2[ctr_var];
		} else {
			x[j] = x2[ctr_var];
		}
		ctr_var++;
	}
	*value = value2 + obj_diff;
}
