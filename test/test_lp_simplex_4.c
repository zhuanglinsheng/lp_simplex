/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#include <lp_simplex/lp_simplex.h>
#include <stdio.h>
#include <assert.h>

/* LP Example
 * (from: https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.linprog.html)
 *
 *         min      - x0 + 4 * x1
 *         s.t.  -3 * x0 +     x1 <=  6
 *                   -x0 - 2 * x1 >= -4
 *                             x1 >= -3
 *
 * The solution is (10, -3) and the optimal value is -22
 */
#define m 2        /* number of constraints */
#define n 2        /* number of variables   */

double obj[] = {-1., 4.};
double constraint_1_coef[] = {-3., 1.};
double constraint_2_coef[] = {-1., -2.};

struct lp_simplex_LinearConstraint constraints[] = {
	{ "", constraint_1_coef,  6., lp_simplex_CONS_T_LE },
	{ "", constraint_2_coef, -4., lp_simplex_CONS_T_GE },
};
struct lp_simplex_VariableBound bounds[] = {
	{ "x0", __lp_simplex_NINF__, __lp_simplex_INF__, lp_simplex_BOUND_T_FR, lp_simplex_VAR_T_REAL },
	{ "x1",            -3, __lp_simplex_INF__, lp_simplex_BOUND_T_LO, lp_simplex_VAR_T_REAL },
};

int main(void)
{
	/* call simplex subroutine */
	double x[n], value;
	int code;
	int state = lp_simplex_lp_simplex(obj, constraints, bounds, m, n, "bland", 1000, x, &value, &code);

	printf("Error code = %u\n", code);
	assert(state == lp_simplex_EXIT_SUCCESS);
	assert(__lp_simplex_ABS__(value + 22) < 1e-8);
	assert(__lp_simplex_ABS__(x[0] - 10.) < 1e-8);
	assert(__lp_simplex_ABS__(x[1] + 3.) < 1e-8);
	return 0;
}
