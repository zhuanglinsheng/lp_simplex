/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#include <lp_simplex/lp_simplex.h>
#include <stdio.h>
#include <assert.h>

/* LP Example
 * (from: https://online-optimizer.appspot.com/?model=builtin:default.mod)
 *
 *         max    3 * x1 + 2 * x2
 *         s.t.       x1 + x2 <=  9
 *                3 * x1 - x2 <= 18
 *                    x1      <=  7
 *                         x2 <=  6
 *                x1, x2 >= 0
 *
 * The solution is (4.5, 4.5) and the optimal value is 22.5
 */
#define m 4        /* number of constraints */
#define n 2        /* number of variables   */

double obj[] = {-3., -2.};        /* transform "max" into "min" */
double constraint_1_coef[] = {1., 1.};
double constraint_2_coef[] = {3., 1.};
double constraint_3_coef[] = {1., 0.};
double constraint_4_coef[] = {0., 1.};

struct lp_simplex_LinearConstraint constraints[] = {
	{ "", constraint_1_coef,  9., lp_simplex_CONS_T_LE },
	{ "", constraint_2_coef, 18., lp_simplex_CONS_T_LE },
	{ "", constraint_3_coef,  7., lp_simplex_CONS_T_LE },
	{ "", constraint_4_coef,  6., lp_simplex_CONS_T_LE }
};

int main(void)
{
	/* call simplex subroutine */
	double x[n], value;
	int code;
	int state = lp_simplex_lp_simplex(obj, constraints, NULL, m, n, "bland", 1000, x, &value, &code);

	printf("Error code = %u\n", code);
	assert(state == lp_simplex_EXIT_SUCCESS);
	assert(__lp_simplex_ABS__(value + 22.5) < 1e-8);
	assert(__lp_simplex_ABS__(x[0] - 4.5) < 1e-8);
	assert(__lp_simplex_ABS__(x[1] - 4.5) < 1e-8);
	return 0;
}
