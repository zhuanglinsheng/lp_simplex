/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#include <lp_simplex/lp_simplex.h>
#include <stdio.h>
#include <assert.h>

/* LP Example
 * (from: https://developers.google.com/optimization/lp/lp_example)
 *
 *         max    3 * x + 4 * y
 *         s.t.       x + 2 * y <= 14
 *                3 * x -     y >=  0
 *                    x -     y <=  2
 *                x, y >= 0
 *
 * The solution is (6., 4.) and the optimal value is 34
 */
#define m 3        /* number of constraints */
#define n 2        /* number of variables   */

double obj[] = {-3., -4.};     /* transform "max" into "min" */
double constraint_1_coef[] = {1., 2.};
double constraint_2_coef[] = {3., -1.};
double constraint_3_coef[] = {1., -1.};

struct lp_simplex_LinearConstraint constraints[] = {
	{ "", constraint_1_coef, 14., lp_simplex_CONS_T_LE },
	{ "", constraint_2_coef,  0., lp_simplex_CONS_T_GE },
	{ "", constraint_3_coef,  2., lp_simplex_CONS_T_LE }
};

int main(void)
{
	/* call simplex subroutine */
	double x[n], value;
	int code;
	int state = lp_simplex_lp_simplex(obj, constraints, NULL, m, n, "bland", 1000, x, &value, &code);

	printf("Error code = %u\n", code);
	assert(state == lp_simplex_EXIT_SUCCESS);
	assert(__lp_simplex_ABS__(value + 34.) < 1e-8);
	assert(__lp_simplex_ABS__(x[0] - 6.) < 1e-8);
	assert(__lp_simplex_ABS__(x[1] - 4.) < 1e-8);
	return 0;
}
