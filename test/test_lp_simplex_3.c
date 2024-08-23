/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#include <lp_simplex/lp_simplex.h>
#include <stdio.h>
#include <assert.h>

/* LP Example: Degeneracy
 *
 *         max    -3/4 * x1 + 20 * x2 - 0.5 * x3 + 6 * x4
 *         s.t.   0.25 * x1 -  8 * x2 -       x3 + 9 * x4 <= 0
 *                0.5  * x1 - 12 * x2 - 0.5 * x3 + 3 * x4 <= 0
 *                                            x3          <= 0
 *                x, y >= 0
 */
#define m 3        /* number of constraints */
#define n 4        /* number of variables   */

double obj[] = {-3./4., 20., -0.5, 6.};
double constraint_1_coef[] = {0.25, -8., -1., 9.};
double constraint_2_coef[] = {0.5, -12., -0.5, 3.};
double constraint_3_coef[] = {0., 0., 1., 0.};

struct optm_LinearConstraint constraints[] = {
	{ "", constraint_1_coef,  0., lp_CONS_T_LE },
	{ "", constraint_2_coef,  0., lp_CONS_T_LE },
	{ "", constraint_3_coef,  1., lp_CONS_T_LE }
};

int main(void)
{
	/* call simplex subroutine that should be degenerated */
	double x[n], value;
	int code;
	int state = lp_simplex(obj, constraints, NULL, m, n, "bland", 1000, x, &value, &code);

	printf("Error code = %u\n", code);
	assert(state == lp_simplex_EXIT_SUCCESS);
	printf("value = %f\nSolution = ", value);
	printf("\n");
	return 0;
}
