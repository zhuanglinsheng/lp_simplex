/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#include <lp_simplex/lp_simplex.h>
#include <stdio.h>
#include <assert.h>

/* LP Example
 * (from: https://sma.epfl.ch/~niemeier/opt09/opt09_ch06.pdf)
 *
 *	min  x1 + x2 + x3
 *	s.t. x1 + 2 * x2 + 3 * x3	=  3
 *	-1 * x1 + 2 * x2 + 6 * x3	=  2
 *		 -4 * x2 - 9 * x3	= -5
 *			   3 * x3 + x4	=  1
 *	x1, x2, x3, x4 >= 0
 */
#define m 4	/* number of constraints */
#define n 4	/* number of variables   */

double objective[] = {1., 1., 1., 0.};
double constraint_1_coef[] = { 1.,  2.,  3., 0.};
double constraint_2_coef[] = {-1.,  2.,  6., 0.};
double constraint_3_coef[] = { 0., -4., -9., 0.};
double constraint_4_coef[] = { 0.,  0.,  3., 1.};

struct optm_LinearConstraint constraints[] = {
	{ "", constraint_1_coef,  3., opt_CONS_T_EQ },
	{ "", constraint_2_coef,  2., opt_CONS_T_EQ },
	{ "", constraint_3_coef, -5., opt_CONS_T_EQ },
	{ "", constraint_4_coef,  1., opt_CONS_T_EQ }
};

int main(void)
{
	/* call simplex subroutine */
	double x[n], value;
	int code;
	int state = lp_simplex(objective, constraints, NULL, m, n, "bland", 1000, x, &value, &code);

	printf("Error code = %u\n", code);
	printf("value = %e\n", value);
	assert(state == lp_simplex_EXIT_SUCCESS);
	assert(__lp_simplex_ABS__(value - 7. / 4.) < 1e-8);
	assert(__lp_simplex_ABS__(x[0] - 0.5) < 1e-8);
	assert(__lp_simplex_ABS__(x[1] - 1.25) < 1e-8);
	assert(__lp_simplex_ABS__(x[2]) < 1e-8);
	assert(__lp_simplex_ABS__(x[3] - 1.) < 1e-8);
	return 0;
}
