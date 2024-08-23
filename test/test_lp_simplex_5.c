/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#include <lp_simplex/lp_simplex.h>
#include <stdio.h>
#include <assert.h>

/* LP Example
 * (from: https://www.mathworks.com/help/optim/ug/linprog.html)
 *
 *         min   -1 * x - 1/3 *  y
 *         s.t.       x +        y <=  2
 *                    x + 0.25 * y <=  1
 *                    x -        y <=  2
 *            -0.25 * x -        y <=  1
 *                   -x -        y <= -1
 *                   -x +        y <=  2
 *
 * The solution is (2/3, 4/3) and the value is -10/9
 */
#define m 6        /* number of constraints */
#define n 2        /* number of variables   */

double obj[] = {-1., -1./3., 0.5};
double constraint_1_coef[] = {  1.,  1.   };
double constraint_2_coef[] = {  1.,  0.25 };
double constraint_3_coef[] = {  1., -1.   };
double constraint_4_coef[] = { -0.25, -1.   };
double constraint_5_coef[] = { -1., -1.   };
double constraint_6_coef[] = { -1.,  1.   };

struct lp_simplex_LinearConstraint constraints[] = {
	{ "", constraint_1_coef,  2., lp_simplex_CONS_T_LE },
	{ "", constraint_2_coef,  1., lp_simplex_CONS_T_LE },
	{ "", constraint_3_coef,  2., lp_simplex_CONS_T_LE },
	{ "", constraint_4_coef,  1., lp_simplex_CONS_T_LE },
	{ "", constraint_5_coef, -1., lp_simplex_CONS_T_LE },
	{ "", constraint_6_coef,  2., lp_simplex_CONS_T_LE }
};
struct lp_simplex_VariableBound bounds[] = {
	{ "x", __lp_simplex_NINF__, __lp_simplex_INF__, lp_simplex_BOUND_T_FR, lp_simplex_VAR_T_REAL },
	{ "y", __lp_simplex_NINF__, __lp_simplex_INF__, lp_simplex_BOUND_T_FR, lp_simplex_VAR_T_REAL }
};

int main(void)
{
	/* call simplex subroutine */
	double x[n], value;
	int code;
	int state = lp_simplex(obj, constraints, bounds, m, n, "bland", 1000, x, &value, &code);

	printf("error = %i\n", code);
	assert(state == lp_simplex_EXIT_SUCCESS);
	assert(__lp_simplex_ABS__(value + 10. / 9.) < 1e-8);
	assert(__lp_simplex_ABS__(x[0] - 2./3.) < 1e-8);
	assert(__lp_simplex_ABS__(x[1] - 4./3.) < 1e-8);
	return 0;
}
