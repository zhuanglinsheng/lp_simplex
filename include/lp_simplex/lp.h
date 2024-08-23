/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#ifndef __LP_H__
#define __LP_H__

#ifdef __cpluscplus
extern "C" {
#endif /* __cplusplus */

#define lp_VAR_T_REAL	0
#define lp_VAR_T_INT	1
#define lp_VAR_T_BIN	2

#define lp_BOUND_T_FR	0	/* free */
#define lp_BOUND_T_UP	1	/* upper bounded */
#define lp_BOUND_T_LO	2	/* lower bounded */
#define lp_BOUND_T_BS	3	/* bounded from both sides */

#define lp_CONS_T_EQ	0	/* equal to */
#define lp_CONS_T_GE	1	/* greater or equal to */
#define lp_CONS_T_LE	2	/* less or equal to */


struct optm_VariableBound {
	char name[8];
	double lb;
	double ub;
	int b_type;
	int v_type;
};

struct optm_LinearConstraint {
	char name[8];
	double * coef;
	double rhs;
	int type;
};

/* Linear programming model */
struct lp_Model {
	int m;			/* number of constraints */
	int n;			/* number of variables */
	double *objective;
	double *coefficients;	/* row major */
	struct optm_LinearConstraint *constraints;
	struct optm_VariableBound *bounds;
};


#ifdef __cpluscplus
}
#endif /* __cpluscplus */

#endif
