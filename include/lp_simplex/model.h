/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#ifndef LP_SIMPLEX_MODEL_H
#define LP_SIMPLEX_MODEL_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#define optm_VAR_T_REAL	0	/* variable type: real */
#define optm_VAR_T_INT	1	/* variable type: integer */
#define optm_VAR_T_BIN	2	/* variable type: binary */

#define optm_BOUND_T_FR	0	/* bound type: free */
#define optm_BOUND_T_UP	1	/* bound type: upper bounded */
#define optm_BOUND_T_LO	2	/* bound type: lower bounded */
#define optm_BOUND_T_BS	3	/* bound type: bounded from both sides */

#define optm_CONS_T_EQ	0	/* constraint type: equal to */
#define optm_CONS_T_GE	1	/* constraint type: greater or equal to */
#define optm_CONS_T_LE	2	/* constraint type: less or equal to */


struct optm_VariableBound {
	char name[9];
	double lb;
	double ub;
	int b_type;
	int v_type;
};

struct optm_LinearConstraint {
	char name[9];
	double *coef;
	double rhs;
	int type;
};

/**
 * Dense linear-programming input model.
 *
 * `coefficients` owns an m-by-n row-major matrix.  Every `constraints[i].coef`
 * points to row i of that matrix and must not be allocated or freed separately.
 */
struct lp_Model {
	int                           m;             /* number of constraints */
	int                           n;             /* number of variables */
	double                       *objective;
	double                       *coefficients;  /* row major */
	struct optm_LinearConstraint *constraints;
	struct optm_VariableBound    *bounds;
};


/** Allocate a zero-initialized dense LP model with nonnegative real variables. */
struct lp_Model *lp_model_create(int m, int n);

/** Release a model and all arrays owned by it. */
void lp_model_free(struct lp_Model *model);


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif
