/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

#ifndef __lp_simplex_LP_H__
#define __lp_simplex_LP_H__

#ifdef __cpluscplus
extern "C" {
#endif /* __cplusplus */

#ifdef USE_BLAS
extern void dscal_(int *n, double *alpha, double *x, int *incx);
extern void daxpy_(int *n, double *alpha, double *x, int *incx, double *y, int *incy);
#endif

#define __lp_simplex_ABS__(x) ((x) >= 0 ? (x) : (-(x)))
#define __lp_simplex_MAX__(x, y) ((x) >= (y) ? (x) : (y))
#define __lp_simplex_MIN__(x, y) ((x) <= (y) ? (x) : (y))

#define __lp_simplex_INF__ (1. / 0.)
#define __lp_simplex_NINF__ (-1. / 0.)

#define lp_simplex_Success			0
#define lp_simplex_MemoryAllocError		1
#define lp_simplex_CondUnsatisfied		2
#define lp_simplex_ExceedIterLimit		3
#define lp_simplex_Singularity			4
#define lp_simplex_OverDetermination		5
#define lp_simplex_Unboundedness		6
#define lp_simplex_Infeasibility		7
#define lp_simplex_Degeneracy			8
#define lp_simplex_PrecisionError		9

#define lp_simplex_EXIT_FAILURE -1
#define lp_simplex_EXIT_SUCCESS 0

#define lp_simplex_VAR_T_REAL	0
#define lp_simplex_VAR_T_INT	1
#define lp_simplex_VAR_T_BIN	2

#define lp_simplex_BOUND_T_FR	0	/* free */
#define lp_simplex_BOUND_T_UP	1	/* upper bounded */
#define lp_simplex_BOUND_T_LO	2	/* lower bounded */
#define lp_simplex_BOUND_T_BS	3	/* bounded from both sides */

#define lp_simplex_CONS_T_EQ	0
#define lp_simplex_CONS_T_GE	1
#define lp_simplex_CONS_T_LE	2

struct lp_simplex_VariableBound {
	char name[16];
	double lb;
	double ub;
	int b_type;
	int v_type;
};

struct lp_simplex_LinearConstraint {
	char name[16];
	double * coef;
	double rhs;
	int type;
};

/* Linear programming model */
struct lp_simplex_Model_LP {
	int m;			/* number of constraints */
	int n;			/* number of variables */
	double *objective;
	double *coefficients;	/* row major */
	struct lp_simplex_LinearConstraint *constraints;
	struct lp_simplex_VariableBound *bounds;
};

/* Importing MPS file and get a `model`
 *
 * Note:
 * 	1. only "strict" MPS format is recognized by this function, which is an
 *		old format with a line of at most 61 columns. See:
 *		https://lpsolve.sourceforge.net/5.5/mps-format.htm
 *	2. the return of this function should be released by `lp_simplex_lp_free`
 *	3. return `NULL` on failure
 */
struct lp_simplex_Model_LP* lp_simplex_lp_readmps(const char *file);

/* Release the LP model */
void lp_simplex_lp_free(struct lp_simplex_Model_LP *model);

/* Simplex algorithm for solving LP of general form
 *
 *	min  c'x
 *	s.t. Ai x =[, >=, <=] bi, i = 1, ..., m
 *		lb <= x <= ub
 *
 * where x is n-dimensional vector and b is m-dimensional vector
 *
 * Parameters
 *	objective	coefficients of objective function (length = n)
 *	constraints	linear constraint array (length = m)
 *	bounds		could be either:
 *				0) `NULL` pointer indicating "x >= 0"
 *				1) variable bound array (length = n)
 *	m		number of linear constraints
 *	n		number of variables
 *	criteria	pivot criteria, including:
 *				0) ""		default
 *				1) "dantzig"	Dantzig's original rule
 *				2) "bland"	Bland's rule
 *				3) "pan97"	Pan (1997)
 *	niter		iteration limit
 *	x		array of solutions (length = n)
 *	value		optimal value of the objective
 *	code		error code (see basic.h)
 *
 * Return: `EXIT_SUCCESS` or `EXIT_FAILURE`
 */
int lp_simplex_lp_simplex(const double *objective, const struct lp_simplex_LinearConstraint *constraints,
			const struct lp_simplex_VariableBound *bounds,
			const int m, const int n, const char *criteria, const int niter,
			double *x, double *value, int *code);

/* Simplex algorithm for solving LP of general form
 * (Wrapper of `lp_simplex_fmin_lp_simplex_full` by taking `lp_simplex_Model_LP` as input)
 *
 * Note:
 * 	1. to use this method, users can either create `lp_simplex_Model_LP` manually
 *		or get a model from `lp_simplex_lp_readmps`
	2. to manually create model, TAKE CARE of the inner relation of struct
		`lp_simplex_Model_LP`, where you should map coefs to constraints
 *
 * Return: `EXIT_SUCCESS` or `EXIT_FAILURE`
 */
int lp_simplex_lp_simplex_wrp(const struct lp_simplex_Model_LP *model, const char *criteria, const int niter,
			double *x, double *value, int *code);

/* Simplex algorithm for solving LP of standard form
 *
 *	min  c'x
 *	s.t. Ai x =(, >=, <=) bi, i = 1, ..., m
 *		x >= 0
 *
 * where x is n-dimensional vector and b is m-dimensional vector
 *
 * Parameters:
 *	objective	coefficients of objective function (length = n)
 *	constraints	linear constraint array (length = m)
 *	m		number of linear constraints
 *	n		number of variables
 *	criteria	pivot criteria, including:
 *				0) ""		default
 *				1) "dantzig"	Dantzig's original rule
 *				2) "bland"	Bland's rule
 *				3) "pan97"	Pan (1997)
 *	x		array of solutions (length = n)
 *	value		optimal value of the objective
 *	code		error code (see basic.h)
 *	niter		iteration limit
 *
 * Return: `EXIT_SUCCESS` or `EXIT_FAILURE`
 */
int lp_simplex_lp_simplex_std(const double *objective, const struct lp_simplex_LinearConstraint *constraints,
			const int m, const int n, const char *criteria, const int niter,
			double *x, double *value, int *code);

/* Simplex algorithm for solving LP of basic representation
 *
 * Return
 *	0: current BFS is NOT optimal (stop before converged)
 *	1: current BSF is optimal
 *	2: LP is unbounded
 *	3: LP is circled more than accepted times (indicating for degeneracy)
 *	9: numerical precision error
 */
int simplex_pivot_bsc(int *epoch, double *table, const int ldtable, int *basis,
			const int m, const int n, const int nreal,
			const char *criteria, const int niter);

void simplex_pivot_core(double *table, const int ldtable,
			const int m, const int n, const int p, const int q,
			const int rule1, const int rule2, const int rule3);

#ifdef __cpluscplus
}
#endif /* __cpluscplus */

#endif
