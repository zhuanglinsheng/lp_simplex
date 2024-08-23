/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

#ifndef __IMPF_LP_H__
#define __IMPF_LP_H__

#ifdef __cpluscplus
extern "C" {
#endif /* __cplusplus */

/* Linear programming model */
struct impf_Model_LP {
	int m;			/* number of constraints */
	int n;			/* number of variables */
	double *objective;
	double *coefficients;	/* row major */
	struct impf_LinearConstraint *constraints;
	struct impf_VariableBound *bounds;
};

/* Importing MPS file and get a `model`
 *
 * Note:
 * 	1. only "strict" MPS format is recognized by this function, which is an
 *		old format with a line of at most 61 columns. See:
 *		https://lpsolve.sourceforge.net/5.5/mps-format.htm
 *	2. the return of this function should be released by `impf_lp_free`
 *	3. return `NULL` on failure
 */
struct impf_Model_LP* lp_readmps(const char *file);

/* Release the LP model */
void impf_lp_free(struct impf_Model_LP *model);

/*******************************************************************************
 * Optimization of "lp-simplex-family"
 ******************************************************************************/

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
int impf_lp_simplex(const double *objective, const struct impf_LinearConstraint *constraints,
		    const struct impf_VariableBound *bounds,
		    const int m, const int n, const char *criteria, const int niter,
		    double *x, double *value, int *code);

/* Simplex algorithm for solving LP of general form
 * (Wrapper of `impf_fmin_lp_simplex_full` by taking `impf_Model_LP` as input)
 *
 * Note:
 * 	1. to use this method, users can either create `impf_Model_LP` manually
 *		or get a model from `impf_lp_readmps`
	2. to manually create model, TAKE CARE of the inner relation of struct
		`impf_Model_LP`, where you should map coefs to constraints
 *
 * Return: `EXIT_SUCCESS` or `EXIT_FAILURE`
 */
int impf_lp_simplex_wrp(const struct impf_Model_LP *model, const char *criteria, const int niter,
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
int impf_lp_simplex_std(const double *objective, const struct impf_LinearConstraint *constraints,
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

/*******************************************************************************
 * Optimization of "lp-interior-point-family"
 ******************************************************************************/

#ifdef __cpluscplus
}
#endif /* __cpluscplus */

#endif
