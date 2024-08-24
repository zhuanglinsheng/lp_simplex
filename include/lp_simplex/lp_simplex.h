/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#ifndef __lp_simplex_LP_H__
#define __lp_simplex_LP_H__

#include "lp.h"
#include "lp_simplex_utils.h"

#ifdef __cpluscplus
extern "C" {
#endif /* __cplusplus */

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

#define lp_simplex_EXIT_FAILURE			-1
#define lp_simplex_EXIT_SUCCESS			0


/* Importing MPS file and get a `model`
 *
 * Note:
 * 	1. only "strict" MPS format is recognized by this function, which is an
 *		old format with a line of at most 61 columns. See:
 *		https://lpsolve.sourceforge.net/5.5/mps-format.htm
 *	2. the return of this function should be released by `lp_simplex_free`
 *	3. return `NULL` on failure
 */
struct lp_Model* lp_simplex_read_mps(const char *file);

/* Release the LP model
 */
void lp_simplex_model_free(struct lp_Model *model);

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
int lp_simplex(const double *objective, const struct optm_LinearConstraint *constraints,
		const struct optm_VariableBound *bounds,
		const int m, const int n, const char *criteria, const int niter,
		double *x, double *value, int *code);

/* Simplex algorithm for solving LP of general form
 * (Wrapper of `lp_simplex_fmin_lp_simplex_full` by taking `lp_simplex_Model_LP` as input)
 *
 * Note:
 * 	1. to use this method, users can either create `lp_simplex_Model_LP` manually
 *		or get a model from `lp_simplex_lp_readmps`
 *	2. to manually create model, TAKE CARE of the inner relation of struct
 *		`lp_simplex_Model_LP`, where you should map coefs to constraints
 *
 * Return: `EXIT_SUCCESS` or `EXIT_FAILURE`
 */
int lp_simplex_wrp(const struct lp_Model *model, const char *criteria, const int niter,
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
int lp_simplex_std(const double *objective, const struct optm_LinearConstraint *constraints,
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
int lp_simplex_bsc(int *epoch, double *table, const int ldtable, int *basis,
			const int m, const int n, const int nreal,
			const char *criteria, const int niter);

/* Key subroutine of pivoting
 *
 * Parameter:
 *	p	idx of variable to leave basis
 *	q	idx of variable to enter basis
 *
 * Work:
 *	rule 1. row_p normalized by dividing y_p_q
 *	rule 2. row_i -= row_p * y_i_q
 *	rule 3. row_0 -= row_p * beta_q
 */
void lp_simplex_pivot_core(double *table, const int ldtable,
			const int m, const int n, const int p, const int q,
			const int rule1, const int rule2, const int rule3);

#ifdef __cpluscplus
}
#endif /* __cpluscplus */

#endif
