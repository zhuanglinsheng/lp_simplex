/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_certificate.c
 * @brief Objective-independent separation of row and variable intervals.
 */

#include "simplex_problem.h"
#include "utils.h"

#include <math.h>
#include <float.h>

/* Project into the row-sign cone before recomputing A^T y. Projection
 * repairs multiplier roundoff; it is not permission to clip column products. */
static long double certificate_multiplier(const struct simplex_Problem *p,
	const double *y, int i, int sign, int nonnegative_equality)
{
	long double value = (long double)(sign * y[i]);
	if (nonnegative_equality)
		return value;
	if ((p->row_type[i] == optm_CONS_T_GE && value < 0.) ||
		(p->row_type[i] == optm_CONS_T_LE && value > 0.))
		return 0.;
	return value;
}

/* Error-free expansion addition. The last component carries the exact
 * sign of the sum. FMA supplies the low part of every binary product;
 * thus cancellation on an unbounded column is checked, not thresholded. */
static void certificate_add(double *expansion, int *length, double value)
{
	int i, count = 0;
	for (i = 0; i < *length; i++) {
		double previous = expansion[i];
		double sum = value + previous;
		double tail = sum - value;
		double error = (value - (sum - tail)) + (previous - tail);
		if (error != 0.)
			expansion[count++] = error;
		value = sum;
	}
	if (value != 0.)
		expansion[count++] = value;
	*length = count;
}

static int certificate_separates(const struct simplex_Problem *p,
	const double *y, int sign, double *expansion, int nonnegative_equality)
{
	long double gap = 0.;
	int i, j, k;

	for (i = 0; i < p->rows; i++) {
		long double weight = certificate_multiplier(
			p, y, i, sign, nonnegative_equality);
		long double term = weight * p->rhs[i];
		if (!isfinite(term))
			return 0;
		if (weight != 0.)
			gap = nextafterl(
				gap + nextafterl(term, -HUGE_VALL), -HUGE_VALL);
	}
	for (j = 0; j < p->columns; j++) {
		const struct optm_VariableBound *bound =
			nonnegative_equality ? NULL : p->bounds + j;
		long double lower = 0., upper = 0., term;
		int length = 0, q;
		double endpoint;

		for (k = p->matrix.column_start[j];
			k < p->matrix.column_start[j + 1]; k++) {
			double weight = (double)certificate_multiplier(p, y,
				p->matrix.row_index[k], sign,
				nonnegative_equality);
			double value = p->matrix.value[k];
			double product = value * weight;
			if (!isfinite(product) ||
				(value != 0. && weight != 0. &&
					fabs(product) < DBL_MIN / DBL_EPSILON))
				return 0;
			certificate_add(expansion, &length,
				fma(value, weight, -product));
			certificate_add(expansion, &length, product);
		}
		if (length == 0)
			continue;
		for (q = 0; q < length; q++) {
			lower = nextafterl(lower + expansion[q], -HUGE_VALL);
			upper = nextafterl(upper + expansion[q], HUGE_VALL);
		}
		if (!isfinite(lower) || !isfinite(upper))
			return 0;
		if (nonnegative_equality) {
			/* z >= 0: a positive column product has no finite
			 * upper support; a nonpositive product has support zero. */
			if (expansion[length - 1] > 0.)
				return 0;
			endpoint = 0.;
		} else if (expansion[length - 1] > 0.)
			endpoint = (bound->b_type == optm_BOUND_T_UP ||
					   bound->b_type == optm_BOUND_T_BS)
					   ? bound->ub
					   : HUGE_VAL;
		else
			endpoint = (bound->b_type == optm_BOUND_T_LO ||
					   bound->b_type == optm_BOUND_T_BS)
					   ? bound->lb
					   : -HUGE_VAL;
		if (!isfinite(endpoint))
			return 0;
		term = endpoint == 0.
			       ? 0.
			       : nextafterl((endpoint > 0. ? upper : lower) *
						    endpoint,
					 HUGE_VALL);
		if (!isfinite(term))
			return 0;
		gap = nextafterl(gap - term, -HUGE_VALL);
	}
	return isfinite(gap) && gap > 0.;
}

static int certificate_verify(const struct simplex_Problem *p, const double *y,
	int nonnegative_equality)
{
	double *expansion;
	int i, maximum_degree = 0, valid;
	for (i = 0; i < p->rows; i++)
		if (!isfinite(y[i]))
			return 0;
	for (i = 0; i < p->columns; i++) {
		int degree = p->matrix.column_start[i + 1] -
			     p->matrix.column_start[i];
		if (degree > maximum_degree)
			maximum_degree = degree;
	}
	expansion = (double *)lp_simplex_malloc(
		(2 * (size_t)maximum_degree + 1) * sizeof(double));
	if (expansion == NULL)
		return 0;
	valid = certificate_separates(
			p, y, 1, expansion, nonnegative_equality) ||
		certificate_separates(
			p, y, -1, expansion, nonnegative_equality);
	lp_simplex_free(expansion);
	return valid;
}

int simplex_problem_verify_infeasibility(
	const struct simplex_Problem *p, const double *y)
{
	return certificate_verify(p, y, 0);
}

int simplex_problem_verify_nonnegative_infeasibility(
	const struct simplex_CscMatrix *matrix, const double *rhs,
	const double *y)
{
	struct simplex_Problem equality;

	/* A borrowed equality/nonnegative view, not a new algorithm or model
	 * transformation. Only the verifier's explicit mode supplies row/bound
	 * semantics; ordinary Problem verification still requires both arrays. */
	lp_simplex_memset(&equality, 0, sizeof(equality));
	equality.rows = matrix->rows;
	equality.columns = matrix->columns;
	equality.matrix = *matrix;
	equality.rhs = (double *)rhs;
	return certificate_verify(&equality, y, 1);
}
