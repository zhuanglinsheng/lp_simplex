/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_basis_certificate.c
 * @brief Optional exact verification of a terminal basis-row certificate.
 */

#include "policy/resources.h"
#include "simplex_problem.h"
#include "utils.h"

#ifdef LP_SIMPLEX_HAVE_GMP
#include <gmp.h>
#include <math.h>

/* Bound the optional cubic solve and rational storage. Exceeding either
 * limit means "unverified", never infeasible or feasible. Logical columns
 * are eliminated first, so these limits apply to structural basic columns. */

static int certificate_size_ok(const mpq_t value)
{
	return mpz_sizeinbase(mpq_numref(value), 2) <= CERTIFICATE_MAX_BITS &&
	       mpz_sizeinbase(mpq_denref(value), 2) <= CERTIFICATE_MAX_BITS;
}

static int certificate_solve(mpq_t *a, int n)
{
	mpq_t pivot, product;
	int i, j, k, p, valid = 0;
	mpq_init(pivot);
	mpq_init(product);
	for (k = 0; k < n; k++) {
		p = -1;
		for (i = k; i < n; i++)
			if (mpq_sgn(a[i * (n + 1) + k]) != 0 &&
				(p < 0 || fabs(mpq_get_d(a[i * (n + 1) + k])) >
						  fabs(mpq_get_d(
							  a[p * (n + 1) + k]))))
				p = i;
		if (p < 0)
			goto done;
		for (j = k; j <= n; j++)
			mpq_swap(a[p * (n + 1) + j], a[k * (n + 1) + j]);
		mpq_set(pivot, a[k * (n + 1) + k]);
		for (j = k + 1; j <= n; j++) {
			mpq_div(a[k * (n + 1) + j], a[k * (n + 1) + j], pivot);
			if (!certificate_size_ok(a[k * (n + 1) + j]))
				goto done;
		}
		mpq_set_ui(a[k * (n + 1) + k], 1, 1);
		for (i = 0; i < n; i++) {
			if (i == k)
				continue;
			mpq_set(pivot, a[i * (n + 1) + k]);
			for (j = k + 1; j <= n; j++) {
				mpq_mul(product, pivot, a[k * (n + 1) + j]);
				mpq_sub(a[i * (n + 1) + j], a[i * (n + 1) + j],
					product);
				if (!certificate_size_ok(a[i * (n + 1) + j]))
					goto done;
			}
			mpq_set_ui(a[i * (n + 1) + k], 0, 1);
		}
	}
	valid = 1;
done:
	mpq_clear(pivot);
	mpq_clear(product);
	return valid;
}

static int certificate_separates(
	const struct simplex_Problem *p, mpq_t *y, int sign)
{
	mpq_t gap, coefficient, value, product, weight;
	int i, j, k, valid = 0;
	mpq_init(gap);
	mpq_init(coefficient);
	mpq_init(value);
	mpq_init(product);
	mpq_init(weight);
	for (i = 0; i < p->rows; i++) {
		mpq_set(weight, y[i]);
		if (sign < 0)
			mpq_neg(weight, weight);
		if ((p->row_type[i] == optm_CONS_T_GE && mpq_sgn(weight) < 0) ||
			(p->row_type[i] == optm_CONS_T_LE &&
				mpq_sgn(weight) > 0))
			goto done;
		mpq_set_d(value, p->rhs[i]);
		mpq_mul(product, weight, value);
		mpq_add(gap, gap, product);
	}
	for (j = 0; j < p->columns; j++) {
		const struct optm_VariableBound *bound = p->bounds + j;
		double endpoint;
		mpq_set_ui(coefficient, 0, 1);
		for (k = p->matrix.column_start[j];
			k < p->matrix.column_start[j + 1]; k++) {
			mpq_set_d(value, p->matrix.value[k]);
			mpq_mul(product, value, y[p->matrix.row_index[k]]);
			mpq_add(coefficient, coefficient, product);
		}
		if (sign < 0)
			mpq_neg(coefficient, coefficient);
		if (mpq_sgn(coefficient) == 0)
			continue;
		if (mpq_sgn(coefficient) > 0)
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
			goto done;
		mpq_set_d(value, endpoint);
		mpq_mul(product, value, coefficient);
		mpq_sub(gap, gap, product);
	}
	valid = mpq_sgn(gap) > 0;
done:
	mpq_clear(gap);
	mpq_clear(coefficient);
	mpq_clear(value);
	mpq_clear(product);
	mpq_clear(weight);
	return valid;
}
#endif

int simplex_problem_verify_basis_infeasibility(
	const struct simplex_Problem *p, const int *basis, int leaving)
{
#ifdef LP_SIMPLEX_HAVE_GMP
	mpq_t *a = NULL, *y = NULL;
	mpq_t value, product;
	int *index = NULL, *position = NULL;
	int i, j, k, n = 0, row = 0, valid = 0, initialized = 0;
	if (leaving < 0 || leaving >= p->rows)
		return 0;
	for (i = 0; i < p->rows; i++)
		if (!isfinite(p->rhs[i]))
			return 0;
	for (i = 0; i < p->matrix.nonzeros; i++)
		if (!isfinite(p->matrix.value[i]))
			return 0;
	index = (int *)lp_simplex_malloc(((size_t)p->rows + 1) * sizeof(int));
	position = (int *)lp_simplex_malloc(
		((size_t)p->rows + p->columns + 1) * sizeof(int));
	if (!index || !position)
		goto done;
	for (i = 0; i < p->rows + p->columns; i++)
		position[i] = -1;
	for (i = 0; i < p->rows; i++) {
		int variable = basis[i];
		if (variable < 0 || variable >= p->rows + p->columns ||
			position[variable] >= 0)
			goto done;
		position[variable] = i;
		if (variable < p->columns)
			n++;
	}
	if (n > CERTIFICATE_MAX_RANK)
		goto done;
	for (i = 0; i < p->rows; i++)
		index[i] = position[p->columns + i] < 0 ? row++ : -1;
	if (row != n)
		goto done;
	a = (mpq_t *)lp_simplex_malloc(
		((size_t)n * (n + 1) + 1) * sizeof(mpq_t));
	y = (mpq_t *)lp_simplex_malloc(((size_t)p->rows + 1) * sizeof(mpq_t));
	if (!a || !y)
		goto done;
	for (i = 0; i < n * (n + 1); i++)
		mpq_init(a[i]);
	for (i = 0; i < p->rows; i++) {
		mpq_init(y[i]);
		if (position[p->columns + i] == leaving)
			mpq_set_si(y[i], -1, 1);
	}
	mpq_init(value);
	mpq_init(product);
	initialized = 1;
	row = 0;
	for (i = 0; i < p->rows; i++)
		if (basis[i] < p->columns) {
			j = basis[i];
			mpq_set_ui(
				a[row * (n + 1) + n], i == leaving ? 1 : 0, 1);
			for (k = p->matrix.column_start[j];
				k < p->matrix.column_start[j + 1]; k++) {
				int r = p->matrix.row_index[k], col = index[r];
				mpq_set_d(value, p->matrix.value[k]);
				if (col >= 0)
					mpq_add(a[row * (n + 1) + col],
						a[row * (n + 1) + col], value);
				else {
					mpq_mul(product, value, y[r]);
					mpq_sub(a[row * (n + 1) + n],
						a[row * (n + 1) + n], product);
				}
			}
			row++;
		}
	if (!certificate_solve(a, n))
		goto done;
	for (i = 0; i < p->rows; i++)
		if (index[i] >= 0)
			mpq_set(y[i], a[index[i] * (n + 1) + n]);
	valid = certificate_separates(p, y, 1) ||
		certificate_separates(p, y, -1);
done:
	if (initialized) {
		for (i = 0; i < n * (n + 1); i++)
			mpq_clear(a[i]);
		for (i = 0; i < p->rows; i++)
			mpq_clear(y[i]);
		mpq_clear(value);
		mpq_clear(product);
	}
	lp_simplex_free(a);
	lp_simplex_free(y);
	lp_simplex_free(index);
	lp_simplex_free(position);
	return valid;
#else
	(void)p;
	(void)basis;
	(void)leaving;
	return 0;
#endif
}
