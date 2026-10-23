/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_pan_bounded.h
 * @brief Native-bound Pan/BDA solving, feasibility recovery and KKT checks.
 */

#ifndef LP_SIMPLEX_PAN_BOUNDED_INTERNAL_H
#define LP_SIMPLEX_PAN_BOUNDED_INTERNAL_H
#include "simplex_pan_standard.h"

#include <lp_simplex/solve.h>


/* Original structural columns and inequality slacks only. No bound rows,
 * free-variable splits, shifts or fixed-variable eliminations. */
struct simplex_PanBounded {
	struct simplex_PanStandard a;
	double *lower, *upper;
};
int simplex_pan_bounded_create(
	struct simplex_PanBounded *b, const struct simplex_Problem *p);
void simplex_pan_bounded_destroy(struct simplex_PanBounded *b);
int simplex_pan_bounded_encode(const struct simplex_PanBounded *b,
	const double *x, double *point, double *residual, double *work,
	double tolerance);
double simplex_pan_bounded_primal(const struct simplex_PanBounded *b,
	const double *x, double *residual, double *work);
int simplex_pan_bounded_check_dual(const struct simplex_PanBounded *b,
	const double *x, const double *y, double tolerance, double *violation);
/* Transactional multiplier search; the audit above remains read-only.
 * limit and iterations belong to the same cumulative solver budget. */
int simplex_pan_bounded_repair_dual(const struct simplex_PanBounded *b,
	const double *x, double *y, double tolerance, int limit,
	int *iterations, double *violation);
/* Exports only audited feasible improvements; 1 means full native KKT. */
int simplex_pan_bounded_pdhg(const struct simplex_PanBounded *b, double *point,
	double *multiplier, const struct lp_simplex_Options *o, int *iterations,
	int budget);
int simplex_pan_bounded_restore(const struct simplex_PanBounded *b,
	double *point, const struct lp_simplex_Options *o, int *iterations,
	long double ceiling);
/* Diagnostic export on success; multipliers use c + A^T y. */
int simplex_pan_bounded_solve_certified(const struct simplex_Problem *p,
	const struct lp_simplex_Options *o, const double *start, double *out,
	double *row_dual, struct lp_simplex_Result *result);
int simplex_pan_bounded_solve(const struct simplex_Problem *p,
	const struct lp_simplex_Options *options, const double *start,
	double *x, struct lp_simplex_Result *result);
#endif
