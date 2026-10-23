/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_pan.h
 * @brief Pan basis-deficiency-allowing generalized simplex solver.
 */

#ifndef LP_SIMPLEX_PAN_INTERNAL_H
#define LP_SIMPLEX_PAN_INTERNAL_H

#include "simplex_problem.h"

#include <lp_simplex/solve.h>


int simplex_pan_solve_problem(const struct simplex_Problem *problem,
	const struct lp_simplex_Options *options, double *x,
	struct lp_simplex_Result *result);

/* A feasible basic point may be supplied to bypass NNLS Phase I. The active
 * support must admit a full-column-rank factor; no optimality is presumed. */
int simplex_pan_solve_from_feasible(const struct simplex_Problem *problem,
	const struct lp_simplex_Options *options, const double *start,
	double *x, struct lp_simplex_Result *result);

#endif
