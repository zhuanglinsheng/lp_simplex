/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#ifndef LP_SIMPLEX_DUAL_INTERNAL_H
#define LP_SIMPLEX_DUAL_INTERNAL_H

#include "simplex_problem.h"

#include <lp_simplex/model.h>
#include <lp_simplex/solve.h>


int simplex_dual_solve(
		const struct lp_Model *model,
		const struct lp_simplex_Options *options,
		double *x, double *row_dual, struct lp_simplex_Result *result,
		int propagate_bounds);

int simplex_dual_solve_problem(
		const struct simplex_Problem *problem,
		const struct lp_simplex_Options *options,
		double *x, double *row_dual, struct lp_simplex_Result *result,
		int propagate_bounds);

#endif
