/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#ifndef LP_SIMPLEX_PHASE_H
#define LP_SIMPLEX_PHASE_H

#include <lp_simplex/model.h>
#include <lp_simplex/solve.h>


int simplex_solve_standard(
		const double *objective,
		const struct optm_LinearConstraint *constraints,
		int m, int n, const struct lp_simplex_Options *options,
		double *x, double *value, int *status, int *iterations);

#endif
