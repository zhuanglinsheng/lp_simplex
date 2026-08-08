/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#ifndef LP_SIMPLEX_SINGLETON_DUAL_INTERNAL_H
#define LP_SIMPLEX_SINGLETON_DUAL_INTERNAL_H

#include <lp_simplex/model.h>
#include <lp_simplex/solve.h>


int simplex_singleton_dual_solve(
		const struct lp_Model *model,
		const struct lp_simplex_Options *options,
		double *x, struct lp_simplex_Result *result, int *applicable);

#endif
