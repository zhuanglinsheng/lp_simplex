/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#ifndef LP_SIMPLEX_PIVOT_H
#define LP_SIMPLEX_PIVOT_H

#include <lp_simplex/solve.h>

void simplex_apply_pivot(
		double *table,
		int ld, int m, int n,
		int leaving_row, int entering_column,
		int normalize, int eliminate_rows, int eliminate_objective);

int simplex_run_pivots(
		int *iteration, double *table, int ld, int *basis,
		int m, int n, int real_columns,
		const struct lp_simplex_Options *options);

#endif
