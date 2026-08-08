/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#ifndef LP_SIMPLEX_TABLEAU_SOLVER_INTERNAL_H
#define LP_SIMPLEX_TABLEAU_SOLVER_INTERNAL_H

#include <lp_simplex/model.h>


int simplex_tableau_solve_model(
		const struct lp_Model *model,
		const char *criteria, int iteration_limit,
		double *x, double *objective, int *status);

#endif
