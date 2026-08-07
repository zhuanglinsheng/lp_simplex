/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#ifndef LP_SIMPLEX_SOLVE_H
#define LP_SIMPLEX_SOLVE_H

#include "model.h"
#include "status.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Solve a continuous LP model with the two-phase tableau simplex method.
 *
 * `x` must provide model->n elements.  `criteria` accepts "bland" (the safe
 * default), "dantzig", or NULL.  The function returns EXIT_SUCCESS only when
 * an optimum has been produced; `status` contains the detailed termination
 * reason.
 */
int lp_simplex_solve(const struct lp_Model *model, const char *criteria,
		int iteration_limit, double *x, double *objective, int *status);

#ifdef __cplusplus
}
#endif

#endif
