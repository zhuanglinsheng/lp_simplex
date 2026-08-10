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

enum lp_simplex_Algorithm {
	lp_simplex_ALGORITHM_TABLEAU = 0,
	lp_simplex_ALGORITHM_DUAL_REVISED = 1,
	/** Pan's generalized simplex with a dynamically deficient basis. */
	lp_simplex_ALGORITHM_PAN_BDA = 2
};

enum lp_simplex_Pricing {
	lp_simplex_PRICING_BLAND = 0,
	lp_simplex_PRICING_DANTZIG = 1,
	lp_simplex_PRICING_DUAL_STEEPEST_EDGE = 2,
	/** Euclidean-normalized violated-constraint selection for Pan/BDA. */
	lp_simplex_PRICING_PAN_NORMALIZED = 3
};

struct lp_simplex_Options {
	int algorithm;
	int pricing;
	int iteration_limit;
	double primal_tolerance;
	double dual_tolerance;
	double pivot_tolerance;
	/** Nonzero enables presolve; defaults to 1. */
	int presolve;
};

struct lp_simplex_Result {
	int status;
	int iterations;
	double objective;
	double primal_infeasibility;
	double dual_infeasibility;
};

/** Fill an option structure with the defaults for the selected algorithm. */
void lp_simplex_default_options(struct lp_simplex_Options *options,
		int algorithm);

/**
 * Solve a continuous LP with the configured simplex implementation.
 *
 * `x` must provide `model->n` elements.  A successful call returns
 * `lp_simplex_EXIT_SUCCESS` and sets `result->status` to
 * `lp_simplex_Success`.
 */
int lp_simplex_solve(
		const struct lp_Model *model,
		const struct lp_simplex_Options *options,
		double *x, struct lp_simplex_Result *result);

#ifdef __cplusplus
}
#endif

#endif
