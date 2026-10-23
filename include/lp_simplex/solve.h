/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file solve.h
 * @brief Public solver options, result diagnostics and solve entry points.
 */

#ifndef LP_SIMPLEX_SOLVE_H
#define LP_SIMPLEX_SOLVE_H

#include "model.h"
#include "status.h"


#ifdef __cplusplus
extern "C" {
#endif

#define LP_SIMPLEX_DEFAULT_AUTO_ITERATION_LIMIT 300000

enum lp_simplex_Algorithm {
	lp_simplex_ALGORITHM_TABLEAU = 0,
	lp_simplex_ALGORITHM_DUAL_REVISED = 1,
	/** Pan's generalized simplex with a dynamically deficient basis. */
	lp_simplex_ALGORITHM_PAN_BDA = 2,
	/** Bounded dual probe, native primal-dual/BDA, numerical fallback. */
	lp_simplex_ALGORITHM_AUTO = 3
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
	/** Shared work budget across all AUTO stages. */
	int iteration_limit;
	double primal_tolerance;
	double dual_tolerance;
	double pivot_tolerance;
	/** Nonzero enables presolve; defaults to 1. */
	int presolve;
};

/** Model space of the reported dual diagnostic; never inferred from status. */
enum lp_simplex_ResidualSpace {
	lp_simplex_RESIDUAL_UNAVAILABLE = 0,
	lp_simplex_RESIDUAL_ORIGINAL = 1,
	lp_simplex_RESIDUAL_PRESOLVED = 2,
	lp_simplex_RESIDUAL_STANDARD = 3,
	lp_simplex_RESIDUAL_SCALED = 4,
	lp_simplex_RESIDUAL_AUXILIARY = 5
};

struct lp_simplex_Result {
	int status;
	int iterations;
	/** Original objective (including offset/sense) at returned x; NaN if absent. */
	double objective;
	/** Maximum absolute original row/bound violation; infinity if x is absent. */
	double primal_infeasibility;
	/** Dual diagnostic in dual_residual_space; NaN means not available.
	 * This value alone is not a primal-dual gap or an optimality certificate. */
	double dual_infeasibility;
	int dual_residual_space;
};

/** Fill an option structure with the defaults for the selected algorithm. */
void lp_simplex_default_options(
	struct lp_simplex_Options *options, int algorithm);

/**
 * Solve a continuous LP with the configured simplex implementation.
 *
 * `x` must provide `model->n` elements.  A successful call returns
 * `lp_simplex_EXIT_SUCCESS` and sets `result->status` to
 * `lp_simplex_Success`.
 */
int lp_simplex_solve(const struct lp_Model *model,
	const struct lp_simplex_Options *options, double *x,
	struct lp_simplex_Result *result);

#ifdef __cplusplus
}
#endif

#endif
