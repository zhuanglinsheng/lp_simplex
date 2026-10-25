/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_dual_backend.h
 * @brief Dual numerical backend recovery and cumulative attempt accounting.
 */

#ifndef LP_SIMPLEX_DUAL_BACKEND_INTERNAL_H
#define LP_SIMPLEX_DUAL_BACKEND_INTERNAL_H

struct simplex_Problem;
struct lp_simplex_Options;
struct lp_simplex_Result;
struct simplex_DualBasisSeed;
struct simplex_DualRetry;

/* Run the default factorization and its eligible numerical fallback within
 * one Dual grant. Dimensions describe this borrowed kernel-space problem;
 * outputs and optional basis storage remain caller-owned. The retry request
 * is already bound to the immutable input and its representation by the entry.
 * Retain/audit candidates, charge both attempts and preserve causal sources.
 * No cross-algorithm selection or AUTO runtime mutation is performed here. */
int simplex_dual_backend_solve(const struct simplex_Problem *problem,
	const struct lp_simplex_Options *options, double *x, double *row_dual,
	struct lp_simplex_Result *result, int propagate_bounds, int columns,
	int rows, int feasibility_only, struct simplex_DualBasisSeed *seed,
	const struct simplex_DualRetry *retry);

#endif
