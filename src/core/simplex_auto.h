/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_auto.h
 * @brief Problem-independent orchestration; no global environment mutation.
 */

#ifndef LP_SIMPLEX_AUTO_INTERNAL_H
#define LP_SIMPLEX_AUTO_INTERNAL_H

#include "simplex_problem.h"
#include "simplex_degeneracy.h"

#include <lp_simplex/solve.h>

#define SIMPLEX_AUTO_PROBE_CAP SIMPLEX_DUAL_EXPLORE_INITIAL
#define SIMPLEX_AUTO_PROBE_MAX SIMPLEX_DUAL_EXPLORE_MAX
#define SIMPLEX_AUTO_DENSE_BYTES (32UL * 1024UL * 1024UL)
#define SIMPLEX_AUTO_CONVERSION_BYTES (64UL * 1024UL * 1024UL)

enum simplex_AutoRoute {
	SIMPLEX_AUTO_HYBRID = 0,
	SIMPLEX_AUTO_EXPLICIT,
	SIMPLEX_AUTO_CLASSIC,
	SIMPLEX_AUTO_STANDARD_IMPLICIT,
	SIMPLEX_AUTO_STANDARD_EXPLICIT,
	SIMPLEX_AUTO_SINGLETON,
	SIMPLEX_AUTO_TABLEAU,
	SIMPLEX_AUTO_ROUTE_COUNT
};

int simplex_auto_route_solve(int route, const struct simplex_Problem *problem,
	const struct lp_simplex_Options *options, const double *seed, double *x,
	struct lp_simplex_Result *result);

/* Resource eligibility, not a mathematical terminal classification. */
int simplex_auto_route_allowed(
	const struct simplex_Problem *problem, int route);

int simplex_auto_solve_problem(const struct simplex_Problem *problem,
	const struct lp_simplex_Options *options, double *x,
	struct lp_simplex_Result *result, int propagate_bounds);

/* Injectable internal kernel boundary, also used for deterministic failure
 * tests. No callback/global state escapes a solve. */
struct simplex_AutoKernels {
	int (*dual)(const struct simplex_Problem *,
		const struct lp_simplex_Options *, double *, double *,
		struct lp_simplex_Result *, int);

	int (*hybrid)(const struct simplex_Problem *,
		const struct lp_simplex_Options *, const double *, double *,
		struct lp_simplex_Result *);
};

/* Optional recovery routes are supplied separately so existing injectable
 * two-kernel users need not initialize additional function pointers. */
struct simplex_AutoRoutes {
	int (*solve)(int, const struct simplex_Problem *,
		const struct lp_simplex_Options *, const double *, double *,
		struct lp_simplex_Result *);
};

int simplex_auto_run_routes(const struct simplex_Problem *problem,
	const struct lp_simplex_Options *options, double *x,
	struct lp_simplex_Result *result, int propagate_bounds,
	const struct simplex_AutoKernels *kernels,
	const struct simplex_AutoRoutes *routes);

int simplex_auto_run(const struct simplex_Problem *problem,
	const struct lp_simplex_Options *options, double *x,
	struct lp_simplex_Result *result, int propagate_bounds,
	const struct simplex_AutoKernels *kernels);

#endif
