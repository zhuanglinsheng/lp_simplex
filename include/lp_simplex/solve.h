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

/** Default cumulative work budget for AUTO; other algorithms use 100000. */
#define LP_SIMPLEX_DEFAULT_AUTO_ITERATION_LIMIT 300000

/**
 * @brief Algorithm choices accepted by lp_simplex_Options.algorithm.
 */
enum lp_simplex_Algorithm {
	/** Two-phase dense tableau; may materialise sparse input as a dense
	 * matrix. */
	lp_simplex_ALGORITHM_TABLEAU = 0,

	/** Sparse dual revised simplex. */
	lp_simplex_ALGORITHM_DUAL_REVISED = 1,

	/** Pan's generalized simplex with a dynamically deficient basis. */
	lp_simplex_ALGORITHM_PAN_BDA = 2,

	/** Automatic staging: bounded dual probe, native BDA and numerical
	 * fallback. */
	lp_simplex_ALGORITHM_AUTO = 3
};

/**
 * Pricing choices. Legal combinations (first listed is the default):
 * - TABLEAU: BLAND, DANTZIG.
 * - DUAL_REVISED and AUTO: DUAL_STEEPEST_EDGE only.
 * - PAN_BDA: DANTZIG, PAN_NORMALIZED.
 * Other algorithm/pricing combinations are rejected, not silently converted.
 */
enum lp_simplex_Pricing {
	/** Bland's index-based selection for tableau. */
	lp_simplex_PRICING_BLAND = 0,

	/** Dantzig selection for tableau or unnormalised selection for Pan/BDA. */
	lp_simplex_PRICING_DANTZIG = 1,

	/** Dual steepest-edge pricing for dual revised and AUTO. */
	lp_simplex_PRICING_DUAL_STEEPEST_EDGE = 2,

	/** Euclidean-normalized violated-constraint selection for Pan/BDA. */
	lp_simplex_PRICING_PAN_NORMALIZED = 3
};

/** Original-model progress. Bounds bracket the optimal objective (including
 * offset and sense); unavailable bounds are -/+infinity. A missing gap is
 * infinity, never zero. Feasibility uses primal_tolerance. Reporting does not
 * change the kernel's optimality tests.
 *
 * point is non-NULL only for a newly improved, independently checked feasible
 * point, in original variable order. Copy it during the callback if needed;
 * all pointers expire on return. stage_iterations is local to the named
 * kernel and may restart when AUTO changes route. */
struct lp_simplex_Progress {
	const char *stage;
	int stage_iterations;
	double cpu_seconds;
	double lower_bound;
	double upper_bound;
	double absolute_gap;
	/** absolute_gap / max(1, |lower_bound|, |upper_bound|). */
	double relative_gap;
	const double *point;
	double primal_infeasibility;
};

typedef void (*lp_simplex_ProgressCallback)(
	void *context, const struct lp_simplex_Progress *progress);

/**
 * @brief Per-call solver configuration.
 *
 * Initialise with lp_simplex_default_options(); zero-initialisation alone is
 * invalid. All tolerances must be finite and strictly positive. Reinitialise
 * defaults when switching algorithms so that pricing remains compatible.
 * No wall-clock timeout is provided by this API.
 * Work budgets include all AUTO stages and retries. Primal tolerance bounds
 * the final original-model maximum absolute row/bound violation. Dual
 * tolerance follows kernel scaling and is not an absolute threshold for every
 * public dual diagnostic.
 */
struct lp_simplex_Options {
	/**< An lp_simplex_Algorithm value. */
	int algorithm;

	/**< Compatible lp_simplex_Pricing value; see the pricing enum. */
	int pricing;

	/**< Positive work budget; AUTO 300000, others 100000. */
	int iteration_limit;

	/**< Original primal feasibility tolerance; default 1e-7. */
	double primal_tolerance;

	/**< Kernel dual search/optimality tolerance; default 1e-7. */
	double dual_tolerance;

	/**< Numerical pivot tolerance; default 1e-9. */
	double pivot_tolerance;

	/**< Exactly 0 disables presolve; exactly 1 enables it (default). */
	int presolve;

	/** Optional synchronous callback; NULL (default) disables progress work.
	 * Sparse dual/native Pan and PDHG emit intermediate events. All public
	 * algorithms emit initial/final events. No optimality assertion is made
	 * solely from a small reported gap. Do not mutate the model or solver
	 * storage from the callback. Rebuild callers when upgrading this struct. */
	lp_simplex_ProgressCallback progress_callback;
	void *progress_context;
	/** Positive CPU seconds between intermediate reports; default 5.
	 * First feasible points and final reports bypass this interval. */
	double progress_interval;
};

/**
 * @brief Model space of the dual diagnostic; never infer it from solve status.
 */
enum lp_simplex_ResidualSpace {
	/** No public dual diagnostic; dual_infeasibility is NaN, not zero. */
	lp_simplex_RESIDUAL_UNAVAILABLE = 0,

	/** Original caller model. */
	lp_simplex_RESIDUAL_ORIGINAL = 1,

	/** Presolved model; original row multipliers have not been recovered. */
	lp_simplex_RESIDUAL_PRESOLVED = 2,

	/** Standard-form model, before recovering original variables/rows. */
	lp_simplex_RESIDUAL_STANDARD = 3,

	/** Scaled model. */
	lp_simplex_RESIDUAL_SCALED = 4,

	/** Auxiliary problem rather than the original objective. */
	lp_simplex_RESIDUAL_AUXILIARY = 5
};

/**
 * @brief Caller-owned solve diagnostics; no fields require separate freeing.
 *
 * A finite candidate point on failure is not an optimal solution. Objective
 * and primal diagnostics describe the returned candidate when available.
 * Missing dual diagnostics do not by themselves invalidate a successful solve.
 * Objective values are recomputed at the returned point with offset and sense
 * restored. Primal residuals use original units; unavailable or nonfinite
 * points/diagnostics give positive infinity. The dual diagnostic is neither a
 * primal-dual gap nor an independent optimality proof; do not replace kernel
 * acceptance with a comparison against dual_tolerance.
 */
struct lp_simplex_Result {
	/**< lp_Status code; lp_simplex_Success means optimal. */
	int status;

	/**< Cumulative work, including stages/retries; not time. */
	int iterations;

	/**< Original objective including offset/sense; NaN if absent. */
	double objective;

	/** Original max row/bound violation; infinity if unavailable. */
	double primal_infeasibility;

	/** KKT sign diagnostic in dual_residual_space; NaN if absent. */
	double dual_infeasibility;

	/**< lp_simplex_ResidualSpace for the dual diagnostic only. */
	int dual_residual_space;
};

/**
 * @brief Fill every option field with defaults for the chosen algorithm.
 *
 * Presolve is enabled; primal/dual tolerances are 1e-7 and pivot tolerance is
 * 1e-9. Pricing defaults are documented by lp_simplex_Pricing; iteration limits
 * by lp_simplex_Options.iteration_limit. Invalid algorithm values are copied
 * into options and subsequently rejected by lp_simplex_solve().
 *
 * @param options Writable output structure; NULL is a no-op.
 * @param algorithm An lp_simplex_Algorithm value (AUTO for automatic selection).
 */
void lp_simplex_default_options(
	struct lp_simplex_Options *options, int algorithm);

/**
 * @brief Solve a continuous LP without taking ownership of caller storage.
 *
 * Initialise options with lp_simplex_default_options(), then apply overrides.
 * Input dimensions must be positive; row types, bounds and sparse storage must
 * be valid. Integer/binary variables are rejected rather than relaxed.
 * Use finite matrix entries, costs, objective offset and RHS values. When both
 * dense and sparse storage exist, keep all representations consistent.
 *
 * x is output only, not a warm start. After input validation it is initialised
 * to NaN; a failed solve may overwrite it with a candidate point. An invalid
 * input leaves x unchanged and initialises result to CondUnsatisfied, zero
 * iterations, NaN objective/dual and infinite primal infeasibility.
 *
 * Success requires kernel optimality plus a finite objective and an original
 * primal residual <= primal_tolerance. A provided dual diagnostic must be
 * finite, but its space may differ from the original model. On any failure,
 * inspect result.status; never treat a finite x or objective alone as success.
 * No primal ray, row multipliers or infeasibility certificate is exported.
 *
 * @param model Non-NULL input model; retained by the caller and not modified.
 * @param options Non-NULL valid solver configuration; not modified.
 * @param x Non-NULL writable array of at least model->n doubles, in original
 * variable order. Must not alias the input model or other output storage.
 * @param result Non-NULL writable diagnostics. If NULL, returns failure
 * immediately without writing outputs.
 * @return lp_simplex_EXIT_SUCCESS (0) exactly when result->status is
 * lp_simplex_Success; otherwise lp_simplex_EXIT_FAILURE (-1).
 */
int lp_simplex_solve(const struct lp_Model *model,
	const struct lp_simplex_Options *options, double *x,
	struct lp_simplex_Result *result);

#ifdef __cplusplus
}
#endif

#endif
