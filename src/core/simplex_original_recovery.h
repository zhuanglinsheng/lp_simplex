/**
 * @file simplex_original_recovery.h
 * @brief Original-model result audit and a bounded cold retry.
 */

#ifndef LP_SIMPLEX_ORIGINAL_RECOVERY_INTERNAL_H
#define LP_SIMPLEX_ORIGINAL_RECOVERY_INTERNAL_H
#include "simplex_execution.h"

struct lp_Model;
struct simplex_ProgressSink;
/* A fresh original-model solve, not a postsolve call or a warm-state import.
 * The adapter must preserve algorithm/configuration and return local work.
 * The model and options remain immutable for this synchronous phase. */
typedef int (*simplex_OriginalSolve)(const struct lp_Model *,
	const struct lp_simplex_Options *, double *, struct lp_simplex_Result *,
	const struct simplex_ProgressSink *,
	const struct simplex_ExecutionScope *);
/* Only a fresh original-model residual may authorize feasibility recovery. */
int simplex_original_recovery_reason(
	int status, int work, double original_primal, double tolerance);
/* The point is in original coordinates; incoming diagnostics can explicitly
 * refer to presolved coordinates. Only fresh original audits authorize a
 * feasibility-based retry. One retry at most, under limit minus prior work.
 * No point, multiplier, basis or factors from the reduced model are imported
 * into the fresh solve. Trace capacity
 * does not decide eligibility, accounting or selection. Terminal classification
 * is independent of point selection; tolerance-feasible counterexamples
 * invalidate infeasibility claims even if their objectives are unavailable. */
int simplex_original_recovery_run(const struct lp_Model *,
	const struct lp_simplex_Options *, double *, struct lp_simplex_Result *,
	int, const struct simplex_ProgressSink *,
	const struct simplex_ExecutionScope *,
	const struct simplex_ExecutionRef *reduced_actor,
	simplex_OriginalSolve);
#endif
