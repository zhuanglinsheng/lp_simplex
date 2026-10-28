/**
 * @file search.h
 * @brief Internal policy defaults; units, rationale and limits:
 *
 * docs/numerical-policy.md. Search/resource defaults are not certificates.
 */

#ifndef LP_SIMPLEX_POLICY_SEARCH_H
#define LP_SIMPLEX_POLICY_SEARCH_H

#include <limits.h>

/* src/core/simplex_auto.h */
#define SIMPLEX_AUTO_PROBE_CAP SIMPLEX_DUAL_EXPLORE_INITIAL
#define SIMPLEX_AUTO_COVER_COLUMN_MARGIN 16
#define SIMPLEX_AUTO_COVER_MIN_SWEEPS 8
#define SIMPLEX_AUTO_INITIALIZER_OBSERVATIONS 8

/* src/dual/simplex_degeneracy.h */
#define SIMPLEX_DUAL_PROGRESS_WINDOW 128
#define SIMPLEX_DUAL_OBSERVATION_WINDOWS 8
#define SIMPLEX_DUAL_EXPLORE_INITIAL \
	(SIMPLEX_DUAL_OBSERVATION_WINDOWS * SIMPLEX_DUAL_PROGRESS_WINDOW)
#define SIMPLEX_DUAL_MIN_LEASE_REDUCTION .01
#define SIMPLEX_DUAL_RECOVERY_ALLOWANCE (16 * SIMPLEX_DUAL_EXPLORE_INITIAL)

/* src/dual/simplex_recovery.h */
#define SIMPLEX_PRIMAL_STALL_LIMIT 64

/* src/dual/simplex_dual_kernel.c */
#define DUAL_PERTURB_SHIFT_CAP 64.
/* Initial basis exchanges have a separate attempt cap, not a new pivot
 * grant. Dimensions are validated nonnegative ints; arguments have no side
 * effects. Saturation preserves a finite cap without signed overflow. */
#define SIMPLEX_DUAL_CRASH_ATTEMPT_LIMIT(rows, variables) \
	((variables) > (INT_MAX - (rows)) / 4 ? INT_MAX \
					      : 4 * (variables) + (rows))

/* Presolve optional search work. */
#define SIMPLEX_PRESOLVE_DOMINANCE_WORK_FACTOR 64
/* Index visits and exact/modular operations; independent of solve work. */
#define SIMPLEX_PRESOLVE_DEPENDENCY_WORK_MULTIPLIER 8
/* Non-pivot restarts need a separate finite bound; they consume no pivots. */
#define SIMPLEX_DUAL_RESTART_ROW_FACTOR 4
#define SIMPLEX_DUAL_RESTART_ALLOWANCE 64
#define SIMPLEX_DUAL_REPORT_PERIOD 256
/* A bound perturbation is limited relative to its per-variable tolerance. */
#define SIMPLEX_BOUND_PERTURB_TOLERANCE_UNITS 64.

#endif
