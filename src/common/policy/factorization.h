/**
 * @file factorization.h
 * @brief Internal policy defaults; units, rationale and limits:
 *
 * docs/numerical-policy.md. Search/resource defaults are not certificates.
 */

#ifndef LP_SIMPLEX_POLICY_FACTORIZATION_H
#define LP_SIMPLEX_POLICY_FACTORIZATION_H

/* src/basis/simplex_basis.c */
#define SIMPLEX_BASIS_UPDATE_LIMIT 512
#define SIMPLEX_BASIS_STABILITY_PIVOT 1e-7

/* src/basis/simplex_sparse_lu.c */
#define SPARSE_LU_PIVOT_TOLERANCE 1e-13
#define SPARSE_LU_MARKOWITZ_THRESHOLD 1e-1
#define SPARSE_LU_GLOBAL_MARKOWITZ_LIMIT 1024
#define SPARSE_LU_FT_UPDATE_LIMIT 512

/* src/dual/simplex_dual_kernel.c */
#define DUAL_RELATIVE_PIVOT_TOLERANCE 1e-8

/* Absolute zero-pivot rejection in basis solve/update coordinates. These
 * legacy cutoffs are not the user's rank/pivot tolerance. Do not merge them
 * with normalized residual or relative-pivot controls. */
#define SIMPLEX_BASIS_ABSOLUTE_PIVOT_FLOOR 1e-14
#define SIMPLEX_BASIS_COMPACT_RESIDUAL_TARGET 1e-11
#define SIMPLEX_BASIS_COMPACT_REFINEMENTS 2
#define SIMPLEX_BASIS_TERMINAL_REFINEMENTS 3
#define SIMPLEX_BASIS_VALIDATION_PERIOD 32
/* Legacy size/cost gates, pending independent timing evidence. */
#define SIMPLEX_BASIS_LARGE_ROWS 4096
#define SIMPLEX_DUAL_STRUCTURAL_PREFERENCE_ROWS 4096
/* Relative comparison/audit floors used by pricing, not feasibility tests. */
#define SIMPLEX_PRICING_TIE_RELATIVE 1e-12
#define SIMPLEX_PRICING_AUDIT_RELATIVE 1e-11
#define SIMPLEX_DUAL_EDGE_WEIGHT_FLOOR 1e-12
/* Updated squared DSE norm must be at least a quarter of its BTRAN norm.
 * This bounds optimistic pricing error, not feasibility or optimality error. */
#define SIMPLEX_DUAL_DSE_ACCEPT_RATIO .25
#define SIMPLEX_PRESOLVE_SUBSTITUTION_RELATIVE 1e-10
#define SIMPLEX_PRESOLVE_BOUND_NOISE_FRACTION 1e-4
#define SIMPLEX_PRESOLVE_BOUND_NOISE_FLOOR 1e-12
#define SIMPLEX_PAN_NEWTON_RANK_TOLERANCE 1e-12
#define SIMPLEX_PAN_CG_ABSOLUTE_TARGET 1e-15
#define SIMPLEX_PAN_CG_RELATIVE_TARGET 1e-10

/* KLU pivot preference; distinct from internal LU zero rejection. */
#define SIMPLEX_KLU_PIVOT_PREFERENCE .1
/* Maximum ratio rounding uncertainty admitted to Harris test. */
#define SIMPLEX_DUAL_RATIO_NOISE_LIMIT .25
/* Relative FTRAN residual stagnation gate. */
#define SIMPLEX_DUAL_FTRAN_PROGRESS_RATIO .99
/* Heap row rescore trigger relative to its stored score. */
#define SIMPLEX_DUAL_FEASIBILITY_REPRICE_RATIO .8
/* Reset Devex reference when norm becomes too small. */
#define SIMPLEX_PRIMAL_DEVEX_RESET_RATIO .1
/* Half an integer violation count; progress signal noise, not feasibility. */
#define SIMPLEX_DUAL_COUNT_NOISE .5
/* Cover certification cancellation guard. */
#define SIMPLEX_COVER_ROUNDING_LIMIT .25
/* Full-column cover certificate refinement factor. */
#define SIMPLEX_COVER_RETRY_TOLERANCE_FACTOR .1

/* Separate performance gates retain distinct names even at equal defaults. */
#define SIMPLEX_BASIS_COMPACT_MIN_ROWS 4096
#define SIMPLEX_BASIS_COMPACT_ROW_RATIO 3
#define SIMPLEX_DUAL_STABLE_ORDER_ROWS 2048
#define SIMPLEX_DUAL_DENSE_MEMBERSHIP_DIVISOR 20

#endif
