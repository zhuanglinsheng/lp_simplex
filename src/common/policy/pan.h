/**
 * @file pan.h
 * @brief Pan search policy. Empirical choices, not LP terminal criteria.
 *
 * Units and failure contracts: docs/numerical-policy.md.
 */

#ifndef LP_SIMPLEX_POLICY_PAN_H
#define LP_SIMPLEX_POLICY_PAN_H

/* Work grants; all remain bounded by the caller's cumulative budget. */
#define SIMPLEX_PAN_PDHG_WORK_CAP 500000
#define SIMPLEX_PAN_DUAL_POLISH_WORK_CAP 8192
/* Matrix-free recovery block; continuation requires residual halving. */
#define SIMPLEX_PAN_DUAL_KRYLOV_WORK_BLOCK 256
#define SIMPLEX_PAN_RESTORE_WORK_CAP 2048
#define SIMPLEX_PAN_NNLS_BLOCK_CAP 32
#define SIMPLEX_PAN_RESTORE_PASSES 64
#define SIMPLEX_PAN_TRIAL_PERIOD 32
#define SIMPLEX_PAN_PREDICT_WORK_CAP 512
#define SIMPLEX_PAN_LOW_PROGRESS_STEPS 3

/* PDHG interventions require a full restart window without 1% progress.
 * Productive epochs retain their metric; the stage budget bounds total work. */
#define SIMPLEX_PDHG_AUDIT_PERIOD 128
#define SIMPLEX_PDHG_REPORT_PERIOD 256
#define SIMPLEX_PDHG_OBSERVE_WORK (2 * SIMPLEX_PDHG_AUDIT_PERIOD)
#define SIMPLEX_PDHG_STAGNANT_WORK SIMPLEX_PDHG_RESTART_MAX_WORK
#define SIMPLEX_PDHG_MERIT_REDUCTION .01
#define SIMPLEX_PDHG_RESTART_REDUCTION .5
#define SIMPLEX_PDHG_RESTART_MIN_WORK 256
#define SIMPLEX_PDHG_RESTART_MAX_WORK 8192
/* Schur-bound contraction; any value strictly between zero and one works. */
#define SIMPLEX_PDHG_STEP_CONTRACTION .99
/* Dimensionless balance controls; the step product is unchanged. */
#define SIMPLEX_PDHG_WEIGHT_MIN 1.
#define SIMPLEX_PDHG_WEIGHT_MAX 1e6
#define SIMPLEX_PDHG_WEIGHT_CHANGE 4.
#define SIMPLEX_PDHG_MANUAL_WEIGHT_MIN 1e-12
#define SIMPLEX_PDHG_MANUAL_WEIGHT_MAX 1e12

/* Inner accuracy is deliberately tighter than the outer acceptance test. */
#define SIMPLEX_PAN_INNER_TOLERANCE_FRACTION 1e-3
#define SIMPLEX_PAN_RESTORE_TARGET_FRACTION .25
#define SIMPLEX_PAN_FEASIBLE_MARGIN .5
#define SIMPLEX_PAN_REPAIR_REDUCTION .5
/* Objective progress in units of dual_tolerance * sum |c_j x_j|.
 * These are search policies; each accepted point is audited separately. */
#define SIMPLEX_PAN_SMALL_GAIN_UNITS 100.
#define SIMPLEX_PAN_FACE_SMALL_GAIN_UNITS 10.
#define SIMPLEX_PAN_FACE_LARGE_GAIN_UNITS 1000.
#define SIMPLEX_PAN_FACE_MAX_FRACTION .01
#define SIMPLEX_PAN_FACE_EXPANSION 4.
#define SIMPLEX_PAN_FACE_CONTRACTION .5

/* Relative NNLS block candidate quality. */
#define SIMPLEX_PAN_NNLS_BLOCK_SCORE_FRACTION .25
/* Fraction of remaining feasibility margin usable by a trial. */
#define SIMPLEX_PAN_TRIAL_MARGIN_FRACTION .25
/* Restore before all of the feasibility allowance is spent. */
#define SIMPLEX_PAN_RESTORE_TRIGGER_FRACTION .75
/* Face repair objective and primal tolerance allowance. */
#define SIMPLEX_PAN_FACE_COST_FRACTION .125
/* Repair subproblem accuracy relative to primal tolerance. */
#define SIMPLEX_PAN_RESTORE_INNER_FRACTION .1
/* Stricter inner accuracy after an unsuccessful audited direction. */
#define SIMPLEX_PAN_INNER_RETRY_FACTOR .01
/* Legacy support compression mass budget; not a zero test. */
#define SIMPLEX_PAN_SUPPORT_MASS_FRACTION .001
/* Legacy native trial gain floor, relative to 1+|cost|. */
#define SIMPLEX_PAN_TRIAL_GAIN_RELATIVE 1e-10
/* Legacy feasible incumbent gain floor, relative to 1+|cost|. */
#define SIMPLEX_PDHG_INCUMBENT_GAIN_RELATIVE 1e-12
/* Sufficient-decrease fraction for the backtracking line search. */
#define SIMPLEX_PAN_FIRST_ORDER_ARMIJO 1e-4
/* Backtracking step contraction. */
#define SIMPLEX_PAN_FIRST_ORDER_BACKTRACK .5
/* Legacy first-order plateau test against checkpoint loss. */
#define SIMPLEX_PAN_FIRST_ORDER_PROGRESS_RATIO .99999
/* Safety margin on first-order Lipschitz estimate. */
#define SIMPLEX_PAN_LIPSCHITZ_MARGIN 1.01
/* Maximum backtracking trials per direction. */
#define SIMPLEX_PAN_FIRST_ORDER_LINE_LIMIT 20
/* Maximum work in a first-order Newton CG proposal. */
#define SIMPLEX_PAN_CG_WORK_CAP 256
/* First independent-column repair cap. */
#define SIMPLEX_PAN_RESTORE_RANK_INITIAL 64
/* Second independent-column repair cap. */
#define SIMPLEX_PAN_RESTORE_RANK_SECOND 256
/* Last independent-column repair cap. */
#define SIMPLEX_PAN_RESTORE_RANK_MAX 2048
/* Finite point-fingerprint history for implicit NNLS cycling. */
#define SIMPLEX_PAN_IMPLICIT_HISTORY 64

/* Legacy zero-column floor for face ordering only; it cannot certify a ray. */
#define SIMPLEX_PAN_FACE_NORM_FLOOR 1e-100
#define SIMPLEX_PAN_RATIO_TIE_RELATIVE 1e-12

/* Bounded finite repair attempts and optional proposal windows. */
#define SIMPLEX_PAN_RESTORE_ATTEMPTS 16
#define SIMPLEX_PAN_TRIAL_SPARSE_PERIOD 16
#define SIMPLEX_PAN_FEASIBILITY_PROBE_WORK 1024
#define SIMPLEX_PAN_BULK_COLUMNS 256
#define SIMPLEX_PAN_STALLED_STEPS 16
#define SIMPLEX_PAN_INNER_RETRIES 2
#define SIMPLEX_PAN_SUPPORT_ROW_FRACTION 4
#define SIMPLEX_PAN_SUPPORT_ROW_SLACK 32
#define SIMPLEX_PAN_NEWTON_FIRST_WORK 512
#define SIMPLEX_PAN_NEWTON_TRIAL_WORK 32
#define SIMPLEX_PAN_NEWTON_PERIOD 128
#define SIMPLEX_PAN_FIRST_ORDER_AUDIT_PERIOD 16
#define SIMPLEX_PAN_FIRST_ORDER_CHECK_PERIOD 128
#define SIMPLEX_PAN_PACKED_REBUILD_COLUMNS 32

#endif
