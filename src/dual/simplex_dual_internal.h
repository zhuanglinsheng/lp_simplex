/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_dual_internal.h
 * @brief Internal dual revised-simplex state and kernel interfaces.
 */

#ifndef LP_SIMPLEX_DUAL_STATE_INTERNAL_H
#define LP_SIMPLEX_DUAL_STATE_INTERNAL_H

#include "simplex_basis.h"
#include "simplex_csc.h"
#include "simplex_degeneracy.h"
#include "simplex_recovery.h"
#include "simplex_dual_feasibility.h"
#include "simplex_sparse_vector.h"

#include <lp_simplex/solve.h>


#define DUAL_STATUS_BASIC 0
#define DUAL_STATUS_LOWER 1
#define DUAL_STATUS_UPPER 2
#define DUAL_STATUS_FIXED 3
#define DUAL_STATUS_FREE 4


struct simplex_DualProfile {
	int enabled;
	double ratio_seconds;
	double total_seconds;
	long rejected_relative;
	long rejected_ftran;
	long bound_flips;
	long flip_batches;
	long lex_solves;
	long lex_choices;
	long main_pivots;
	long zero_dual_steps;
	long flat_objective_steps;
};


/* Owners for the dense solver workspaces.  The state keeps typed aliases for
 * hot-path readability; destruction only releases these four roots. */
struct simplex_DualWorkspace {
	int *index;
	unsigned char *status;
	double *column;
	double *row;
};


/* Shared only by the dual orchestration modules; not part of the public API. */
struct simplex_DualState {
	int rows;
	int structural;
	const struct simplex_Problem *progress_problem;
	int variables;
	int iterations;
	struct simplex_DualExplore exploration;
	int reinversions;
	int structural_basic;
	const struct lp_simplex_Options *options;
	struct simplex_CscMatrix matrix;
	struct simplex_Basis factor;
	struct simplex_DualWorkspace workspace;
	int *basis;
	int *position;
	unsigned char *status;
	double *lower;
	double *upper;
	double *cost;
	double *original_cost;
	double *value;
	double *basic_value;
	double *basic_lower;
	double *basic_upper;
	double *feasibility_tolerance;
	double *reduced;
	double *pi;
	double *rho;
	struct simplex_SparseVector packed_rho;
	struct simplex_SparseVector packed_alpha;
	int packed_alpha_valid;
	struct simplex_SparseVector packed_direction;
	struct simplex_SparseVector flip_rhs;
	double flip_leaving_change;
	double *alpha;
	double *breakpoint;
	double *direction;
	double *work;
	double *edge_weight;
	struct simplex_DualFeasibility feasibility;
	int *candidate_sign;
	int *candidate_index;
	int candidate_count;
	int ratio_active_count;
	int *nonbasic_index;
	int *nonbasic_slot;
	int nonbasic_count;
	struct simplex_DegeneracyControl degeneracy;
	struct simplex_DegeneracyControl crash_history;
	struct simplex_Recovery recovery;
	int original_dual_recompute;
	int crash_tracking_active;
	int pan_enabled;
	int stable_candidate_order;
	int regular_columns;
	int numerically_stressed;
	int modal_column_degree;
	struct simplex_DualProfile profile;
	double ratio_minimum;
	int ratio_minimum_valid;
	int ratio_ordered;
	int ratio_heap_valid;
	int pricing_validation_countdown;
	int cost_shift_active;
	int cost_shift_count;
	int stable_phase_one;
	/* Number of deterministic face-perturbation epochs applied to the
	 * working objective.  Perturbation is repeatable: each epoch shifts only
	 * variables whose accumulated shift is still below the per-variable cap. */
	long perturb_epochs;
	/* Diagnostic gate (environment): LP_SIMPLEX_SINGLE_PERTURB_EPOCH restores
	 * the legacy one-shot face perturbation. */
	int perturb_multi_epoch;
	int primal_cleanup_active;
	int stop_at_primal_feasible;
	/* Symbolic cost perturbations are anchored to the nonbasic bound signs
	 * when lexicographic recovery starts.  Storage is allocated on demand. */
	void *lex_storage;
	double *lex_rho;
	int *lex_candidate;
	int *lex_slot;
	unsigned char *lex_anchor_status;
	int lex_initialized;
	int ratio_failure_status;
	/* A leaving row whose complete Harris set had no stable actual FTRAN
	 * pivot is skipped once so Pan can move along another face direction. */
	int pan_deferred_row;
	/* Consecutive step restarts that committed no pivot.  Restarts do not
	 * count as iterations, so this budget is what stops a cycle between two
	 * alternating deferred rows from spinning forever inside dual_run(). */
	int stalled_restarts;
};

#endif
