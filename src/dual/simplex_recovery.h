/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_recovery.h
 * @brief Coordination of objective cleanup, primal anti-stalling and bound
 * rollback.
 */

#ifndef LP_SIMPLEX_RECOVERY_INTERNAL_H
#define LP_SIMPLEX_RECOVERY_INTERNAL_H

#include "simplex_degeneracy.h"


struct simplex_DualState;

#define SIMPLEX_PRIMAL_STALL_LIMIT 64

enum simplex_RecoveryStage {
	SIMPLEX_RECOVERY_WORKING_DUAL = 0,
	SIMPLEX_RECOVERY_PRIMAL = 1,
	SIMPLEX_RECOVERY_PERTURBED_PRIMAL = 2,
	SIMPLEX_RECOVERY_ORIGINAL_DUAL = 3,
	SIMPLEX_RECOVERY_ORIGINAL = 4
};

enum simplex_PrimalRecoveryAction {
	SIMPLEX_PRIMAL_CONTINUE = 0,
	SIMPLEX_PRIMAL_PERTURB_BOUNDS = 1,
	SIMPLEX_PRIMAL_BLAND = 2,
	SIMPLEX_PRIMAL_REBUILD = 3,
	SIMPLEX_PRIMAL_RECOVERY_FAILED = 4
};

struct simplex_Recovery {
	int stage;
	int original_cost_locked;
	int primal_tracking;
	int flat_streak;
	int bland;
	int bland_steps;
	int cycle_guard;
	int bounds_attempted;
	int bounds_active;
	int disable_bounds;
	int force_bounds;
	double *saved_bounds;
	struct simplex_DegeneracyControl primal_history;
	long primal_stalls;
	long primal_cycles;
	long primal_rebuilds;
	long bound_epochs;
	long bound_restores;
	long dual_repairs;
	long protected_costs;
};

void simplex_recovery_initialize(struct simplex_Recovery *recovery);
void simplex_recovery_destroy(struct simplex_DualState *state);
void simplex_recovery_begin_primal(struct simplex_Recovery *recovery);
void simplex_recovery_seal_objective(struct simplex_Recovery *recovery);
int simplex_recovery_cost_changes_allowed(
	const struct simplex_Recovery *recovery);
int simplex_recovery_progress_certified(double gain, double error);
int simplex_recovery_observe_primal(struct simplex_Recovery *recovery,
	double gain, double error, int repeated_state);

/* A single bounded epoch per solver state, retained across numerical restarts.
 * Only active basic bounds are moved, so
 * the current point stays feasible without moving nonbasic values.
 * Returns -1 on allocation failure, 0 when no bound changed, 1 on success. */
int simplex_recovery_perturb_bounds(struct simplex_DualState *state);

/* Restore bounds exactly, repair fixed nonbasic status/value, and release
 * the snapshot. The caller must recompute basic values and feasibility. */
int simplex_recovery_restore_bounds(struct simplex_DualState *state);

/* Internal orchestration entry: primal-feasible basis, restored original
 * costs. Never certifies optimality; returns to the outer dual loop. */
int simplex_dual_primal_cleanup(struct simplex_DualState *state, int *status);
int simplex_dual_continue(struct simplex_DualState *state, int *status);

#endif
