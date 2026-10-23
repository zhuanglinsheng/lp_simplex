/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_recovery.c
 * @brief Coordinate objective cleanup, primal anti-stalling and bound rollback.
 */

#include "simplex_recovery.h"
#include "simplex_dual_internal.h"
#include "utils.h"

#include <math.h>


void simplex_recovery_initialize(struct simplex_Recovery *recovery)
{
	lp_simplex_memset(recovery, 0, sizeof(*recovery));
}


void simplex_recovery_seal_objective(struct simplex_Recovery *recovery)
{
	recovery->original_cost_locked = 1;
	if (recovery->stage == SIMPLEX_RECOVERY_WORKING_DUAL)
		recovery->stage = SIMPLEX_RECOVERY_ORIGINAL;
}


int simplex_recovery_cost_changes_allowed(
	const struct simplex_Recovery *recovery)
{
	return !recovery->original_cost_locked &&
	       recovery->stage == SIMPLEX_RECOVERY_WORKING_DUAL;
}


void simplex_recovery_begin_primal(struct simplex_Recovery *recovery)
{
	simplex_recovery_seal_objective(recovery);
	recovery->stage = SIMPLEX_RECOVERY_PRIMAL;
	recovery->primal_tracking = 1;
	recovery->flat_streak = 0;
	recovery->bland = recovery->cycle_guard;
	recovery->bland_steps = 0;
	simplex_degeneracy_reset_state(&recovery->primal_history);
}


int simplex_recovery_progress_certified(const double gain, const double error)
{
	return isfinite(gain) && isfinite(error) && error >= 0. && gain > error;
}


int simplex_recovery_observe_primal(struct simplex_Recovery *recovery,
	const double gain, const double error, const int repeated_state)
{
	if (simplex_recovery_progress_certified(gain, error)) {
		recovery->flat_streak = 0;
		recovery->bland = recovery->cycle_guard;
		recovery->bland_steps = 0;
		return recovery->bland ? SIMPLEX_PRIMAL_BLAND
				       : SIMPLEX_PRIMAL_CONTINUE;
	}
	if (recovery->flat_streak < SIMPLEX_PRIMAL_STALL_LIMIT)
		recovery->flat_streak++;

	recovery->primal_stalls++;
	if (repeated_state) {
		recovery->primal_cycles++;
		if (recovery->cycle_guard) {
			if (recovery->primal_rebuilds > 0)
				return SIMPLEX_PRIMAL_RECOVERY_FAILED;
			recovery->primal_rebuilds++;
			return SIMPLEX_PRIMAL_REBUILD;
		}
		recovery->cycle_guard = recovery->bland = 1;
	}
	/* Repetition takes priority over numeric perturbation. A cycle guard is
	 * retained even after progress; numeric microsteps cannot turn it off. */
	if (recovery->cycle_guard)
		return SIMPLEX_PRIMAL_BLAND;
	/* A flat objective is not proof of a cycle. Limit the slower pricing
	 * probe, then let normal pricing try another path. Keep face history:
	 * an actual repeated basis/status still latches Bland permanently. */
	if (recovery->bland &&
		++recovery->bland_steps >= SIMPLEX_PRIMAL_STALL_LIMIT) {
		recovery->bland = 0;
		recovery->bland_steps = 0;
		recovery->flat_streak = 0;
		return SIMPLEX_PRIMAL_CONTINUE;
	}
	if (recovery->flat_streak >= SIMPLEX_PRIMAL_STALL_LIMIT) {
		if (!recovery->bounds_attempted && !recovery->disable_bounds)
			return SIMPLEX_PRIMAL_PERTURB_BOUNDS;
		if (!recovery->bland)
			recovery->bland_steps = 0;
		recovery->bland = 1;
		return SIMPLEX_PRIMAL_BLAND;
	}
	return SIMPLEX_PRIMAL_CONTINUE;
}


static double recovery_margin(
	const int variable, const double bound, const double tolerance)
{
	unsigned int order = (unsigned int)variable + 1U;

	order ^= order >> 16;
	order *= 0x7feb352dU;
	order ^= order >> 15;
	order *= 0x846ca68bU;
	order ^= order >> 16;
	return 5. * tolerance *
	       __lp_simplex_MAX__(1., __lp_simplex_ABS__(bound)) *
	       (1. + (double)order / 4294967296.);
}


int simplex_recovery_perturb_bounds(struct simplex_DualState *state)
{
	struct simplex_Recovery *recovery = &state->recovery;
	const double tolerance = state->options->primal_tolerance;
	int i, changed = 0;

	if (recovery->bounds_attempted || recovery->bounds_active ||
		recovery->cycle_guard || recovery->disable_bounds)
		return 0;
	recovery->bounds_attempted = 1;
	if (state->variables == 0)
		return 0;
	recovery->saved_bounds = (double *)lp_simplex_malloc(
		(size_t)state->variables * 2 * sizeof(double));
	if (recovery->saved_bounds == NULL)
		return -1;
	lp_simplex_memcpy(recovery->saved_bounds, state->lower,
		(size_t)state->variables * sizeof(double));
	lp_simplex_memcpy(recovery->saved_bounds + state->variables,
		state->upper, (size_t)state->variables * sizeof(double));
	for (i = 0; i < state->rows; i++) {
		int j = state->basis[i];
		double value = state->basic_value[i];
		double lower = state->lower[j], upper = state->upper[j];
		double distance, margin, shifted, budget;

		if (!isfinite(value))
			continue;
		if (isfinite(lower)) {
			distance = value - lower;
			margin = recovery_margin(j, lower, tolerance);
			budget = 64. * tolerance *
				 __lp_simplex_MAX__(1., fabs(lower));
			shifted = nextafter(
				__lp_simplex_MIN__(lower, value) - margin,
				-HUGE_VAL);
			if (distance <= tolerance && isfinite(shifted) &&
				lower - shifted > 0. &&
				lower - shifted <= budget) {
				state->lower[j] = shifted;
				changed++;
			}
		}
		if (isfinite(upper)) {
			distance = upper - value;
			margin = recovery_margin(j, upper, tolerance);
			budget = 64. * tolerance *
				 __lp_simplex_MAX__(1., fabs(upper));
			shifted = nextafter(
				__lp_simplex_MAX__(upper, value) + margin,
				HUGE_VAL);
			if (distance <= tolerance && isfinite(shifted) &&
				shifted - upper > 0. &&
				shifted - upper <= budget) {
				state->upper[j] = shifted;
				changed++;
			}
		}
		state->basic_lower[i] = state->lower[j];
		state->basic_upper[i] = state->upper[j];
	}
	if (changed == 0) {
		lp_simplex_free(recovery->saved_bounds);
		recovery->saved_bounds = NULL;
		return 0;
	}
	recovery->bounds_active = 1;
	recovery->stage = SIMPLEX_RECOVERY_PERTURBED_PRIMAL;
	recovery->bound_epochs++;
	recovery->flat_streak = 0;
	recovery->bland = 0;
	/* Basis/status hashes from the unperturbed domain cannot detect cycles
	 * in this new problem. Numeric recovery alone does not clear history. */
	simplex_degeneracy_reset_state(&recovery->primal_history);
	return 1;
}


int simplex_recovery_restore_bounds(struct simplex_DualState *state)
{
	struct simplex_Recovery *recovery = &state->recovery;
	int j, i;

	if (!recovery->bounds_active)
		return 0;
	lp_simplex_memcpy(state->lower, recovery->saved_bounds,
		(size_t)state->variables * sizeof(double));
	lp_simplex_memcpy(state->upper,
		recovery->saved_bounds + state->variables,
		(size_t)state->variables * sizeof(double));
	for (j = 0; j < state->variables; j++) {
		if (state->position[j] >= 0)
			continue;
		if (state->lower[j] == state->upper[j])
			state->status[j] = DUAL_STATUS_FIXED;
		if (state->status[j] == DUAL_STATUS_LOWER ||
			state->status[j] == DUAL_STATUS_FIXED)
			state->value[j] = state->lower[j];
		else if (state->status[j] == DUAL_STATUS_UPPER)
			state->value[j] = state->upper[j];
	}
	for (i = 0; i < state->rows; i++) {
		state->basic_lower[i] = state->lower[state->basis[i]];
		state->basic_upper[i] = state->upper[state->basis[i]];
	}
	lp_simplex_free(recovery->saved_bounds);
	recovery->saved_bounds = NULL;
	recovery->bounds_active = 0;
	recovery->bound_restores++;
	recovery->stage = SIMPLEX_RECOVERY_ORIGINAL_DUAL;
	simplex_degeneracy_reset_state(&recovery->primal_history);
	simplex_degeneracy_reset_state(&state->degeneracy);
	state->lex_initialized = 0;
	return 1;
}


void simplex_recovery_destroy(struct simplex_DualState *state)
{
	simplex_recovery_restore_bounds(state);
	simplex_degeneracy_destroy(&state->recovery.primal_history);
}
