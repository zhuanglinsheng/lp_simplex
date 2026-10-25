/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_dual_kernel.c
 * @brief Dual iteration, numerical safeguards and terminal certification.
 */

#include "simplex_execution.h"
#include "policy/factorization.h"
#include "policy/search.h"
#include "policy/controls.h"
#include "simplex_dual.h"
#include "simplex_dual_kernel.h"
#include "simplex_stage.h"
#include "simplex_progress.h"
#include "simplex_dual_bounds.h"
#include "simplex_dual_internal.h"
#include "simplex_dual_pricing.h"
#include "simplex_dual_primal_pricing.h"
#include "simplex_dual_state.h"
#include "linalg.h"
#include "utils.h"

#include <lp_simplex/status.h>

#include <math.h>
#include <float.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

/* A resource refusal belongs to this owned factor operation. It cannot
 * authorize another factor attempt inside the same kernel invocation. */
static int dual_factorize(struct simplex_DualState *state)
{
	if (simplex_basis_storage_failed(&state->factor))
		return lp_simplex_EXIT_FAILURE;
	return simplex_basis_factorize(&state->factor);
}

static int dual_basis_failure_status(
	const struct simplex_DualState *state, int numerical_status)
{
	return simplex_basis_storage_failed(&state->factor)
		       ? lp_simplex_MemoryAllocError
		       : numerical_status;
}

enum dual_Step {
	DUAL_STEP_FAILED = -1,
	DUAL_STEP_RESTART = 0,
	DUAL_STEP_READY = 1,
	DUAL_STEP_FINISHED = 2,
	DUAL_STEP_RETRY = 3
};

/* Values that belong to one pivot attempt.  Keeping them out of DualState
 * makes the lifetime of transient iteration data explicit. */
struct dual_Iteration {
	int p;
	int q;
	int kappa;
	double target;
	double primal_infeasibility;
	double theta;
	double delta;
	double new_value;
	double merit_before;
};

static int dual_is_finite_lower(double value)
{
	return value > __lp_simplex_NINF__;
}

static int dual_is_finite_upper(double value)
{
	return value < __lp_simplex_INF__;
}

/* Separate the initial dual face before symmetry has selected a fragile
 * basis.  Unlike later feasibility-restoring shifts, this is a bounded,
 * one-shot change to the original objective, removed by the usual cleanup. */
static void dual_perturb_initial_cost(
	struct simplex_DualState *state, const int recovery)
{
	int j, shifted = 0;

	if (!state->controls.initial_cost_perturb && !recovery)
		return;
	for (j = 0; j < state->structural; j++) {
		unsigned int order = (unsigned int)j + 1U;
		double sign, margin, cost;

		if (state->position[j] >= 0 ||
			(state->status[j] != DUAL_STATUS_LOWER &&
				state->status[j] != DUAL_STATUS_UPPER))
			continue;
		order ^= order >> 16;
		order *= 0x7feb352dU;
		order ^= order >> 15;
		order *= 0x846ca68bU;
		order ^= order >> 16;
		sign = state->status[j] == DUAL_STATUS_LOWER ? 1. : -1.;
		margin = 5. * state->options->dual_tolerance *
			 (1. + __lp_simplex_ABS__(state->original_cost[j])) *
			 (1. + (double)order / 4294967296.);
		cost = nextafter(state->cost[j] + sign * margin,
			sign > 0. ? HUGE_VAL : -HUGE_VAL);
		if (!isfinite(cost) || !isfinite(margin))
			continue;
		state->cost[j] = cost;
		shifted++;
	}
	if (shifted > 0) {
		state->costs_version++;
		state->cost_shift_active = 1;
	}
	if (state->profile.enabled)
		fprintf(stderr,
			"dual profile: initial cost perturb columns=%d\n",
			shifted);
}

static int dual_pan_structure_enabled(const struct simplex_DualState *state)
{
	return state->rows >= SIMPLEX_DUAL_STRUCTURAL_PREFERENCE_ROWS &&
	       state->structural_basic < state->rows;
}

static void dual_change_status(struct simplex_DualState *state,
	const int variable, const unsigned char new_status)
{
	simplex_degeneracy_update_status(&state->degeneracy, variable,
		state->status[variable], new_status);

	if (state->recovery.primal_tracking)
		simplex_degeneracy_update_status(
			&state->recovery.primal_history, variable,
			state->status[variable], new_status);
	if (state->crash_tracking_active)
		simplex_degeneracy_update_status(&state->crash_history,
			variable, state->status[variable], new_status);
	state->status[variable] = new_status;
}

static void dual_change_basis(struct simplex_DualState *state,
	const int position, const int new_variable)
{
	simplex_degeneracy_update_basis(&state->degeneracy, position,
		state->basis[position], new_variable);

	if (state->recovery.primal_tracking)
		simplex_degeneracy_update_basis(&state->recovery.primal_history,
			position, state->basis[position], new_variable);
	if (state->crash_tracking_active)
		simplex_degeneracy_update_basis(&state->crash_history, position,
			state->basis[position], new_variable);
	state->basis[position] = new_variable;
}

static unsigned char dual_status_for_cost(
	const struct simplex_DualState *state, const int variable)
{
	if (state->lower[variable] == state->upper[variable])
		return DUAL_STATUS_FIXED;

	if (state->cost[variable] >= 0. &&
		dual_is_finite_lower(state->lower[variable]))
		return DUAL_STATUS_LOWER;
	if (state->cost[variable] < 0. &&
		dual_is_finite_upper(state->upper[variable]))
		return DUAL_STATUS_UPPER;
	if (dual_is_finite_lower(state->lower[variable]))
		return DUAL_STATUS_LOWER;
	if (dual_is_finite_upper(state->upper[variable]))
		return DUAL_STATUS_UPPER;
	return DUAL_STATUS_FREE;
}

static void dual_set_nonbasic_value(
	struct simplex_DualState *state, const int variable)
{
	switch (state->status[variable]) {
	case DUAL_STATUS_LOWER:
	case DUAL_STATUS_FIXED:
		state->value[variable] = state->lower[variable];
		break;
	case DUAL_STATUS_UPPER:
		state->value[variable] = state->upper[variable];
		break;
	case DUAL_STATUS_FREE:
		state->value[variable] = 0.;
		break;
	}
}

static void dual_initialize_basis(
	struct simplex_DualState *state, const struct simplex_Problem *problem)
{
	int i, j;

	state->lex_initialized = 0;
	for (j = 0; j < state->variables; j++)
		state->position[j] = -1;
	for (j = 0; j < problem->columns; j++) {
		state->status[j] = dual_status_for_cost(state, j);
		dual_set_nonbasic_value(state, j);
	}
	for (i = 0; i < problem->rows; i++) {
		int variable = problem->columns + i;

		state->basis[i] = variable;
		state->position[variable] = i;
		state->status[variable] = DUAL_STATUS_BASIC;
		state->value[variable] = 0.;
	}
}

static double dual_column_dot(const struct simplex_DualState *state,
	const int variable, const double *vector)
{
	if (variable >= state->structural)
		return -vector[variable - state->structural];

	return simplex_csc_column_dot(&state->matrix, variable, vector);
}

static int dual_compute_reduced_costs_impl(
	struct simplex_DualState *state, const int refined)
{
	int i, j;

	for (i = 0; i < state->rows; i++)
		state->pi[i] = state->cost[state->basis[i]];
	if ((refined ? simplex_basis_btran_refined(&state->factor, state->pi)
		     : simplex_basis_btran(&state->factor, state->pi)) ==
		lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	for (j = 0; j < state->variables; j++)
		state->reduced[j] =
			state->position[j] >= 0
				? 0.
				: state->cost[j] -
					  dual_column_dot(state, j, state->pi);
	return lp_simplex_EXIT_SUCCESS;
}

static int dual_compute_reduced_costs(struct simplex_DualState *state)
{
	return dual_compute_reduced_costs_impl(state, 0);
}

static double dual_variable_infeasibility(
	const struct simplex_DualState *state, const int variable)
{
	double reduced = state->reduced[variable];

	if (!isfinite(reduced))
		return HUGE_VAL;
	switch (state->status[variable]) {
	case DUAL_STATUS_LOWER:
		return reduced < 0. ? -reduced : 0.;
	case DUAL_STATUS_UPPER:
		return reduced > 0. ? reduced : 0.;
	case DUAL_STATUS_FREE:
		return __lp_simplex_ABS__(reduced);
	default:
		return 0.;
	}
}

static int dual_bound_status_for_reduced(const struct simplex_DualState *state,
	const int variable, const double reduced, unsigned char *status)
{
	if (state->lower[variable] == state->upper[variable]) {
		*status = DUAL_STATUS_FIXED;
		return 1;
	}
	if (reduced >= -state->options->dual_tolerance &&
		dual_is_finite_lower(state->lower[variable])) {
		*status = DUAL_STATUS_LOWER;
		return 1;
	}
	if (reduced <= state->options->dual_tolerance &&
		dual_is_finite_upper(state->upper[variable])) {
		*status = DUAL_STATUS_UPPER;
		return 1;
	}
	if (__lp_simplex_ABS__(reduced) <= state->options->dual_tolerance &&
		!dual_is_finite_lower(state->lower[variable]) &&
		!dual_is_finite_upper(state->upper[variable])) {
		*status = DUAL_STATUS_FREE;
		return 1;
	}
	return 0;
}

static int dual_crash_finish(struct simplex_DualState *state, const int result)
{
	state->crash_tracking_active = 0;
	state->lex_initialized = 0;
	return result;
}

/* Build a dual-feasible crash basis without changing the immutable CSC model. */
static int dual_crash(struct simplex_DualState *state)
{
	int attempt, limit = SIMPLEX_DUAL_CRASH_ATTEMPT_LIMIT(
			     state->rows, state->variables);

	simplex_degeneracy_reset_state(&state->crash_history);
	state->crash_tracking_active = 1;
	for (attempt = 0; attempt < limit; attempt++) {
		int candidate, candidates = 0, i, p = -1, q = -1;
		unsigned char leaving_status = DUAL_STATUS_FREE;

		if (simplex_degeneracy_record_state(&state->crash_history,
			    state->basis, state->rows, state->status,
			    state->variables)) {
			if (state->profile.enabled)
				fprintf(stderr,
					"dual profile: crash repeated state attempt=%d\n",
					attempt);
			return dual_crash_finish(
				state, lp_simplex_EXIT_FAILURE);
		}
		if (dual_compute_reduced_costs(state) ==
			lp_simplex_EXIT_FAILURE)
			return dual_crash_finish(
				state, lp_simplex_EXIT_FAILURE);
		/* Snapshot the infeasible nonbasics once for this basis.  The previous
		 * rejection loop repeatedly rescanned every variable from index zero. */
		for (i = 0; i < state->variables; i++)
			if (state->position[i] < 0 &&
				dual_variable_infeasibility(state, i) >
					state->options->dual_tolerance)
				state->candidate_index[candidates++] = i;
		if (attempt == 0 && state->profile.enabled)
			fprintf(stderr,
				"dual profile: crash initial infeasible variables=%d\n",
				candidates);
		if (candidates == 0)
			return dual_crash_finish(
				state, lp_simplex_EXIT_SUCCESS);

		for (candidate = 0; candidate < candidates; candidate++) {
			double best_pivot = 0.;

			q = state->candidate_index[candidate];
			simplex_csc_column_to_dense(&state->matrix,
				state->structural, q, state->direction);
			if (simplex_basis_ftran(
				    &state->factor, state->direction) ==
				lp_simplex_EXIT_FAILURE)
				return dual_crash_finish(
					state, lp_simplex_EXIT_FAILURE);
			p = -1;
			for (i = 0; i < state->rows; i++) {
				int leaving = state->basis[i];
				double pivot = state->direction[i];
				double leaving_reduced;
				unsigned char proposed;

				if (__lp_simplex_ABS__(pivot) <=
					state->options->pivot_tolerance)
					continue;
				leaving_reduced = -state->reduced[q] / pivot;
				if (!dual_bound_status_for_reduced(state,
					    leaving, leaving_reduced,
					    &proposed))
					continue;
				if (__lp_simplex_ABS__(pivot) > best_pivot) {
					best_pivot = __lp_simplex_ABS__(pivot);
					p = i;
					leaving_status = proposed;
				}
			}
			if (p >= 0)
				break;
		}
		if (p < 0) {
			if (state->profile.enabled)
				fprintf(stderr,
					"dual profile: crash stalled attempt=%d rejected=%d\n",
					attempt, candidates);
			return dual_crash_finish(
				state, lp_simplex_EXIT_FAILURE);
		}
		{
			int leaving = state->basis[p];

			state->structural_basic +=
				(q < state->structural) -
				(leaving < state->structural);
			state->position[leaving] = -1;
			dual_change_status(state, leaving, leaving_status);
			dual_set_nonbasic_value(state, leaving);
		}
		dual_change_basis(state, p, q);
		state->position[q] = p;
		state->feasibility_tolerance[p] = 0.;
		dual_change_status(state, q, DUAL_STATUS_BASIC);
		{
			int update = simplex_basis_update(
				&state->factor, p, state->direction);

			if (update == 1) {
				if (dual_factorize(state) ==
					lp_simplex_EXIT_FAILURE)
					return dual_crash_finish(
						state, lp_simplex_EXIT_FAILURE);
			} else if (update == lp_simplex_EXIT_FAILURE) {
				return dual_crash_finish(
					state, lp_simplex_EXIT_FAILURE);
			}
		}
	}
	if (state->profile.enabled)
		fprintf(stderr, "dual profile: crash iteration limit=%d\n",
			limit);
	return dual_crash_finish(state, lp_simplex_EXIT_FAILURE);
}

static int dual_compute_primal_values_impl(
	struct simplex_DualState *state, const int refined)
{
	int i, j;

	lp_simplex_memset(state->work, 0, (size_t)state->rows * sizeof(double));
	for (j = 0; j < state->variables; j++) {
		if (state->position[j] >= 0)
			continue;
		dual_set_nonbasic_value(state, j);
		if (state->value[j] != 0.)
			simplex_csc_column_axpy(&state->matrix,
				state->structural, j, -state->value[j],
				state->work);
	}
	if ((refined ? simplex_basis_ftran_refined(&state->factor, state->work)
		     : simplex_basis_ftran(&state->factor, state->work)) ==
		lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	for (i = 0; i < state->rows; i++) {
		state->value[state->basis[i]] = state->basic_value[i] =
			state->work[i];
		state->basic_lower[i] = state->lower[state->basis[i]];
		state->basic_upper[i] = state->upper[state->basis[i]];
	}
	return lp_simplex_EXIT_SUCCESS;
}

static int dual_compute_primal_values(struct simplex_DualState *state)
{
	return dual_compute_primal_values_impl(state, 0);
}

static int dual_compute_primal_values_refined(struct simplex_DualState *state)
{
	return dual_compute_primal_values_impl(state, 1);
}

static double dual_max_dual_infeasibility(
	const struct simplex_DualState *state);
static int dual_apply_cost_shifts(struct simplex_DualState *state);

/* Controller-level reinversion transaction.  All exact-state consumers are
 * refreshed together so the factor, primal/dual values and hypersparse
 * infeasibility queues can never describe different bases. */
static int dual_reinvert(struct simplex_DualState *state)
{
	if (dual_factorize(state) == lp_simplex_EXIT_FAILURE ||
		dual_compute_primal_values(state) == lp_simplex_EXIT_FAILURE ||
		dual_compute_reduced_costs(state) == lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;

	if (!state->primal_cleanup_active &&
		dual_max_dual_infeasibility(state) >
			state->options->dual_tolerance) {
		if (dual_compute_reduced_costs_impl(state, 1) ==
			lp_simplex_EXIT_FAILURE)
			return lp_simplex_EXIT_FAILURE;
		simplex_basis_use_product_updates(&state->factor);
		dual_apply_cost_shifts(state);
	}
	state->reinversions++;
	simplex_dual_feasibility_rebuild(state);
	/* Recovery is an event consumed by the exact-state transaction, not by
	 * individual callers.  This keeps every reinversion path consistent and
	 * records recovery only after factor, primal, dual and heap state agree. */
	if (state->pan_enabled &&
		simplex_degeneracy_take_recovery(&state->degeneracy))
		simplex_degeneracy_record_recovery(&state->degeneracy);
	return lp_simplex_EXIT_SUCCESS;
}

static int dual_update_edge_weights(
	struct simplex_DualState *state, const int p)
{
	int k;
	double pivot = state->direction[p];
	double old_p = state->edge_weight[p];

	if (__lp_simplex_ABS__(pivot) <= state->options->pivot_tolerance)
		return lp_simplex_EXIT_FAILURE;
	for (k = 0; k < state->packed_direction.count; k++) {
		int i = state->packed_direction.index[k];
		double ratio;
		double updated;

		if (i == p)
			continue;
		ratio = state->packed_direction.value[k] / pivot;
		updated = state->edge_weight[i] - 2. * ratio * state->work[i] +
			  ratio * ratio * old_p;
		state->edge_weight[i] = __lp_simplex_MAX__(
			updated, SIMPLEX_DUAL_EDGE_WEIGHT_FLOOR);
	}
	state->edge_weight[p] = __lp_simplex_MAX__(
		old_p / (pivot * pivot), SIMPLEX_DUAL_EDGE_WEIGHT_FLOOR);
	return lp_simplex_EXIT_SUCCESS;
}

static int dual_opposite_bound_crossed(const struct simplex_DualState *state,
	const int variable, const int sign, const double new_value,
	double *movement, unsigned char *new_status)
{
	if (sign > 0 && dual_is_finite_upper(state->upper[variable]) &&
		new_value > state->upper[variable] +
				    state->options->primal_tolerance) {
		*movement = state->upper[variable] - state->value[variable];
		*new_status = DUAL_STATUS_UPPER;
		return 1;
	}
	if (sign < 0 && dual_is_finite_lower(state->lower[variable]) &&
		new_value < state->lower[variable] -
				    state->options->primal_tolerance) {
		*movement = state->lower[variable] - state->value[variable];
		*new_status = DUAL_STATUS_LOWER;
		return 1;
	}
	return 0;
}

static int dual_accumulate_bound_flip(struct simplex_DualState *state,
	const int variable, const double movement)
{
	int k;

	if (variable >= state->structural)
		return simplex_sparse_vector_add(&state->flip_rhs,
			variable - state->structural, -movement);
	for (k = state->matrix.column_start[variable];
		k < state->matrix.column_start[variable + 1]; k++)
		if (simplex_sparse_vector_add(&state->flip_rhs,
			    state->matrix.row_index[k],
			    movement * state->matrix.value[k]) ==
			lp_simplex_EXIT_FAILURE)
			return lp_simplex_EXIT_FAILURE;
	return lp_simplex_EXIT_SUCCESS;
}

static int dual_flush_bound_flips(
	struct simplex_DualState *state, const int leaving_position)
{
	int i;

	if (state->flip_rhs.count == 0)
		return lp_simplex_EXIT_SUCCESS;
	lp_simplex_memcpy(state->direction, state->flip_rhs.dense,
		(size_t)state->rows * sizeof(double));
	if (simplex_basis_ftran(&state->factor, state->direction) ==
		lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	/* Replace the accumulated pricing-row estimate by the actual batched
	 * solve, including the leaving row.  Tiny alpha errors times large bound
	 * movements must not become permanent errors in the primal state. */
	state->basic_value[leaving_position] +=
		state->flip_leaving_change - state->direction[leaving_position];
	state->flip_leaving_change = 0.;
	simplex_sparse_vector_pack(
		&state->packed_direction, state->direction, 0.);
	for (i = 0; i < state->packed_direction.count; i++) {
		int row = state->packed_direction.index[i];
		/* The leaving-row value is maintained incrementally from alpha so
		 * subsequent flips can be classified before this batched FTRAN. */
		if (row != leaving_position)
			state->basic_value[row] -=
				state->packed_direction.value[i];
	}
	simplex_dual_feasibility_update_packed(state, &state->packed_direction);
	simplex_dual_feasibility_update(state, leaving_position);
	simplex_sparse_vector_clear(&state->flip_rhs);
	state->profile.flip_batches++;
	return lp_simplex_EXIT_SUCCESS;
}

static double dual_max_dual_infeasibility(const struct simplex_DualState *state)
{
	int j;
	double maximum = 0.;

	for (j = 0; j < state->variables; j++) {
		double violation;

		if (state->position[j] >= 0)
			continue;
		violation = dual_variable_infeasibility(state, j);
		if (violation > maximum)
			maximum = violation;
	}
	return maximum;
}

/* A local reduced-cost tolerance alone does not control objective accuracy:
 * a tiny wrong-sign cost on a wide bounded variable permits a large objective
 * error.  Bound widths convert those residuals into objective units.  An
 * unbounded variable retains the existing local residual test (unit weight). */
static double dual_bound_weight(
	const struct simplex_DualState *state, const int j)
{
	if (dual_is_finite_lower(state->lower[j]) &&
		dual_is_finite_upper(state->upper[j]))
		return state->upper[j] - state->lower[j];

	return 1.;
}

static int dual_bound_error_exceeds_budget(
	const struct simplex_DualState *state)
{
	long double error = 0., objective = 0.;
	int j;

	for (j = 0; j < state->variables; j++) {
		objective +=
			(long double)state->original_cost[j] * state->value[j];
		if (state->position[j] < 0)
			error += (long double)dual_variable_infeasibility(
					 state, j) *
				 dual_bound_weight(state, j);
	}
	return !isfinite(error) || !isfinite(objective) ||
	       error > state->options->dual_tolerance * (1. + fabsl(objective));
}

/* Restore dual feasibility after an exact reinversion by changing only the
 * objective coefficients of infeasible nonbasic variables.  Since basic
 * costs are untouched, pi remains valid and cost[j]-=reduced[j] places each
 * shifted reduced cost exactly at zero.  Shifts are retained across pivots and
 * removed transactionally before accepting optimality for the original LP. */
static int dual_apply_cost_shifts(struct simplex_DualState *state)
{
	double maximum_relative_shift = 0.;
	int j, shifted = 0;

	if (!simplex_recovery_cost_changes_allowed(&state->recovery))
		return 0;
	for (j = 0; j < state->variables; j++) {
		double violation;

		if (state->position[j] >= 0)
			continue;
		violation = dual_variable_infeasibility(state, j);
		if (violation <= state->options->dual_tolerance)
			continue;
		state->cost[j] -= state->reduced[j];
		state->reduced[j] = 0.;
		shifted++;
	}
	if (shifted > 0) {
		int recovery_pending = state->degeneracy.recovery_pending;

		state->lex_initialized = 0;
		simplex_degeneracy_reset_state(&state->degeneracy);
		state->degeneracy.recovery_pending = recovery_pending;
		state->costs_version++;
		state->cost_shift_active = 1;
		state->cost_shift_count += shifted;
		if (state->profile.enabled) {
			for (j = 0; j < state->variables; j++)
				maximum_relative_shift = __lp_simplex_MAX__(
					maximum_relative_shift,
					__lp_simplex_ABS__(
						state->cost[j] -
						state->original_cost[j]) /
						__lp_simplex_MAX__(1.,
							__lp_simplex_ABS__(
								state->original_cost
									[j])));
			fprintf(stderr,
				"dual profile: phase-I cost shifts=%d total=%d max-relative=%.17g\n",
				shifted, state->cost_shift_count,
				maximum_relative_shift);
		}
	}
	return shifted;
}

/* Separate a flat dual face using margins determined by the existing
 * feasibility/error budget.  No instance-size gate, random seed, or tunable
 * perturbation amplitude is needed.  Basic costs (and pi) remain unchanged;
 * the original objective is restored before primal cleanup and certification.
 *
 * Perturbation is repeatable across stagnant faces: every epoch shifts the
 * currently tied nonbasic variables whose accumulated shift is still below
 * DUAL_PERTURB_SHIFT_CAP times their base margin, and the epoch number is
 * folded into the decorrelating hash so repeated shifts of one variable stay
 * independent.  When every tied variable is saturated, tie-breaking escalates
 * to lexicographic mode rather than grinding on the face. */
static void dual_perturb_face(struct simplex_DualState *state)
{
	int j, perturbed = 0;
	/* Numerical recovery and a detected cycle take priority over a new
	 * working objective. Saturation of one candidate set alone is weaker:
	 * moving pi may expose other variables with unused perturbation budget. */
	if (!simplex_recovery_cost_changes_allowed(&state->recovery) ||
		state->degeneracy.repeated_states > 0 ||
		state->degeneracy.recovery_pending)
		return;
	if (!state->perturb_multi_epoch && state->perturb_epochs > 0)
		return;
	for (j = 0; j < state->variables; j++) {
		double scale, base, margin, shifted, sign, accumulated;
		int k;
		unsigned int order = (unsigned int)j + 1U;

		if (state->status[j] != DUAL_STATUS_LOWER &&
			state->status[j] != DUAL_STATUS_UPPER)
			continue;
		if (__lp_simplex_ABS__(state->reduced[j]) >
			state->options->dual_tolerance)
			continue;
		sign = state->status[j] == DUAL_STATUS_LOWER ? 1. : -1.;
		scale = __lp_simplex_ABS__(state->cost[j]);
		if (j < state->structural)
			for (k = state->matrix.column_start[j];
				k < state->matrix.column_start[j + 1]; k++)
				scale += __lp_simplex_ABS__(
					state->matrix.value[k] *
					state->pi[state->matrix.row_index[k]]);
		else
			scale += __lp_simplex_ABS__(
				state->pi[j - state->structural]);
		base = __lp_simplex_MAX__(
			state->options->dual_tolerance, DBL_EPSILON * scale);
		/* Each variable accumulates at most a fixed multiple of its own
		 * margin across epochs, keeping the working objective within a
		 * tolerance-derived distance of the original one. */
		accumulated = __lp_simplex_ABS__(
			state->cost[j] - state->original_cost[j]);
		if (accumulated >= DUAL_PERTURB_SHIFT_CAP * base)
			continue;
		/* A bijective integer mixer decorrelates the ordering from sparse
		 * row/column numbering (linear ranks can preserve model symmetry).
		 * The epoch keeps repeated shifts of one variable independent. */
		order ^=
			(unsigned int)(state->perturb_epochs + 1) * 0x9e3779b9U;
		order ^= order >> 16;
		order *= 0x7feb352dU;
		order ^= order >> 15;
		order *= 0x846ca68bU;
		order ^= order >> 16;
		margin = base * (1. + (double)order / 4294967296.);
		shifted = nextafter(state->cost[j] + sign * margin,
			sign > 0. ? HUGE_VAL : -HUGE_VAL);
		if (!isfinite(shifted) || !isfinite(margin))
			continue;
		state->reduced[j] += shifted - state->cost[j];
		state->cost[j] = shifted;
		perturbed++;
	}
	state->lex_initialized = 0;
	/* This is a new working objective; old face fingerprints no longer
	 * describe zero objective progress for it. */
	simplex_degeneracy_reset_state(&state->degeneracy);
	if (perturbed > 0) {
		state->perturb_epochs++;
		state->costs_version++;
		state->cost_shift_active = 1;
		state->degeneracy.phase = SIMPLEX_PAN_NORMAL;
	} else
		simplex_degeneracy_escalate_lexicographic(&state->degeneracy);
}

static int dual_restore_original_cost(struct simplex_DualState *state)
{
	simplex_recovery_seal_objective(&state->recovery);
	if (!state->cost_shift_active)
		return lp_simplex_EXIT_SUCCESS;
	lp_simplex_memcpy(state->cost, state->original_cost,
		(size_t)state->variables * sizeof(double));
	state->costs_version++;
	state->cost_shift_active = 0;
	state->lex_initialized = 0;
	simplex_degeneracy_reset_state(&state->degeneracy);
	return dual_compute_reduced_costs_impl(state, 1);
}

/* Scale a terminal Farkas gap by the cancellation in rho^T(-N x_N).
 * This distinguishes a genuine infeasibility certificate from roundoff in a
 * redundant zero-RHS row with very large primal activities. */
static double dual_certificate_scale(struct simplex_DualState *state)
{
	long double scale = 1.;
	int i, j;

	lp_simplex_memset(state->work, 0, (size_t)state->rows * sizeof(double));
	for (j = 0; j < state->variables; j++) {
		if (state->position[j] >= 0 || state->value[j] == 0.)
			continue;
		simplex_csc_column_axpy(&state->matrix, state->structural, j,
			-state->value[j], state->work);
	}
	for (i = 0; i < state->rows; i++)
		scale += fabsl((long double)state->rho[i] * state->work[i]);
	return (double)scale;
}

static void dual_update_dual_values(struct simplex_DualState *state,
	const int kappa, const double theta, const int leaving,
	const int entering)
{
	int i;
	double tau = -kappa * theta;
	clock_t started = 0;

	if (state->profile.enabled)
		started = clock();

	if (tau != 0.) {
		for (i = 0; i < state->rows; i++)
			state->pi[i] += tau * state->rho[i];
		if (state->packed_alpha_valid) {
			if (state->profile.enabled) {
				state->profile.reduced_sparse_calls++;
				state->profile.reduced_sparse_touched +=
					state->packed_alpha.count;
			}
			for (i = 0; i < state->packed_alpha.count; i++) {
				int j = state->packed_alpha.index[i];
				/* Reduced costs are maintained only for the active
				 * nonbasic set.  Basic entries must remain exactly zero. */
				if (state->nonbasic_slot[j] >= 0)
					state->reduced[j] -=
						tau *
						state->packed_alpha.value[i];
			}
		} else {
			if (state->profile.enabled) {
				state->profile.reduced_dense_calls++;
				state->profile.reduced_dense_columns +=
					state->nonbasic_count;
			}
			if (state->nonbasic_count >=
				state->variables -
					state->variables /
						SIMPLEX_DUAL_DENSE_MEMBERSHIP_DIVISOR) {
				for (i = 0; i < state->variables; i++)
					if (state->nonbasic_slot[i] >= 0)
						state->reduced[i] -=
							tau * state->alpha[i];
			} else
				for (i = 0; i < state->nonbasic_count; i++) {
					int j = state->nonbasic_index[i];

					state->reduced[j] -=
						tau * state->alpha[j];
				}
		}
	}
	state->reduced[entering] = 0.;
	state->reduced[leaving] = -tau;
	if (state->profile.enabled)
		state->profile.reduced_update_seconds +=
			(double)(clock() - started) / CLOCKS_PER_SEC;
}

static void dual_rebuild_iteration_indexes(struct simplex_DualState *state)
{
	int i;

	for (i = 0; i < state->rows; i++)
		state->edge_weight[i] = 1.;
	simplex_dual_feasibility_rebuild(state);
	lp_simplex_memset(
		state->alpha, 0, (size_t)state->variables * sizeof(double));
	lp_simplex_memset(state->candidate_sign, 0,
		(size_t)state->variables * sizeof(int));
	state->candidate_count = 0;
	simplex_dual_nonbasic_initialize(state);
	state->structural_basic = 0;
	for (i = 0; i < state->rows; i++)
		if (state->basis[i] < state->structural)
			state->structural_basic++;
}

/* Use the actual restored objective, rather than treating every positive
 * floating-point step as progress. The absolute activity bounds cancellation
 * in the two objective sums; local primal/dual tolerances cover the step. */
static long double dual_cleanup_objective(
	const struct simplex_DualState *state, long double *activity)
{
	long double objective = 0.;
	int j;
	*activity = 0.;
	for (j = 0; j < state->variables; j++) {
		long double contribution =
			(long double)state->cost[j] * state->value[j];

		objective += contribution;
		*activity += fabsl(contribution);
	}
	return objective;
}

static int dual_cleanup_perturb_bounds(
	struct simplex_DualState *state, int *status)
{
	int changed = simplex_recovery_perturb_bounds(state);

	if (changed < 0) {
		*status = lp_simplex_MemoryAllocError;
		return lp_simplex_EXIT_FAILURE;
	}
	if (changed > 0) {
		dual_rebuild_iteration_indexes(state);
		simplex_degeneracy_record_state(&state->recovery.primal_history,
			state->basis, state->rows, state->status,
			state->variables);
		if (state->profile.enabled)
			fprintf(stderr,
				"dual profile: primal bound perturb epoch=%ld iteration=%d\n",
				state->recovery.bound_epochs,
				state->iterations);
	} else
		state->recovery.bland = 1;
	return lp_simplex_EXIT_SUCCESS;
}

static int dual_cleanup_observe(struct simplex_DualState *state, int *status,
	const long double before, const long double before_activity,
	const double step, const double reduced)
{
	long double activity, after, error, rounding;
	double gain;
	int repeated, action;

	after = dual_cleanup_objective(state, &activity);
	/* Account for the number of summed terms and the platform's actual long
	 * double precision (which may equal double), rather than a fixed ulp count. */
	rounding = ((long double)state->variables + 2.) * LDBL_EPSILON;
	if (rounding >= SIMPLEX_DUAL_RATIO_NOISE_LIMIT) {
		*status = lp_simplex_PrecisionError;
		return lp_simplex_EXIT_FAILURE;
	}
	error = 2. * rounding / (1. - rounding) *
			(1. + activity + before_activity) +
		state->options->primal_tolerance * fabs(reduced) +
		state->options->dual_tolerance * fabs(step);
	gain = (double)(before - after);
	if (!isfinite(before) || !isfinite(after) || !isfinite(error) ||
		!isfinite(gain) || !isfinite((double)error)) {
		*status = lp_simplex_PrecisionError;
		return lp_simplex_EXIT_FAILURE;
	}
	if (simplex_recovery_progress_certified(gain, (double)error))
		simplex_degeneracy_reset_state(&state->recovery.primal_history);
	repeated = simplex_degeneracy_record_state(
		&state->recovery.primal_history, state->basis, state->rows,
		state->status, state->variables);
	action = simplex_recovery_observe_primal(
		&state->recovery, gain, (double)error, repeated);
	if (action == SIMPLEX_PRIMAL_PERTURB_BOUNDS)
		return dual_cleanup_perturb_bounds(state, status);
	if (action == SIMPLEX_PRIMAL_RECOVERY_FAILED) {
		*status = lp_simplex_PrecisionError;
		return lp_simplex_EXIT_FAILURE;
	}
	if (action == SIMPLEX_PRIMAL_REBUILD) {
		/* Bland can still revisit a floating-point basis. Rebuild once,
		 * retain its guard, and reject persistent repetition as numerical
		 * failure instead of consuming the remaining budget in a cycle. */
		if (dual_reinvert(state) == lp_simplex_EXIT_FAILURE) {
			*status = dual_basis_failure_status(
				state, lp_simplex_Singularity);
			return lp_simplex_EXIT_FAILURE;
		}
		if (simplex_dual_feasibility_merit(state) > 0.) {
			*status = lp_simplex_PrecisionError;
			return lp_simplex_EXIT_FAILURE;
		}
		simplex_degeneracy_reset_state(&state->recovery.primal_history);
		simplex_degeneracy_record_state(&state->recovery.primal_history,
			state->basis, state->rows, state->status,
			state->variables);
	}
	return lp_simplex_EXIT_SUCCESS;
}

static int dual_auto_checkpoint(struct simplex_DualState *state,
	double residual, int phase, double noise, long cycles);

/* Optimize the restored original objective from a primal-feasible basis.
 * Dual Phase-I cost shifts deliberately end in exactly this state: primal
 * feasibility is already available, while some original reduced costs have
 * the wrong sign.  A revised-primal cleanup is the mathematically natural
 * transition; trying to manufacture another dual-feasible crash basis throws
 * away the feasible point and can cycle on highly degenerate models. */
static int dual_primal_cleanup_impl(struct simplex_DualState *state,
	int *status, double *weight, double *ref, long *stats)
{
	int cleanup_iterations = 0;
	/* Real progress, cycle detection, one bound epoch and Bland are ordered
	 * by the shared recovery controller. Exact zeros remain a separate metric. */
	int degenerate_steps = 0;
	long cleanup_zero_steps = 0, cleanup_flips = 0;
	int use_devex = weight != NULL &&
			(state->seeded_basis || state->controls.cleanup_devex);
	int pricing_checkpoint = 0;
	long pricing_zeros = 0;
	double pricing_residual = dual_max_dual_infeasibility(state);

	stats[4] = use_devex;

	while (state->iterations < state->options->iteration_limit) {
		long double objective_before, objective_activity;
		double entering_reduced;
		int bland_mode = state->recovery.bland;
		double step = HUGE_VAL, opposite_step = HUGE_VAL;
		double relaxed_step = HUGE_VAL, best_pivot = 0.;
		double entering_violation = 0.;
		int entering = -1, entering_sign = 0;
		int polish_bounds = dual_bound_error_exceeds_budget(state);
		int i, leaving_position = -1;
		unsigned char leaving_status = DUAL_STATUS_FREE;

		objective_before =
			dual_cleanup_objective(state, &objective_activity);
		if (state->profile.enabled && cleanup_iterations > 0 &&
			(cleanup_iterations & (cleanup_iterations - 1)) == 0) {
			double objective = 0.;

			for (i = 0; i < state->variables; i++)
				objective += state->cost[i] * state->value[i];
			fprintf(stderr,
				"dual profile: primal cleanup progress=%d dual=%.6g objective=%.12g zero=%ld flips=%ld bland=%d\n",
				cleanup_iterations,
				dual_max_dual_infeasibility(state), objective,
				cleanup_zero_steps, cleanup_flips, bland_mode);
		}

		/* Devex prices by violation per approximate edge length, including
		 * polishing.  The bound-error budget still controls which violations
		 * must be resolved and when cleanup may finish.  The legacy path
		 * prices material aggregate bound errors by objective contribution.
		 * Bland mode takes the lowest-index violated variable instead of the
		 * largest violator. */
		for (i = 0; i < state->variables; i++) {
			double violation, score;

			if (state->position[i] >= 0)
				continue;
			violation = dual_variable_infeasibility(state, i);
			score = polish_bounds ? violation * dual_bound_weight(
								    state, i)
					      : violation;
			if (use_devex && !bland_mode)
				score = violation / sqrt(weight[i]);
			if (violation <= 0. || score <= entering_violation ||
				(!polish_bounds &&
					violation <=
						state->options->dual_tolerance))
				continue;
			entering = i;
			entering_violation = score;
			if (state->status[i] == DUAL_STATUS_LOWER)
				entering_sign = 1;
			else if (state->status[i] == DUAL_STATUS_UPPER)
				entering_sign = -1;
			else
				entering_sign = state->reduced[i] < 0. ? 1 : -1;
			if (bland_mode)
				break;
		}
		if (entering < 0) {
			if (state->profile.enabled)
				fprintf(stderr,
					"dual profile: primal cleanup iterations=%d\n",
					cleanup_iterations);
			return lp_simplex_EXIT_SUCCESS;
		}

		/* Zero pivots alone do not imply failure. Trial Devex only when a
		 * complete observation window also fails to reduce dual error. The
		 * basis and primal point stay live; the normal lease bounds the trial. */
		if (!use_devex && weight &&
			cleanup_iterations - pricing_checkpoint >=
				SIMPLEX_DUAL_PROGRESS_WINDOW) {
			double error = dual_max_dual_infeasibility(state);

			if (error > state->options->dual_tolerance &&
				error >= SIMPLEX_DUAL_FTRAN_PROGRESS_RATIO *
						 pricing_residual &&
				2 * (cleanup_zero_steps - pricing_zeros) >=
					cleanup_iterations -
						pricing_checkpoint) {
				use_devex = 1;
				stats[4] = 1;
				simplex_primal_devex_initialize(
					state, weight, ref);
				if (state->profile.enabled)
					fprintf(stderr,
						"auto: cleanup pricing=devex reason=stagnant-dual-and-zero-pivots iterations=%d\n",
						state->iterations);
			}
			pricing_checkpoint = cleanup_iterations;
			pricing_zeros = cleanup_zero_steps;
			pricing_residual = error;
		}

		/* Original-objective cleanup uses a dual residual, not the
		 * Phase-I feasibility merit. Sample only at completed pivots. */
		if (state->iterations - state->exploration.checkpoint >=
				SIMPLEX_DUAL_PROGRESS_WINDOW &&
			dual_auto_checkpoint(state,
				dual_max_dual_infeasibility(state),
				10 + polish_bounds + (use_devex ? 2 : 0),
				state->options->dual_tolerance,
				state->recovery.primal_history
					.repeated_states)) {
			*status = state->exploration.reason ==
						  SIMPLEX_DUAL_EXPLORE_NUMERICAL
					  ? lp_simplex_PrecisionError
					  : lp_simplex_ExceedIterLimit;
			return lp_simplex_EXIT_FAILURE;
		}

		entering_reduced = state->reduced[entering];
		simplex_csc_column_to_dense(&state->matrix, state->structural,
			entering, state->direction);
		if (simplex_basis_ftran_refined(&state->factor,
			    state->direction) == lp_simplex_EXIT_FAILURE) {
			*status = dual_basis_failure_status(
				state, lp_simplex_Singularity);
			return lp_simplex_EXIT_FAILURE;
		}
		/* Harris first pass bounds the permitted primal error.  The second
		 * pass chooses a stable pivot inside that interval, rather than a
		 * nearly singular row winning by a few rounding bits. */
		for (i = 0; i < state->rows; i++) {
			double change = -entering_sign * state->direction[i];
			double distance;

			if (__lp_simplex_ABS__(change) <=
				state->options->pivot_tolerance)
				continue;
			if (change < 0. &&
				dual_is_finite_lower(state->basic_lower[i]))
				distance = state->basic_value[i] -
					   state->basic_lower[i];
			else if (change > 0. &&
				 dual_is_finite_upper(state->basic_upper[i]))
				distance = state->basic_upper[i] -
					   state->basic_value[i];
			else
				continue;
			relaxed_step = __lp_simplex_MIN__(relaxed_step,
				(__lp_simplex_MAX__(0., distance) +
					state->options->primal_tolerance) /
					__lp_simplex_ABS__(change));
		}
		for (i = 0; i < state->rows; i++) {
			double change = -entering_sign * state->direction[i];
			double candidate;
			unsigned char candidate_status;

			if (__lp_simplex_ABS__(state->direction[i]) <=
				state->options->pivot_tolerance)
				continue;
			if (change < 0. &&
				dual_is_finite_lower(state->basic_lower[i])) {
				candidate = (state->basic_value[i] -
						    state->basic_lower[i]) /
					    -change;
				candidate_status = DUAL_STATUS_LOWER;
			} else if (change > 0. &&
				   dual_is_finite_upper(
					   state->basic_upper[i])) {
				candidate = (state->basic_upper[i] -
						    state->basic_value[i]) /
					    change;
				candidate_status = DUAL_STATUS_UPPER;
			} else
				continue;
			candidate = __lp_simplex_MAX__(0., candidate);
			if (bland_mode) {
				/* Exact minimum ratio; ties go to the lowest-index basic
				 * variable.  No Harris window and no pivot-strength
				 * preference while Bland's rule governs the cleanup. */
				if (candidate < step ||
					(candidate == step &&
						(leaving_position < 0 ||
							state->basis[i] <
								state->basis
									[leaving_position]))) {
					step = candidate;
					leaving_position = i;
					leaving_status = candidate_status;
				}
			} else if (
				candidate <= relaxed_step &&
				(__lp_simplex_ABS__(state->direction[i]) >
						best_pivot ||
					(__lp_simplex_ABS__(
						 state->direction[i]) ==
							best_pivot &&
						(leaving_position < 0 ||
							state->basis[i] <
								state->basis
									[leaving_position])))) {
				best_pivot =
					__lp_simplex_ABS__(state->direction[i]);
				step = candidate;
				leaving_position = i;
				leaving_status = candidate_status;
			}
		}
		if (entering_sign > 0 &&
			dual_is_finite_upper(state->upper[entering]))
			opposite_step = __lp_simplex_MAX__(
				0., state->upper[entering] -
					    state->value[entering]);
		else if (entering_sign < 0 &&
			 dual_is_finite_lower(state->lower[entering]))
			opposite_step = __lp_simplex_MAX__(
				0., state->value[entering] -
					    state->lower[entering]);

		if (opposite_step < step ||
			(bland_mode && isfinite(opposite_step) &&
				opposite_step == step &&
				(leaving_position < 0 ||
					entering <
						state->basis
							[leaving_position]))) {
			dual_change_status(state, entering,
				entering_sign > 0 ? DUAL_STATUS_UPPER
						  : DUAL_STATUS_LOWER);
			dual_set_nonbasic_value(state, entering);
			if (dual_compute_primal_values_refined(state) ==
				lp_simplex_EXIT_FAILURE) {
				*status = lp_simplex_PrecisionError;
				return lp_simplex_EXIT_FAILURE;
			}
			simplex_dual_feasibility_rebuild(state);
			degenerate_steps =
				opposite_step <= 0. ? degenerate_steps + 1 : 0;
			cleanup_zero_steps += opposite_step <= 0.;
			cleanup_flips++;
			stats[0]++;
			stats[1] += opposite_step <= 0.;
			stats[2] += bland_mode;
			stats[3] =
				__lp_simplex_MAX__(stats[3], degenerate_steps);
			state->iterations++;
			cleanup_iterations++;
			if (dual_cleanup_observe(state, status,
				    objective_before, objective_activity,
				    opposite_step, entering_reduced) ==
				lp_simplex_EXIT_FAILURE)
				return lp_simplex_EXIT_FAILURE;
			continue;
		}
		if (leaving_position < 0) {
			*status = lp_simplex_Unboundedness;
			return lp_simplex_EXIT_FAILURE;
		}
		if (use_devex && simplex_primal_devex_update(state, weight, ref,
					 entering, leaving_position) ==
					 lp_simplex_EXIT_FAILURE) {
			*status = lp_simplex_PrecisionError;
			return lp_simplex_EXIT_FAILURE;
		}
		{
			int leaving = state->basis[leaving_position];
			int update;

			simplex_dual_nonbasic_remove(state, entering);
			state->structural_basic +=
				(entering < state->structural) -
				(leaving < state->structural);
			state->position[leaving] = -1;
			dual_change_status(state, leaving, leaving_status);
			if (state->lower[leaving] == state->upper[leaving])
				dual_change_status(
					state, leaving, DUAL_STATUS_FIXED);
			else
				simplex_dual_nonbasic_add(state, leaving);
			dual_change_basis(state, leaving_position, entering);
			state->position[entering] = leaving_position;
			dual_change_status(state, entering, DUAL_STATUS_BASIC);
			update = simplex_basis_update(&state->factor,
				leaving_position, state->direction);
			if (update == lp_simplex_EXIT_FAILURE) {
				*status = simplex_basis_storage_failed(
						  &state->factor)
						  ? lp_simplex_MemoryAllocError
						  : lp_simplex_Singularity;
				return lp_simplex_EXIT_FAILURE;
			}
			if (update == 1) {
				if (dual_reinvert(state) ==
					lp_simplex_EXIT_FAILURE) {
					*status = dual_basis_failure_status(
						state, lp_simplex_Singularity);
					return lp_simplex_EXIT_FAILURE;
				}
			} else if (dual_compute_primal_values_refined(state) ==
					   lp_simplex_EXIT_FAILURE ||
				   dual_compute_reduced_costs(state) ==
					   lp_simplex_EXIT_FAILURE) {
				*status = lp_simplex_PrecisionError;
				return lp_simplex_EXIT_FAILURE;
			}
			/* Reinversion already rebuilds feasibility. Other exchanges
			 * need only the row queues: nonbasic membership was updated
			 * above, and dual pricing scratch is unused during cleanup. */
			if (update != 1)
				simplex_dual_feasibility_rebuild(state);
		}
		degenerate_steps = step <= 0. ? degenerate_steps + 1 : 0;
		cleanup_zero_steps += step <= 0.;
		stats[0]++;
		stats[1] += step <= 0.;
		stats[2] += bland_mode;
		stats[3] = __lp_simplex_MAX__(stats[3], degenerate_steps);
		state->iterations++;
		cleanup_iterations++;
		if (dual_cleanup_observe(state, status, objective_before,
			    objective_activity, step,
			    entering_reduced) == lp_simplex_EXIT_FAILURE)
			return lp_simplex_EXIT_FAILURE;
	}
	*status = lp_simplex_ExceedIterLimit;
	return lp_simplex_EXIT_FAILURE;
}

int simplex_dual_primal_cleanup(struct simplex_DualState *state, int *status)
{
	int result, restored;
	long stats[5] = {0, 0, 0, 0, 0};
	double *storage = NULL;
	double *weight = NULL, *ref = NULL;

	/* The stagnation-triggered Devex trial below also applies to an explicit
	 * Dual solve. Without optional scratch that recovery was unreachable and
	 * a long degenerate Dantzig cleanup could exhaust numerical stability. */
	if (state->seeded_basis || state->controls.cleanup_devex ||
		state->options->algorithm == lp_simplex_ALGORITHM_AUTO ||
		state->options->algorithm ==
			lp_simplex_ALGORITHM_DUAL_REVISED) {
		storage = (double *)lp_simplex_malloc(
			(size_t)state->variables * 2 * sizeof(double));
		if (storage == NULL) {
			if (state->controls.cleanup_devex) {
				*status = lp_simplex_MemoryAllocError;
				return lp_simplex_EXIT_FAILURE;
			}
			/* Adaptive pricing is optional; retain the live legacy path
			 * if its extra workspace cannot be allocated. */
		} else {
			weight = storage;
			ref = storage + state->variables;
			simplex_primal_devex_initialize(state, weight, ref);
		}
	}
	simplex_recovery_begin_primal(&state->recovery);
	simplex_degeneracy_reset_state(&state->degeneracy);
	state->lex_initialized = 0;
	simplex_degeneracy_record_state(&state->recovery.primal_history,
		state->basis, state->rows, state->status, state->variables);
	state->primal_cleanup_active = 1;
	/* Initialize once at the phase boundary, then maintain membership at
	 * each exchange instead of rescanning all columns after every pivot. */
	dual_rebuild_iteration_indexes(state);
	result = lp_simplex_EXIT_FAILURE;
	if (state->recovery.force_bounds && !state->recovery.bounds_attempted &&
		dual_cleanup_perturb_bounds(state, status) ==
			lp_simplex_EXIT_FAILURE)
		goto RESTORE;
	result = dual_primal_cleanup_impl(state, status, weight, ref, stats);
RESTORE:
	/* Every exit, including a pivot limit or failed factor update, restores
	 * the original domain. Only the outer dual loop may certify the result. */
	restored = simplex_recovery_restore_bounds(state);
	state->recovery.primal_tracking = 0;
	if (restored) {
		lp_simplex_memset(state->feasibility_tolerance, 0,
			(size_t)state->rows * sizeof(double));
		if (dual_compute_primal_values_refined(state) ==
				lp_simplex_EXIT_FAILURE ||
			dual_compute_reduced_costs_impl(state, 1) ==
				lp_simplex_EXIT_FAILURE) {
			*status = lp_simplex_PrecisionError;
			result = lp_simplex_EXIT_FAILURE;
		}
		simplex_dual_feasibility_rebuild(state);
		if (result == lp_simplex_EXIT_FAILURE &&
			(*status == lp_simplex_Infeasibility ||
				*status == lp_simplex_Unboundedness))
			*status = lp_simplex_PrecisionError;
		if (result == lp_simplex_EXIT_SUCCESS &&
			simplex_dual_feasibility_merit(state) > 0.)
			state->recovery.dual_repairs++;
		if (state->profile.enabled)
			fprintf(stderr,
				"dual profile: original bounds restored merit=%.12g dual=%.6g\n",
				simplex_dual_feasibility_merit(state),
				dual_max_dual_infeasibility(state));
	} else
		state->recovery.stage = SIMPLEX_RECOVERY_ORIGINAL;
	state->primal_cleanup_active = 0;
	/* Dual edge weights and symbolic anchors described a basis from before
	 * primal exchanges. Reinitialize them before the outer dual repair. */
	simplex_degeneracy_reset_state(&state->degeneracy);
	state->lex_initialized = 0;
	dual_rebuild_iteration_indexes(state);
	if (state->profile.enabled)
		fprintf(stderr,
			"dual profile: cleanup summary pricing=%s iterations=%ld "
			"zero=%ld bland=%ld max_zero_streak=%ld result=%d status=%d "
			"stalls=%ld cycles=%ld bounds=%ld/%ld repairs=%ld\n",
			stats[4] ? "devex" : "dantzig", stats[0], stats[1],
			stats[2], stats[3], result,
			result == lp_simplex_EXIT_SUCCESS ? lp_simplex_Success
							  : *status,
			state->recovery.primal_stalls,
			state->recovery.primal_cycles,
			state->recovery.bound_epochs,
			state->recovery.bound_restores,
			state->recovery.dual_repairs);
	lp_simplex_free(storage);
	return result;
}

static int dual_select_leaving(struct simplex_DualState *state, int *status,
	struct dual_Iteration *iteration)
{
	int prefer_structural;
	double dual_error;

	prefer_structural =
		state->pan_enabled &&
		simplex_degeneracy_should_probe(&state->degeneracy) &&
		dual_pan_structure_enabled(state);
	/* A lexicographic dual ratio preserves symbolic feasibility for any
	 * infeasible leaving row.  Retain steepest-edge pricing. */
	iteration->merit_before = simplex_dual_feasibility_merit(state);
	iteration->p =
		simplex_dual_feasibility_choose(state, &iteration->target,
			&iteration->kappa, &iteration->primal_infeasibility,
			prefer_structural, 0, state->pan_deferred_row);
	if (iteration->p >= 0)
		return DUAL_STEP_READY;
	if (dual_reinvert(state) == lp_simplex_EXIT_FAILURE) {
		*status = dual_basis_failure_status(
			state, lp_simplex_Singularity);
		return DUAL_STEP_FAILED;
	}
	iteration->merit_before = simplex_dual_feasibility_merit(state);
	iteration->p =
		simplex_dual_feasibility_choose(state, &iteration->target,
			&iteration->kappa, &iteration->primal_infeasibility,
			prefer_structural, 0, state->pan_deferred_row);
	if (iteration->p >= 0)
		return DUAL_STEP_READY;
	/* Certify an apparent primal optimum against the current complete basis.
	 * Iterative refinement is deliberately terminal-only: it prevents an
	 * inaccurate success result without perturbing the normal pivot path. */
	if (dual_compute_primal_values_refined(state) ==
		lp_simplex_EXIT_FAILURE) {
		*status = lp_simplex_PrecisionError;
		return DUAL_STEP_FAILED;
	}
	simplex_dual_feasibility_rebuild(state);
	iteration->merit_before = simplex_dual_feasibility_merit(state);
	iteration->p =
		simplex_dual_feasibility_choose(state, &iteration->target,
			&iteration->kappa, &iteration->primal_infeasibility,
			prefer_structural, 0, state->pan_deferred_row);
	if (iteration->p >= 0)
		return DUAL_STEP_READY;
	if (state->stop_at_primal_feasible) {
		*status = lp_simplex_CondUnsatisfied;
		return DUAL_STEP_FINISHED;
	}
	dual_error = dual_max_dual_infeasibility(state);
	if (!isfinite(dual_error)) {
		*status = lp_simplex_PrecisionError;
		return DUAL_STEP_FAILED;
	}
	if (state->cost_shift_active &&
		dual_error <= state->options->dual_tolerance) {
		if (dual_restore_original_cost(state) ==
			lp_simplex_EXIT_FAILURE) {
			*status = lp_simplex_PrecisionError;
			return DUAL_STEP_FAILED;
		}
		dual_error = dual_max_dual_infeasibility(state);
	}
	if (dual_error <= state->options->dual_tolerance &&
		!dual_bound_error_exceeds_budget(state)) {
		if (state->recovery.bounds_active) {
			*status = lp_simplex_PrecisionError;
			return DUAL_STEP_FAILED;
		}
		simplex_recovery_seal_objective(&state->recovery);
		state->recovery.stage = SIMPLEX_RECOVERY_ORIGINAL;
		*status = lp_simplex_Success;
		return DUAL_STEP_FINISHED;
	}
	/* At this point the complete basis has passed primal certification.  Any
	 * remaining wrong-sign reduced costs belong to a primal-simplex cleanup;
	 * a dual crash would destructively discard the feasible basis on failure. */
	if (simplex_dual_primal_cleanup(state, status) ==
		lp_simplex_EXIT_FAILURE)
		return DUAL_STEP_FAILED;
	dual_rebuild_iteration_indexes(state);
	return DUAL_STEP_RESTART;
}

static int dual_prepare_leaving_row(struct simplex_DualState *state,
	int *status, const struct dual_Iteration *iteration)
{
	double certificate_scale;
	double exact_weight = 0.;
	int i, candidates;

	lp_simplex_memset(state->rho, 0, (size_t)state->rows * sizeof(double));
	state->rho[iteration->p] = 1.;
	if (simplex_basis_btran(&state->factor, state->rho) ==
		lp_simplex_EXIT_FAILURE) {
		*status = dual_basis_failure_status(
			state, lp_simplex_Singularity);
		return DUAL_STEP_FAILED;
	}
	for (i = 0; i < state->rows; i++)
		exact_weight += state->rho[i] * state->rho[i];
	state->edge_weight[iteration->p] = __lp_simplex_MAX__(
		exact_weight, SIMPLEX_DUAL_EDGE_WEIGHT_FLOOR);
	candidates = simplex_dual_prepare_ratio_test(state, iteration->kappa);
	if (candidates < 0) {
		*status = state->ratio_failure_status;
		return DUAL_STEP_FAILED;
	}
	if (candidates > 0)
		return DUAL_STEP_READY;
	if (!simplex_dual_nonbasic_is_consistent(state)) {
		if (state->profile.enabled)
			fprintf(stderr,
				"dual profile: repairing nonbasic index set\n");
		simplex_dual_nonbasic_initialize(state);
		candidates = simplex_dual_prepare_ratio_test(
			state, iteration->kappa);
		if (candidates < 0) {
			*status = state->ratio_failure_status;
			return DUAL_STEP_FAILED;
		}
		if (candidates > 0)
			return DUAL_STEP_RESTART;
	}
	if (simplex_basis_update_count(&state->factor) != 0) {
		if (dual_reinvert(state) == lp_simplex_EXIT_SUCCESS)
			return DUAL_STEP_RESTART;
		if (simplex_basis_storage_failed(&state->factor)) {
			*status = lp_simplex_MemoryAllocError;
			return DUAL_STEP_FAILED;
		}
	}
	certificate_scale = dual_certificate_scale(state);
	if (state->basis[iteration->p] >= state->structural &&
		iteration->primal_infeasibility <=
			state->options->primal_tolerance * certificate_scale) {
		state->feasibility_tolerance[iteration->p] =
			iteration->primal_infeasibility;
		simplex_dual_feasibility_update(state, iteration->p);
		return DUAL_STEP_RESTART;
	}
	if (state->profile.enabled)
		fprintf(stderr,
			"dual profile: infeasible certificate row=%d basic=%d value=%.17g lower=%.17g upper=%.17g target=%.17g violation=%.3g updates=%d\n",
			iteration->p, state->basis[iteration->p],
			state->basic_value[iteration->p],
			state->basic_lower[iteration->p],
			state->basic_upper[iteration->p], iteration->target,
			iteration->primal_infeasibility,
			simplex_basis_update_count(&state->factor));
	state->certificate_row = iteration->p;
	*status = lp_simplex_Infeasibility;
	return DUAL_STEP_FAILED;
}

/* Observe progress only after a complete pivot transaction: factor, primal,
 * dual, and bound-status state must describe the same basis. */
static void dual_observe_progress(struct simplex_DualState *state,
	const struct dual_Iteration *iteration, const double theta)
{
	double pan_step, pan_error, pan_merit_error;

	state->pan_deferred_row = -1;
	state->stalled_restarts = 0;
	pan_step = theta * iteration->primal_infeasibility;
	/* First-order error propagation for theta*v, using the solver's
	 * feasibility tolerances rather than an activation-count threshold. */
	pan_error =
		state->options->dual_tolerance *
			iteration->primal_infeasibility +
		state->options->primal_tolerance * __lp_simplex_ABS__(theta) +
		state->options->dual_tolerance *
			state->options->primal_tolerance;
	pan_merit_error =
		(double)state->rows * state->options->primal_tolerance;
	state->profile.main_pivots++;
	state->profile.zero_dual_steps += theta == 0.;
	state->profile.flat_objective_steps +=
		__lp_simplex_ABS__(pan_step) <= pan_error;
	if (state->profile.enabled &&
		(state->profile.main_pivots &
			(state->profile.main_pivots - 1)) == 0)
		fprintf(stderr,
			"dual profile: main progress=%ld zero_theta=%ld flat=%ld "
			"merit=%.12g perturb_epochs=%ld cycles=%ld\n",
			state->profile.main_pivots,
			state->profile.zero_dual_steps,
			state->profile.flat_objective_steps,
			simplex_dual_feasibility_merit(state),
			state->perturb_epochs,
			state->degeneracy.repeated_states);
	if (!state->pan_enabled)
		return;
	simplex_degeneracy_observe(&state->degeneracy, pan_step, pan_error,
		iteration->merit_before, simplex_dual_feasibility_merit(state),
		pan_merit_error);
	if (!simplex_degeneracy_should_probe(&state->degeneracy))
		state->lex_initialized = 0;
	if (simplex_degeneracy_should_probe(&state->degeneracy)) {
		simplex_basis_request_compact(&state->factor);
		simplex_degeneracy_record_state(&state->degeneracy,
			state->basis, state->rows, state->status,
			state->variables);
	}
	/* Detect repetition before a cost perturbation can erase this face's
	 * history. A detected cycle retains the symbolic ratio ordering. */
	if (state->degeneracy.stagnant_face)
		dual_perturb_face(state);
}

static int dual_handle_exhausted_ratio(struct simplex_DualState *state,
	int *status, const struct dual_Iteration *iteration,
	const int bound_flips, const int rejected_relative)
{
	double certificate_scale;

	if (bound_flips > 0) {
		/* Every usable pivot was rejected after tentative flips.  Do not
		 * advance pi without a basis exchange: B^T*pi=c_B would be broken.
		 * Flush the primal state, then use the common numerical recovery. */
		if (dual_flush_bound_flips(state, iteration->p) ==
			lp_simplex_EXIT_FAILURE) {
			*status = dual_basis_failure_status(
				state, lp_simplex_Singularity);
			return DUAL_STEP_FAILED;
		}
		*status = lp_simplex_PrecisionError;
		return DUAL_STEP_FAILED;
	}
	if (rejected_relative > 0) {
		/* Candidate exhaustion is direct evidence that the approximate ratio
		 * row and the actual FTRAN pivots disagree.  Enter Pan face mode and
		 * perform one exact-pivot pass; do not tune a tolerance repeatedly. */
		if (state->pan_enabled)
			simplex_degeneracy_request_recovery(&state->degeneracy);
		if (simplex_basis_update_count(&state->factor) != 0) {
			if (dual_reinvert(state) == lp_simplex_EXIT_FAILURE) {
				*status = lp_simplex_PrecisionError;
				return DUAL_STEP_FAILED;
			}
		}
		/* Never accept an arbitrarily weak pivot.  Defer this leaving row
		 * once and let the feasibility queue expose another direction on the
		 * same face.  Returning to an already deferred row means no alternate
		 * direction exists, so normal certification below must decide. */
		if (state->pan_deferred_row != iteration->p) {
			state->pan_deferred_row = iteration->p;
			return DUAL_STEP_RESTART;
		}
		*status = lp_simplex_PrecisionError;
		return DUAL_STEP_FAILED;
	}
	if (simplex_basis_update_count(&state->factor) != 0) {
		if (dual_reinvert(state) == lp_simplex_EXIT_SUCCESS)
			return DUAL_STEP_RESTART;
		if (simplex_basis_storage_failed(&state->factor)) {
			*status = lp_simplex_MemoryAllocError;
			return DUAL_STEP_FAILED;
		}
	}
	certificate_scale = dual_certificate_scale(state);
	if (state->basis[iteration->p] >= state->structural &&
		iteration->primal_infeasibility <=
			state->options->primal_tolerance * certificate_scale) {
		state->feasibility_tolerance[iteration->p] =
			iteration->primal_infeasibility;
		simplex_dual_feasibility_update(state, iteration->p);
		return DUAL_STEP_RESTART;
	}
	if (state->profile.enabled)
		fprintf(stderr,
			"dual profile: exhausted ratio row=%d basic=%d value=%.17g lower=%.17g upper=%.17g target=%.17g violation=%.3g rejected=%d\n",
			iteration->p, state->basis[iteration->p],
			state->basic_value[iteration->p],
			state->basic_lower[iteration->p],
			state->basic_upper[iteration->p], iteration->target,
			iteration->primal_infeasibility, rejected_relative);
	state->certificate_row = iteration->p;
	*status = lp_simplex_Infeasibility;
	return DUAL_STEP_FAILED;
}

/* Ratio tests can clip a tiny wrong-sign reduced cost to zero.  Make that
 * change explicit in the working objective before a flip or basis exchange;
 * otherwise c_B and the maintained pi describe different problems, and the
 * next BTRAN can amplify a tolerance-sized discrepancy into a large error. */
static void dual_align_entering_cost(
	struct simplex_DualState *state, const struct dual_Iteration *iteration)
{
	int q = iteration->q;
	double target = -iteration->kappa * iteration->theta * state->alpha[q];
	double correction = target - state->reduced[q];

	if (correction == 0.)
		return;
	if (!simplex_recovery_cost_changes_allowed(&state->recovery)) {
		/* Clipping a primary ratio must not silently unseal the original
		 * objective. Rebuild real reduced costs after the basis transaction. */
		state->original_dual_recompute = 1;
		state->recovery.protected_costs++;
		return;
	}
	state->cost[q] += correction;
	state->reduced[q] = target;
	state->costs_version++;
	state->cost_shift_active = 1;
	state->cost_shift_count++;
}

static int dual_try_batched_bound_flip(struct simplex_DualState *state,
	int *status, const struct dual_Iteration *iteration)
{
	double tentative_delta, tentative_value, movement = 0.;
	unsigned char flipped_status = DUAL_STATUS_FREE;
	int sign;

	if (state->ratio_active_count <= 1 || state->numerically_stressed ||
		state->regular_columns ||
		__lp_simplex_ABS__(state->alpha[iteration->q]) <=
			state->options->pivot_tolerance)
		return DUAL_STEP_READY;
	tentative_delta =
		(state->basic_value[iteration->p] - iteration->target) /
		state->alpha[iteration->q];
	tentative_value = state->value[iteration->q] + tentative_delta;
	sign = state->candidate_sign[iteration->q];
	if (sign * tentative_delta < -state->options->primal_tolerance ||
		state->status[iteration->q] == DUAL_STATUS_FREE ||
		!dual_opposite_bound_crossed(state, iteration->q, sign,
			tentative_value, &movement, &flipped_status) ||
		movement == 0.)
		return DUAL_STEP_READY;
	if (dual_accumulate_bound_flip(state, iteration->q, movement) ==
		lp_simplex_EXIT_FAILURE) {
		*status = lp_simplex_PrecisionError;
		return DUAL_STEP_FAILED;
	}
	dual_align_entering_cost(state, iteration);
	state->flip_leaving_change += state->alpha[iteration->q] * movement;
	state->basic_value[iteration->p] -=
		state->alpha[iteration->q] * movement;
	state->value[iteration->q] += movement;
	dual_change_status(state, iteration->q, flipped_status);
	simplex_dual_ratio_remove(state, iteration->q);
	state->profile.bound_flips++;
	return DUAL_STEP_RETRY;
}

static int dual_build_entering_direction(struct simplex_DualState *state,
	int *status, struct dual_Iteration *iteration, int *weight_work_ready,
	int *rejected_relative, int *bound_flipped)
{
	double direction_maximum = 0.;
	double movement = 0.;
	unsigned char flipped_status = DUAL_STATUS_FREE;
	int i, sign;
	*bound_flipped = 0;
	simplex_csc_column_to_dense(&state->matrix, state->structural,
		iteration->q, state->direction);
	if (!*weight_work_ready) {
		lp_simplex_memcpy(state->work, state->rho,
			(size_t)state->rows * sizeof(double));
		if (simplex_basis_ftran_pair(&state->factor, state->direction,
			    state->work) == lp_simplex_EXIT_FAILURE) {
			state->profile.rejected_ftran++;
			simplex_dual_ratio_remove(state, iteration->q);
			return DUAL_STEP_RETRY;
		}
		*weight_work_ready = 1;
	} else if (simplex_basis_ftran(&state->factor, state->direction) ==
		   lp_simplex_EXIT_FAILURE) {
		state->profile.rejected_ftran++;
		simplex_dual_ratio_remove(state, iteration->q);
		return DUAL_STEP_RETRY;
	}
	simplex_sparse_vector_pack_maximum(&state->packed_direction,
		state->direction, 0., &direction_maximum);
	if (!isfinite(direction_maximum)) {
		*status = lp_simplex_PrecisionError;
		return DUAL_STEP_FAILED;
	}
	if (__lp_simplex_ABS__(state->direction[iteration->p]) <=
			state->options->pivot_tolerance ||
		__lp_simplex_ABS__(state->direction[iteration->p]) <
			DUAL_RELATIVE_PIVOT_TOLERANCE * direction_maximum) {
		simplex_dual_ratio_remove(state, iteration->q);
		(*rejected_relative)++;
		state->profile.rejected_relative++;
		return DUAL_STEP_RETRY;
	}
	if (__lp_simplex_ABS__(state->direction[iteration->p] -
			       state->alpha[iteration->q]) >
		100. * state->options->pivot_tolerance *
			__lp_simplex_MAX__(
				1., __lp_simplex_ABS__(
					    state->alpha[iteration->q]))) {
		if (state->pan_enabled)
			simplex_degeneracy_request_recovery(&state->degeneracy);
		if (simplex_basis_update_count(&state->factor) != 0 &&
			dual_reinvert(state) == lp_simplex_EXIT_SUCCESS) {
			return DUAL_STEP_RESTART;
		}
		*status = lp_simplex_PrecisionError;
		return DUAL_STEP_FAILED;
	}
	iteration->delta =
		(state->basic_value[iteration->p] - iteration->target) /
		state->direction[iteration->p];
	sign = state->candidate_sign[iteration->q];
	iteration->new_value = state->value[iteration->q] + iteration->delta;
	/* Keep a final pivot: a new basic variable may be outside its bounds in
	 * dual simplex.  Flipping the last candidate instead would leave no
	 * basis exchange to make the dual step consistent with c_B. */
	if (sign * iteration->delta >= -state->options->primal_tolerance &&
		(state->ratio_active_count == 1 ||
			state->status[iteration->q] == DUAL_STATUS_FREE ||
			!dual_opposite_bound_crossed(state, iteration->q, sign,
				iteration->new_value, &movement,
				&flipped_status)))
		return DUAL_STEP_READY;
	if (state->status[iteration->q] == DUAL_STATUS_FREE || movement == 0.) {
		simplex_dual_ratio_remove(state, iteration->q);
		return DUAL_STEP_RETRY;
	}
	dual_align_entering_cost(state, iteration);
	for (i = 0; i < state->packed_direction.count; i++) {
		int row = state->packed_direction.index[i];

		state->basic_value[row] -=
			state->packed_direction.value[i] * movement;
	}
	simplex_dual_feasibility_update_packed(state, &state->packed_direction);
	state->value[iteration->q] += movement;
	dual_change_status(state, iteration->q, flipped_status);
	simplex_dual_ratio_remove(state, iteration->q);
	state->profile.bound_flips++;
	*bound_flipped = 1;
	return DUAL_STEP_RETRY;
}

static int dual_candidate_would_flip(const struct simplex_DualState *state,
	const struct dual_Iteration *iteration)
{
	int q = iteration->q;
	double movement;
	double delta = (state->basic_value[iteration->p] - iteration->target) /
		       state->alpha[q];
	unsigned char status;

	return state->candidate_sign[q] * delta >=
		       -state->options->primal_tolerance &&
	       dual_opposite_bound_crossed(state, q, state->candidate_sign[q],
		       state->value[q] + delta, &movement, &status);
}

static int dual_select_entering(struct simplex_DualState *state, int *status,
	struct dual_Iteration *iteration)
{
	int bound_flips = 0;
	int rejected_relative = 0;
	int weight_work_ready = 0;

	state->flip_leaving_change = 0.;
	for (;;) {
		int bound_flipped, step;

		iteration->q = simplex_dual_next_ratio_candidate(state,
			&iteration->theta,
			state->pan_enabled && simplex_degeneracy_should_probe(
						      &state->degeneracy),
			state->basis[iteration->p], iteration->kappa);
		state->ratio_minimum_valid = 0;
		if (iteration->q == SIMPLEX_DUAL_RATIO_ERROR) {
			*status = state->ratio_failure_status;
			return DUAL_STEP_FAILED;
		}
		if (iteration->q < 0) {
			step = dual_handle_exhausted_ratio(state, status,
				iteration, bound_flips, rejected_relative);
			if (step == DUAL_STEP_RETRY)
				continue;
			return step;
		}
		/* A Harris pivot may lie beyond an earlier breakpoint.  Flipping it
		 * first and later pivoting at a smaller theta leaves its reduced cost
		 * on the wrong side of the new bound.  Switch to monotone order before
		 * committing any flip, including the non-batched FTRAN path. */
		if (!state->ratio_ordered &&
			dual_candidate_would_flip(state, iteration)) {
			state->ratio_ordered = 1;
			continue;
		}
		step = dual_try_batched_bound_flip(state, status, iteration);
		if (step == DUAL_STEP_FAILED)
			return step;
		if (step == DUAL_STEP_RETRY) {
			bound_flips++;
			continue;
		}
		if (dual_flush_bound_flips(state, iteration->p) ==
			lp_simplex_EXIT_FAILURE) {
			*status = dual_basis_failure_status(
				state, lp_simplex_Singularity);
			return DUAL_STEP_FAILED;
		}
		step = dual_build_entering_direction(state, status, iteration,
			&weight_work_ready, &rejected_relative, &bound_flipped);
		if (step == DUAL_STEP_RETRY) {
			state->ratio_ordered = 1;
			if (bound_flipped)
				bound_flips++;
			continue;
		}
		return step;
	}
}

static void dual_exchange_basis(
	struct simplex_DualState *state, const struct dual_Iteration *iteration)
{
	int leaving = state->basis[iteration->p];

	dual_align_entering_cost(state, iteration);
	dual_update_dual_values(state, iteration->kappa, iteration->theta,
		leaving, iteration->q);
	simplex_dual_nonbasic_remove(state, iteration->q);
	state->structural_basic += (iteration->q < state->structural) -
				   (leaving < state->structural);
	state->value[leaving] = iteration->target;
	dual_change_status(state, leaving,
		iteration->kappa > 0 ? DUAL_STATUS_LOWER : DUAL_STATUS_UPPER);
	if (state->lower[leaving] == state->upper[leaving])
		dual_change_status(state, leaving, DUAL_STATUS_FIXED);
	else
		simplex_dual_nonbasic_add(state, leaving);
	state->position[leaving] = -1;
	dual_change_basis(state, iteration->p, iteration->q);
	state->basic_value[iteration->p] = iteration->new_value;
	state->basic_lower[iteration->p] = state->lower[iteration->q];
	state->basic_upper[iteration->p] = state->upper[iteration->q];
	state->position[iteration->q] = iteration->p;
	state->feasibility_tolerance[iteration->p] = 0.;
	dual_change_status(state, iteration->q, DUAL_STATUS_BASIC);
	state->alpha[iteration->q] = 0.;
	state->value[iteration->q] = iteration->new_value;
}

static int dual_update_factorization(struct simplex_DualState *state,
	const int p, const int compact_was_active, int *recomputed_dual)
{
	int compact_desired, update;

	compact_desired =
		simplex_basis_compact_requested(&state->factor) &&
		state->rows >= SIMPLEX_DUAL_STRUCTURAL_PREFERENCE_ROWS &&
		state->structural_basic * 3 <= state->rows;
	if (!compact_was_active && compact_desired)
		update = 1;
	else
		update = simplex_basis_update_packed(&state->factor, p,
			state->direction, state->packed_direction.index,
			state->packed_direction.count);
	if (update == 1) {
		update = dual_reinvert(state);
		if (update == lp_simplex_EXIT_SUCCESS)
			*recomputed_dual = 1;
	}
	return update;
}

static int dual_validate_recomputed_dual(
	struct simplex_DualState *state, int *status, const int recomputed_dual)
{
	if (recomputed_dual && dual_max_dual_infeasibility(state) >
				       10. * state->options->dual_tolerance) {
		*status = lp_simplex_PrecisionError;
		return lp_simplex_EXIT_FAILURE;
	}
	return lp_simplex_EXIT_SUCCESS;
}

static int dual_commit_iteration(struct simplex_DualState *state, int *status,
	const struct dual_Iteration *iteration)
{
	int compact_was_active, i, recomputed_dual = 0;

	for (i = 0; i < state->packed_direction.count; i++) {
		int row = state->packed_direction.index[i];

		if (row != iteration->p)
			state->basic_value[row] -=
				state->packed_direction.value[i] *
				iteration->delta;
	}
	if (dual_update_edge_weights(state, iteration->p) ==
		lp_simplex_EXIT_FAILURE) {
		*status = lp_simplex_PrecisionError;
		return lp_simplex_EXIT_FAILURE;
	}
	dual_exchange_basis(state, iteration);
	simplex_dual_feasibility_update_packed(state, &state->packed_direction);
	state->iterations++;
	compact_was_active = simplex_basis_compact_active(&state->factor);
	if (dual_update_factorization(state, iteration->p, compact_was_active,
		    &recomputed_dual) == lp_simplex_EXIT_FAILURE) {
		*status = simplex_basis_storage_failed(&state->factor)
				  ? lp_simplex_MemoryAllocError
				  : lp_simplex_Singularity;
		return lp_simplex_EXIT_FAILURE;
	}
	if (state->original_dual_recompute) {
		state->original_dual_recompute = 0;
		if (dual_compute_reduced_costs_impl(state, 1) ==
			lp_simplex_EXIT_FAILURE) {
			*status = lp_simplex_PrecisionError;
			return lp_simplex_EXIT_FAILURE;
		}
		recomputed_dual = 1;
	}
	if (dual_validate_recomputed_dual(state, status, recomputed_dual) ==
		lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	dual_observe_progress(state, iteration, iteration->theta);
	return lp_simplex_EXIT_SUCCESS;
}

/* Scheduling audits immutable kernel input at lease boundaries, independent
 * of display callbacks. A basic value is authoritative; value[] may still
 * contain its last nonbasic value between pivot commits. */
static unsigned int dual_observe_input(struct simplex_DualState *s, int phase)
{
	struct simplex_StageReport *r = &s->input_observation;
	struct simplex_Observation samples[5];
	unsigned int improvements = 0;
	int i, slot = phase == 0 ? 0 : phase - 9;
	double noise[5];
	static const unsigned int flag[5] = {SIMPLEX_INPUT_PRIMAL_PROGRESS,
		SIMPLEX_INPUT_COUNT_PROGRESS, SIMPLEX_INPUT_OBJECTIVE_PROGRESS,
		SIMPLEX_INPUT_BOUND_PROGRESS, SIMPLEX_INPUT_DUAL_PROGRESS};

	if (!s->progress_problem)
		return 0;
	if (slot < 0 || slot > 4)
		slot = 0;
	for (i = 0; i < s->rows; i++)
		s->value[s->basis[i]] = s->basic_value[i];
	simplex_stage_observe(s->progress_problem, s->options, s->value, s->pi,
		-1, s->iterations, r);
	samples[0] = r->feasibility_merit;
	samples[1] = r->violation_count;
	samples[2] = r->feasible_objective;
	samples[3] = r->lower_bound;
	samples[4] = r->dual;
	noise[0] = ((double)s->rows + s->structural) *
		   s->options->primal_tolerance;
	noise[1] = SIMPLEX_DUAL_COUNT_NOISE;
	noise[2] = s->options->dual_tolerance * fabs(samples[2].value);
	noise[3] = s->options->dual_tolerance * fabs(samples[3].value);
	noise[4] = s->options->dual_tolerance;
	for (i = 0; i < 5; i++) {
		int change = simplex_stage_progress_observe(
			&s->input_progress[slot][i], samples + i, s->iterations,
			SIMPLEX_DUAL_MIN_LEASE_REDUCTION, noise[i]);
		if (change == SIMPLEX_OBSERVATION_IMPROVED)
			improvements |= flag[i];
		else if (change == SIMPLEX_OBSERVATION_BASELINE && i == 2)
			improvements |= SIMPLEX_INPUT_FEASIBLE_AVAILABLE;
		else if (change == SIMPLEX_OBSERVATION_BASELINE && i == 3)
			improvements |= SIMPLEX_INPUT_BOUND_AVAILABLE;
	}
	r->input_improvements = improvements;
	if (improvements & (SIMPLEX_INPUT_FEASIBLE_AVAILABLE |
				   SIMPLEX_INPUT_OBJECTIVE_PROGRESS |
				   SIMPLEX_INPUT_BOUND_AVAILABLE |
				   SIMPLEX_INPUT_BOUND_PROGRESS))
		simplex_progress_observe(s->progress_problem, s->value, s->pi,
			-1., "dual-input-audit", s->iterations);
	return improvements;
}

/* Leases are renewed inside this live state. No basis or factor is lost
 * at a checkpoint. The caller's iteration_limit is the hard exploration cap. */
static int dual_auto_checkpoint(struct simplex_DualState *s, double residual,
	int phase, double noise, long cycles)
{
	int old_limit, yield, j;
	unsigned int input_progress;
	double secondary, secondary_noise;

	if (s->options->algorithm != lp_simplex_ALGORITHM_AUTO ||
		s->stop_at_primal_feasible)
		return 0;
	/* Material Phase-I cost repairs temporarily invalidate monotone merit
	 * comparisons. Give that recovery a finite initial allowance; symbolic
	 * tie perturbations alone do not qualify, and epochs never reset it. */
	if (phase == 0 && s->cost_shift_count > 0 &&
		!s->exploration.recovery_until)
		s->exploration.recovery_until = SIMPLEX_DUAL_RECOVERY_ALLOWANCE;
	if (s->exploration.initialized && s->exploration.phase == phase &&
		s->iterations < s->exploration.limit)
		return 0;
	old_limit = s->exploration.limit;
	input_progress = dual_observe_input(s, phase);
	secondary = (double)s->feasibility.count;
	secondary_noise = SIMPLEX_DUAL_COUNT_NOISE;
	if (phase != 0) {
		long double objective = 0.;

		for (j = 0; j < s->variables; j++)
			objective +=
				(long double)s->original_cost[j] * s->value[j];
		secondary = (double)objective;
		secondary_noise =
			s->options->dual_tolerance * (1. + fabs(secondary));
	}
	yield = simplex_dual_explore_yield_pair(&s->exploration, s->iterations,
		s->options->iteration_limit, phase,
		phase == 0 ? s->bounds_version : s->costs_version, residual,
		noise, cycles, secondary, secondary_noise);
	/* Changed working data cannot buy progress by subtraction against an
	 * older problem. An independent immutable-input improvement can justify
	 * renewal; otherwise let the caller reclassify the state explicitly. */
	if (yield &&
		simplex_dual_explore_reassess(&s->exploration,
			s->options->iteration_limit, input_progress, cycles))
		yield = 0;
	if (s->profile.enabled && (yield || old_limit != s->exploration.limit))
		fprintf(stderr,
			"auto: dual-checkpoint action=%s iterations=%d lease=%d hard=%d "
			"phase=%d residual=%.12g secondary=%.12g input-progress=%u reason=%d\n",
			yield ? "yield" : "continue", s->iterations,
			s->exploration.limit, s->options->iteration_limit,
			phase, residual, secondary, input_progress,
			s->exploration.reason);
	return yield;
}

static int dual_auto_yield(struct simplex_DualState *s)
{
	double merit = simplex_dual_feasibility_merit(s);

	/* Incremental subtraction can leave a negative aggregate although every
	 * row merit is nonnegative. Rebuild from physical row violations before
	 * making a scheduling decision; this does not relax feasibility. */
	if (merit < 0.) {
		simplex_dual_feasibility_rebuild(s);
		merit = simplex_dual_feasibility_merit(s);
	}
	return dual_auto_checkpoint(s, merit, 0,
		(double)s->rows * s->options->primal_tolerance,
		s->degeneracy.repeated_states);
}

static void dual_report_progress(struct simplex_DualState *s, int force)
{
	const struct simplex_Problem *p = s->progress_problem;
	const struct simplex_ProgressSink *sink;
	double *point;
	int i;
	if (!p || !(sink = p->progress) ||
		(!force && !sink->due(sink->context)))
		return;
	point = lp_simplex_malloc((size_t)s->variables * sizeof(double));
	if (!point)
		return;
	lp_simplex_memcpy(
		point, s->value, (size_t)s->variables * sizeof(double));
	for (i = 0; i < s->rows; i++)
		point[s->basis[i]] = s->basic_value[i];
	simplex_progress_emit(p, point, s->pi, -1., "dual", s->iterations, 1);
	lp_simplex_free(point);
}

int simplex_dual_continue(struct simplex_DualState *state, int *status)
{
	/*
	 * A restart commits no pivot, so `iterations` -- and therefore the
	 * caller's limit -- cannot bound it.  The only legitimate producer is the
	 * "defer this leaving row once" rule, which fires at most once per row
	 * before the same row returns and the precision path decides.  With two
	 * or more violated rows the feasibility queue can hand the two deferred
	 * rows back and forth, so an explicit budget is required to keep the loop
	 * finite.  Exceeding it is reported, never silently absorbed.
	 */
	int restart_budget =
		state->rows > (INT_MAX - SIMPLEX_DUAL_RESTART_ALLOWANCE) /
					SIMPLEX_DUAL_RESTART_ROW_FACTOR
			? INT_MAX
			: SIMPLEX_DUAL_RESTART_ROW_FACTOR * state->rows +
				  SIMPLEX_DUAL_RESTART_ALLOWANCE;

	while (state->iterations < state->options->iteration_limit) {
		struct dual_Iteration iteration;
		int step;
		if (simplex_basis_storage_failed(&state->factor)) {
			*status = lp_simplex_MemoryAllocError;
			return lp_simplex_EXIT_FAILURE;
		}

		if ((state->iterations % SIMPLEX_DUAL_REPORT_PERIOD) == 0)
			dual_report_progress(state, 0);
		if (dual_auto_yield(state)) {
			*status = state->exploration.reason ==
						  SIMPLEX_DUAL_EXPLORE_NUMERICAL
					  ? lp_simplex_PrecisionError
					  : lp_simplex_ExceedIterLimit;
			return lp_simplex_EXIT_FAILURE;
		}
		if (state->pan_enabled && state->degeneracy.recovery_pending) {
			if (dual_reinvert(state) == lp_simplex_EXIT_FAILURE ||
				dual_max_dual_infeasibility(state) >
					10. * state->options->dual_tolerance) {
				*status = lp_simplex_PrecisionError;
				return lp_simplex_EXIT_FAILURE;
			}
		}
		simplex_sparse_vector_clear(&state->flip_rhs);
		step = dual_select_leaving(state, status, &iteration);
		if (step == DUAL_STEP_FINISHED)
			return lp_simplex_EXIT_SUCCESS;
		if (step == DUAL_STEP_FAILED)
			return lp_simplex_EXIT_FAILURE;
		if (step == DUAL_STEP_RESTART)
			goto RESTART;
		step = dual_prepare_leaving_row(state, status, &iteration);
		if (step == DUAL_STEP_FAILED)
			return lp_simplex_EXIT_FAILURE;
		if (step == DUAL_STEP_RESTART)
			goto RESTART;
		step = dual_select_entering(state, status, &iteration);
		if (step == DUAL_STEP_FAILED)
			return lp_simplex_EXIT_FAILURE;
		if (step == DUAL_STEP_RESTART)
			goto RESTART;
		if (dual_commit_iteration(state, status, &iteration) ==
			lp_simplex_EXIT_FAILURE)
			return lp_simplex_EXIT_FAILURE;
		continue;
RESTART:
		if (++state->stalled_restarts > restart_budget) {
			*status = lp_simplex_PrecisionError;
			return lp_simplex_EXIT_FAILURE;
		}
	}
	*status = lp_simplex_ExceedIterLimit;
	return lp_simplex_EXIT_FAILURE;
}

static double dual_primal_infeasibility(const struct simplex_DualState *state)
{
	int i, j;
	double maximum = 0.;

	for (j = 0; j < state->variables; j++) {
		double violation = 0.;

		if (state->value[j] < state->lower[j])
			violation = state->lower[j] - state->value[j];
		else if (state->value[j] > state->upper[j])
			violation = state->value[j] - state->upper[j];
		if (violation > maximum)
			maximum = violation;
	}
	lp_simplex_memset(state->work, 0, (size_t)state->rows * sizeof(double));
	for (j = 0; j < state->structural; j++) {
		if (state->value[j] != 0.)
			simplex_csc_column_axpy(&state->matrix,
				state->structural, j, state->value[j],
				state->work);
	}
	for (i = 0; i < state->rows; i++) {
		double residual = __lp_simplex_ABS__(
			state->work[i] - state->value[state->structural + i]);

		if (residual > maximum)
			maximum = residual;
	}
	return maximum;
}

/* After numerical failure, retain the exact B=-I/product-update recovery,
 * then allow one separated cold face. Both retain the remaining pivot
 * grant and restore the original objective before success. */
static int dual_initialize_stable_phase_one(struct simplex_DualState *state,
	const struct simplex_Problem *problem, const int separate_costs)
{
	/* A representation restart is allowed only after the primal transaction
	 * has restored every bound. Its one-epoch budget survives this restart. */
	if (state->recovery.bounds_active)
		return lp_simplex_EXIT_FAILURE;

	state->recovery.stage = SIMPLEX_RECOVERY_WORKING_DUAL;
	state->recovery.original_cost_locked = 0;
	state->recovery.primal_tracking = 0;
	state->original_dual_recompute = 0;
	simplex_degeneracy_reset_state(&state->recovery.primal_history);
	state->stable_phase_one = !separate_costs;
	state->perturb_epochs = 0;
	state->costs_version++;
	state->cost_shift_active = 0;
	state->pan_deferred_row = -1;
	state->stalled_restarts = 0;
	/* Row allowances belonged to the discarded basis, not the logical one. */
	lp_simplex_memset(state->feasibility_tolerance, 0,
		(size_t)state->rows * sizeof(double));
	lp_simplex_memcpy(state->cost, state->original_cost,
		(size_t)state->variables * sizeof(double));
	simplex_degeneracy_reset_state(&state->degeneracy);
	if (separate_costs)
		simplex_basis_use_default_updates(&state->factor);
	else
		simplex_basis_use_product_updates(&state->factor);
	dual_initialize_basis(state, problem);
	/* The failed path supplies evidence for one deterministic face change.
	 * Its costs are restored and checked before success; the restart retains
	 * the original pivot grant and cannot restart itself. */
	dual_perturb_initial_cost(state, separate_costs);
	if (dual_factorize(state) != lp_simplex_EXIT_SUCCESS ||
		(separate_costs &&
			(dual_crash(state) != lp_simplex_EXIT_SUCCESS ||
				dual_factorize(state) !=
					lp_simplex_EXIT_SUCCESS)) ||
		dual_compute_primal_values(state) != lp_simplex_EXIT_SUCCESS ||
		dual_compute_reduced_costs(state) != lp_simplex_EXIT_SUCCESS)
		return lp_simplex_EXIT_FAILURE;
	dual_apply_cost_shifts(state);
	dual_rebuild_iteration_indexes(state);
	if (state->profile.enabled)
		fprintf(stderr, "dual profile: %s Phase I at iteration=%d\n",
			separate_costs ? "separated cost" : "stable logical",
			state->iterations);
	return lp_simplex_EXIT_SUCCESS;
}

static int dual_start_stable_phase_one(struct simplex_DualState *state,
	const struct simplex_Problem *problem, const int separate_costs)
{
	struct simplex_StageReport report;
	struct lp_simplex_Result result;
	int id = simplex_execution_begin(&problem->execution,
		SIMPLEX_EXECUTION_STABLE_INITIALIZE, separate_costs,
		SIMPLEX_REASON_NUMERICAL_FAILURE,
		simplex_execution_parent(&problem->execution));
	int status = dual_initialize_stable_phase_one(
		state, problem, separate_costs);

	simplex_stage_reset(&report);
	lp_simplex_memset(&result, 0, sizeof(result));
	result.status = status ? dual_basis_failure_status(
					 state, lp_simplex_PrecisionError)
			       : lp_simplex_CondUnsatisfied;
	result.iterations = 0; /* No pivots are performed by initialization. */
	report.bounds_epoch = state->bounds_version;
	report.cost_epoch = state->costs_version;
	report.stop = status ? simplex_stage_stop_from_status(result.status)
			     : SIMPLEX_STOP_HANDOFF;
	simplex_execution_end(problem->execution.history, id, &result, &report);
	return status;
}

static int dual_import_basis(
	struct simplex_DualState *state, struct simplex_DualBasisSeed *seed)
{
	int i, j, variable, align = seed->prefer_dual;
	double primal_error, dual_error;

	seed->import_attempted = 1;
	seed->rejection = SIMPLEX_BASIS_LAYOUT_REJECTED;
	for (j = 0; j < state->variables; j++)
		state->position[j] = -1;
	for (i = 0; i < state->rows; i++) {
		variable = seed->index[i];
		if (variable < 0 || variable >= state->variables ||
			state->position[variable] >= 0)
			return lp_simplex_EXIT_FAILURE;
		state->basis[i] = variable;
		state->position[variable] = i;
	}
	for (j = 0; j < state->variables; j++) {
		unsigned char status = seed->status[j];

		if (state->position[j] >= 0) {
			if (status != DUAL_STATUS_BASIC)
				return lp_simplex_EXIT_FAILURE;
			state->status[j] = DUAL_STATUS_BASIC;
			continue;
		}
		if (!((status == DUAL_STATUS_FIXED &&
			      state->lower[j] == state->upper[j]) ||
			    (status == DUAL_STATUS_LOWER &&
				    dual_is_finite_lower(state->lower[j])) ||
			    (status == DUAL_STATUS_UPPER &&
				    dual_is_finite_upper(state->upper[j])) ||
			    (status == DUAL_STATUS_FREE &&
				    !dual_is_finite_lower(state->lower[j]) &&
				    !dual_is_finite_upper(state->upper[j]))))
			return lp_simplex_EXIT_FAILURE;
		state->status[j] = state->lower[j] == state->upper[j]
					   ? DUAL_STATUS_FIXED
					   : status;
		dual_set_nonbasic_value(state, j);
	}
	seed->rejection = SIMPLEX_BASIS_FACTORIZATION_REJECTED;
	if (dual_factorize(state)) {
		if (simplex_basis_storage_failed(&state->factor))
			seed->rejection = SIMPLEX_BASIS_RESOURCE_REJECTED;
		return lp_simplex_EXIT_FAILURE;
	}
	seed->rejection = SIMPLEX_BASIS_RESIDUAL_REJECTED;
	if (dual_compute_reduced_costs_impl(state, 1))
		return lp_simplex_EXIT_FAILURE;
	for (j = 0; j < state->variables; j++)
		if (!isfinite(state->reduced[j]))
			return lp_simplex_EXIT_FAILURE;
	/* Bounded new columns with negative reduced cost can start at their upper
	 * bound. This preserves the old dual basis and lets dual simplex repair the
	 * primal state. Check every required bound before changing any status: if an
	 * unbounded column blocks this transaction, keep the feasible primal seed. */
	if (align)
		for (j = 0; j < state->variables; j++) {
			unsigned char status;

			if (state->position[j] < 0 &&
				dual_variable_infeasibility(state, j) >
					state->options->dual_tolerance &&
				!dual_bound_status_for_reduced(
					state, j, state->reduced[j], &status)) {
				align = 0;
				break;
			}
		}
	if (align)
		for (j = 0; j < state->variables; j++) {
			unsigned char status;

			if (state->position[j] < 0 &&
				dual_variable_infeasibility(state, j) >
					state->options->dual_tolerance) {
				(void)dual_bound_status_for_reduced(
					state, j, state->reduced[j], &status);
				state->status[j] = status;
				dual_set_nonbasic_value(state, j);
			}
		}
	if (dual_compute_primal_values_refined(state))
		return lp_simplex_EXIT_FAILURE;
	for (i = 0; i < state->rows; i++)
		if (!isfinite(state->basic_value[i]))
			return lp_simplex_EXIT_FAILURE;
	dual_rebuild_iteration_indexes(state);
	/* A feasible basis admits primal cleanup when newly added columns have
	 * negative reduced costs. A seed that lost feasibility is only a hint. */
	primal_error = simplex_dual_feasibility_merit(state);
	dual_error = dual_max_dual_infeasibility(state);
	if (!isfinite(primal_error) || !isfinite(dual_error))
		return lp_simplex_EXIT_FAILURE;
	if (primal_error != 0. &&
		!(align && dual_error <= state->options->dual_tolerance)) {
		seed->rejection = SIMPLEX_BASIS_FEASIBILITY_REJECTED;
		return lp_simplex_EXIT_FAILURE;
	}
	seed->rejection = SIMPLEX_BASIS_NOT_REJECTED;
	return lp_simplex_EXIT_SUCCESS;
}

static int dual_initial_resource_exit(struct simplex_DualState *state,
	const struct simplex_Problem *problem, struct lp_simplex_Result *result)
{
	result->status = lp_simplex_MemoryAllocError;
	result->iterations = state->iterations;
	if (problem->report)
		problem->report->stop = SIMPLEX_STOP_RESOURCE;
	simplex_dual_state_destroy(state);
	return lp_simplex_EXIT_FAILURE;
}

int simplex_dual_kernel_solve(const struct simplex_Problem *problem,
	const struct lp_simplex_Options *options, double *x, double *row_dual,
	struct lp_simplex_Result *result, const int propagate_bounds,
	const int basis_backend, const int feasibility_only,
	struct simplex_DualBasisSeed *seed)
{
	struct simplex_DualState state;
	int j, initialized, status = lp_simplex_CondUnsatisfied;
	int solve_state, stable_restart_attempted = 0,
			 cost_restart_attempted = 0;
	int seeded = seed && seed->valid;
	double objective = 0.;
	clock_t profile_started = 0;

	simplex_dual_seed_begin(seed);

	profile_started = clock();
	if (simplex_dual_state_create(&state, problem, options,
		    basis_backend) == lp_simplex_EXIT_FAILURE) {
		result->status = lp_simplex_MemoryAllocError;
		if (problem->report)
			problem->report->stop = SIMPLEX_STOP_RESOURCE;
		return lp_simplex_EXIT_FAILURE;
	}
	cost_restart_attempted = state.controls.initial_cost_perturb;
	/* No seed arrays were touched by failed state allocation. Once import
	 * begins, a rejected/numerically failed attempt must export afresh. */
	if (seed)
		seed->valid = 0;
	state.progress_problem = problem;
	state.stop_at_primal_feasible = feasibility_only;
	simplex_dual_initialize_bounds(&state, problem, propagate_bounds);
	lp_simplex_memcpy(state.original_cost, state.cost,
		(size_t)state.variables * sizeof(double));
	if (seeded) {
		int imported = dual_import_basis(&state, seed);
		if (problem->report) {
			problem->report->basis_import_attempted =
				seed->import_attempted;
			problem->report->basis_rejection = seed->rejection;
		}
		seed->accepted = imported == lp_simplex_EXIT_SUCCESS;
		if (imported && simplex_basis_storage_failed(&state.factor)) {
			/* Resource refusal does not disprove the caller's retained hint. */
			seed->valid = 1;
			return dual_initial_resource_exit(
				&state, problem, result);
		}
	}
	if (seeded && seed->accepted) {
		state.seeded_basis = 1;
	} else {
		dual_initialize_basis(&state, problem);
		dual_perturb_initial_cost(&state, 0);
		if (dual_factorize(&state) != lp_simplex_EXIT_SUCCESS ||
			dual_crash(&state) != lp_simplex_EXIT_SUCCESS ||
			dual_factorize(&state) != lp_simplex_EXIT_SUCCESS ||
			dual_compute_primal_values(&state) !=
				lp_simplex_EXIT_SUCCESS ||
			dual_compute_reduced_costs(&state) !=
				lp_simplex_EXIT_SUCCESS) {
			if (simplex_basis_storage_failed(&state.factor))
				return dual_initial_resource_exit(
					&state, problem, result);
			stable_restart_attempted = 1;
			initialized =
				dual_start_stable_phase_one(&state, problem, 0);
			if (initialized != lp_simplex_EXIT_SUCCESS &&
				!cost_restart_attempted &&
				!simplex_basis_storage_failed(&state.factor)) {
				cost_restart_attempted = 1;
				initialized = dual_start_stable_phase_one(
					&state, problem, 1);
			}
			if (initialized != lp_simplex_EXIT_SUCCESS) {
				if (simplex_basis_storage_failed(&state.factor))
					return dual_initial_resource_exit(
						&state, problem, result);
				result->status = lp_simplex_PrecisionError;
				if (problem->report)
					problem->report->stop =
						SIMPLEX_STOP_NUMERICAL;
				simplex_dual_state_destroy(&state);
				return lp_simplex_EXIT_FAILURE;
			}
		} else
			dual_rebuild_iteration_indexes(&state);
	}
	if (seed && state.profile.enabled)
		fprintf(stderr, "dual profile: basis seed accepted=%d\n",
			seed->accepted);

	for (;;) {
		solve_state = simplex_dual_continue(&state, &status);
		if (simplex_basis_storage_failed(&state.factor)) {
			status = lp_simplex_MemoryAllocError;
			solve_state = lp_simplex_EXIT_FAILURE;
		}
		/* Candidate exhaustion proposes a row combination. Only a fresh
		 * check against immutable constraints and bounds may terminate.
		 * The objective, shifted or original, is not part of this proof. */
		if (status == lp_simplex_Infeasibility) {
			int verified = simplex_problem_verify_infeasibility(
				problem, state.rho);
			const char *method =
				verified ? "row-combination" : "unverified";
			if (!verified &&
				simplex_problem_verify_basis_infeasibility(
					problem, state.basis,
					state.certificate_row)) {
				verified = 1;
				method = "exact-basis-row";
			}
			if (state.profile.enabled)
				fprintf(stderr,
					"dual profile: terminal proof=%s\n",
					method);
			if (!verified) {
				status = lp_simplex_PrecisionError;
				solve_state = lp_simplex_EXIT_FAILURE;
			}
		}
		/* Each recovery representation has one grant within the caller budget. */
		if (solve_state == lp_simplex_EXIT_SUCCESS ||
			state.iterations >= options->iteration_limit ||
			(seed && seed->accepted &&
				simplex_basis_resolved_backend(&state.factor) ==
					SIMPLEX_BASIS_BACKEND_KLU) ||
			(status != lp_simplex_PrecisionError &&
				status != lp_simplex_Singularity))
			break;
		if (!stable_restart_attempted) {
			stable_restart_attempted = 1;
			if (dual_start_stable_phase_one(&state, problem, 0) ==
				lp_simplex_EXIT_SUCCESS) {
				status = lp_simplex_CondUnsatisfied;
				continue;
			}
			if (simplex_basis_storage_failed(&state.factor)) {
				status = lp_simplex_MemoryAllocError;
				break;
			}
		}
		/* A failed degenerate trajectory authorizes one cold face change.
		 * Do not repeat it or spend pivots after the caller's grant expires. */
		if (!cost_restart_attempted) {
			cost_restart_attempted = 1;
			if (dual_start_stable_phase_one(&state, problem, 1) ==
				lp_simplex_EXIT_SUCCESS) {
				status = lp_simplex_CondUnsatisfied;
				continue;
			}
			if (simplex_basis_storage_failed(&state.factor)) {
				status = lp_simplex_MemoryAllocError;
				break;
			}
		}
		status = dual_basis_failure_status(&state, status);
		break;
	}

	dual_report_progress(&state, 1);
	for (j = 0; j < state.rows; j++)
		state.value[state.basis[j]] = state.basic_value[j];
	for (j = 0; j < problem->columns; j++) {
		x[j] = state.value[j];
		objective += problem->objective[j] * x[j];
	}
	if (problem->report) {
		struct simplex_StageReport *r = problem->report;
		struct simplex_BasisProfile work;

		r->feasibility_merit =
			state.input_observation.feasibility_merit;
		r->violation_count = state.input_observation.violation_count;
		r->feasible_objective =
			state.input_observation.feasible_objective;
		r->input_improvements =
			state.input_observation.input_improvements;
		r->input_audits = state.input_observation.input_audits;
		r->input_progress_space = SIMPLEX_SPACE_INPUT;
		r->stop = simplex_stage_stop_from_status(status);
		if (status == lp_simplex_ExceedIterLimit &&
			state.iterations < options->iteration_limit) {
			if (state.exploration.reason ==
				SIMPLEX_DUAL_EXPLORE_STAGNATION)
				r->stop = SIMPLEX_STOP_STAGNATION;
			else if (state.exploration.reason ==
				 SIMPLEX_DUAL_EXPLORE_CYCLE)
				r->stop = SIMPLEX_STOP_CYCLE;
			else if (state.exploration.reason ==
				 SIMPLEX_DUAL_EXPLORE_MEASUREMENT)
				r->stop = SIMPLEX_STOP_REASSESS;
			else
				r->stop = SIMPLEX_STOP_UNKNOWN;
		}
		if (status == lp_simplex_Singularity)
			r->failure = SIMPLEX_FAILURE_FACTORIZATION;
		r->phase = state.primal_cleanup_active ? 1 : 0;
		r->bounds_epoch = state.bounds_version;
		r->cost_epoch = state.costs_version;
		r->perturbed_bounds = state.recovery.bounds_active;
		r->perturbed_costs = state.cost_shift_active;
		r->working_primal.value = dual_primal_infeasibility(&state);
		r->working_dual.value = dual_max_dual_infeasibility(&state);
		r->working_primal.quality = r->working_dual.quality =
			SIMPLEX_EVIDENCE_ESTIMATED;
		r->working_primal.space = r->working_dual.space =
			SIMPLEX_SPACE_WORKING;
		r->working_primal.work = r->working_dual.work =
			state.iterations;
		r->working_primal.epoch = r->bounds_epoch;
		r->working_dual.epoch = r->cost_epoch;
		r->repeated_states = state.degeneracy.repeated_states;
		r->rebuilds = state.reinversions;
		r->point_available = 1;
		/* No factor survives destruction. A basis hint is not a claim of
		 * either feasibility; its importer must refactor and classify. */
		r->factor_available = 0;
		r->work.cpu_seconds =
			(double)(clock() - profile_started) / CLOCKS_PER_SEC;
		simplex_basis_get_profile(&state.factor, &work);
		r->work.basis_counts_available = 1;
		r->work.factorizations = work.factor_calls;
		r->work.ftran = work.ftran_calls;
		r->work.btran = work.btran_calls;
	}
	result->status = status;
	result->iterations = state.iterations;
	result->objective = objective;
	result->primal_infeasibility = dual_primal_infeasibility(&state);
	result->dual_infeasibility = dual_max_dual_infeasibility(&state);
	if (row_dual != NULL)
		for (j = 0; j < problem->rows; j++)
			row_dual[j] = state.pi[j];
	if (seed && (status == lp_simplex_Success ||
			    status == lp_simplex_ExceedIterLimit)) {
		lp_simplex_memcpy(seed->index, state.basis,
			(size_t)state.rows * sizeof(int));
		lp_simplex_memcpy(seed->status, state.status,
			(size_t)state.variables * sizeof(unsigned char));
		simplex_dual_seed_export(seed, &problem->execution);
		if (problem->report)
			problem->report->basis_available = 1;
	}
	if (state.profile.enabled) {
		struct simplex_BasisProfile basis_profile;
		int active_rank = 0;

		simplex_basis_get_profile(&state.factor, &basis_profile);
		for (j = 0; j < state.rows; j++) {
			int variable = state.basis[j];

			if (state.basic_value[j] >
					state.lower[variable] +
						state.options
							->primal_tolerance &&
				state.basic_value[j] <
					state.upper[variable] -
						state.options->primal_tolerance)
				active_rank++;
		}
		state.profile.total_seconds =
			(double)(clock() - profile_started) /
			(double)CLOCKS_PER_SEC;
		fprintf(stderr,
			"dual profile: pricing detail setup=%.6f sparse=%.6f/%ld "
			"dense=%.6f/%ld audit=%.6f/%ld certify=%.6f/%ld "
			"scatter_entries=%.0f sparse_touched=%.0f dense_columns=%.0f "
			"audit_columns=%.0f reduced_update=%.6f sparse_updates=%ld/%.0f "
			"dense_updates=%ld/%.0f\n",
			state.profile.ratio_setup_seconds,
			state.profile.ratio_sparse_seconds,
			state.profile.ratio_sparse_calls,
			state.profile.ratio_dense_seconds,
			state.profile.ratio_dense_calls,
			state.profile.ratio_audit_seconds,
			state.profile.ratio_audit_calls,
			state.profile.ratio_certify_seconds,
			state.profile.ratio_certify_calls,
			state.profile.ratio_scatter_entries,
			state.profile.ratio_sparse_touched,
			state.profile.ratio_dense_columns,
			state.profile.ratio_audit_columns,
			state.profile.reduced_update_seconds,
			state.profile.reduced_sparse_calls,
			state.profile.reduced_sparse_touched,
			state.profile.reduced_dense_calls,
			state.profile.reduced_dense_columns);
		fprintf(stderr,
			"dual profile: main summary pivots=%ld zero_theta=%ld flat=%ld\n",
			state.profile.main_pivots,
			state.profile.zero_dual_steps,
			state.profile.flat_objective_steps);
		fprintf(stderr,
			"dual profile: recovery stage=%d locked=%d primal_stalls=%ld "
			"primal_cycles=%ld rebuilds=%ld bounds=%ld/%ld repairs=%ld protected_costs=%ld\n",
			state.recovery.stage,
			state.recovery.original_cost_locked,
			state.recovery.primal_stalls,
			state.recovery.primal_cycles,
			state.recovery.primal_rebuilds,
			state.recovery.bound_epochs,
			state.recovery.bound_restores,
			state.recovery.dual_repairs,
			state.recovery.protected_costs);
		fprintf(stderr,
			"dual profile: total=%.6f factor=%.6f/%ld ftran=%.6f/%ld "
			"btran=%.6f/%ld ratio=%.6f residual=%.6f "
			"pan_degenerate=%ld pan_activations=%ld "
			"pan_cycles=%ld pan_recovery=%ld/%ld "
			"pan_rank=%d/%d factor_core=%d/%d compact=%ld[%d,%d] "
			"eta_nnz=%ld/%ld reinvert=%ld/%ld cert=%ld/%ld refine=%ld/%ld "
			"reject=%ld/%ld flips=%ld/%ld lex=%ld/%ld perturb_epochs=%ld "
			"lex_escal=%ld\n",
			state.profile.total_seconds,
			basis_profile.factor_seconds,
			basis_profile.factor_calls, basis_profile.ftran_seconds,
			basis_profile.ftran_calls, basis_profile.btran_seconds,
			basis_profile.btran_calls, state.profile.ratio_seconds,
			state.profile.total_seconds -
				basis_profile.factor_seconds -
				basis_profile.ftran_seconds -
				basis_profile.btran_seconds -
				state.profile.ratio_seconds,
			state.degeneracy.degenerate_pivots,
			state.degeneracy.activations,
			state.degeneracy.repeated_states,
			state.degeneracy.recoveries,
			state.degeneracy.recovery_requests, active_rank,
			state.rows, basis_profile.factor_size, state.rows,
			basis_profile.compact_calls, basis_profile.compact_min,
			basis_profile.compact_max, basis_profile.eta_nonzeros,
			basis_profile.eta_slots,
			basis_profile.fill_reinversions,
			basis_profile.stability_reinversions,
			basis_profile.compact_ftran_validations,
			basis_profile.compact_btran_validations,
			basis_profile.compact_ftran_refinements,
			basis_profile.compact_btran_refinements,
			state.profile.rejected_relative,
			state.profile.rejected_ftran, state.profile.bound_flips,
			state.profile.flip_batches, state.profile.lex_choices,
			state.profile.lex_solves, state.perturb_epochs,
			state.degeneracy.lex_escalations);
	}
	simplex_dual_state_destroy(&state);
	return solve_state;
}
