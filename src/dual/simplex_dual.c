/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
/* Bounded dual revised simplex over an immutable CSC constraint matrix. */
#include "simplex_dual.h"
#include "simplex_scaling.h"
#include "simplex_dual_bounds.h"
#include "simplex_dual_internal.h"
#include "simplex_dual_pricing.h"
#include "simplex_dual_state.h"
#include "linalg.h"
#include "utils.h"

#include <lp_simplex/status.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>


#define DUAL_RELATIVE_PIVOT_TOLERANCE 1e-8

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


static int dual_pan_structure_enabled(const struct simplex_DualState *state)
{
	return state->rows >= 4096 && state->structural_basic < state->rows;
}


static void dual_change_status(
		struct simplex_DualState *state, const int variable,
		const unsigned char new_status)
{
	simplex_degeneracy_update_status(&state->degeneracy, variable,
		state->status[variable], new_status);
	if (state->crash_tracking_active)
		simplex_degeneracy_update_status(&state->crash_history, variable,
			state->status[variable], new_status);
	state->status[variable] = new_status;
}


static void dual_change_basis(
		struct simplex_DualState *state, const int position,
		const int new_variable)
{
	simplex_degeneracy_update_basis(&state->degeneracy, position,
		state->basis[position], new_variable);
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


static double dual_column_dot(
		const struct simplex_DualState *state,
		const int variable, const double *vector)
{
	if (variable >= state->structural)
		return -vector[variable - state->structural];
	return simplex_csc_column_dot(&state->matrix, variable, vector);
}


static int dual_compute_reduced_costs(struct simplex_DualState *state)
{
	int i, j;
	for (i = 0; i < state->rows; i++)
		state->pi[i] = state->cost[state->basis[i]];
	if (simplex_basis_btran(&state->factor, state->pi) == lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	for (j = 0; j < state->variables; j++)
		state->reduced[j] = state->cost[j] -
			dual_column_dot(state, j, state->pi);
	for (i = 0; i < state->rows; i++)
		state->reduced[state->basis[i]] = 0.;
	return lp_simplex_EXIT_SUCCESS;
}


static double dual_variable_infeasibility(
		const struct simplex_DualState *state, const int variable)
{
	double reduced = state->reduced[variable];
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


static int dual_bound_status_for_reduced(
		const struct simplex_DualState *state, const int variable,
		const double reduced, unsigned char *status)
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


static int dual_crash_finish(
		struct simplex_DualState *state, const int result)
{
	state->crash_tracking_active = 0;
	return result;
}


/* Build a dual-feasible crash basis without changing the immutable CSC model. */
static int dual_crash(
		struct simplex_DualState *state, const int largest_violation)
{
	int attempt, limit = 4 * state->variables + state->rows;
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
			return dual_crash_finish(state, lp_simplex_EXIT_FAILURE);
		}
		if (dual_compute_reduced_costs(state) == lp_simplex_EXIT_FAILURE)
			return dual_crash_finish(state, lp_simplex_EXIT_FAILURE);
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
			return dual_crash_finish(state, lp_simplex_EXIT_SUCCESS);
		/* Repair the largest dual violation first.  Index-order crash selection
		 * can exchange a long sequence of nearly feasible variables and revisit
		 * the same degenerate bases (MAROS/SHARE1B).  Moving only the best item to
		 * the front preserves the stable fallback order when that column has no
		 * admissible leaving variable. */
		if (largest_violation) {
			int best = 0;
			double best_violation = dual_variable_infeasibility(
				state, state->candidate_index[0]);
			for (i = 1; i < candidates; i++) {
				double violation = dual_variable_infeasibility(
					state, state->candidate_index[i]);
				if (violation > best_violation) {
					best = i;
					best_violation = violation;
				}
			}
			if (best != 0) {
				int variable = state->candidate_index[0];
				state->candidate_index[0] = state->candidate_index[best];
				state->candidate_index[best] = variable;
			}
		}
		for (candidate = 0; candidate < candidates; candidate++) {
			double best_pivot = 0.;
			q = state->candidate_index[candidate];
			simplex_csc_column_to_dense(&state->matrix, state->structural,
						    q, state->direction);
			if (simplex_basis_ftran(&state->factor, state->direction) ==
			    lp_simplex_EXIT_FAILURE)
				return dual_crash_finish(state, lp_simplex_EXIT_FAILURE);
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
				if (!dual_bound_status_for_reduced(state, leaving,
							   leaving_reduced, &proposed))
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
			return dual_crash_finish(state, lp_simplex_EXIT_FAILURE);
		}
		{
			int leaving = state->basis[p];
			state->structural_basic +=
				(q < state->structural) - (leaving < state->structural);
			state->position[leaving] = -1;
			dual_change_status(state, leaving, leaving_status);
			dual_set_nonbasic_value(state, leaving);
		}
		dual_change_basis(state, p, q);
		state->position[q] = p;
		state->feasibility_tolerance[p] = 0.;
		dual_change_status(state, q, DUAL_STATUS_BASIC);
		{
			int update = simplex_basis_update(&state->factor, p,
						  state->direction);
			if (update == 1) {
				if (simplex_basis_factorize(&state->factor) ==
				    lp_simplex_EXIT_FAILURE)
					return dual_crash_finish(
						state, lp_simplex_EXIT_FAILURE);
			} else if (update == lp_simplex_EXIT_FAILURE) {
				return dual_crash_finish(state, lp_simplex_EXIT_FAILURE);
			}
		}
	}
	if (state->profile.enabled)
		fprintf(stderr, "dual profile: crash iteration limit=%d\n", limit);
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
			simplex_csc_column_axpy(&state->matrix, state->structural,
						  j, -state->value[j], state->work);
	}
	if ((refined
	     ? simplex_basis_ftran_refined(&state->factor, state->work)
	     : simplex_basis_ftran(&state->factor, state->work)) ==
	    lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	for (i = 0; i < state->rows; i++) {
		state->value[state->basis[i]] = state->basic_value[i] = state->work[i];
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


/* Controller-level reinversion transaction.  All exact-state consumers are
 * refreshed together so the factor, primal/dual values and hypersparse
 * infeasibility queues can never describe different bases. */
static int dual_reinvert(struct simplex_DualState *state)
{
	if (simplex_basis_factorize(&state->factor) == lp_simplex_EXIT_FAILURE ||
	    dual_compute_primal_values(state) == lp_simplex_EXIT_FAILURE ||
	    dual_compute_reduced_costs(state) == lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
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
		state->edge_weight[i] = __lp_simplex_MAX__(updated, 1e-12);
	}
	state->edge_weight[p] = __lp_simplex_MAX__(old_p / (pivot * pivot), 1e-12);
	return lp_simplex_EXIT_SUCCESS;
}


static int dual_opposite_bound_crossed(
		const struct simplex_DualState *state, const int variable,
		const int sign, const double new_value, double *movement,
		unsigned char *new_status)
{
	if (sign > 0 && dual_is_finite_upper(state->upper[variable]) &&
	    new_value > state->upper[variable] + state->options->primal_tolerance) {
		*movement = state->upper[variable] - state->value[variable];
		*new_status = DUAL_STATUS_UPPER;
		return 1;
	}
	if (sign < 0 && dual_is_finite_lower(state->lower[variable]) &&
	    new_value < state->lower[variable] - state->options->primal_tolerance) {
		*movement = state->lower[variable] - state->value[variable];
		*new_status = DUAL_STATUS_LOWER;
		return 1;
	}
	return 0;
}


static int dual_accumulate_bound_flip(
		struct simplex_DualState *state, const int variable,
		const double movement)
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
	simplex_dual_feasibility_update_packed(
		state, &state->packed_direction);
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


/* Restore dual feasibility after an exact reinversion by changing only the
 * objective coefficients of infeasible nonbasic variables.  Since basic
 * costs are untouched, pi remains valid and cost[j]-=reduced[j] places each
 * shifted reduced cost exactly at zero.  Shifts are retained across pivots and
 * removed transactionally before accepting optimality for the original LP. */
static int dual_apply_cost_shifts(struct simplex_DualState *state)
{
	double maximum_relative_shift = 0.;
	int j, shifted = 0;
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
		state->cost_shift_active = 1;
		state->cost_shift_count += shifted;
		if (state->profile.enabled) {
			for (j = 0; j < state->variables; j++)
				maximum_relative_shift = __lp_simplex_MAX__(
					maximum_relative_shift,
					__lp_simplex_ABS__(state->cost[j] -
						state->original_cost[j]) /
					__lp_simplex_MAX__(1.,
						__lp_simplex_ABS__(state->original_cost[j])));
			fprintf(stderr,
				"dual profile: phase-I cost shifts=%d total=%d max-relative=%.17g\n",
				shifted, state->cost_shift_count, maximum_relative_shift);
		}
	}
	return shifted;
}


static int dual_restore_original_cost(struct simplex_DualState *state)
{
	if (!state->cost_shift_active)
		return lp_simplex_EXIT_SUCCESS;
	lp_simplex_memcpy(state->cost, state->original_cost,
		(size_t)state->variables * sizeof(double));
	state->cost_shift_active = 0;
	return dual_compute_reduced_costs(state);
}


/* Scale a terminal Farkas gap by the cancellation in rho^T(-N x_N).
 * This distinguishes a genuine infeasibility certificate from roundoff in a
 * redundant zero-RHS row with very large primal activities. */
static double dual_certificate_scale(struct simplex_DualState *state)
{
	long double scale = 1.;
	int i, j;
	lp_simplex_memset(state->work, 0,
		(size_t)state->rows * sizeof(double));
	for (j = 0; j < state->variables; j++) {
		if (state->position[j] >= 0 || state->value[j] == 0.)
			continue;
		simplex_csc_column_axpy(&state->matrix, state->structural,
			j, -state->value[j], state->work);
	}
	for (i = 0; i < state->rows; i++)
		scale += fabsl((long double)state->rho[i] * state->work[i]);
	return (double)scale;
}


static void dual_update_dual_values(
		struct simplex_DualState *state, const int kappa,
		const double theta, const int leaving, const int entering)
{
	int i;
	double tau = -kappa * theta;
	if (tau != 0.) {
		for (i = 0; i < state->rows; i++)
			state->pi[i] += tau * state->rho[i];
		if (state->packed_alpha_valid) {
			for (i = 0; i < state->packed_alpha.count; i++) {
				int j = state->packed_alpha.index[i];
				/* Reduced costs are maintained only for the active
				 * nonbasic set.  Basic entries must remain exactly zero. */
				if (state->nonbasic_slot[j] >= 0)
					state->reduced[j] -=
						tau * state->packed_alpha.value[i];
			}
		} else {
			for (i = 0; i < state->nonbasic_count; i++) {
				int j = state->nonbasic_index[i];
				state->reduced[j] -= tau * state->alpha[j];
			}
		}
	}
	state->reduced[entering] = 0.;
	state->reduced[leaving] = -tau;
}


static void dual_rebuild_iteration_indexes(struct simplex_DualState *state)
{
	int i;
	for (i = 0; i < state->rows; i++)
		state->edge_weight[i] = 1.;
	simplex_dual_feasibility_rebuild(state);
	lp_simplex_memset(state->alpha, 0,
		(size_t)state->variables * sizeof(double));
	lp_simplex_memset(state->candidate_sign, 0,
		(size_t)state->variables * sizeof(int));
	state->candidate_count = 0;
	simplex_dual_nonbasic_initialize(state);
	state->structural_basic = 0;
	for (i = 0; i < state->rows; i++)
		if (state->basis[i] < state->structural)
			state->structural_basic++;
}


/* Optimize the restored original objective from a primal-feasible basis.
 * Dual Phase-I cost shifts deliberately end in exactly this state: primal
 * feasibility is already available, while some original reduced costs have
 * the wrong sign.  A revised-primal cleanup is the mathematically natural
 * transition; trying to manufacture another dual-feasible crash basis throws
 * away the feasible point and can cycle on highly degenerate models. */
static int dual_primal_cleanup(
		struct simplex_DualState *state, int *status)
{
	int cleanup_iterations = 0;
	while (state->iterations < state->options->iteration_limit) {
		double step = HUGE_VAL, opposite_step = HUGE_VAL;
		double entering_violation = 0.;
		int entering = -1, entering_sign = 0;
		int i, leaving_position = -1;
		unsigned char leaving_status = DUAL_STATUS_FREE;
		if (state->profile.enabled && cleanup_iterations > 0 &&
		    (cleanup_iterations & (cleanup_iterations - 1)) == 0)
			fprintf(stderr,
				"dual profile: primal cleanup progress=%d dual=%.6g\n",
				cleanup_iterations, dual_max_dual_infeasibility(state));

		/* Price by original reduced-cost violation.  This is the bounded-primal
		 * analogue of Dantzig pricing; index order resolves exact ties. */
		for (i = 0; i < state->variables; i++) {
			double violation;
			if (state->position[i] >= 0)
				continue;
			violation = dual_variable_infeasibility(state, i);
			if (violation <= state->options->dual_tolerance ||
			    violation <= entering_violation)
				continue;
			entering = i;
			entering_violation = violation;
			if (state->status[i] == DUAL_STATUS_LOWER)
				entering_sign = 1;
			else if (state->status[i] == DUAL_STATUS_UPPER)
				entering_sign = -1;
			else
				entering_sign = state->reduced[i] < 0. ? 1 : -1;
		}
		if (entering < 0) {
			if (state->profile.enabled)
				fprintf(stderr,
					"dual profile: primal cleanup iterations=%d\n",
					cleanup_iterations);
			return lp_simplex_EXIT_SUCCESS;
		}

		simplex_csc_column_to_dense(&state->matrix, state->structural,
			entering, state->direction);
		if (simplex_basis_ftran(&state->factor, state->direction) ==
		    lp_simplex_EXIT_FAILURE) {
			*status = lp_simplex_Singularity;
			return lp_simplex_EXIT_FAILURE;
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
					state->basic_lower[i]) / -change;
				candidate_status = DUAL_STATUS_LOWER;
			} else if (change > 0. &&
				   dual_is_finite_upper(state->basic_upper[i])) {
				candidate = (state->basic_upper[i] -
					state->basic_value[i]) / change;
				candidate_status = DUAL_STATUS_UPPER;
			} else
				continue;
			candidate = __lp_simplex_MAX__(0., candidate);
			if (candidate < step ||
			    (candidate == step && leaving_position >= 0 &&
			     state->basis[i] < state->basis[leaving_position])) {
				step = candidate;
				leaving_position = i;
				leaving_status = candidate_status;
			}
		}
		if (entering_sign > 0 &&
		    dual_is_finite_upper(state->upper[entering]))
			opposite_step = __lp_simplex_MAX__(0.,
				state->upper[entering] - state->value[entering]);
		else if (entering_sign < 0 &&
			 dual_is_finite_lower(state->lower[entering]))
			opposite_step = __lp_simplex_MAX__(0.,
				state->value[entering] - state->lower[entering]);

		if (opposite_step < step) {
			dual_change_status(state, entering, entering_sign > 0
				? DUAL_STATUS_UPPER : DUAL_STATUS_LOWER);
			dual_set_nonbasic_value(state, entering);
			if (dual_compute_primal_values(state) == lp_simplex_EXIT_FAILURE) {
				*status = lp_simplex_PrecisionError;
				return lp_simplex_EXIT_FAILURE;
			}
			simplex_dual_feasibility_rebuild(state);
			state->iterations++;
			cleanup_iterations++;
			continue;
		}
		if (leaving_position < 0) {
			*status = lp_simplex_Unboundedness;
			return lp_simplex_EXIT_FAILURE;
		}
		{
			int leaving = state->basis[leaving_position];
			int update;
			state->position[leaving] = -1;
			dual_change_status(state, leaving, leaving_status);
			if (state->lower[leaving] == state->upper[leaving])
				dual_change_status(state, leaving, DUAL_STATUS_FIXED);
			dual_change_basis(state, leaving_position, entering);
			state->position[entering] = leaving_position;
			dual_change_status(state, entering, DUAL_STATUS_BASIC);
			update = simplex_basis_update(&state->factor,
				leaving_position, state->direction);
			if (update == lp_simplex_EXIT_FAILURE) {
				*status = lp_simplex_Singularity;
				return lp_simplex_EXIT_FAILURE;
			}
			if (update == 1) {
				if (dual_reinvert(state) == lp_simplex_EXIT_FAILURE) {
					*status = lp_simplex_Singularity;
					return lp_simplex_EXIT_FAILURE;
				}
			} else if (dual_compute_primal_values(state) ==
					lp_simplex_EXIT_FAILURE ||
				  dual_compute_reduced_costs(state) ==
					lp_simplex_EXIT_FAILURE) {
				*status = lp_simplex_PrecisionError;
				return lp_simplex_EXIT_FAILURE;
			}
		}
		state->iterations++;
		cleanup_iterations++;
		dual_rebuild_iteration_indexes(state);
	}
	*status = lp_simplex_ExceedIterLimit;
	return lp_simplex_EXIT_FAILURE;
}


static int dual_select_leaving(
		struct simplex_DualState *state, int *status,
		struct dual_Iteration *iteration)
{
	int prefer_structural;
	int lexicographic;
	double dual_error;
	prefer_structural = state->pan_enabled &&
		simplex_degeneracy_should_probe(&state->degeneracy) &&
		dual_pan_structure_enabled(state);
	lexicographic = state->pan_enabled &&
		simplex_degeneracy_is_lexicographic(&state->degeneracy);
	iteration->merit_before = simplex_dual_feasibility_merit(state);
	iteration->p = simplex_dual_feasibility_choose(state,
		&iteration->target, &iteration->kappa,
		&iteration->primal_infeasibility, prefer_structural,
		lexicographic, state->pan_deferred_row);
	if (iteration->p >= 0)
		return DUAL_STEP_READY;
	if (dual_reinvert(state) == lp_simplex_EXIT_FAILURE) {
		*status = lp_simplex_Singularity;
		return DUAL_STEP_FAILED;
	}
	iteration->merit_before = simplex_dual_feasibility_merit(state);
	iteration->p = simplex_dual_feasibility_choose(state,
		&iteration->target, &iteration->kappa,
		&iteration->primal_infeasibility, prefer_structural,
		lexicographic, state->pan_deferred_row);
	if (iteration->p >= 0)
		return DUAL_STEP_READY;
	/* Certify an apparent primal optimum against the current complete basis.
	 * Iterative refinement is deliberately terminal-only: it prevents an
	 * inaccurate success result without perturbing the normal pivot path. */
	if (dual_compute_primal_values_refined(state) == lp_simplex_EXIT_FAILURE) {
		*status = lp_simplex_PrecisionError;
		return DUAL_STEP_FAILED;
	}
	simplex_dual_feasibility_rebuild(state);
	iteration->merit_before = simplex_dual_feasibility_merit(state);
	iteration->p = simplex_dual_feasibility_choose(state,
		&iteration->target, &iteration->kappa,
		&iteration->primal_infeasibility, prefer_structural,
		lexicographic, state->pan_deferred_row);
	if (iteration->p >= 0)
		return DUAL_STEP_READY;
	dual_error = dual_max_dual_infeasibility(state);
	if (state->cost_shift_active &&
	    dual_error <= state->options->dual_tolerance) {
		if (dual_restore_original_cost(state) == lp_simplex_EXIT_FAILURE) {
			*status = lp_simplex_PrecisionError;
			return DUAL_STEP_FAILED;
		}
		dual_error = dual_max_dual_infeasibility(state);
	}
	if (dual_error <= state->options->dual_tolerance) {
		*status = lp_simplex_Success;
		return DUAL_STEP_FINISHED;
	}
	/* At this point the complete basis has passed primal certification.  Any
	 * remaining wrong-sign reduced costs belong to a primal-simplex cleanup;
	 * a dual crash would destructively discard the feasible basis on failure. */
	if (dual_primal_cleanup(state, status) == lp_simplex_EXIT_FAILURE)
		return DUAL_STEP_FAILED;
	dual_rebuild_iteration_indexes(state);
	return DUAL_STEP_RESTART;
}


static int dual_prepare_leaving_row(
		struct simplex_DualState *state, int *status,
		const struct dual_Iteration *iteration)
{
	double certificate_scale;
	double exact_weight = 0.;
	int i;
	lp_simplex_memset(state->rho, 0,
		(size_t)state->rows * sizeof(double));
	state->rho[iteration->p] = 1.;
	if (simplex_basis_btran(&state->factor, state->rho) ==
	    lp_simplex_EXIT_FAILURE) {
		*status = lp_simplex_Singularity;
		return DUAL_STEP_FAILED;
	}
	for (i = 0; i < state->rows; i++)
		exact_weight += state->rho[i] * state->rho[i];
	state->edge_weight[iteration->p] =
		__lp_simplex_MAX__(exact_weight, 1e-12);
	if (simplex_dual_prepare_ratio_test(state, iteration->kappa) != 0)
		return DUAL_STEP_READY;
	if (!simplex_dual_nonbasic_is_consistent(state)) {
		if (state->profile.enabled)
			fprintf(stderr,
				"dual profile: repairing nonbasic index set\n");
		simplex_dual_nonbasic_initialize(state);
		if (simplex_dual_prepare_ratio_test(state, iteration->kappa) != 0)
			return DUAL_STEP_RESTART;
	}
	if (simplex_basis_update_count(&state->factor) != 0 &&
	    dual_reinvert(state) == lp_simplex_EXIT_SUCCESS)
		return DUAL_STEP_RESTART;
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
	*status = lp_simplex_Infeasibility;
	return DUAL_STEP_FAILED;
}


static int dual_handle_exhausted_ratio(
		struct simplex_DualState *state, int *status,
		const struct dual_Iteration *iteration, const int bound_flips,
		int *rejected_relative, double *relative_pivot_tolerance)
{
	double certificate_scale;
	if (bound_flips > 0) {
		if (dual_flush_bound_flips(state, iteration->p) ==
		    lp_simplex_EXIT_FAILURE) {
			*status = lp_simplex_Singularity;
			return DUAL_STEP_FAILED;
		}
		return DUAL_STEP_RESTART;
	}
	if (*rejected_relative > 0 && *relative_pivot_tolerance > 0.) {
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
	if (simplex_basis_update_count(&state->factor) != 0 &&
	    dual_reinvert(state) == lp_simplex_EXIT_SUCCESS)
		return DUAL_STEP_RESTART;
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
			iteration->primal_infeasibility, *rejected_relative);
	*status = lp_simplex_Infeasibility;
	return DUAL_STEP_FAILED;
}


static int dual_try_batched_bound_flip(
		struct simplex_DualState *state, int *status,
		const struct dual_Iteration *iteration)
{
	double tentative_delta, tentative_value, movement = 0.;
	unsigned char flipped_status = DUAL_STATUS_FREE;
	int sign;
	if (state->numerically_stressed || state->regular_columns ||
	    __lp_simplex_ABS__(state->alpha[iteration->q]) <=
	    state->options->pivot_tolerance)
		return DUAL_STEP_READY;
	tentative_delta = (state->basic_value[iteration->p] -
		iteration->target) / state->alpha[iteration->q];
	tentative_value = state->value[iteration->q] + tentative_delta;
	sign = state->candidate_sign[iteration->q];
	if (sign * tentative_delta < -state->options->primal_tolerance ||
	    state->status[iteration->q] == DUAL_STATUS_FREE ||
	    !dual_opposite_bound_crossed(state, iteration->q, sign,
		tentative_value, &movement, &flipped_status) || movement == 0.)
		return DUAL_STEP_READY;
	if (dual_accumulate_bound_flip(state, iteration->q, movement) ==
	    lp_simplex_EXIT_FAILURE) {
		*status = lp_simplex_PrecisionError;
		return DUAL_STEP_FAILED;
	}
	state->basic_value[iteration->p] -=
		state->alpha[iteration->q] * movement;
	state->value[iteration->q] += movement;
	dual_change_status(state, iteration->q, flipped_status);
	state->candidate_sign[iteration->q] = 0;
	state->profile.bound_flips++;
	return DUAL_STEP_RETRY;
}


static int dual_build_entering_direction(
		struct simplex_DualState *state, int *status,
		struct dual_Iteration *iteration, int *weight_work_ready,
		int *rejected_relative, int *bound_flipped,
		const double relative_pivot_tolerance)
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
			state->candidate_sign[iteration->q] = 0;
			return DUAL_STEP_RETRY;
		}
		*weight_work_ready = 1;
	} else if (simplex_basis_ftran(&state->factor, state->direction) ==
		   lp_simplex_EXIT_FAILURE) {
		state->profile.rejected_ftran++;
		state->candidate_sign[iteration->q] = 0;
		return DUAL_STEP_RETRY;
	}
	simplex_sparse_vector_pack(&state->packed_direction,
		state->direction, 0.);
	for (i = 0; i < state->packed_direction.count; i++)
		direction_maximum = __lp_simplex_MAX__(direction_maximum,
			__lp_simplex_ABS__(state->packed_direction.value[i]));
	if (__lp_simplex_ABS__(state->direction[iteration->p]) <=
	    state->options->pivot_tolerance ||
	    __lp_simplex_ABS__(state->direction[iteration->p]) <
	    relative_pivot_tolerance * direction_maximum) {
		state->candidate_sign[iteration->q] = 0;
		(*rejected_relative)++;
		state->profile.rejected_relative++;
		return DUAL_STEP_RETRY;
	}
	if (__lp_simplex_ABS__(state->direction[iteration->p] -
		state->alpha[iteration->q]) >
	    100. * state->options->pivot_tolerance *
	    __lp_simplex_MAX__(1., __lp_simplex_ABS__(state->alpha[iteration->q]))) {
		if (state->pan_enabled)
			simplex_degeneracy_request_recovery(&state->degeneracy);
		if (simplex_basis_update_count(&state->factor) != 0 &&
		    dual_reinvert(state) == lp_simplex_EXIT_SUCCESS) {
			return DUAL_STEP_RESTART;
		}
		*status = lp_simplex_PrecisionError;
		return DUAL_STEP_FAILED;
	}
	iteration->delta = (state->basic_value[iteration->p] -
		iteration->target) / state->direction[iteration->p];
	sign = state->candidate_sign[iteration->q];
	iteration->new_value = state->value[iteration->q] + iteration->delta;
	if (sign * iteration->delta >= -state->options->primal_tolerance &&
	    (state->status[iteration->q] == DUAL_STATUS_FREE ||
	     !dual_opposite_bound_crossed(state, iteration->q, sign,
		iteration->new_value, &movement, &flipped_status)))
		return DUAL_STEP_READY;
	if (state->status[iteration->q] == DUAL_STATUS_FREE || movement == 0.) {
		state->candidate_sign[iteration->q] = 0;
		return DUAL_STEP_RETRY;
	}
	for (i = 0; i < state->packed_direction.count; i++) {
		int row = state->packed_direction.index[i];
		state->basic_value[row] -=
			state->packed_direction.value[i] * movement;
	}
	simplex_dual_feasibility_update_packed(state, &state->packed_direction);
	state->value[iteration->q] += movement;
	dual_change_status(state, iteration->q, flipped_status);
	state->candidate_sign[iteration->q] = 0;
	state->profile.bound_flips++;
	*bound_flipped = 1;
	return DUAL_STEP_RETRY;
}


static int dual_select_entering(
		struct simplex_DualState *state, int *status,
		struct dual_Iteration *iteration)
{
	int bound_flips = 0;
	int rejected_relative = 0;
	int weight_work_ready = 0;
	double relative_pivot_tolerance = DUAL_RELATIVE_PIVOT_TOLERANCE;
	for (;;) {
		int bound_flipped, step;
		iteration->q = simplex_dual_next_ratio_candidate(state,
			&iteration->theta,
			state->pan_enabled &&
				simplex_degeneracy_should_probe(&state->degeneracy),
			state->basis[iteration->p]);
		state->ratio_minimum_valid = 0;
		if (iteration->q < 0) {
			step = dual_handle_exhausted_ratio(state, status, iteration,
				bound_flips, &rejected_relative,
				&relative_pivot_tolerance);
			if (step == DUAL_STEP_RETRY)
				continue;
			return step;
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
			*status = lp_simplex_Singularity;
			return DUAL_STEP_FAILED;
		}
		step = dual_build_entering_direction(state, status, iteration,
			&weight_work_ready, &rejected_relative, &bound_flipped,
			relative_pivot_tolerance);
		if (step == DUAL_STEP_RETRY) {
			if (bound_flipped)
				bound_flips++;
			continue;
		}
		return step;
	}
}


static void dual_exchange_basis(
		struct simplex_DualState *state,
		const struct dual_Iteration *iteration)
{
	int leaving = state->basis[iteration->p];
	dual_update_dual_values(state, iteration->kappa, iteration->theta,
		leaving, iteration->q);
	simplex_dual_nonbasic_remove(state, iteration->q);
	state->structural_basic += (iteration->q < state->structural) -
		(leaving < state->structural);
	state->value[leaving] = iteration->target;
	dual_change_status(state, leaving, iteration->kappa > 0
		? DUAL_STATUS_LOWER : DUAL_STATUS_UPPER);
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


static int dual_update_factorization(
		struct simplex_DualState *state, const int p,
		const int compact_was_active, int *recomputed_dual)
{
	int compact_desired, update;
	compact_desired = simplex_basis_compact_requested(&state->factor) &&
		state->rows >= 4096 && state->structural_basic * 3 <= state->rows;
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
		struct simplex_DualState *state, int *status,
		const int recomputed_dual)
{
	if (!recomputed_dual || dual_max_dual_infeasibility(state) <=
	    10. * state->options->dual_tolerance)
		return lp_simplex_EXIT_SUCCESS;
	if (dual_reinvert(state) == lp_simplex_EXIT_FAILURE) {
		*status = lp_simplex_PrecisionError;
		return lp_simplex_EXIT_FAILURE;
	}
	if (dual_max_dual_infeasibility(state) <=
	    10. * state->options->dual_tolerance)
		return lp_simplex_EXIT_SUCCESS;
	if (state->cost_shift_allowed && dual_apply_cost_shifts(state) > 0 &&
	    dual_max_dual_infeasibility(state) <=
	    state->options->dual_tolerance)
		return lp_simplex_EXIT_SUCCESS;
	if (dual_crash(state, 0) == lp_simplex_EXIT_FAILURE ||
	    dual_reinvert(state) == lp_simplex_EXIT_FAILURE ||
	    dual_max_dual_infeasibility(state) >
	    10. * state->options->dual_tolerance) {
		*status = lp_simplex_PrecisionError;
		return lp_simplex_EXIT_FAILURE;
	}
	dual_rebuild_iteration_indexes(state);
	return lp_simplex_EXIT_SUCCESS;
}


static int dual_commit_iteration(
		struct simplex_DualState *state, int *status,
		const struct dual_Iteration *iteration)
{
	int compact_was_active, i, recomputed_dual = 0;
	double pan_step, pan_error, pan_merit_error;
	for (i = 0; i < state->packed_direction.count; i++) {
		int row = state->packed_direction.index[i];
		if (row != iteration->p)
			state->basic_value[row] -=
				state->packed_direction.value[i] * iteration->delta;
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
		*status = lp_simplex_Singularity;
		return lp_simplex_EXIT_FAILURE;
	}
	if (dual_validate_recomputed_dual(state, status, recomputed_dual) ==
	    lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	state->pan_deferred_row = -1;
	if (!state->pan_enabled)
		return lp_simplex_EXIT_SUCCESS;
	pan_step = iteration->theta * iteration->primal_infeasibility;
	/* First-order error propagation for theta*v, using the solver's
	 * feasibility tolerances rather than an activation-count threshold. */
	pan_error = state->options->dual_tolerance *
		iteration->primal_infeasibility +
		state->options->primal_tolerance *
		__lp_simplex_ABS__(iteration->theta) +
		state->options->dual_tolerance *
		state->options->primal_tolerance;
	pan_merit_error = (double)state->rows *
		state->options->primal_tolerance;
	simplex_degeneracy_observe(&state->degeneracy, pan_step, pan_error,
		iteration->merit_before,
		simplex_dual_feasibility_merit(state), pan_merit_error);
	if (simplex_degeneracy_should_probe(&state->degeneracy)) {
		simplex_basis_request_compact(&state->factor);
		simplex_degeneracy_record_state(&state->degeneracy,
			state->basis, state->rows, state->status, state->variables);
	}
	return lp_simplex_EXIT_SUCCESS;
}


static int dual_run(struct simplex_DualState *state, int *status)
{
	while (state->iterations < state->options->iteration_limit) {
		struct dual_Iteration iteration;
		int step;
		simplex_sparse_vector_clear(&state->flip_rhs);
		step = dual_select_leaving(state, status, &iteration);
		if (step == DUAL_STEP_FINISHED)
			return lp_simplex_EXIT_SUCCESS;
		if (step == DUAL_STEP_FAILED)
			return lp_simplex_EXIT_FAILURE;
		if (step == DUAL_STEP_RESTART)
			continue;
		step = dual_prepare_leaving_row(state, status, &iteration);
		if (step == DUAL_STEP_FAILED)
			return lp_simplex_EXIT_FAILURE;
		if (step == DUAL_STEP_RESTART)
			continue;
		step = dual_select_entering(state, status, &iteration);
		if (step == DUAL_STEP_FAILED)
			return lp_simplex_EXIT_FAILURE;
		if (step == DUAL_STEP_RESTART)
			continue;
		if (dual_commit_iteration(state, status, &iteration) ==
		    lp_simplex_EXIT_FAILURE)
			return lp_simplex_EXIT_FAILURE;
	}
	*status = lp_simplex_ExceedIterLimit;
	return lp_simplex_EXIT_FAILURE;
}

static double dual_primal_infeasibility(
		const struct simplex_DualState *state)
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
	lp_simplex_memset(state->work, 0,
		(size_t)state->rows * sizeof(double));
	for (j = 0; j < state->structural; j++) {
		if (state->value[j] != 0.)
			simplex_csc_column_axpy(&state->matrix, state->structural,
						  j, state->value[j], state->work);
	}
	for (i = 0; i < state->rows; i++) {
		double residual = __lp_simplex_ABS__(
			state->work[i] - state->value[state->structural + i]);
		if (residual > maximum)
			maximum = residual;
	}
	return maximum;
}


static int dual_solve_problem_kernel(
		const struct simplex_Problem *problem,
		const struct lp_simplex_Options *options,
		double *x, double *row_dual, struct lp_simplex_Result *result,
		const int propagate_bounds)
{
	struct simplex_DualState state;
	int j, status = lp_simplex_CondUnsatisfied;
	int solve_state;
	double objective = 0.;
	clock_t profile_started = 0;
	if (getenv("LP_SIMPLEX_PROFILE") != NULL)
		profile_started = clock();
	if (simplex_dual_state_create(&state, problem, options) ==
	    lp_simplex_EXIT_FAILURE) {
		result->status = lp_simplex_MemoryAllocError;
		return lp_simplex_EXIT_FAILURE;
	}
	simplex_dual_initialize_bounds(&state, problem, propagate_bounds);
	state.cost_shift_allowed = !propagate_bounds;
	lp_simplex_memcpy(state.original_cost, state.cost,
		(size_t)state.variables * sizeof(double));
	dual_initialize_basis(&state, problem);
	{
		const char *stage = NULL;
		int fallback_crash = 0;
		if (simplex_basis_factorize(&state.factor) == lp_simplex_EXIT_FAILURE)
			stage = "initial factorization";
		else if (dual_crash(&state, 0) == lp_simplex_EXIT_FAILURE) {
			/* Preserve the established index-ordered crash for normal
			 * models.  Only after it demonstrably cycles, reset the complete
			 * logical basis.  Small models first retain the targeted
			 * largest-violation crash; large raw-model certification uses the
			 * exact logical basis with a reversible Phase-I objective. */
			if (state.rows <= 1024) {
				dual_initialize_basis(&state, problem);
				if (simplex_basis_factorize(&state.factor) ==
					lp_simplex_EXIT_FAILURE ||
				    dual_crash(&state, 1) == lp_simplex_EXIT_FAILURE)
					stage = "dual crash";
				else
					fallback_crash = 1;
			} else if (state.cost_shift_allowed) {
				dual_initialize_basis(&state, problem);
				if (simplex_basis_factorize(&state.factor) ==
						lp_simplex_EXIT_FAILURE ||
				    dual_compute_reduced_costs(&state) ==
						lp_simplex_EXIT_FAILURE)
					stage = "phase-I logical basis";
				else {
					dual_apply_cost_shifts(&state);
					fallback_crash = 1;
				}
			} else
				stage = "dual crash";
		}
		if (stage == NULL && !fallback_crash &&
		    simplex_basis_factorize(&state.factor) ==
			 lp_simplex_EXIT_FAILURE)
			stage = "crash-basis factorization";
		if (stage == NULL &&
		    dual_compute_primal_values(&state) == lp_simplex_EXIT_FAILURE)
			stage = "initial primal values";
		if (stage == NULL &&
		    dual_compute_reduced_costs(&state) == lp_simplex_EXIT_FAILURE)
			stage = "initial reduced costs";
		if (stage == NULL)
			for (j = 0; j < state.rows; j++)
				state.edge_weight[j] = 1.;
		if (stage != NULL) {
			if (state.profile.enabled)
				fprintf(stderr, "dual profile: initialization failed at %s\n",
					stage);
			result->status = lp_simplex_PrecisionError;
			simplex_dual_state_destroy(&state);
			return lp_simplex_EXIT_FAILURE;
		}
	}
	simplex_dual_feasibility_rebuild(&state);
	lp_simplex_memset(state.alpha, 0,
		(size_t)state.variables * sizeof(double));
	lp_simplex_memset(state.candidate_sign, 0,
		(size_t)state.variables * sizeof(int));
	state.candidate_count = 0;
	simplex_dual_nonbasic_initialize(&state);
	state.structural_basic = 0;
	for (j = 0; j < state.rows; j++)
		if (state.basis[j] < state.structural)
			state.structural_basic++;
	solve_state = dual_run(&state, &status);
	/* A shifted auxiliary objective can establish primal infeasibility only
	 * after the original costs have been restored and certified.  Until then,
	 * treat a numerical certificate failure honestly as a precision error. */
	if (state.cost_shift_active && status == lp_simplex_Infeasibility) {
		dual_restore_original_cost(&state);
		status = lp_simplex_PrecisionError;
		solve_state = lp_simplex_EXIT_FAILURE;
	}
	for (j = 0; j < state.rows; j++)
		state.value[state.basis[j]] = state.basic_value[j];
	for (j = 0; j < problem->columns; j++) {
		x[j] = state.value[j];
		objective += problem->objective[j] * x[j];
	}
	result->status = status;
	result->iterations = state.iterations;
	result->objective = objective;
	result->primal_infeasibility = dual_primal_infeasibility(&state);
	result->dual_infeasibility = dual_max_dual_infeasibility(&state);
	if (row_dual != NULL)
		for (j = 0; j < problem->rows; j++)
			row_dual[j] = state.pi[j];
	if (state.profile.enabled) {
		struct simplex_BasisProfile basis_profile;
		int active_rank = 0;
		simplex_basis_get_profile(&state.factor, &basis_profile);
		for (j = 0; j < state.rows; j++) {
			int variable = state.basis[j];
			if (state.basic_value[j] > state.lower[variable] +
			    state.options->primal_tolerance &&
			    state.basic_value[j] < state.upper[variable] -
			    state.options->primal_tolerance)
				active_rank++;
		}
		state.profile.total_seconds =
			(double)(clock() - profile_started) / (double)CLOCKS_PER_SEC;
		fprintf(stderr,
			"dual profile: total=%.6f factor=%.6f/%ld ftran=%.6f/%ld "
			"btran=%.6f/%ld ratio=%.6f residual=%.6f "
			"pan_degenerate=%ld pan_activations=%ld pan_probes=%ld "
			"pan_cycles=%ld pan_recovery=%ld/%ld "
			"pan_rank=%d/%d factor_core=%d/%d compact=%ld[%d,%d] "
			"eta_nnz=%ld/%ld reinvert=%ld/%ld cert=%ld/%ld refine=%ld/%ld "
			"reject=%ld/%ld flips=%ld/%ld\n",
			state.profile.total_seconds,
			basis_profile.factor_seconds,
			basis_profile.factor_calls,
			basis_profile.ftran_seconds,
			basis_profile.ftran_calls,
			basis_profile.btran_seconds,
			basis_profile.btran_calls,
			state.profile.ratio_seconds,
			state.profile.total_seconds -
			basis_profile.factor_seconds -
			basis_profile.ftran_seconds -
			basis_profile.btran_seconds -
			state.profile.ratio_seconds,
			state.degeneracy.degenerate_pivots,
			state.degeneracy.activations,
			state.degeneracy.probes,
			state.degeneracy.repeated_states,
			state.degeneracy.recoveries,
			state.degeneracy.recovery_requests,
			active_rank, state.rows,
			basis_profile.factor_size, state.rows,
			basis_profile.compact_calls,
			basis_profile.compact_min,
			basis_profile.compact_max,
			basis_profile.eta_nonzeros,
			basis_profile.eta_slots,
			basis_profile.fill_reinversions,
			basis_profile.stability_reinversions,
			basis_profile.compact_ftran_validations,
			basis_profile.compact_btran_validations,
			basis_profile.compact_ftran_refinements,
			basis_profile.compact_btran_refinements,
			state.profile.rejected_relative,
			state.profile.rejected_ftran,
			state.profile.bound_flips,
			state.profile.flip_batches);
	}
	simplex_dual_state_destroy(&state);
	return solve_state;
}


static double dual_problem_primal_infeasibility(
		const struct simplex_Problem *problem, const double *x)
{
	double maximum = 0.;
	int i, j, k;
	for (j = 0; j < problem->columns; j++) {
		double violation = 0.;
		int type = problem->bounds[j].b_type;
		if ((type == optm_BOUND_T_LO || type == optm_BOUND_T_BS) &&
		    x[j] < problem->bounds[j].lb)
			violation = problem->bounds[j].lb - x[j];
		if ((type == optm_BOUND_T_UP || type == optm_BOUND_T_BS) &&
		    x[j] > problem->bounds[j].ub)
			violation = __lp_simplex_MAX__(violation,
				x[j] - problem->bounds[j].ub);
		maximum = __lp_simplex_MAX__(maximum, violation);
	}
	for (i = 0; i < problem->rows; i++) {
		double activity = 0.;
		double violation;
		for (k = problem->matrix.row_start[i];
		     k < problem->matrix.row_start[i + 1]; k++)
			activity += problem->matrix.row_value[k] *
				x[problem->matrix.column_index[k]];
		violation = activity - problem->rhs[i];
		if (problem->row_type[i] == optm_CONS_T_EQ)
			violation = __lp_simplex_ABS__(violation);
		else if (problem->row_type[i] == optm_CONS_T_GE)
			violation = violation < 0. ? -violation : 0.;
		else
			violation = violation > 0. ? violation : 0.;
		maximum = __lp_simplex_MAX__(maximum, violation);
	}
	return maximum;
}


int simplex_dual_solve_problem(
		const struct simplex_Problem *problem,
		const struct lp_simplex_Options *options,
		double *x, double *row_dual, struct lp_simplex_Result *result,
		const int propagate_bounds)
{
	struct simplex_Scaling scaling;
	struct simplex_Problem *scaled;
	double *dual = row_dual;
	int i, j, state;
	if (getenv("LP_SIMPLEX_DISABLE_SCALING") != NULL)
		return dual_solve_problem_kernel(problem, options, x, row_dual,
			result, propagate_bounds);
	scaled = simplex_scaling_create_problem(problem, &scaling);
	if (scaled == NULL) {
		result->status = lp_simplex_MemoryAllocError;
		return lp_simplex_EXIT_FAILURE;
	}
	if (dual == NULL && problem->rows > 0)
		dual = (double *)lp_simplex_malloc(
			(size_t)problem->rows * sizeof(double));
	if (dual == NULL && problem->rows > 0) {
		simplex_problem_free(scaled);
		simplex_scaling_destroy(&scaling);
		result->status = lp_simplex_MemoryAllocError;
		return lp_simplex_EXIT_FAILURE;
	}
	/* The kernel may reject an initial factorization before producing a
	 * solution.  Initialize the recovery buffers so that failure reporting does
	 * not read indeterminate caller storage. */
	for (j = 0; j < problem->columns; j++)
		x[j] = 0.;
	for (i = 0; i < problem->rows; i++)
		dual[i] = 0.;
	state = dual_solve_problem_kernel(scaled, options, x, dual, result,
		propagate_bounds);
	simplex_scaling_recover_primal(&scaling, x);
	simplex_scaling_recover_dual(&scaling, dual);
	result->objective = 0.;
	for (j = 0; j < problem->columns; j++)
		result->objective += problem->objective[j] * x[j];
	result->primal_infeasibility =
		dual_problem_primal_infeasibility(problem, x);
	if (result->status == lp_simplex_Success &&
	    result->primal_infeasibility > 10. * options->primal_tolerance) {
		result->status = lp_simplex_PrecisionError;
		state = lp_simplex_EXIT_FAILURE;
	}
	if (row_dual == NULL)
		lp_simplex_free(dual);
	simplex_problem_free(scaled);
	simplex_scaling_destroy(&scaling);
	return state;
}


int simplex_dual_solve(
		const struct lp_Model *model,
		const struct lp_simplex_Options *options,
		double *x, double *row_dual, struct lp_simplex_Result *result,
		const int propagate_bounds)
{
	struct simplex_Problem problem;
	int state;
	if (simplex_problem_from_model(&problem, model) ==
	    lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	state = simplex_dual_solve_problem(&problem, options, x, row_dual,
		result, propagate_bounds);
	simplex_problem_destroy(&problem);
	return state;
}
