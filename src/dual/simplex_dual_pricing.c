/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
/* Nonbasic active set, hypersparse pricing, and dual ratio candidate policy. */
#include "simplex_dual_pricing.h"
#include "utils.h"

#include <lp_simplex/status.h>

#include <time.h>
#include <stdio.h>


int simplex_dual_pricing_create(struct simplex_DualState *state)
{
	int variables = state->variables;
	if (simplex_sparse_vector_create(&state->packed_rho, state->rows) ==
	    lp_simplex_EXIT_FAILURE ||
	    simplex_sparse_vector_create(&state->packed_alpha, variables) ==
	    lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	state->candidate_sign = (int *)lp_simplex_malloc(
		(size_t)4 * variables * sizeof(int));
	state->candidate_index = state->candidate_sign != NULL
		? state->candidate_sign + variables : NULL;
	state->nonbasic_index = state->candidate_sign != NULL
		? state->candidate_sign + 2 * variables : NULL;
	state->nonbasic_slot = state->candidate_sign != NULL
		? state->candidate_sign + 3 * variables : NULL;
	if (state->candidate_sign == NULL || state->candidate_index == NULL ||
	    state->nonbasic_index == NULL || state->nonbasic_slot == NULL)
		return lp_simplex_EXIT_FAILURE;
	lp_simplex_memset(state->candidate_sign, 0,
		(size_t)variables * sizeof(int));
	state->candidate_count = 0;
	state->nonbasic_count = 0;
	return lp_simplex_EXIT_SUCCESS;
}


void simplex_dual_pricing_destroy(struct simplex_DualState *state)
{
	simplex_sparse_vector_destroy(&state->packed_rho);
	simplex_sparse_vector_destroy(&state->packed_alpha);
	lp_simplex_free(state->candidate_sign);
	state->candidate_sign = NULL;
	state->candidate_index = NULL;
	state->nonbasic_index = NULL;
	state->nonbasic_slot = NULL;
	state->candidate_count = 0;
	state->nonbasic_count = 0;
}


void simplex_dual_nonbasic_add(
		struct simplex_DualState *state, const int variable)
{
	if (state->nonbasic_slot[variable] >= 0)
		return;
	state->nonbasic_slot[variable] = state->nonbasic_count;
	state->nonbasic_index[state->nonbasic_count++] = variable;
}


void simplex_dual_nonbasic_remove(
		struct simplex_DualState *state, const int variable)
{
	int slot = state->nonbasic_slot[variable];
	int last;
	if (slot < 0)
		return;
	last = state->nonbasic_index[--state->nonbasic_count];
	state->nonbasic_index[slot] = last;
	state->nonbasic_slot[last] = slot;
	state->nonbasic_slot[variable] = -1;
}


void simplex_dual_nonbasic_initialize(struct simplex_DualState *state)
{
	int j;
	state->nonbasic_count = 0;
	for (j = 0; j < state->variables; j++)
		state->nonbasic_slot[j] = -1;
	for (j = 0; j < state->variables; j++)
		if (state->position[j] < 0 &&
		    state->status[j] != DUAL_STATUS_FIXED)
			simplex_dual_nonbasic_add(state, j);
}


int simplex_dual_nonbasic_is_consistent(
		const struct simplex_DualState *state)
{
	int count = 0, j;
	for (j = 0; j < state->variables; j++) {
		int expected = state->position[j] < 0 &&
			state->status[j] != DUAL_STATUS_FIXED;
		if (expected)
			count++;
		if (expected != (state->nonbasic_slot[j] >= 0))
			return 0;
	}
	return count == state->nonbasic_count;
}


static double pricing_column_dot(
		const struct simplex_DualState *state, const int variable)
{
	if (variable >= state->structural)
		return -state->rho[variable - state->structural];
	return simplex_csc_column_dot(&state->matrix, variable, state->rho);
}


static int pricing_candidate_sign(
		const struct simplex_DualState *state,
		const int variable, const int kappa, const double alpha)
{
	if (state->status[variable] == DUAL_STATUS_LOWER)
		return 1;
	if (state->status[variable] == DUAL_STATUS_UPPER)
		return -1;
	if (state->status[variable] == DUAL_STATUS_FREE)
		return kappa * alpha > 0. ? -1 : 1;
	return 0;
}


static int pricing_significantly_greater(
		const double value, const double incumbent)
{
	return value > incumbent + 1e-12 *
		(1. + __lp_simplex_MAX__(__lp_simplex_ABS__(value),
			__lp_simplex_ABS__(incumbent)));
}


static int pricing_nearly_equal(const double left, const double right)
{
	return __lp_simplex_ABS__(left - right) <= 1e-12 *
		(1. + __lp_simplex_MAX__(__lp_simplex_ABS__(left),
			__lp_simplex_ABS__(right)));
}


static void pricing_add_candidate(
		struct simplex_DualState *state, const int j, const int kappa,
		const int alpha_ready, const int harris_enabled,
		int *candidates, double *minimum)
{
	int sign;
	double denominator, signed_reduced, theta, relaxed_theta;
	state->candidate_sign[j] = 0;
	if (!alpha_ready)
		state->alpha[j] = pricing_column_dot(state, j);
	sign = pricing_candidate_sign(state, j, kappa, state->alpha[j]);
	denominator = -kappa * sign * state->alpha[j];
	if (denominator <= state->options->pivot_tolerance)
		return;
	signed_reduced = sign * state->reduced[j];
	if (signed_reduced < 0. &&
	    signed_reduced >= -state->options->dual_tolerance)
		signed_reduced = 0.;
	if (signed_reduced < 0.)
		return;
	theta = signed_reduced / denominator;
	/* Harris first pass: permit the reduced cost to move to the edge of the
	 * dual feasibility tolerance.  The second pass can then prefer a much
	 * stronger pivot anywhere inside this safe step window instead of being
	 * tied to the exact minimum ratio. */
	relaxed_theta = (signed_reduced + state->options->dual_tolerance) /
		denominator;
	state->candidate_sign[j] = sign;
	state->breakpoint[j] = theta;
	state->candidate_index[*candidates] = j;
	if ((harris_enabled ? relaxed_theta : theta) < *minimum)
		*minimum = harris_enabled ? relaxed_theta : theta;
	(*candidates)++;
}


int simplex_dual_prepare_ratio_test(
		struct simplex_DualState *state, const int kappa)
{
	int candidate, candidates = 0, i, j;
	int alpha_ready = 0;
	int harris_enabled = !simplex_degeneracy_is_stressed(
		&state->degeneracy);
	long row_entries = 0;
	double minimum = __lp_simplex_INF__;
	clock_t started = 0;
	if (state->profile.enabled)
		started = clock();
	for (candidate = 0; candidate < state->candidate_count; candidate++)
		state->candidate_sign[state->candidate_index[candidate]] = 0;
	simplex_sparse_vector_pack(&state->packed_rho, state->rho, 0.);
	simplex_sparse_vector_clear(&state->packed_alpha);
	for (i = 0; i < state->packed_rho.count; i++) {
		int row = state->packed_rho.index[i];
		row_entries += state->matrix.row_start[row + 1] -
			state->matrix.row_start[row];
	}
	if (row_entries * 4 < state->matrix.nonzeros) {
		for (i = 0; i < state->packed_rho.count; i++) {
			int k;
			int row = state->packed_rho.index[i];
			double value = state->packed_rho.value[i];
			for (k = state->matrix.row_start[row];
			     k < state->matrix.row_start[row + 1]; k++)
				simplex_sparse_vector_add(&state->packed_alpha,
					state->matrix.column_index[k],
					value * state->matrix.row_value[k]);
			simplex_sparse_vector_set(&state->packed_alpha,
				state->structural + row, -value);
		}
		for (i = 0; i < state->packed_alpha.count; i++)
			state->alpha[state->packed_alpha.index[i]] =
				state->packed_alpha.value[i];
		alpha_ready = 1;
	}
	state->packed_alpha_valid = alpha_ready;
	if (alpha_ready) {
		if (state->pricing_validation_countdown-- <= 0) {
			double maximum_error = 0.;
			int missing = 0;
			for (i = 0; i < state->nonbasic_count; i++) {
				double exact, packed, error;
				j = state->nonbasic_index[i];
				exact = pricing_column_dot(state, j);
				packed = simplex_sparse_vector_get(
					&state->packed_alpha, j);
				error = __lp_simplex_ABS__(exact - packed);
				maximum_error = __lp_simplex_MAX__(maximum_error, error);
				if (error > 1e-11 * (1. + __lp_simplex_ABS__(exact))) {
					simplex_sparse_vector_set(
						&state->packed_alpha, j, exact);
					missing++;
				}
			}
			if (state->profile.enabled)
				fprintf(stderr,
					"dual profile: hypersparse alpha audit error=%.3g repaired=%d touched=%d/%d\n",
					maximum_error, missing, state->packed_alpha.count,
					state->nonbasic_count);
			state->pricing_validation_countdown = 127;
		}
		for (i = 0; i < state->packed_alpha.count; i++) {
			j = state->packed_alpha.index[i];
			if (state->nonbasic_slot[j] < 0)
				continue;
			state->candidate_sign[j] = 0;
			pricing_add_candidate(state, j, kappa, 1, harris_enabled,
				&candidates, &minimum);
		}
	} else for (i = 0; i < state->nonbasic_count; i++) {
		j = state->nonbasic_index[i];
		state->candidate_sign[j] = 0;
		pricing_add_candidate(state, j, kappa, 0, harris_enabled,
			&candidates, &minimum);
	}
	/* A sparse pricing pass is only a fast candidate generator, never an
	 * infeasibility certificate.  Cancellation and stale sparsity can omit an
	 * eligible column; certify an empty set with a full CSC dot-product scan. */
	if (candidates == 0 && alpha_ready) {
		minimum = __lp_simplex_INF__;
		for (i = 0; i < state->nonbasic_count; i++) {
			j = state->nonbasic_index[i];
			pricing_add_candidate(state, j, kappa, 0, harris_enabled,
				&candidates, &minimum);
		}
		/* The exact certification pass supersedes the hypersparse alpha.
		 * Keeping packed_alpha_valid here would update reduced costs for only
		 * the old sparse support after choosing a column found by the full scan. */
		state->packed_alpha_valid = 0;
	}
	state->ratio_minimum = minimum;
	state->ratio_minimum_valid = 1;
	state->candidate_count = candidates;
	if (state->profile.enabled)
		state->profile.ratio_seconds +=
			(double)(clock() - started) / (double)CLOCKS_PER_SEC;
	return candidates;
}


int simplex_dual_next_ratio_candidate(
		struct simplex_DualState *state, double *chosen_theta,
		const int pan_perturbation, const int leaving_variable)
{
	int candidate, j, q = -1;
	int stable_order = state->stable_candidate_order || pan_perturbation;
	double minimum = state->ratio_minimum_valid
		? state->ratio_minimum : __lp_simplex_INF__;
	double best_pivot = 0.;
	if (!state->ratio_minimum_valid)
		for (candidate = 0; candidate < state->candidate_count; candidate++) {
			j = state->candidate_index[candidate];
			if (state->candidate_sign[j] != 0 &&
			    state->breakpoint[j] < minimum)
				minimum = state->breakpoint[j];
		}
	if (minimum == __lp_simplex_INF__)
		return -1;
	if (pan_perturbation && minimum <= 10. * state->options->dual_tolerance) {
		int pan_candidate = -1;
		int best_rank_delta = 2;
		int structure_enabled =
			simplex_basis_compact_ever_active(&state->factor) ||
			(state->rows >= 4096 && state->structural_basic < state->rows);
		double best_structural_score = -1.;
		for (candidate = 0; candidate < state->candidate_count; candidate++) {
			int nonzeros, rank_delta;
			double structural_score;
			j = state->candidate_index[candidate];
			if (state->candidate_sign[j] == 0 ||
			    state->breakpoint[j] > minimum +
			    1e-12 * (1. + __lp_simplex_ABS__(minimum)))
				continue;
			rank_delta = (j < state->structural) -
				(leaving_variable < state->structural);
			nonzeros = j < state->structural
				? state->matrix.column_start[j + 1] -
				  state->matrix.column_start[j] : 1;
			structural_score = __lp_simplex_ABS__(state->alpha[j]) /
				(1. + nonzeros);
			if (pan_candidate < 0 ||
			    (!structure_enabled && j < pan_candidate) ||
			    (structure_enabled &&
			     (rank_delta < best_rank_delta ||
			      (rank_delta == best_rank_delta &&
				       (pricing_significantly_greater(structural_score,
					best_structural_score) ||
				        (pricing_nearly_equal(structural_score,
					best_structural_score) &&
			         j < pan_candidate)))))) {
				pan_candidate = j;
				best_rank_delta = rank_delta;
				best_structural_score = structural_score;
			}
		}
		if (pan_candidate >= 0) {
			*chosen_theta = state->breakpoint[pan_candidate];
			simplex_degeneracy_record_probe(&state->degeneracy);
			return pan_candidate;
		}
	}
	for (candidate = 0; candidate < state->candidate_count; candidate++) {
		double theta, relaxed;
		j = state->candidate_index[candidate];
		if (state->candidate_sign[j] == 0)
			continue;
		theta = state->breakpoint[j];
		relaxed = minimum + 1e-12 * (1. + __lp_simplex_ABS__(minimum));
		if (theta <= relaxed &&
		    ((!stable_order &&
		      (__lp_simplex_ABS__(state->alpha[j]) > best_pivot ||
		       (__lp_simplex_ABS__(state->alpha[j]) == best_pivot &&
		        (q < 0 || j < q)))) ||
		     (stable_order &&
		      (pricing_significantly_greater(
			__lp_simplex_ABS__(state->alpha[j]), best_pivot) ||
		       (pricing_nearly_equal(
			__lp_simplex_ABS__(state->alpha[j]), best_pivot) &&
		        (q < 0 || j < q)))))) {
			best_pivot = __lp_simplex_ABS__(state->alpha[j]);
			q = j;
			*chosen_theta = theta;
		}
	}
	return q;
}
