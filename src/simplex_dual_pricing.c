/* Nonbasic active set, hypersparse pricing, and dual ratio candidate policy. */
#include "simplex_dual_pricing.h"
#include "utils.h"
#include <lp_simplex/status.h>
#include <time.h>


int simplex_dual_pricing_create(struct simplex_DualState *state)
{
	int variables = state->variables;
	if (simplex_sparse_vector_create(&state->packed_rho, state->rows) ==
	    lp_simplex_EXIT_FAILURE ||
	    simplex_sparse_vector_create(&state->packed_alpha, variables) ==
	    lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	state->candidate_sign = (int *)lp_simplex_malloc(
		(size_t)variables * sizeof(int));
	state->candidate_index = (int *)lp_simplex_malloc(
		(size_t)variables * sizeof(int));
	state->nonbasic_index = (int *)lp_simplex_malloc(
		(size_t)variables * sizeof(int));
	state->nonbasic_slot = (int *)lp_simplex_malloc(
		(size_t)variables * sizeof(int));
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
	lp_simplex_free(state->candidate_index);
	lp_simplex_free(state->nonbasic_index);
	lp_simplex_free(state->nonbasic_slot);
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


static void pricing_add_candidate(
		struct simplex_DualState *state, const int j, const int kappa,
		const int alpha_ready, int *candidates, double *minimum)
{
	int sign;
	double denominator, signed_reduced, theta;
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
	state->candidate_sign[j] = sign;
	state->breakpoint[j] = theta;
	state->candidate_index[*candidates] = j;
	if (theta < *minimum)
		*minimum = theta;
	(*candidates)++;
}


int simplex_dual_prepare_ratio_test(
		struct simplex_DualState *state, const int kappa)
{
	int candidate, candidates = 0, i, j;
	int alpha_ready = 0;
	long row_entries = 0;
	double minimum = __lp_simplex_INF__;
	clock_t started = 0;
	if (state->profile_enabled)
		started = clock();
	for (candidate = 0; candidate < state->candidate_count; candidate++)
		state->candidate_sign[state->candidate_index[candidate]] = 0;
	simplex_sparse_vector_pack(&state->packed_rho, state->rho, 0.);
	for (i = 0; i < state->packed_rho.count; i++) {
		int row = state->packed_rho.index[i];
		row_entries += state->matrix.row_start[row + 1] -
			state->matrix.row_start[row];
	}
	if (row_entries * 4 < state->matrix.nonzeros) {
		lp_simplex_memset(state->alpha, 0,
			(size_t)state->variables * sizeof(double));
		for (i = 0; i < state->packed_rho.count; i++) {
			int k;
			int row = state->packed_rho.index[i];
			double value = state->packed_rho.value[i];
			for (k = state->matrix.row_start[row];
			     k < state->matrix.row_start[row + 1]; k++)
				state->alpha[state->matrix.column_index[k]] +=
					value * state->matrix.row_value[k];
			state->alpha[state->structural + row] = -value;
		}
		alpha_ready = 1;
	}
	simplex_sparse_vector_clear(&state->packed_alpha);
	state->packed_alpha_valid = alpha_ready;
	for (i = 0; i < state->nonbasic_count; i++) {
		j = state->nonbasic_index[i];
		state->candidate_sign[j] = 0;
		pricing_add_candidate(state, j, kappa, alpha_ready,
			&candidates, &minimum);
		if (alpha_ready && state->alpha[j] != 0.)
			simplex_sparse_vector_set(
				&state->packed_alpha, j, state->alpha[j]);
	}
	state->ratio_minimum = minimum;
	state->ratio_minimum_valid = 1;
	state->candidate_count = candidates;
	if (state->profile_enabled)
		state->profile_ratio_seconds +=
			(double)(clock() - started) / (double)CLOCKS_PER_SEC;
	return candidates;
}


int simplex_dual_next_ratio_candidate(
		struct simplex_DualState *state, double *chosen_theta,
		const int pan_perturbation, const int leaving_variable)
{
	int candidate, j, q = -1;
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
		int structure_enabled = state->factor.compact_ever_active ||
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
			       (structural_score > best_structural_score ||
			        (structural_score == best_structural_score &&
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
		    (__lp_simplex_ABS__(state->alpha[j]) > best_pivot ||
		     (__lp_simplex_ABS__(state->alpha[j]) == best_pivot &&
		      (q < 0 || j < q)))) {
			best_pivot = __lp_simplex_ABS__(state->alpha[j]);
			q = j;
			*chosen_theta = theta;
		}
	}
	return q;
}
