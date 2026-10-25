/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_dual_pricing.c
 * @brief Nonbasic active set, hypersparse pricing, and dual ratio candidate
 * policy.
 */

#include "policy/factorization.h"
#include "simplex_dual_pricing.h"
#include "utils.h"
#include "simplex_storage.h"

#include <lp_simplex/status.h>

#include <time.h>
#include <stdio.h>
#include <math.h>

int simplex_dual_pricing_create(struct simplex_DualState *state)
{
	int variables = state->variables;
	size_t bytes;

	if (variables < 0 || !simplex_storage_bytes((size_t)variables,
				     4 * sizeof(int), (size_t)-1, &bytes))
		return lp_simplex_EXIT_FAILURE;
	if (simplex_sparse_vector_create(&state->packed_rho, state->rows) ==
			lp_simplex_EXIT_FAILURE ||
		simplex_sparse_vector_create(&state->packed_alpha, variables) ==
			lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	state->candidate_sign = (int *)lp_simplex_malloc(bytes);
	state->candidate_index = state->candidate_sign != NULL
					 ? state->candidate_sign + variables
					 : NULL;
	state->nonbasic_index =
		state->candidate_sign != NULL
			? state->candidate_sign + (size_t)2 * variables
			: NULL;
	state->nonbasic_slot =
		state->candidate_sign != NULL
			? state->candidate_sign + (size_t)3 * variables
			: NULL;
	if (state->candidate_sign == NULL || state->candidate_index == NULL ||
		state->nonbasic_index == NULL || state->nonbasic_slot == NULL)
		return lp_simplex_EXIT_FAILURE;
	lp_simplex_memset(
		state->candidate_sign, 0, (size_t)variables * sizeof(int));
	state->candidate_count = 0;
	state->nonbasic_count = 0;
	return lp_simplex_EXIT_SUCCESS;
}

void simplex_dual_pricing_destroy(struct simplex_DualState *state)
{
	simplex_sparse_vector_destroy(&state->packed_rho);
	simplex_sparse_vector_destroy(&state->packed_alpha);
	lp_simplex_free(state->candidate_sign);
	lp_simplex_free(state->lex_storage);
	state->lex_storage = NULL;
	state->lex_rho = NULL;
	state->lex_candidate = NULL;
	state->lex_slot = NULL;
	state->lex_anchor_status = NULL;
	state->lex_initialized = 0;
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

int simplex_dual_nonbasic_is_consistent(const struct simplex_DualState *state)
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
	int k;
	double result = 0.;
	const struct simplex_CscMatrix *matrix = &state->matrix;

	if (variable >= state->structural)
		return -state->rho[variable - state->structural];

	/* Keep the short CSC reduction in the pricing translation unit so the
	 * compiler can fuse it into the candidate loop without per-column calls.
	 * Preserve the CSC summation order used by the certification pass. */
	for (k = matrix->column_start[variable];
		k < matrix->column_start[variable + 1]; k++)
		result += state->rho[matrix->row_index[k]] * matrix->value[k];
	return result;
}

static int pricing_candidate_sign(const struct simplex_DualState *state,
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
	return value >
	       incumbent +
		       SIMPLEX_PRICING_TIE_RELATIVE *
			       (1. + __lp_simplex_MAX__(
					     __lp_simplex_ABS__(value),
					     __lp_simplex_ABS__(incumbent)));
}

static int pricing_nearly_equal(const double left, const double right)
{
	return __lp_simplex_ABS__(left - right) <=
	       SIMPLEX_PRICING_TIE_RELATIVE *
		       (1. + __lp_simplex_MAX__(__lp_simplex_ABS__(left),
				     __lp_simplex_ABS__(right)));
}

static void pricing_add_candidate(struct simplex_DualState *state, const int j,
	const int kappa, const int alpha_ready, const int harris_enabled,
	int *candidates, double *minimum)
{
	int sign;
	double denominator, signed_reduced, theta, relaxed_theta;

	state->candidate_sign[j] = 0;
	if (!alpha_ready)
		state->alpha[j] = pricing_column_dot(state, j);
	/* Non-finite arithmetic cannot certify an empty admissible set. */
	if (!isfinite(state->alpha[j]) || !isfinite(state->reduced[j])) {
		state->ratio_failure_status = lp_simplex_PrecisionError;
		return;
	}
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
	relaxed_theta =
		(signed_reduced + state->options->dual_tolerance) / denominator;
	if (!isfinite(theta) || !isfinite(relaxed_theta)) {
		state->ratio_failure_status = lp_simplex_PrecisionError;
		return;
	}
	state->candidate_sign[j] = sign;
	state->breakpoint[j] = theta;
	state->candidate_index[*candidates] = j;
	if ((harris_enabled ? relaxed_theta : theta) < *minimum)
		*minimum = harris_enabled ? relaxed_theta : theta;
	(*candidates)++;
}

void simplex_dual_ratio_remove(
	struct simplex_DualState *state, const int variable)
{
	if (state->candidate_sign[variable] != 0) {
		state->candidate_sign[variable] = 0;
		state->ratio_active_count--;
	}
}

/* CSR scatter visits the same column many times. Accumulate only in dense
 * storage, record its first touch, then materialize packed values once. This
 * avoids a checked add/set call and two numeric stores for every matrix entry.
 * Matrix indices are validated when the problem is built. Keep row order and
 * exact zeros: cancellation must not change the support used by the audit. */
static void pricing_scatter_alpha(struct simplex_DualState *state)
{
	struct simplex_SparseVector *alpha = &state->packed_alpha;
	const struct simplex_CscMatrix *matrix = &state->matrix;
	int i, k;

	for (i = 0; i < state->packed_rho.count; i++) {
		int row = state->packed_rho.index[i];
		double rho = state->packed_rho.value[i];

		for (k = matrix->row_start[row]; k < matrix->row_start[row + 1];
			k++) {
			int column = matrix->column_index[k];
			double value = rho * matrix->row_value[k];

			if (alpha->mark[column] != alpha->generation) {
				if (value == 0.)
					continue;
				alpha->mark[column] = alpha->generation;
				alpha->slot[column] = alpha->count;
				alpha->index[alpha->count++] = column;
				alpha->dense[column] = value;
			} else
				alpha->dense[column] += value;
		}
		/* Logical columns have disjoint support and are touched only once. */
		k = state->structural + row;
		alpha->mark[k] = alpha->generation;
		alpha->slot[k] = alpha->count;
		alpha->index[alpha->count++] = k;
		alpha->dense[k] = -rho;
	}
	for (i = 0; i < alpha->count; i++)
		alpha->value[i] = alpha->dense[alpha->index[i]];
}

/* With at least 95% nonbasics, sequential CSC access avoids following a
 * swap-deleted index list at a membership-scan overhead of at most 5.3%.
 * Less-wide models retain the compact list.
 * This changes enumeration order only, not the admissibility or ratio tests. */
static void pricing_dense_candidates(struct simplex_DualState *state,
	const int kappa, const int harris_enabled, int *candidates,
	double *minimum)
{
	int i, j;

	if (state->nonbasic_count >=
		state->variables -
			state->variables /
				SIMPLEX_DUAL_DENSE_MEMBERSHIP_DIVISOR) {
		for (j = 0; j < state->variables; j++) {
			if (state->nonbasic_slot[j] < 0)
				continue;
			pricing_add_candidate(state, j, kappa, 0,
				harris_enabled, candidates, minimum);
		}
	} else
		for (i = 0; i < state->nonbasic_count; i++) {
			j = state->nonbasic_index[i];
			pricing_add_candidate(state, j, kappa, 0,
				harris_enabled, candidates, minimum);
		}
}

int simplex_dual_prepare_ratio_test(
	struct simplex_DualState *state, const int kappa)
{
	int candidate, candidates = 0, i, j;
	int alpha_ready = 0;
	int harris_enabled =
		!simplex_degeneracy_is_stressed(&state->degeneracy);
	long row_entries = 0;
	double minimum = __lp_simplex_INF__;
	double audit_seconds = 0.;
	clock_t started = 0, path_started = 0;

	if (state->profile.enabled)
		started = clock();
	state->ratio_failure_status = lp_simplex_Success;
	for (candidate = 0; candidate < state->candidate_count; candidate++)
		state->candidate_sign[state->candidate_index[candidate]] = 0;
	simplex_sparse_vector_pack(&state->packed_rho, state->rho, 0.);
	simplex_sparse_vector_clear(&state->packed_alpha);
	for (i = 0; i < state->packed_rho.count; i++) {
		int row = state->packed_rho.index[i];

		row_entries += state->matrix.row_start[row + 1] -
			       state->matrix.row_start[row];
	}
	if (state->profile.enabled) {
		path_started = clock();
		state->profile.ratio_setup_seconds +=
			(double)(path_started - started) / CLOCKS_PER_SEC;
	}
	if (row_entries * 4 < state->matrix.nonzeros) {
		pricing_scatter_alpha(state);
		alpha_ready = 1;
	}
	state->packed_alpha_valid = alpha_ready;
	if (alpha_ready) {
		if (state->pricing_validation_countdown-- <= 0) {
			double maximum_error = 0.;
			int missing = 0;
			clock_t audit_started = 0;

			if (state->profile.enabled)
				audit_started = clock();

			for (i = 0; i < state->nonbasic_count; i++) {
				double exact, packed, error;

				j = state->nonbasic_index[i];
				exact = pricing_column_dot(state, j);
				packed = simplex_sparse_vector_get(
					&state->packed_alpha, j);
				error = __lp_simplex_ABS__(exact - packed);
				maximum_error = __lp_simplex_MAX__(
					maximum_error, error);
				if (error > SIMPLEX_PRICING_AUDIT_RELATIVE *
						    (1. + __lp_simplex_ABS__(
								  exact))) {
					simplex_sparse_vector_set(
						&state->packed_alpha, j, exact);
					missing++;
				}
			}
			if (state->profile.enabled) {
				audit_seconds =
					(double)(clock() - audit_started) /
					CLOCKS_PER_SEC;
				state->profile.ratio_audit_seconds +=
					audit_seconds;
				state->profile.ratio_audit_calls++;
				state->profile.ratio_audit_columns +=
					state->nonbasic_count;
				fprintf(stderr,
					"dual profile: hypersparse alpha audit error=%.3g repaired=%d touched=%d/%d\n",
					maximum_error, missing,
					state->packed_alpha.count,
					state->nonbasic_count);
			}
			state->pricing_validation_countdown = 127;
		}
		/* Audit repairs update the shared dense ratio row directly. */
		for (i = 0; i < state->packed_alpha.count; i++) {
			j = state->packed_alpha.index[i];
			if (state->nonbasic_slot[j] < 0)
				continue;
			state->candidate_sign[j] = 0;
			pricing_add_candidate(state, j, kappa, 1,
				harris_enabled, &candidates, &minimum);
		}
	} else
		pricing_dense_candidates(
			state, kappa, harris_enabled, &candidates, &minimum);
	if (state->profile.enabled) {
		double elapsed =
			(double)(clock() - path_started) / CLOCKS_PER_SEC;

		if (alpha_ready) {
			state->profile.ratio_sparse_calls++;
			state->profile.ratio_scatter_entries += row_entries;
			state->profile.ratio_sparse_touched +=
				state->packed_alpha.count;
			state->profile.ratio_sparse_seconds +=
				elapsed - audit_seconds;
		} else {
			state->profile.ratio_dense_calls++;
			state->profile.ratio_dense_columns +=
				state->nonbasic_count;
			state->profile.ratio_dense_seconds += elapsed;
		}
	}
	/* A sparse pricing pass is only a fast candidate generator, never an
	 * infeasibility certificate.  Cancellation and stale sparsity can omit an
	 * eligible column; certify an empty set with a full CSC dot-product scan. */
	if (candidates == 0 && alpha_ready) {
		if (state->profile.enabled) {
			path_started = clock();
			state->profile.ratio_certify_calls++;
		}
		minimum = __lp_simplex_INF__;
		pricing_dense_candidates(
			state, kappa, harris_enabled, &candidates, &minimum);
		/* The exact certification pass supersedes the hypersparse alpha.
		 * Keeping packed_alpha_valid here would update reduced costs for only
		 * the old sparse support after choosing a column found by the full scan. */
		state->packed_alpha_valid = 0;
		if (state->profile.enabled)
			state->profile.ratio_certify_seconds +=
				(double)(clock() - path_started) /
				CLOCKS_PER_SEC;
	}
	state->ratio_ordered = 0;
	state->ratio_heap_valid = 0;
	state->ratio_minimum = minimum;
	state->ratio_minimum_valid = 1;
	state->candidate_count = candidates;
	state->ratio_active_count = candidates;
	if (state->profile.enabled)
		state->profile.ratio_seconds +=
			(double)(clock() - started) / (double)CLOCKS_PER_SEC;
	return state->ratio_failure_status == lp_simplex_Success
		       ? candidates
		       : SIMPLEX_DUAL_RATIO_ERROR;
}

struct pricing_LexStorage {
	size_t numeric, indices, total;
};

/* This owner checks its combined allocation independently of the dimensions
 * admitted by the separate startup workspaces. No total-memory quota here. */
static int pricing_lex_storage(int rows, int variables, size_t limit,
	struct pricing_LexStorage *layout)
{
	struct pricing_LexStorage next = {0};
	*layout = next;
	if (rows < 0 || variables < rows ||
		!simplex_storage_bytes(
			(size_t)rows, sizeof(double), limit, &next.numeric) ||
		!simplex_storage_bytes((size_t)variables, 2 * sizeof(int),
			limit, &next.indices) ||
		next.numeric > limit - next.indices)
		return 0;
	next.total = next.numeric + next.indices;
	if ((size_t)variables > limit - next.total)
		return 0;
	next.total += (size_t)variables;
	*layout = next;
	return 1;
}

static int pricing_lex_workspace(struct simplex_DualState *state)
{
	struct pricing_LexStorage layout;
	unsigned char *storage;
	if (state->lex_storage != NULL)
		return lp_simplex_EXIT_SUCCESS;
	if (!pricing_lex_storage(
		    state->rows, state->variables, (size_t)-1, &layout))
		return lp_simplex_EXIT_FAILURE;
	storage = (unsigned char *)lp_simplex_malloc(layout.total);
	if (storage == NULL)
		return lp_simplex_EXIT_FAILURE;
	state->lex_storage = storage;
	state->lex_rho = (double *)storage;
	state->lex_candidate = (int *)(storage + layout.numeric);
	state->lex_slot = state->lex_candidate + state->variables;
	state->lex_anchor_status = storage + layout.numeric + layout.indices;
	return lp_simplex_EXIT_SUCCESS;
}

/* Compare symbolic reduced-cost ratios without choosing a numeric epsilon.
 * At activation, perturb each nonbasic lower/upper cost by +/- eps^j and
 * leave basic/free/fixed costs unchanged.  Coordinate j of a reduced cost is
 * sign_j * (delta_qj - e_position(j)^T B^-1 a_q).  Evaluate only coordinates
 * still needed to distinguish tied candidates, using at most one BTRAN per
 * basic coordinate, rather than one FTRAN for every tied column. */
static int pricing_lexicographic_candidate(struct simplex_DualState *state,
	const double minimum, const int leaving_variable, const int kappa)
{
	int count = 0, i, j, coordinate;
	int leaving_position = state->position[leaving_variable];

	if (pricing_lex_workspace(state) == lp_simplex_EXIT_FAILURE) {
		state->ratio_failure_status = lp_simplex_MemoryAllocError;
		return SIMPLEX_DUAL_RATIO_ERROR;
	}
	if (!state->lex_initialized) {
		lp_simplex_memcpy(state->lex_anchor_status, state->status,
			(size_t)state->variables);
		state->lex_initialized = 1;
	}
	for (j = 0; j < state->variables; j++)
		state->lex_slot[j] = -1;
	for (i = 0; i < state->candidate_count; i++) {
		j = state->candidate_index[i];
		/* Symbolic coordinates break exact primary ties.  Harris' numerical
		 * window must not override the primary objective in this mode. */
		if (state->candidate_sign[j] != 0 &&
			state->breakpoint[j] == minimum) {
			state->lex_slot[j] = count;
			state->lex_candidate[count++] = j;
		}
	}
	for (coordinate = 0; coordinate < state->variables && count > 1;
		coordinate++) {
		int position = state->position[coordinate];
		int anchor = state->lex_anchor_status[coordinate];
		double sign;

		if (anchor != DUAL_STATUS_LOWER && anchor != DUAL_STATUS_UPPER)
			continue;
		sign = anchor == DUAL_STATUS_LOWER ? 1. : -1.;
		if (position < 0) {
			int slot = state->lex_slot[coordinate];

			if (slot >= 0) {
				double denominator =
					-kappa * state->alpha[coordinate];

				if (sign / denominator < 0.) {
					state->lex_candidate[0] = coordinate;
					count = 1;
					break;
				}
				state->lex_candidate[slot] =
					state->lex_candidate[--count];
				state->lex_slot[state->lex_candidate[slot]] =
					slot;
				state->lex_slot[coordinate] = -1;
			}
		} else {
			const double *rho = state->rho;
			long double best = 0.;
			int old_count = count;

			if (position != leaving_position) {
				lp_simplex_memset(state->lex_rho, 0,
					(size_t)state->rows * sizeof(double));
				state->lex_rho[position] = 1.;
				if (simplex_basis_btran(
					    &state->factor, state->lex_rho) !=
					lp_simplex_EXIT_SUCCESS) {
					state->ratio_failure_status =
						lp_simplex_Singularity;
					return SIMPLEX_DUAL_RATIO_ERROR;
				}
				state->profile.lex_solves++;
				rho = state->lex_rho;
			}
			for (i = 0; i < old_count; i++)
				state->lex_slot[state->lex_candidate[i]] = -1;
			count = 0;
			for (i = 0; i < old_count; i++) {
				double coefficient, denominator;
				long double key;

				j = state->lex_candidate[i];
				coefficient =
					j < state->structural
						? simplex_csc_column_dot(
							  &state->matrix, j,
							  rho)
						: -rho[j - state->structural];
				denominator = -kappa * state->alpha[j];
				key = -(long double)sign * coefficient /
				      denominator;
				if (!isfinite(key)) {
					state->ratio_failure_status =
						lp_simplex_PrecisionError;
					return SIMPLEX_DUAL_RATIO_ERROR;
				}
				if (count == 0 || key < best) {
					best = key;
					count = 0;
				}
				if (key == best)
					state->lex_candidate[count++] = j;
			}
			for (i = 0; i < count; i++)
				state->lex_slot[state->lex_candidate[i]] = i;
		}
	}
	if (count == 0)
		return -1;
	j = state->lex_candidate[0];
	for (i = 1; i < count; i++)
		if (state->lex_candidate[i] < j)
			j = state->lex_candidate[i];
	state->profile.lex_choices++;
	return j;
}

/* Long-step flips consume breakpoints in increasing order.  Build the heap
 * only on a retry: the usual single-pivot Harris path keeps its linear scan.
 * Reuse the candidate array; consumed entries already have candidate_sign=0. */
static int pricing_ratio_before(
	const struct simplex_DualState *state, const int left, const int right)
{
	if (state->breakpoint[left] != state->breakpoint[right])
		return state->breakpoint[left] < state->breakpoint[right];

	if (__lp_simplex_ABS__(state->alpha[left]) !=
		__lp_simplex_ABS__(state->alpha[right]))
		return __lp_simplex_ABS__(state->alpha[left]) >
		       __lp_simplex_ABS__(state->alpha[right]);
	return left < right;
}

static void pricing_ratio_sift(struct simplex_DualState *state, int position)
{
	int variable;
	int child;

	if (position >= state->candidate_count / 2)
		return;
	variable = state->candidate_index[position];
	while (position < state->candidate_count / 2) {
		child = 2 * position + 1;
		if (child + 1 < state->candidate_count &&
			pricing_ratio_before(state,
				state->candidate_index[child + 1],
				state->candidate_index[child]))
			child++;
		if (!pricing_ratio_before(
			    state, state->candidate_index[child], variable))
			break;
		state->candidate_index[position] =
			state->candidate_index[child];
		position = child;
	}
	state->candidate_index[position] = variable;
}

static int pricing_next_breakpoint(struct simplex_DualState *state)
{
	int i;

	if (!state->ratio_heap_valid) {
		for (i = state->candidate_count / 2; i > 0; i--)
			pricing_ratio_sift(state, i - 1);
		state->ratio_heap_valid = 1;
	}
	while (state->candidate_count > 0) {
		int variable = state->candidate_index[0];

		if (state->candidate_sign[variable] != 0)
			return variable;
		state->candidate_index[0] =
			state->candidate_index[--state->candidate_count];
		if (state->candidate_count > 0)
			pricing_ratio_sift(state, 0);
	}
	return -1;
}

int simplex_dual_next_ratio_candidate(struct simplex_DualState *state,
	double *chosen_theta, const int pan_perturbation,
	const int leaving_variable, const int kappa)
{
	int candidate, j, q = -1;
	int stable_order = state->stable_candidate_order ||
			   simplex_degeneracy_is_stressed(&state->degeneracy);
	double minimum = state->ratio_minimum_valid ? state->ratio_minimum
						    : __lp_simplex_INF__;
	double best_pivot = 0.;

	if (state->ratio_failure_status != lp_simplex_Success)
		return SIMPLEX_DUAL_RATIO_ERROR;
	if (state->ratio_ordered) {
		q = pricing_next_breakpoint(state);
		if (q < 0)
			return -1;
		minimum = state->breakpoint[q];
		if (pan_perturbation &&
			simplex_degeneracy_is_lexicographic(&state->degeneracy))
			q = pricing_lexicographic_candidate(
				state, minimum, leaving_variable, kappa);
		if (q >= 0)
			*chosen_theta = state->breakpoint[q];
		return q;
	}
	if (!state->ratio_minimum_valid)
		for (candidate = 0; candidate < state->candidate_count;
			candidate++) {
			j = state->candidate_index[candidate];
			if (state->candidate_sign[j] != 0 &&
				state->breakpoint[j] < minimum)
				minimum = state->breakpoint[j];
		}
	if (minimum == __lp_simplex_INF__)
		return -1;
	if (pan_perturbation &&
		simplex_degeneracy_is_lexicographic(&state->degeneracy)) {
		q = pricing_lexicographic_candidate(
			state, minimum, leaving_variable, kappa);
		if (q >= 0)
			*chosen_theta = state->breakpoint[q];
		return q;
	}

	for (candidate = 0; candidate < state->candidate_count; candidate++) {
		double theta, relaxed;

		j = state->candidate_index[candidate];
		if (state->candidate_sign[j] == 0)
			continue;
		theta = state->breakpoint[j];
		relaxed = minimum + SIMPLEX_PRICING_TIE_RELATIVE *
					    (1. + __lp_simplex_ABS__(minimum));
		if (theta <= relaxed &&
			((!stable_order &&
				 (__lp_simplex_ABS__(state->alpha[j]) >
						 best_pivot ||
					 (__lp_simplex_ABS__(state->alpha[j]) ==
							 best_pivot &&
						 (q < 0 || j < q)))) ||
				(stable_order &&
					(pricing_significantly_greater(
						 __lp_simplex_ABS__(
							 state->alpha[j]),
						 best_pivot) ||
						(pricing_nearly_equal(
							 __lp_simplex_ABS__(
								 state->alpha
									 [j]),
							 best_pivot) &&
							(q < 0 || j < q)))))) {
			best_pivot = __lp_simplex_ABS__(state->alpha[j]);
			q = j;
			*chosen_theta = theta;
		}
	}
	return q;
}
