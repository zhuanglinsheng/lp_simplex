/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_degeneracy.c
 * @brief Event-driven Pan face-state controller.
 */

#include "simplex_degeneracy.h"
#include "utils.h"

#include <lp_simplex/status.h>

#include <math.h>


static unsigned int simplex_degeneracy_mix(unsigned int value)
{
	value ^= value >> 16;
	value *= 0x7feb352dU;
	value ^= value >> 15;
	value *= 0x846ca68bU;
	value ^= value >> 16;
	return value;
}


static unsigned int simplex_degeneracy_basis_item(
	const int position, const int variable)
{
	return simplex_degeneracy_mix((unsigned int)position * 0x9e3779b9U ^
				      (unsigned int)variable + 0x85ebca6bU);
}


static unsigned int simplex_degeneracy_status_item(
	const int variable, const unsigned char status)
{
	return simplex_degeneracy_mix((unsigned int)variable * 0xc2b2ae35U ^
				      (unsigned int)status + 0x27d4eb2fU);
}


static unsigned int simplex_degeneracy_pair_hash(
	const unsigned int basis, const unsigned int status)
{
	return simplex_degeneracy_mix(
		basis ^ simplex_degeneracy_mix(status + 0x9e3779b9U));
}


static int simplex_degeneracy_table_reserve(
	struct simplex_DegeneracyControl *control, const int capacity)
{
	unsigned int *basis;
	unsigned int *status;
	unsigned char *used;
	int next = 16;
	int i;

	while (next < capacity)
		next *= 2;
	basis = (unsigned int *)lp_simplex_malloc(
		(size_t)next * sizeof(unsigned int));
	status = (unsigned int *)lp_simplex_malloc(
		(size_t)next * sizeof(unsigned int));
	used = (unsigned char *)lp_simplex_malloc((size_t)next);
	if (basis == NULL || status == NULL || used == NULL) {
		lp_simplex_free(basis);
		lp_simplex_free(status);
		lp_simplex_free(used);
		return lp_simplex_EXIT_FAILURE;
	}
	lp_simplex_memset(used, 0, (size_t)next);
	for (i = 0; i < control->fingerprint_capacity; i++)
		if (control->fingerprint_used[i]) {
			unsigned int hash = simplex_degeneracy_pair_hash(
				control->fingerprint_basis[i],
				control->fingerprint_status[i]);
			int slot = (int)(hash & (unsigned int)(next - 1));

			while (used[slot])
				slot = (slot + 1) & (next - 1);
			used[slot] = 1;
			basis[slot] = control->fingerprint_basis[i];
			status[slot] = control->fingerprint_status[i];
		}
	lp_simplex_free(control->fingerprint_basis);
	lp_simplex_free(control->fingerprint_status);
	lp_simplex_free(control->fingerprint_used);
	control->fingerprint_basis = basis;
	control->fingerprint_status = status;
	control->fingerprint_used = used;
	control->fingerprint_capacity = next;
	return lp_simplex_EXIT_SUCCESS;
}


static void simplex_degeneracy_clear_face(
	struct simplex_DegeneracyControl *control)
{
	control->fingerprint_count = 0;

	if (control->fingerprint_used != NULL)
		lp_simplex_memset(control->fingerprint_used, 0,
			(size_t)control->fingerprint_capacity);
}


void simplex_degeneracy_initialize(struct simplex_DegeneracyControl *control)
{
	lp_simplex_memset(control, 0, sizeof(*control));
}


void simplex_degeneracy_destroy(struct simplex_DegeneracyControl *control)
{
	if (control == NULL)
		return;
	lp_simplex_free(control->fingerprint_basis);
	lp_simplex_free(control->fingerprint_status);
	lp_simplex_free(control->fingerprint_used);
	lp_simplex_memset(control, 0, sizeof(*control));
}


void simplex_degeneracy_reset_state(struct simplex_DegeneracyControl *control)
{
	control->phase = SIMPLEX_PAN_NORMAL;
	control->stagnant_face = 0;
	control->recovery_pending = 0;
	control->fingerprint_count = 0;
	control->fingerprint_initialized = 0;
	control->current_basis_hash = 0;
	control->current_status_hash = 0;
	if (control->fingerprint_used != NULL)
		lp_simplex_memset(control->fingerprint_used, 0,
			(size_t)control->fingerprint_capacity);
}


void simplex_degeneracy_observe(struct simplex_DegeneracyControl *control,
	const double dual_step, const double dual_error_bound,
	const double merit_before, const double merit_after,
	const double merit_error_bound)
{
	int zero_objective = __lp_simplex_ABS__(dual_step) <= dual_error_bound;

	if (!isfinite(dual_step) || !isfinite(dual_error_bound) ||
		!isfinite(merit_before) || !isfinite(merit_after) ||
		!isfinite(merit_error_bound) || dual_error_bound < 0. ||
		merit_error_bound < 0.) {
		simplex_degeneracy_request_recovery(control);
		return;
	}
	/* Keep history for the entire dual-objective face.  Progress in primal
	 * merit is useful, but it can regress without leaving that face. */
	if (zero_objective)
		control->degenerate_pivots++;
	if (zero_objective) {
		if (control->phase == SIMPLEX_PAN_NORMAL) {
			control->phase = SIMPLEX_PAN_FACE;
			control->activations++;
			control->face_best_merit = merit_before;
			simplex_degeneracy_clear_face(control);
		}
		if (merit_after + merit_error_bound < control->face_best_merit)
			control->face_best_merit = merit_after;
		else
			control->stagnant_face = 1;
		return;
	}
	control->phase = SIMPLEX_PAN_NORMAL;
	control->stagnant_face = 0;
	/* Progress is not evidence that an independently requested numerical
	 * rebuild happened. Only take_recovery() consumes that event. */
	simplex_degeneracy_clear_face(control);
}


int simplex_degeneracy_should_probe(
	const struct simplex_DegeneracyControl *control)
{
	return control->phase != SIMPLEX_PAN_NORMAL;
}


int simplex_degeneracy_is_stressed(
	const struct simplex_DegeneracyControl *control)
{
	/* A newly detected flat face is not numerical stress: retain Harris'
	 * stability envelope while Pan changes only the tie policy.  Widen the
	 * candidate set only after an actual cycle or consistency failure. */
	return control->phase == SIMPLEX_PAN_LEXICOGRAPHIC ||
	       control->recovery_pending;
}


int simplex_degeneracy_is_lexicographic(
	const struct simplex_DegeneracyControl *control)
{
	return control->phase == SIMPLEX_PAN_LEXICOGRAPHIC;
}


/* Cost perturbation can no longer separate the current face (every tied
 * variable is saturated).  Escalate to symbolic lexicographic tie-breaking
 * directly instead of waiting for an exact basis/status repetition. */
void simplex_degeneracy_escalate_lexicographic(
	struct simplex_DegeneracyControl *control)
{
	if (control->phase != SIMPLEX_PAN_LEXICOGRAPHIC) {
		control->phase = SIMPLEX_PAN_LEXICOGRAPHIC;
		control->lex_escalations++;
	}
}


void simplex_degeneracy_request_recovery(
	struct simplex_DegeneracyControl *control)
{
	if (control->phase == SIMPLEX_PAN_NORMAL) {
		control->phase = SIMPLEX_PAN_FACE;
		control->face_best_merit = __lp_simplex_INF__;
		control->activations++;
	}
	if (!control->recovery_pending) {
		control->recovery_pending = 1;
		control->recovery_requests++;
	}
}


int simplex_degeneracy_take_recovery(struct simplex_DegeneracyControl *control)
{
	if (!control->recovery_pending)
		return 0;

	control->recovery_pending = 0;
	return 1;
}


void simplex_degeneracy_record_recovery(
	struct simplex_DegeneracyControl *control)
{
	control->recoveries++;
}


void simplex_degeneracy_update_basis(struct simplex_DegeneracyControl *control,
	const int position, const int old_variable, const int new_variable)
{
	if (!control->fingerprint_initialized || old_variable == new_variable)
		return;

	control->current_basis_hash ^=
		simplex_degeneracy_basis_item(position, old_variable) ^
		simplex_degeneracy_basis_item(position, new_variable);
}


void simplex_degeneracy_update_status(struct simplex_DegeneracyControl *control,
	const int variable, const unsigned char old_status,
	const unsigned char new_status)
{
	if (!control->fingerprint_initialized || old_status == new_status)
		return;

	control->current_status_hash ^=
		simplex_degeneracy_status_item(variable, old_status) ^
		simplex_degeneracy_status_item(variable, new_status);
}


int simplex_degeneracy_record_state(struct simplex_DegeneracyControl *control,
	const int *basis, const int rows, const unsigned char *status,
	const int variables)
{
	unsigned int hash;
	int slot;
	int i;

	if (!control->fingerprint_initialized) {
		for (i = 0; i < rows; i++)
			control->current_basis_hash ^=
				simplex_degeneracy_basis_item(i, basis[i]);
		for (i = 0; i < variables; i++)
			control->current_status_hash ^=
				simplex_degeneracy_status_item(i, status[i]);
		control->fingerprint_initialized = 1;
	}
	if (control->fingerprint_capacity == 0 ||
		(control->fingerprint_count + 1) * 2 >=
			control->fingerprint_capacity)
		if (simplex_degeneracy_table_reserve(
			    control, control->fingerprint_capacity > 0
					     ? 2 * control->fingerprint_capacity
					     : 16) == lp_simplex_EXIT_FAILURE)
			return 0;
	hash = simplex_degeneracy_pair_hash(
		control->current_basis_hash, control->current_status_hash);
	slot = (int)(hash & (unsigned int)(control->fingerprint_capacity - 1));
	while (control->fingerprint_used[slot]) {
		if (control->fingerprint_basis[slot] ==
				control->current_basis_hash &&
			control->fingerprint_status[slot] ==
				control->current_status_hash) {
			control->phase = SIMPLEX_PAN_LEXICOGRAPHIC;
			control->repeated_states++;
			return 1;
		}
		slot = (slot + 1) & (control->fingerprint_capacity - 1);
	}
	control->fingerprint_used[slot] = 1;
	control->fingerprint_basis[slot] = control->current_basis_hash;
	control->fingerprint_status[slot] = control->current_status_hash;
	control->fingerprint_count++;
	return 0;
}


int simplex_dual_explore_yield_pair(struct simplex_DualExplore *h,
	int iterations, int hard_limit, int phase, long epoch, double residual,
	double noise, long cycles, double secondary, double secondary_noise)
{
	int changed, next, slot = phase == 0 ? 0 : phase - 9;
	int secondary_progress;
	unsigned int bit;

	if (slot < 0 || slot > 4)
		slot = 0;
	bit = 1U << slot;

	h->reason = SIMPLEX_DUAL_EXPLORE_CONTINUE;
	if (!isfinite(residual) || residual < 0.) {
		h->reason = SIMPLEX_DUAL_EXPLORE_NUMERICAL;
		return 1;
	}
	if (!h->limit)
		h->limit = hard_limit < SIMPLEX_DUAL_EXPLORE_INITIAL
				   ? hard_limit
				   : SIMPLEX_DUAL_EXPLORE_INITIAL;
	if (!(h->seen & bit)) {
		h->baseline[slot] = residual;
		h->secondary[slot] = secondary;
		h->seen |= bit;
		/* The initial measurement starts the first observation lease. */
		if (!h->initialized)
			h->lease_seen |= bit;
	}
	changed = !h->initialized || h->phase != phase || h->epoch != epoch;
	if (changed) {
		h->initialized = 1;
		h->phase = phase;
		h->epoch = epoch;
		h->checkpoint = iterations;
	}
	if (iterations - h->checkpoint >= SIMPLEX_DUAL_PROGRESS_WINDOW)
		h->checkpoint = iterations;
	if (iterations < h->limit || iterations >= hard_limit)
		return 0;
	/* A phase receives one observation lease. Subsequent renewals must beat
	 * its best value at the previous lease boundary, even after bound epochs
	 * or phase oscillations. A return from a worse peak buys no progress. */
	secondary_progress = isfinite(secondary) &&
			     isfinite(h->secondary[slot]) &&
			     h->secondary[slot] - secondary > secondary_noise &&
			     h->secondary[slot] - secondary >
				     .01 * (1. + fabs(h->secondary[slot]));
	if (cycles - h->cycles >= 2) {
		h->reason = SIMPLEX_DUAL_EXPLORE_CYCLE;
		return 1;
	}
	if ((h->lease_seen & bit) && iterations >= h->recovery_until &&
		residual > noise && !secondary_progress &&
		(h->baseline[slot] - residual <= noise ||
			residual >= .99 * h->baseline[slot])) {
		h->reason = SIMPLEX_DUAL_EXPLORE_STAGNATION;
		return 1;
	}
	next = h->limit > hard_limit / 2 ? hard_limit : 2 * h->limit;
	h->limit = next;
	h->cycles = cycles;
	h->lease_seen = h->seen;
	if (secondary < h->secondary[slot])
		h->secondary[slot] = secondary;
	if (residual < h->baseline[slot])
		h->baseline[slot] = residual;
	return 0;
}


int simplex_dual_explore_yield(struct simplex_DualExplore *h, int iterations,
	int hard_limit, int phase, long epoch, double residual, double noise,
	long cycles)
{
	return simplex_dual_explore_yield_pair(h, iterations, hard_limit, phase,
		epoch, residual, noise, cycles, HUGE_VAL, 0.);
}
