/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 *
 * Degeneracy detection and activation policy for Pan-style face pivots.
 */
#include "simplex_degeneracy.h"
#include "utils.h"


#define PAN_HISTORY_LENGTH 16
#define PAN_HISTORY_MASK 0xffffU
#define PAN_WINDOW_TRIGGER 15
#define PAN_STREAK_TRIGGER 12
#define PAN_EXIT_PROGRESS 1
#define PAN_COOLDOWN 32


void simplex_degeneracy_initialize(struct simplex_DegeneracyControl *control)
{
	lp_simplex_memset(control, 0, sizeof(*control));
}


static int simplex_degeneracy_history_count(unsigned int history)
{
	int count = 0;
	while (history != 0U) {
		count += (int)(history & 1U);
		history >>= 1;
	}
	return count;
}


void simplex_degeneracy_observe(
		struct simplex_DegeneracyControl *control,
		const double dual_step, const double tolerance,
		const int allow_reactivation)
{
	int degenerate = __lp_simplex_ABS__(dual_step) <= tolerance;
	control->history = ((control->history << 1) |
		(degenerate ? 1U : 0U)) & PAN_HISTORY_MASK;
	if (control->history_count < PAN_HISTORY_LENGTH)
		control->history_count++;
	if (control->cooldown > 0)
		control->cooldown--;
	if (degenerate) {
		control->consecutive++;
		control->improving = 0;
		control->degenerate_pivots++;
	} else {
		control->consecutive = 0;
		control->improving++;
	}
	if (control->active) {
		control->active_iterations++;
		if (control->improving >= PAN_EXIT_PROGRESS) {
			control->active = 0;
			control->active_iterations = 0;
			control->cooldown = PAN_COOLDOWN;
		} else if (control->active_iterations >= 64) {
			control->active = 0;
			control->active_iterations = 0;
			if (allow_reactivation)
				control->cooldown = PAN_COOLDOWN;
			else
				control->suppressed = 1;
		}
		return;
	}
	if (control->cooldown == 0 &&
	    !control->suppressed &&
	    (control->consecutive >= PAN_STREAK_TRIGGER ||
	     (control->history_count == PAN_HISTORY_LENGTH &&
	      simplex_degeneracy_history_count(control->history) >=
	      PAN_WINDOW_TRIGGER))) {
		control->active = 1;
		control->active_iterations = 0;
		control->improving = 0;
		control->activations++;
	}
}


int simplex_degeneracy_should_probe(
		const struct simplex_DegeneracyControl *control)
{
	return control->active && !control->suppressed;
}


int simplex_degeneracy_is_stressed(
		const struct simplex_DegeneracyControl *control)
{
	/* Switch ratio policy before full Pan activation.  A half-degenerate
	 * recent window is enough evidence that a wider Harris set can change the
	 * face-walking trajectory more than it improves pivot stability. */
	return control->active || control->consecutive >= 4 ||
		(control->history_count == PAN_HISTORY_LENGTH &&
		 simplex_degeneracy_history_count(control->history) >= 8);
}


void simplex_degeneracy_record_probe(struct simplex_DegeneracyControl *control)
{
	control->probes++;
}


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


void simplex_degeneracy_update_basis(
		struct simplex_DegeneracyControl *control,
		const int position, const int old_variable, const int new_variable)
{
	if (!control->fingerprint_initialized || old_variable == new_variable)
		return;
	control->current_basis_hash ^=
		simplex_degeneracy_basis_item(position, old_variable) ^
		simplex_degeneracy_basis_item(position, new_variable);
}


void simplex_degeneracy_update_status(
		struct simplex_DegeneracyControl *control,
		const int variable, const unsigned char old_status,
		const unsigned char new_status)
{
	if (!control->fingerprint_initialized || old_status == new_status)
		return;
	control->current_status_hash ^=
		simplex_degeneracy_status_item(variable, old_status) ^
		simplex_degeneracy_status_item(variable, new_status);
}


int simplex_degeneracy_record_state(
		struct simplex_DegeneracyControl *control,
		const int *basis, const int rows,
		const unsigned char *status, const int variables)
{
	int i;
	unsigned int basis_hash;
	unsigned int status_hash;
	if (!control->fingerprint_initialized) {
		control->current_basis_hash = 0U;
		control->current_status_hash = 0U;
		for (i = 0; i < rows; i++)
			control->current_basis_hash ^=
				simplex_degeneracy_basis_item(i, basis[i]);
		for (i = 0; i < variables; i++)
			control->current_status_hash ^=
				simplex_degeneracy_status_item(i, status[i]);
		control->fingerprint_initialized = 1;
	}
	basis_hash = control->current_basis_hash;
	status_hash = control->current_status_hash;
	for (i = 0; i < control->fingerprint_count; i++) {
		if (control->fingerprint_basis[i] == basis_hash &&
		    control->fingerprint_status[i] == status_hash) {
			control->active = 0;
			control->suppressed = 1;
			return 1;
		}
	}
	control->fingerprint_basis[control->fingerprint_next] = basis_hash;
	control->fingerprint_status[control->fingerprint_next] = status_hash;
	control->fingerprint_next = (control->fingerprint_next + 1) & 63;
	if (control->fingerprint_count < 64)
		control->fingerprint_count++;
	return 0;
}
