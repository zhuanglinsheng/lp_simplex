/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#ifndef LP_SIMPLEX_DEGENERACY_INTERNAL_H
#define LP_SIMPLEX_DEGENERACY_INTERNAL_H


/* Controller for the bounded-deficiency (Pan) anti-stalling mode. */
struct simplex_DegeneracyControl {
	unsigned int history;
	int history_count;
	int consecutive;
	int improving;
	int cooldown;
	int active;
	int active_iterations;
	int suppressed;
	int fingerprint_count;
	int fingerprint_next;
	int fingerprint_initialized;
	unsigned int current_basis_hash;
	unsigned int current_status_hash;
	unsigned int fingerprint_basis[64];
	unsigned int fingerprint_status[64];
	long degenerate_pivots;
	long activations;
	long probes;
};


void simplex_degeneracy_initialize(struct simplex_DegeneracyControl *control);

void simplex_degeneracy_observe(
		struct simplex_DegeneracyControl *control,
		double dual_step, double tolerance, int allow_reactivation);

int simplex_degeneracy_should_probe(
		const struct simplex_DegeneracyControl *control);

int simplex_degeneracy_is_stressed(
		const struct simplex_DegeneracyControl *control);

void simplex_degeneracy_record_probe(
		struct simplex_DegeneracyControl *control);

void simplex_degeneracy_update_basis(
		struct simplex_DegeneracyControl *control,
		int position, int old_variable, int new_variable);

void simplex_degeneracy_update_status(
		struct simplex_DegeneracyControl *control, int variable,
		unsigned char old_status, unsigned char new_status);

int simplex_degeneracy_record_state(
		struct simplex_DegeneracyControl *control,
		const int *basis, int rows,
		const unsigned char *status, int variables);

#endif
