/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_degeneracy.h
 * @brief Event-driven controller for Pan face pivots in dual revised simplex.
 */

#ifndef LP_SIMPLEX_DEGENERACY_INTERNAL_H
#define LP_SIMPLEX_DEGENERACY_INTERNAL_H


enum simplex_PanPhase {
	SIMPLEX_PAN_NORMAL = 0,
	SIMPLEX_PAN_FACE = 1,
	SIMPLEX_PAN_LEXICOGRAPHIC = 2
};


struct simplex_DegeneracyControl {
	int phase;
	int stagnant_face;
	int recovery_pending;
	int fingerprint_count;
	int fingerprint_capacity;
	int fingerprint_initialized;
	double face_best_merit;
	unsigned int current_basis_hash;
	unsigned int current_status_hash;
	unsigned int *fingerprint_basis;
	unsigned int *fingerprint_status;
	unsigned char *fingerprint_used;
	long degenerate_pivots;
	long activations;
	long repeated_states;
	long recovery_requests;
	long recoveries;
	/* Stagnant faces where every tied variable had already reached its
	 * perturbation cap, so tie-breaking escalated to lexicographic mode
	 * without waiting for an exact state repetition. */
	long lex_escalations;
};


void simplex_degeneracy_initialize(struct simplex_DegeneracyControl *control);

void simplex_degeneracy_destroy(struct simplex_DegeneracyControl *control);

/* Clear basis-dependent fingerprints while retaining the allocated hash table
 * and lifetime counters. */
void simplex_degeneracy_reset_state(struct simplex_DegeneracyControl *control);

/* Retain the whole zero-objective face, even when primal merit improves.
 * Signal stagnation when merit fails to beat the best seen on that face.
 * A certified objective step leaves the face.  The caller supplies the same
 * feasibility-scaled error bound used by the simplex acceptance test. */
void simplex_degeneracy_observe(struct simplex_DegeneracyControl *control,
	double dual_step, double dual_error_bound, double merit_before,
	double merit_after, double merit_error_bound);

int simplex_degeneracy_should_probe(
	const struct simplex_DegeneracyControl *control);

int simplex_degeneracy_is_stressed(
	const struct simplex_DegeneracyControl *control);

int simplex_degeneracy_is_lexicographic(
	const struct simplex_DegeneracyControl *control);

/* Escalate to lexicographic tie-breaking when perturbation is exhausted. */
void simplex_degeneracy_escalate_lexicographic(
	struct simplex_DegeneracyControl *control);

void simplex_degeneracy_request_recovery(
	struct simplex_DegeneracyControl *control);

int simplex_degeneracy_take_recovery(struct simplex_DegeneracyControl *control);

void simplex_degeneracy_record_recovery(
	struct simplex_DegeneracyControl *control);

void simplex_degeneracy_update_basis(struct simplex_DegeneracyControl *control,
	int position, int old_variable, int new_variable);

void simplex_degeneracy_update_status(struct simplex_DegeneracyControl *control,
	int variable, unsigned char old_status, unsigned char new_status);

/* Returns one when a state already seen on the current zero-progress face is
 * revisited.  Repetition escalates to lexicographic Pan mode; it never
 * suppresses Pan. */
int simplex_degeneracy_record_state(struct simplex_DegeneracyControl *control,
	const int *basis, int rows, const unsigned char *status, int variables);

#define SIMPLEX_DUAL_PROGRESS_WINDOW 128
#define SIMPLEX_DUAL_EXPLORE_INITIAL 1024
#define SIMPLEX_DUAL_RECOVERY_ALLOWANCE (16 * SIMPLEX_DUAL_EXPLORE_INITIAL)

/* Bounded leases inside one live solve, never a restart. A lease can renew
 * on local residual reduction or once after a change of measurement space.
 * The hard budget belongs to the caller, so oscillations cannot renew forever. */
enum simplex_DualExploreReason {
	SIMPLEX_DUAL_EXPLORE_CONTINUE = 0,
	SIMPLEX_DUAL_EXPLORE_NUMERICAL,
	SIMPLEX_DUAL_EXPLORE_CYCLE,
	SIMPLEX_DUAL_EXPLORE_STAGNATION
};
struct simplex_DualExplore {
	int checkpoint, limit, phase, initialized;
	int recovery_until, reason;
	long epoch, cycles;
	/* Separate measurement spaces survive working-bound recovery epochs. */
	double baseline[5], secondary[5];
	unsigned int seen, lease_seen;
};
int simplex_dual_explore_yield(struct simplex_DualExplore *history,
	int iterations, int hard_limit, int phase, long epoch, double residual,
	double noise, long cycles);

/* A second, independently scaled progress signal (violation count or original
 * objective) can justify continuation even while the maximum residual rises. */
int simplex_dual_explore_yield_pair(struct simplex_DualExplore *history,
	int iterations, int hard_limit, int phase, long epoch, double residual,
	double noise, long cycles, double secondary, double secondary_noise);

#endif
