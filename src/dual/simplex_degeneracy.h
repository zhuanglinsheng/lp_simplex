/* Event-driven controller for Pan face pivots in dual revised simplex. */
#ifndef LP_SIMPLEX_DEGENERACY_INTERNAL_H
#define LP_SIMPLEX_DEGENERACY_INTERNAL_H


enum simplex_PanPhase {
	SIMPLEX_PAN_NORMAL = 0,
	SIMPLEX_PAN_FACE = 1,
	SIMPLEX_PAN_LEXICOGRAPHIC = 2
};


struct simplex_DegeneracyControl {
	int phase;
	int recovery_pending;
	int fingerprint_count;
	int fingerprint_capacity;
	int fingerprint_initialized;
	unsigned int current_basis_hash;
	unsigned int current_status_hash;
	unsigned int *fingerprint_basis;
	unsigned int *fingerprint_status;
	unsigned char *fingerprint_used;
	long degenerate_pivots;
	long activations;
	long probes;
	long repeated_states;
	long recovery_requests;
	long recoveries;
};


void simplex_degeneracy_initialize(struct simplex_DegeneracyControl *control);

void simplex_degeneracy_destroy(struct simplex_DegeneracyControl *control);

/* A pivot enters the Pan face on zero numerical progress and leaves it on
 * strict progress.  The caller supplies the same feasibility-scaled error
 * bound used by the simplex acceptance test. */
void simplex_degeneracy_observe(
		struct simplex_DegeneracyControl *control,
		double dual_step, double dual_error_bound,
		double merit_before, double merit_after,
		double merit_error_bound);

int simplex_degeneracy_should_probe(
		const struct simplex_DegeneracyControl *control);

int simplex_degeneracy_is_stressed(
		const struct simplex_DegeneracyControl *control);

int simplex_degeneracy_is_lexicographic(
		const struct simplex_DegeneracyControl *control);

void simplex_degeneracy_request_recovery(
		struct simplex_DegeneracyControl *control);

int simplex_degeneracy_take_recovery(
		struct simplex_DegeneracyControl *control);

void simplex_degeneracy_record_recovery(
		struct simplex_DegeneracyControl *control);

void simplex_degeneracy_record_probe(
		struct simplex_DegeneracyControl *control);

void simplex_degeneracy_update_basis(
		struct simplex_DegeneracyControl *control,
		int position, int old_variable, int new_variable);

void simplex_degeneracy_update_status(
		struct simplex_DegeneracyControl *control, int variable,
		unsigned char old_status, unsigned char new_status);

/* Returns one when a state already seen on the current zero-progress face is
 * revisited.  Repetition escalates to lexicographic Pan mode; it never
 * suppresses Pan. */
int simplex_degeneracy_record_state(
		struct simplex_DegeneracyControl *control,
		const int *basis, int rows,
		const unsigned char *status, int variables);

#endif
