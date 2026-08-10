#include "simplex_degeneracy.h"
#include <stdio.h>


static int test_face_entry_and_progress_exit(void)
{
	struct simplex_DegeneracyControl control;
	int valid;
	simplex_degeneracy_initialize(&control);
	simplex_degeneracy_observe(&control, 0., 1e-8, 2., 1., 1e-8);
	if (simplex_degeneracy_should_probe(&control)) {
		simplex_degeneracy_destroy(&control);
		return 0;
	}
	simplex_degeneracy_observe(&control, 0., 1e-8, 2., 2., 1e-8);
	valid = simplex_degeneracy_should_probe(&control) &&
		control.activations == 1;
	simplex_degeneracy_observe(&control, 1., 1e-8, 2., 2., 1e-8);
	valid = valid && !simplex_degeneracy_should_probe(&control);
	simplex_degeneracy_destroy(&control);
	return valid;
}


static int test_repeated_state_escalates(void)
{
	int basis[3] = {3, 1, 4};
	unsigned char status[6] = {0, 1, 2, 0, 3, 1};
	struct simplex_DegeneracyControl control;
	int valid;
	simplex_degeneracy_initialize(&control);
	simplex_degeneracy_observe(&control, 0., 1e-8, 2., 2., 1e-8);
	valid = !simplex_degeneracy_record_state(
		&control, basis, 3, status, 6) &&
		simplex_degeneracy_record_state(
			&control, basis, 3, status, 6) &&
		simplex_degeneracy_is_lexicographic(&control) &&
		simplex_degeneracy_should_probe(&control) &&
		control.repeated_states == 1;
	simplex_degeneracy_destroy(&control);
	return valid;
}


static int test_incremental_fingerprint_and_growth(void)
{
	int basis[3] = {3, 1, 4};
	unsigned char status[64];
	struct simplex_DegeneracyControl control;
	int i;
	int valid = 1;
	for (i = 0; i < 64; i++)
		status[i] = 1;
	status[3] = status[1] = status[4] = 0;
	simplex_degeneracy_initialize(&control);
	simplex_degeneracy_observe(&control, 0., 1e-8, 2., 2., 1e-8);
	if (simplex_degeneracy_record_state(&control, basis, 3, status, 64))
		valid = 0;
	for (i = 5; valid && i < 40; i++) {
		int old = basis[1];
		simplex_degeneracy_update_basis(&control, 1, old, i);
		simplex_degeneracy_update_status(&control, old, status[old], 1);
		simplex_degeneracy_update_status(&control, i, status[i], 0);
		status[old] = 1;
		status[i] = 0;
		basis[1] = i;
		if (simplex_degeneracy_record_state(
			&control, basis, 3, status, 64))
			valid = 0;
	}
	valid = valid && control.fingerprint_capacity >= 64;
	simplex_degeneracy_destroy(&control);
	return valid;
}


static int test_recovery_event(void)
{
	struct simplex_DegeneracyControl control;
	int valid;
	simplex_degeneracy_initialize(&control);
	simplex_degeneracy_request_recovery(&control);
	simplex_degeneracy_request_recovery(&control);
	valid = simplex_degeneracy_should_probe(&control) &&
		control.recovery_requests == 1 &&
		simplex_degeneracy_take_recovery(&control) &&
		!simplex_degeneracy_take_recovery(&control);
	simplex_degeneracy_record_recovery(&control);
	valid = valid && control.recoveries == 1;
	simplex_degeneracy_destroy(&control);
	return valid;
}


int main(void)
{
	if (!test_face_entry_and_progress_exit() ||
	    !test_repeated_state_escalates() ||
	    !test_incremental_fingerprint_and_growth() ||
	    !test_recovery_event()) {
		printf("degeneracy controller regression failed\n");
		return 1;
	}
	printf("degeneracy controller regression passed\n");
	return 0;
}
