#include "simplex_degeneracy.h"
#include <stdio.h>


static int test_streak_activation(void)
{
	int i;
	struct simplex_DegeneracyControl control;
	simplex_degeneracy_initialize(&control);
	for (i = 0; i < 11; i++)
		simplex_degeneracy_observe(&control, 0., 1e-8, 1);
	if (simplex_degeneracy_should_probe(&control))
		return 0;
	simplex_degeneracy_observe(&control, 0., 1e-8, 1);
	if (!simplex_degeneracy_should_probe(&control))
		return 0;
	simplex_degeneracy_observe(&control, 1., 1e-8, 1);
	return !simplex_degeneracy_should_probe(&control);
}


static int test_window_activation(void)
{
	int i;
	struct simplex_DegeneracyControl control;
	simplex_degeneracy_initialize(&control);
	for (i = 0; i < 16; i++)
		simplex_degeneracy_observe(&control, i == 7 ? 1. : 0., 1e-8, 1);
	return simplex_degeneracy_should_probe(&control);
}


static int test_repeated_state_suppresses(void)
{
	int basis[3] = {3, 1, 4};
	unsigned char status[6] = {0, 1, 2, 0, 3, 1};
	struct simplex_DegeneracyControl control;
	simplex_degeneracy_initialize(&control);
	if (simplex_degeneracy_record_state(&control, basis, 3, status, 6))
		return 0;
	if (!simplex_degeneracy_record_state(&control, basis, 3, status, 6))
		return 0;
	return control.suppressed;
}


static int test_activation_budget(void)
{
	int i;
	struct simplex_DegeneracyControl control;
	simplex_degeneracy_initialize(&control);
	for (i = 0; i < 12; i++)
		simplex_degeneracy_observe(&control, 0., 1e-8, 1);
	if (!simplex_degeneracy_should_probe(&control))
		return 0;
	for (i = 0; i < 64; i++)
		simplex_degeneracy_observe(&control, 0., 1e-8, 1);
	return !control.suppressed && control.cooldown > 0 &&
		!simplex_degeneracy_should_probe(&control);
}


static int test_noncompact_budget_suppresses(void)
{
	int i;
	struct simplex_DegeneracyControl control;
	simplex_degeneracy_initialize(&control);
	for (i = 0; i < 12; i++)
		simplex_degeneracy_observe(&control, 0., 1e-8, 0);
	if (!simplex_degeneracy_should_probe(&control))
		return 0;
	for (i = 0; i < 64; i++)
		simplex_degeneracy_observe(&control, 0., 1e-8, 0);
	return control.suppressed && !simplex_degeneracy_should_probe(&control);
}


static int test_incremental_fingerprint(void)
{
	int basis[3] = {3, 1, 4};
	unsigned char status[6] = {0, 1, 2, 0, 3, 1};
	struct simplex_DegeneracyControl control;
	simplex_degeneracy_initialize(&control);
	if (simplex_degeneracy_record_state(&control, basis, 3, status, 6))
		return 0;
	simplex_degeneracy_update_basis(&control, 1, 1, 5);
	simplex_degeneracy_update_status(&control, 1, 1, 2);
	simplex_degeneracy_update_status(&control, 5, 1, 0);
	basis[1] = 5;
	status[1] = 2;
	status[5] = 0;
	if (simplex_degeneracy_record_state(&control, basis, 3, status, 6))
		return 0;
	simplex_degeneracy_update_basis(&control, 1, 5, 1);
	simplex_degeneracy_update_status(&control, 1, 2, 1);
	simplex_degeneracy_update_status(&control, 5, 0, 1);
	basis[1] = 1;
	status[1] = 1;
	status[5] = 1;
	return simplex_degeneracy_record_state(&control, basis, 3, status, 6);
}


int main(void)
{
	if (!test_streak_activation() || !test_window_activation() ||
	    !test_repeated_state_suppresses() || !test_activation_budget() ||
	    !test_noncompact_budget_suppresses() ||
	    !test_incremental_fingerprint()) {
		printf("degeneracy controller regression failed\n");
		return 1;
	}
	printf("degeneracy controller regression passed\n");
	return 0;
}
