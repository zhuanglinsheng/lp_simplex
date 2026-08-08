#ifndef LP_SIMPLEX_DUAL_PRICING_INTERNAL_H
#define LP_SIMPLEX_DUAL_PRICING_INTERNAL_H

#include "simplex_dual_internal.h"

int simplex_dual_pricing_create(struct simplex_DualState *state);
void simplex_dual_pricing_destroy(struct simplex_DualState *state);
void simplex_dual_nonbasic_add(struct simplex_DualState *state, int variable);
void simplex_dual_nonbasic_remove(struct simplex_DualState *state, int variable);
void simplex_dual_nonbasic_initialize(struct simplex_DualState *state);
int simplex_dual_nonbasic_is_consistent(const struct simplex_DualState *state);
int simplex_dual_prepare_ratio_test(struct simplex_DualState *state, int kappa);
int simplex_dual_next_ratio_candidate(struct simplex_DualState *state,
		double *chosen_theta, int pan_perturbation, int leaving_variable);

#endif
