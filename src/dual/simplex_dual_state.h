/* Internal lifecycle for the dual revised-simplex state. */
#ifndef LP_SIMPLEX_DUAL_STATE_LIFECYCLE_INTERNAL_H
#define LP_SIMPLEX_DUAL_STATE_LIFECYCLE_INTERNAL_H

struct lp_simplex_Options;
struct simplex_DualState;
struct simplex_Problem;

int simplex_dual_state_create(
		struct simplex_DualState *state,
		const struct simplex_Problem *problem,
		const struct lp_simplex_Options *options);

void simplex_dual_state_destroy(struct simplex_DualState *state);

#endif
