/* Dual crash-bound initialization and event-driven interval propagation. */
#ifndef LP_SIMPLEX_DUAL_BOUNDS_INTERNAL_H
#define LP_SIMPLEX_DUAL_BOUNDS_INTERNAL_H

#include "simplex_dual_internal.h"
#include "simplex_problem.h"


void simplex_dual_initialize_bounds(
		struct simplex_DualState *state,
		const struct simplex_Problem *problem, int propagate);

#endif
