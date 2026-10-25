/**
 * @file test_problem_header.c
 * @brief Compile the canonical model without algorithm or journal headers.
 */

#include "simplex_execution_scope.h"
#include "simplex_problem.h"

/* This target has no algorithm include paths and does not link a solver.
 * Opaque observers and context references must not import their owners. */
#if defined(LP_SIMPLEX_EXECUTION_INTERNAL_H) || \
	defined(LP_SIMPLEX_STAGE_INTERNAL_H) || \
	defined(LP_SIMPLEX_SOLVE_CONTEXT_INTERNAL_H)
#error Canonical model headers must not import journal or solver state
#endif
