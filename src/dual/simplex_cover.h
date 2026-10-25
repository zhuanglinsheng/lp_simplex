/**
 * @file simplex_cover.h
 * @brief Internal restricted-column solver for nonnegative covering LPs.
 */

#ifndef LP_SIMPLEX_COVER_INTERNAL_H
#define LP_SIMPLEX_COVER_INTERNAL_H
#include "simplex_problem.h"
#include <lp_simplex/solve.h>

/* Failure is nonterminal for the full LP. Only full primal/dual certification
 * may return Success; an infeasible restricted master is never a full proof.
 * Work counts master pivots plus one per full pricing pass. */
int simplex_cover_solve(const struct simplex_Problem *problem,
	const struct lp_simplex_Options *options, double *x, double *row_dual,
	struct lp_simplex_Result *result);
#endif
