/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_pan_phase1.h
 * @brief Feasible starting points for the Pan/BDA solver.
 */

#ifndef LP_SIMPLEX_PAN_PHASE1_INTERNAL_H
#define LP_SIMPLEX_PAN_PHASE1_INTERNAL_H

#include "simplex_pan_basis.h"


/* Deterministic normalized-correlation shortlist, at most 32 candidates.
 * implicit_cost adds standard->objective as an unmaterialized final row. */
int simplex_pan_nnls_shortlist(const struct simplex_PanBasis *basis,
	const double *residual, const unsigned char *blocked, double tolerance,
	int limit, int *candidate, int implicit_cost);


int simplex_pan_phase1_run(struct simplex_PanBasis *basis, double *value,
	double primal_tolerance, double dual_tolerance, int iteration_limit,
	int *iterations, int *status);

/* Mixed least squares: free_mask columns may have either sign and remain in
 * the passive set. Returns 0 for an exact fit, 1 for an NNLS stationary
 * residual, -1 for numerical/budget failure. A stationary residual is not an
 * LP infeasibility claim. residual_out may be NULL. */
int simplex_pan_phase1_mixed(struct simplex_PanBasis *basis, double *value,
	const unsigned char *free_mask, double *residual_out,
	double primal_tolerance, double dual_tolerance, int iteration_limit,
	int *iterations, int *status);

/* Optional original-LP trial check after positivity restoration. Return 2
 * only when the caller certifies an actual finite descent step. */
int simplex_pan_phase1_mixed_audited(struct simplex_PanBasis *basis,
	double *value, const unsigned char *free_mask, double *residual_out,
	double primal_tolerance, double dual_tolerance, int iteration_limit,
	int *iterations, int *status, int (*accept)(void *, const double *),
	void *context);

#endif
