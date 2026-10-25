/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_dual_kernel.h
 * @brief One Dual attempt with an explicitly selected numerical backend.
 */

#ifndef LP_SIMPLEX_DUAL_KERNEL_INTERNAL_H
#define LP_SIMPLEX_DUAL_KERNEL_INTERNAL_H

struct simplex_Problem;
struct lp_simplex_Options;
struct lp_simplex_Result;
struct simplex_DualBasisSeed;
struct simplex_DualState;

/* Kernel-owned loops. Cleanup starts from a primal-feasible basis with
 * restored original costs and returns to the outer Dual loop; it does not
 * independently certify optimality. */
int simplex_dual_primal_cleanup(struct simplex_DualState *state, int *status);
int simplex_dual_continue(struct simplex_DualState *state, int *status);

/* Own the iteration state, factors, numerical safeguards and mathematical
 * classification for one attempt. Input/options/scope and optional seed
 * storage are borrowed for the call; output arrays belong to the caller.
 * The caller initializes result and supplies the complete attempt work grant.
 * This layer does not scale models, choose another backend or operate AUTO.
 * NULL row_dual is supported for the feasibility-only service. */
int simplex_dual_kernel_solve(const struct simplex_Problem *problem,
	const struct lp_simplex_Options *options, double *x, double *row_dual,
	struct lp_simplex_Result *result, int propagate_bounds,
	int basis_backend, int feasibility_only,
	struct simplex_DualBasisSeed *seed);

#endif
