/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_dual_primal_pricing.c
 * @brief Primal Devex weight initialization and basis-exchange updates.
 */

#include "policy/factorization.h"
#include "simplex_dual_primal_pricing.h"
#include "simplex_dual_internal.h"
#include "utils.h"

#include <lp_simplex/status.h>

#include <math.h>

/* Devex uses a fixed reference set of nonbasic variables.  The projected
 * entering edge norm is ref[q] + sum ref[basis[i]] * direction[i]^2.
 * One BTRAN supplies all pivot-row coefficients for the weight update. */
void simplex_primal_devex_initialize(
	const struct simplex_DualState *state, double *weight, double *ref)
{
	int j;

	for (j = 0; j < state->variables; j++) {
		weight[j] = 1.;
		ref[j] = state->position[j] < 0 ? 1. : 0.;
	}
}

int simplex_primal_devex_update(struct simplex_DualState *state, double *weight,
	double *ref, const int entering, const int leaving_position)
{
	double norm = ref[entering];
	double pivot = state->direction[leaving_position];
	int i, j;
	int packed = state->work != NULL && state->packed_rho.dense != NULL;
	double *rho_work = packed ? state->work : state->rho;

	for (i = 0; i < state->rows; i++)
		norm += ref[state->basis[i]] * state->direction[i] *
			state->direction[i];
	/* Restart an inaccurate reference framework before the exchange, when
	 * the new entering reference norm is exactly one. */
	if (!isfinite(norm) ||
		norm < SIMPLEX_PRIMAL_DEVEX_RESET_RATIO * weight[entering] ||
		norm > 10. * weight[entering]) {
		simplex_primal_devex_initialize(state, weight, ref);
		norm = 1.;
	}
	/* Cleanup retains its fully refined solve; pack from independent scratch
	 * so the shared rho dense view never contains untracked coordinates. */
	lp_simplex_memset(rho_work, 0, (size_t)state->rows * sizeof(double));
	rho_work[leaving_position] = 1.;
	if (simplex_basis_btran_refined(&state->factor, rho_work) ==
		lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	if (packed) {
		(void)simplex_sparse_vector_pack(
			&state->packed_rho, rho_work, 0.);
		state->rho = state->packed_rho.dense;
	}
	for (j = 0; j < state->variables; j++) {
		double alpha, ratio, estimate;

		if (state->position[j] >= 0 || j == entering)
			continue;
		alpha = j < state->structural
				? simplex_csc_column_dot(
					  &state->matrix, j, state->rho)
				: -state->rho[j - state->structural];
		ratio = alpha / pivot;
		estimate = norm * ratio * ratio;
		if (isfinite(estimate))
			weight[j] = __lp_simplex_MAX__(weight[j], estimate);
	}
	weight[state->basis[leaving_position]] =
		__lp_simplex_MAX__(1., norm / (pivot * pivot));
	if (!isfinite(weight[state->basis[leaving_position]]))
		simplex_primal_devex_initialize(state, weight, ref);
	return lp_simplex_EXIT_SUCCESS;
}
