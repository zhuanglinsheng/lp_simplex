/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_dual_primal_pricing.h
 * @brief Private Devex pricing for the Dual kernel's primal cleanup.
 */

#ifndef LP_SIMPLEX_DUAL_PRIMAL_PRICING_INTERNAL_H
#define LP_SIMPLEX_DUAL_PRIMAL_PRICING_INTERNAL_H

struct simplex_DualState;

/* Caller owns two arrays of variables doubles.  Update before exchanging the
 * basis; direction must contain the FTRAN of the entering column. */
void simplex_primal_devex_initialize(
	const struct simplex_DualState *state, double *weight, double *ref);
int simplex_primal_devex_update(struct simplex_DualState *state, double *weight,
	double *ref, int entering, int leaving_position);

#endif
