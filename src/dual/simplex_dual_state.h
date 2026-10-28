/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_dual_state.h
 * @brief Internal lifecycle for the dual revised-simplex state.
 */

#ifndef LP_SIMPLEX_DUAL_STATE_LIFECYCLE_INTERNAL_H
#define LP_SIMPLEX_DUAL_STATE_LIFECYCLE_INTERNAL_H

struct lp_simplex_Options;
struct simplex_DualState;
struct simplex_Problem;

int simplex_dual_state_create(struct simplex_DualState *state,
	const struct simplex_Problem *problem,
	const struct lp_simplex_Options *options, int basis_backend);

void simplex_dual_state_destroy(struct simplex_DualState *state);
/* Optional, one-attempt allocation; false retains the established dense RHS. */
int simplex_dual_packed_work_prepare(struct simplex_DualState *state);

#endif
