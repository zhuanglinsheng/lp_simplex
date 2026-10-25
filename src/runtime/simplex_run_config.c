/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_run_config.c
 * @brief Per-solve runtime configuration and experimental control snapshots.
 */

#include "simplex_run_config.h"
#include "simplex_solve_context.h"
#include "policy/controls.h"
#include <stdlib.h>
#include <string.h>

void simplex_run_config_init(struct simplex_RunConfig *c, int experimental)
{
	memset(c, 0, sizeof(*c));
	c->dual.disable_bounds = c->dual.pan_enabled = c->dual.multi_epoch = 1;
	simplex_pan_config_clean(&c->pan);
	if (!experimental)
		return;
	c->experimental_environment = 1;
	c->disable_presolve = getenv(SIMPLEX_ENV_DISABLE_PRESOLVE) != NULL;
	c->disable_scaling = getenv(SIMPLEX_ENV_DISABLE_SCALING) != NULL;
	c->dual.force_bounds = getenv(SIMPLEX_ENV_FORCE_BOUND_PERTURB) != NULL;
	c->dual.disable_bounds =
		getenv(SIMPLEX_ENV_DISABLE_BOUND_PERTURB) != NULL ||
		(getenv(SIMPLEX_ENV_PRIMAL_BOUND_PERTURB) == NULL &&
			!c->dual.force_bounds);
	c->dual.pan_enabled = getenv(SIMPLEX_ENV_DISABLE_PAN) == NULL;
	c->dual.multi_epoch = getenv(SIMPLEX_ENV_SINGLE_PERTURB_EPOCH) == NULL;
	c->dual.initial_cost_perturb =
		getenv(SIMPLEX_ENV_INITIAL_COST_PERTURB) != NULL;
	c->dual.cleanup_devex = getenv(SIMPLEX_ENV_CLEANUP_DEVEX) != NULL;
	c->dual.stable_klu = getenv(SIMPLEX_ENV_STABLE_KLU) != NULL;
	simplex_pan_config_init(&c->pan, NULL);
}

const struct simplex_RunConfig *simplex_run_config_get(
	const struct simplex_SolveContext *context,
	struct simplex_RunConfig *local)
{
	if (context && context->configured)
		return &context->controls;
	/* Unbound private compatibility entries retain the legacy experiment API.
	 * Public and standalone AUTO owners always bind a resolved configuration. */
	simplex_run_config_init(local, 1);
	return local;
}

unsigned long simplex_run_dual_signature(const struct simplex_DualControls *c)
{
	return (unsigned long)c->cleanup_devex |
	       ((unsigned long)c->disable_bounds << 1) |
	       ((unsigned long)!c->pan_enabled << 2) |
	       ((unsigned long)c->force_bounds << 3) |
	       ((unsigned long)c->initial_cost_perturb << 4) |
	       ((unsigned long)!c->multi_epoch << 5) |
	       ((unsigned long)c->stable_klu << 6);
}
