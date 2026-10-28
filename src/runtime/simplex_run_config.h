/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_run_config.h
 * @brief Per-solve runtime configuration and experimental control snapshots.
 */

#ifndef LP_SIMPLEX_RUN_CONFIG_INTERNAL_H
#define LP_SIMPLEX_RUN_CONFIG_INTERNAL_H
#include "simplex_pan_config.h"

struct simplex_DualControls {
	int force_bounds, disable_bounds, pan_enabled, multi_epoch;
	int initial_cost_perturb, cleanup_devex, stable_klu;
};

/* Resolved once by a public solve owner. PROFILE is output-only and excluded.
 * Numerical kernels receive values, never an AUTO/experimental mode flag. */
struct simplex_RunConfig {
	int experimental_environment, disable_presolve, disable_scaling;
	int presolve_dependencies;
	struct simplex_DualControls dual;
	struct simplex_PanConfig pan;
};
struct simplex_SolveContext;
void simplex_run_config_init(
	struct simplex_RunConfig *config, int experimental);
const struct simplex_RunConfig *simplex_run_config_get(
	const struct simplex_SolveContext *context,
	struct simplex_RunConfig *local);
unsigned long simplex_run_dual_signature(
	const struct simplex_DualControls *config);
#endif
