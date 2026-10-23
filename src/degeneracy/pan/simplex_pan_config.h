/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_pan_config.h
 * @brief Resolved Pan controls.
 *
 * Environment strings are parsed at kernel entry; QR and NNLS receive numerical
 * controls, never an AUTO-mode flag.
 */

#ifndef LP_SIMPLEX_PAN_CONFIG_INTERNAL_H
#define LP_SIMPLEX_PAN_CONFIG_INTERNAL_H
#include <lp_simplex/solve.h>
/* Internal route tag; public input validation never accepts this value. */
#define SIMPLEX_PAN_AUTO_EXPLICIT (-100)
#define SIMPLEX_PAN_AUTO_STANDARD_IMPLICIT (-101)
#define SIMPLEX_PAN_AUTO_STANDARD_EXPLICIT (-102)
const char *simplex_pan_setting(
	const struct lp_simplex_Options *options, const char *name);

struct simplex_PanBasis;
struct simplex_PanConfig {
	int implicit_cone, implicit_updates, implicit_bulk;
	int pdhg_iterations, restore, adaptive, predict, dual_start;
	int dual_polish_iterations, first_order_iterations, newton_qr,
		dual_center;
	int nnls_block_limit, block_crash, singleton_crash, full_refinement;
	int pdhg_weight_set;
	double pdhg_weight;
};
/* NULL options resolves the legacy/manual environment defaults. */
void simplex_pan_config_init(struct simplex_PanConfig *config,
	const struct lp_simplex_Options *options);
void simplex_pan_basis_configure(
	struct simplex_PanBasis *basis, const struct simplex_PanConfig *config);
#endif
