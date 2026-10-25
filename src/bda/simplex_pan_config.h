/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_pan_config.h
 * @brief Resolved Pan controls.
 *
 * Environment strings are parsed only by explicit configuration resolution.
 * QR and NNLS receive numerical controls, never an AUTO-mode flag.
 */

#ifndef LP_SIMPLEX_PAN_CONFIG_INTERNAL_H
#define LP_SIMPLEX_PAN_CONFIG_INTERNAL_H
#include <lp_simplex/solve.h>

/* Model representation and direction/oracle implementation are independent
 * of lp_simplex_Options.algorithm and of the QR numerical backend. */
enum simplex_PanModelForm {
	SIMPLEX_PAN_STANDARD,
	SIMPLEX_PAN_NATIVE_BOUNDS
};

enum simplex_PanDirectionMethod {
	SIMPLEX_PAN_BASIS_DIRECTION,
	SIMPLEX_PAN_CONE_DIRECTION
};

enum simplex_PanOracle {
	SIMPLEX_PAN_NO_ORACLE,
	SIMPLEX_PAN_EXPLICIT_ORACLE,
	SIMPLEX_PAN_IMPLICIT_ORACLE
};

/* Nonnegative integral work/epoch setting; malformed input disables it,
 * positive overflow saturates at INT_MAX. No process environment writes. */
int simplex_pan_integer_setting_value(const char *value);

struct simplex_PanBasis;

struct simplex_PanConfig {
	enum simplex_PanModelForm model_form;
	enum simplex_PanDirectionMethod direction_method;
	enum simplex_PanOracle oracle;
	int allow_oracle_recovery;
	int certificate_dual_repair;
	int implicit_updates, implicit_bulk;
	int compatible_pricing, stabilize, affine_multiplier, complete_start;
	int face_guard, cone_batch, first_order_enabled, dual_polish_enabled;
	int center_epoch_enabled, center_epoch_limit;
	int pdhg_iterations, restore, adaptive, predict, dual_start;
	int dual_polish_iterations, first_order_iterations, newton_qr,
		dual_center;
	int nnls_block_limit, block_crash, singleton_crash, full_refinement;
	int pdhg_weight_set;
	double pdhg_weight;
};

/* Deterministic production defaults; no environment access. Invalid
 * combinations are rejected before allocating a numerical workspace. */
int simplex_pan_config_default(struct simplex_PanConfig *config, int model_form,
	int direction_method, int oracle);
int simplex_pan_config_valid(const struct simplex_PanConfig *config);
/* Historical no-environment manual defaults, not AUTO route tuning. */
void simplex_pan_config_clean(struct simplex_PanConfig *config);
struct simplex_Problem;
void simplex_pan_config_resolve(struct simplex_PanConfig *config,
	const struct simplex_Problem *problem,
	const struct lp_simplex_Options *options);
/* NULL options resolves the legacy/manual environment defaults. */
void simplex_pan_config_init(struct simplex_PanConfig *config,
	const struct lp_simplex_Options *options);
void simplex_pan_basis_configure(
	struct simplex_PanBasis *basis, const struct simplex_PanConfig *config);
#endif
