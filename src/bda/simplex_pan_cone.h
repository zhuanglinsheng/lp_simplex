/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_pan_cone.h
 * @brief Tangent-cone directions and original-model acceptance checks for
 * Pan/BDA.
 */

#ifndef LP_SIMPLEX_PAN_CONE_INTERNAL_H
#define LP_SIMPLEX_PAN_CONE_INTERNAL_H
#include "simplex_pan_standard.h"

#include <lp_simplex/solve.h>

struct simplex_PanBasis;
int simplex_pan_cone_create(struct simplex_PanStandard *cone,
	const struct simplex_PanStandard *standard, double scale,
	double *row_scale, int implicit);

struct simplex_PanDirectionTrial {
	const struct simplex_PanStandard *a;
	const struct lp_simplex_Options *options;
	const double *point;
	double *trial, *residual, *work;
	int *iterations;
	int calls, checks, accepted;
};

/* Sparse original-model finite-step certification, sampled every 32 calls.
 * Does not modify point; consumes one work unit per actual audit. */
int simplex_pan_cone_trial_step(void *context, const double *direction);

long double simplex_pan_cone_cost(
	const struct simplex_PanStandard *standard, const double *point);
long double simplex_pan_cone_audit_cost(
	const struct simplex_PanStandard *standard, const double *point);
int simplex_pan_cone_check_dual(const struct simplex_PanStandard *standard,
	const double *point, const double *multiplier, double tolerance,
	double *violation);
/* Optional acceptance audits the proposed direction in the caller's model.
 * NULL accepts the raw oracle's own checks; it does not certify the input LP. */
int simplex_pan_cone_implicit_audited(struct simplex_PanBasis *basis,
	double *direction, const unsigned char *free_mask, double *residual,
	double tolerance, int limit, int *iterations, int *status,
	int (*accept)(void *, const double *), void *context);

/* Resolved controls are borrowed and remain unchanged throughout a call.
 * Numerical entries below are the actual implementations, without default-
 * resolving or report-discarding forwarding signatures. */
struct simplex_PanConfig;

/* Compress dependent positive support by feasible null-space exchanges.
 * Returns 1 if support shrank, 0 if already independent, -1 on failure;
 * failure restores the complete input point. */
int simplex_pan_support_compress_configured(const struct simplex_PanStandard *a,
	double *point, const struct lp_simplex_Options *options,
	const struct simplex_PanConfig *cfg);

/* Experimental multi-column feasible-direction variant. Existing positive
 * columns are free in an augmented NNLS problem; inactive columns are
 * nonnegative. Every accepted direction and terminal dual is independently
 * checked on the unaugmented LP. */
int simplex_pan_cone_phase2_configured(
	const struct simplex_PanStandard *standard,
	struct simplex_PanBasis *initial_basis, double *point,
	double *multiplier, const struct lp_simplex_Options *options,
	int *iterations, int *status, double *dual_infeasibility,
	const struct simplex_PanConfig *cfg);

/* 1: checked dual; 2: audited homogeneous descent; 0: prediction only.
 * Shared work and raw stop cause; no input-LP status promotion. */
int simplex_pan_cone_first_order_configured_reported(
	const struct simplex_PanStandard *a, const double *point,
	double *multiplier, const struct lp_simplex_Options *options,
	int *iterations, double *violation, double *direction,
	const unsigned char *free_mask, int warm,
	const struct simplex_PanConfig *cfg, int *service_stop);

/* Matrix-free dual Newton-CG with a caller-owned acceptance audit and raw
 * stop cause. Prediction alone is never a certificate. */
int simplex_pan_cone_predict_audited_configured_reported(
	const struct simplex_PanStandard *a, const double *point,
	double *multiplier, const struct lp_simplex_Options *options,
	int *iterations, double *violation, double *direction,
	const unsigned char *free_mask, int budget,
	int (*accept)(void *, const double *), void *context,
	const struct simplex_PanConfig *cfg, int *service_stop);

#endif
