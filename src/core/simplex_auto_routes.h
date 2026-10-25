/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_auto_routes.h
 * @brief AUTO capability checks, resource admission and kernel adapters.
 */

#ifndef LP_SIMPLEX_AUTO_ROUTES_INTERNAL_H
#define LP_SIMPLEX_AUTO_ROUTES_INTERNAL_H

#include "simplex_auto.h"
#include "simplex_dual.h"
#include "simplex_problem_identity.h"

int simplex_auto_route_solve(int route, const struct simplex_Problem *problem,
	const struct lp_simplex_Options *options, const double *seed, double *x,
	struct lp_simplex_Result *result);

/* Resource eligibility, not a mathematical terminal classification. */
int simplex_auto_route_allowed(
	const struct simplex_Problem *problem, int route);

struct simplex_Presolve;
struct simplex_PanConfig;

/* Borrowed, solve-local binding. Model/configuration/presolve map remain
 * immutable. original is an n-column workspace that must not overlap the
 * reduced output buffer. target is captured before AUTO installs its read-only view. */
struct simplex_AutoOriginalInitializer {
	const struct lp_Model *model;
	struct lp_Model binding;
	const struct simplex_Presolve *presolve;
	const struct simplex_PanConfig *controls;
	double *original;
	struct simplex_ExecutionScope scope;
	struct simplex_ProblemIdentity target;
};

int simplex_auto_routes_initialize_original(void *context,
	const struct simplex_Problem *input,
	const struct lp_simplex_Options *options, double *point, int *work);
int simplex_auto_routes_initialize(const struct simplex_AutoInitializer *,
	const struct simplex_Problem *, const struct lp_simplex_Options *,
	double *, int *);

/* First blocking model/resource check; estimates are not total RSS bounds. */
enum simplex_AutoModelAdmission {
	SIMPLEX_AUTO_MODEL_UNCHECKED = 0,
	SIMPLEX_AUTO_MODEL_READY,
	SIMPLEX_AUTO_MODEL_ROUTE,
	SIMPLEX_AUTO_MODEL_SHAPE,
	SIMPLEX_AUTO_MODEL_INDEX,
	SIMPLEX_AUTO_MODEL_CONVERSION,
	SIMPLEX_AUTO_MODEL_DENSE,
	SIMPLEX_AUTO_MODEL_DENSE_QR,
	SIMPLEX_AUTO_MODEL_SINGLETON_BOUNDS,
	SIMPLEX_AUTO_MODEL_SINGLETON_ROWS,
	SIMPLEX_AUTO_MODEL_SINGLETON_PAIRS
};

int simplex_auto_route_admission(const struct simplex_Problem *, int route);
const char *simplex_auto_route_admission_name(int admission);

/* Binding and capability queries do not invoke kernels or alter AUTO state.
 * Availability checks callable bindings only; model/resource eligibility is
 * separate. Missing optional bindings must not reserve work. */
int simplex_auto_routes_available(int route, const struct simplex_AutoKernels *,
	const struct simplex_AutoRoutes *);
void simplex_auto_routes_defaults(
	struct simplex_AutoKernels *kernels, struct simplex_AutoRoutes *routes);
int simplex_auto_routes_has_checkpoint(
	const struct simplex_AutoKernels *kernels);
/* Execute the controller's selected action. Outputs belong to this attempt;
 * the controller audits, charges and accepts them. A checkpoint may change
 * only through its kernel's export contract, never by route selection here. */
int simplex_auto_routes_invoke(int route, const struct simplex_Problem *problem,
	const struct lp_simplex_Options *options, const double *seed, double *x,
	struct lp_simplex_Result *result, int propagate_bounds,
	const struct simplex_AutoKernels *kernels,
	const struct simplex_AutoRoutes *routes,
	struct simplex_DualCheckpoint *checkpoint);

#endif
