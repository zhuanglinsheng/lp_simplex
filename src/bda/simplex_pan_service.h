/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_pan_service.h
 * @brief BDA subprocess roles, model bindings and bounded work handoffs.
 */

#ifndef LP_SIMPLEX_PAN_SERVICE_H
#define LP_SIMPLEX_PAN_SERVICE_H

#include "simplex_pan_service_types.h"
#include "simplex_pan_standard.h"
#include "simplex_execution.h"
#include "simplex_problem_identity.h"

struct simplex_PanBasis;
struct simplex_PanBounded;
struct simplex_PanConfig;
struct simplex_PanFaceRecovery;

/* Solve-local protection state; no eligibility is reconstructed from journals. */
struct simplex_PanServiceContext {
	const struct simplex_Problem *input;
	struct simplex_ExecutionScope scope;
	struct simplex_ProblemIdentity identity;
	unsigned long representations;
	int parent, failed;
	struct simplex_ExecutionRef parent_actor;
};

struct simplex_PanServiceFrame {
	struct simplex_PanServiceOutcome outcome;
	struct simplex_PanStandard identity;
	struct simplex_ProblemIdentity input;
	int saved_parent;
	struct simplex_ExecutionRef saved_parent_actor, actor;
	struct simplex_ExecutionRef trigger_actor;
	const struct simplex_PanBounded *bounded;
	const struct simplex_Problem *bounded_input;
	const double *lower, *upper;
	struct simplex_PanServiceContext *context;
	const struct simplex_PanBasis *basis;
	struct simplex_StageReport *sink;
};

int simplex_pan_service_context_valid(
	const struct simplex_PanServiceContext *context);
void simplex_pan_service_bind(struct simplex_PanStandard *a,
	struct simplex_PanServiceContext *context);
void simplex_pan_service_derive(struct simplex_PanStandard *a,
	const struct simplex_PanStandard *source);
void simplex_pan_service_context_init(struct simplex_PanServiceContext *context,
	const struct simplex_Problem *input);
/* Grants are local leases against a shared cumulative counter. The frame
 * captures immutable input identity; finish conservatively charges a corrupt
 * report to the lease end and clears all output permissions. */
int simplex_pan_service_begin(struct simplex_PanServiceFrame *frame,
	const struct simplex_PanStandard *a, int service, int begin, int limit,
	int requested);
int simplex_pan_service_finish(struct simplex_PanServiceFrame *frame,
	int *iterations, int code, int status, int stop, int output,
	struct simplex_PanServiceOutcome *outcome);

int simplex_pan_service_nnls(struct simplex_PanBasis *basis, double *point,
	double primal_tolerance, double dual_tolerance, int limit,
	int *iterations, int *status,
	struct simplex_PanServiceOutcome *outcome);
/* Zero mathematical-work lease; raw face adjustments are not usable-point
 * permissions. Source is the actual multiplier search, never the last event. */
int simplex_pan_service_face(struct simplex_PanFaceRecovery *owner,
	const struct simplex_PanBounded *b, double *point,
	const double *multiplier, const struct lp_simplex_Options *options,
	int *iterations, double *saved, double *residual, double *work,
	int opportunity, int steps, int low_progress,
	const struct simplex_PanServiceOutcome *multiplier_search,
	struct simplex_PanServiceOutcome *outcome);
int simplex_pan_service_dual_feasible(struct simplex_PanServiceContext *context,
	const struct simplex_Problem *input,
	const struct lp_simplex_Options *options, double *point,
	int *iterations, struct lp_simplex_Result *result,
	struct simplex_PanServiceOutcome *outcome);
int simplex_pan_service_pdhg(const struct simplex_PanBounded *b, double *point,
	double *multiplier, const struct lp_simplex_Options *options,
	int *iterations, int budget, const struct simplex_PanConfig *config,
	struct simplex_PanServiceOutcome *outcome);
/* Optional producer evidence creates a completed observation cause. Neither
 * its absence nor loss of its journal reference changes repair eligibility. */
int simplex_pan_service_restore(const struct simplex_PanBounded *b,
	double *point, const struct lp_simplex_Options *options,
	int *iterations, long double ceiling,
	const struct simplex_PanConfig *config,
	const struct simplex_PanServiceTrigger *trigger,
	struct simplex_PanServiceOutcome *outcome);
/* Standard-form correction has the same role/lease protocol, with a distinct
 * equality-point output. A no-op returns 0 without new output authority. */
int simplex_pan_service_restore_standard(const struct simplex_PanStandard *a,
	double *point, const struct lp_simplex_Options *options,
	long double ceiling, int *iterations,
	const struct simplex_PanConfig *config,
	const struct simplex_PanServiceTrigger *trigger,
	struct simplex_PanServiceOutcome *outcome);
int simplex_pan_service_first_order(const struct simplex_PanStandard *a,
	const double *point, double *multiplier,
	const struct lp_simplex_Options *options, int *iterations,
	double *violation, double *direction, const unsigned char *free_mask,
	int warm, const struct simplex_PanConfig *config,
	struct simplex_PanServiceOutcome *outcome);
int simplex_pan_service_predict(const struct simplex_PanStandard *a,
	const double *point, double *multiplier,
	const struct lp_simplex_Options *options, int *iterations,
	double *violation, double *direction, const unsigned char *free_mask,
	int budget, int (*accept)(void *, const double *), void *accept_context,
	const struct simplex_PanConfig *config,
	struct simplex_PanServiceOutcome *outcome);
int simplex_pan_service_cone(struct simplex_PanBasis *basis, int implicit,
	double *direction, const unsigned char *free_mask, double *residual,
	double tolerance, int limit, int *iterations, int *status,
	int (*accept)(void *, const double *), void *accept_context,
	struct simplex_PanServiceOutcome *outcome);
#endif
