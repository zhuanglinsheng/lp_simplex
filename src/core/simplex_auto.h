/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_auto.h
 * @brief Problem-independent orchestration; no global environment mutation.
 */

#ifndef LP_SIMPLEX_AUTO_INTERNAL_H
#define LP_SIMPLEX_AUTO_INTERNAL_H

#include "policy/search.h"
#include "simplex_problem.h"
#include "simplex_stage.h"

#include <lp_simplex/solve.h>

enum simplex_AutoRoute {
	SIMPLEX_AUTO_HYBRID = 0,
	SIMPLEX_AUTO_EXPLICIT,
	SIMPLEX_AUTO_CLASSIC,
	SIMPLEX_AUTO_STANDARD_IMPLICIT,
	SIMPLEX_AUTO_STANDARD_EXPLICIT,
	SIMPLEX_AUTO_SINGLETON,
	SIMPLEX_AUTO_TABLEAU,
	SIMPLEX_AUTO_ROUTE_COUNT
};

/* Optional feasibility search follows a dual yield without a feasible point.
 * Input carries this action's scope and private report, with the exact policy
 * grant in options. It cannot claim a terminal LP status, basis or bound. */
struct simplex_AutoInitializer {
	void *context;
	int (*run)(void *, const struct simplex_Problem *,
		const struct lp_simplex_Options *, double *, int *);
};

/* Injectable internal kernel boundary, also used for deterministic failure
 * tests. No callback/global state escapes a solve. */
struct simplex_AutoKernels {
	int (*dual)(const struct simplex_Problem *,
		const struct lp_simplex_Options *, double *, double *,
		struct lp_simplex_Result *, int);

	int (*hybrid)(const struct simplex_Problem *,
		const struct lp_simplex_Options *, const double *, double *,
		struct lp_simplex_Result *);
};

/* dual is required by AUTO; a missing binding stops before numerical work.
 * hybrid is optional and missing means no native BDA recovery capability. */

/* Optional recovery-route capability is independent of the two primary
 * kernel bindings. A missing solve callback admits no recovery routes. */
struct simplex_AutoRoutes {
	int (*solve)(int, const struct simplex_Problem *,
		const struct lp_simplex_Options *, const double *, double *,
		struct lp_simplex_Result *);
};

enum simplex_AutoPrimaryStage {
	SIMPLEX_AUTO_DUAL_EXPLORE = -1,
	SIMPLEX_AUTO_DUAL_FALLBACK = -2,
	SIMPLEX_AUTO_COVER = -3,
	SIMPLEX_AUTO_INITIALIZER = -4
};

enum simplex_AutoPurpose {
	SIMPLEX_AUTO_PRIMARY_SEARCH = 0,
	SIMPLEX_AUTO_SEEK_FEASIBLE,
	SIMPLEX_AUTO_OPTIMIZE_POINT,
	SIMPLEX_AUTO_REBUILD_BASIS,
	SIMPLEX_AUTO_REFORMULATE,
	SIMPLEX_AUTO_REPLACE_ORACLE,
	SIMPLEX_AUTO_RESTRICT_COLUMNS,
	SIMPLEX_AUTO_ALTERNATIVE_ALGORITHM,
	SIMPLEX_AUTO_REPLACE_DIRECTION
};

enum simplex_AutoValidation {
	SIMPLEX_AUTO_VALIDATION_UNKNOWN = 0,
	SIMPLEX_AUTO_POINT_VALID,
	SIMPLEX_AUTO_POINT_INVALID,
	SIMPLEX_AUTO_TERMINAL_VALID,
	SIMPLEX_AUTO_TERMINAL_INVALID
};

/* A bounded per-solve journal. Reports describe the stage's immutable input;
 * they are not interchangeable with original-model observations after presolve.
 * Each route runs at most once; retrying it requires a new solve/checkpoint. */
enum simplex_AutoPrerequisite {
	SIMPLEX_AUTO_UNCERTIFIED = 0,
	SIMPLEX_AUTO_AUDITED_POINT,
	SIMPLEX_AUTO_BASIS_HINT,
	SIMPLEX_AUTO_IMPLICIT_ORACLE_FAILURE
};

enum simplex_AutoFailureAction {
	SIMPLEX_AUTO_NO_RECOVERY = 0,
	SIMPLEX_AUTO_NEXT_ROUTE
};

/* References into two independent diagnostic journals. Neither is an input
 * identity, a saved numerical state, or evidence that an action is applicable. */
struct simplex_AutoOrigin {
	int record, event;
	struct simplex_ExecutionRef actor;
};

struct simplex_AutoAction {
	int stage, grant, prerequisite, purpose, reason;
	struct simplex_ExecutionRef actor, cause_actor, point_actor,
		basis_actor;
	/* Solve-wide causal event; zero when history is absent or truncated. */
	int execution_id;
	/* Applies only to retryable failure, with remaining shared budget. */
	int failure_action;
	/* A failed stage can cause a transition while a different stage supplies
	 * its point or basis. Indices refer only to earlier records in this solve. */
	int cause, point_source, basis_source;
	int cause_event, point_event, basis_event;
};

struct simplex_AutoStageRecord {
	struct simplex_AutoAction action;
	int work, reported_work, accepted_point, basis_accepted;
	int kernel_status, status, validation;
	/* Import acceptance and fresh exported output are independent. */
	int basis_exported, basis_producer_event;
	struct simplex_ExecutionRef basis_producer_actor;
	/* Scheduler audit is separate from the kernel's own diagnostics. */
	struct simplex_Observation primal_audit, objective_audit;
	struct simplex_StageReport report;
};

enum simplex_AutoExit {
	SIMPLEX_AUTO_EXIT_UNKNOWN = 0,
	SIMPLEX_AUTO_EXIT_TERMINAL,
	SIMPLEX_AUTO_EXIT_BUDGET,
	SIMPLEX_AUTO_EXIT_ROUTES,
	SIMPLEX_AUTO_EXIT_INPUT_CHANGED,
	SIMPLEX_AUTO_EXIT_RESOURCE,
	SIMPLEX_AUTO_EXIT_AUDIT,
	SIMPLEX_AUTO_EXIT_CAPABILITY
};

/* Output only. Global termination and the returned point have independent
 * provenance; neither reference replaces a model identity or certificate. */
struct simplex_AutoOutcome {
	int exit, state, status, work, stop, execution_id;
	int point_available, restored;
	struct simplex_AutoOrigin cause, point;
};

/* Admission observations are diagnostic output, never policy/runtime state. */
enum simplex_AutoBinding {
	SIMPLEX_AUTO_BINDING_UNCHECKED = 0,
	SIMPLEX_AUTO_BINDING_MISSING,
	SIMPLEX_AUTO_BINDING_READY
};

enum simplex_AutoCondition {
	SIMPLEX_AUTO_CONDITION_UNCHECKED = 0,
	SIMPLEX_AUTO_CONDITION_READY,
	SIMPLEX_AUTO_CONDITION_ORACLE
};

struct simplex_AutoAdmission {
	int binding, model, condition;
	struct simplex_AutoOrigin prerequisite;
};

/* Output-only journal and summary. AUTO does not recover runtime state from it. */
struct simplex_AutoHistory {
	struct simplex_AutoStageRecord stages[SIMPLEX_AUTO_ROUTE_COUNT + 4];
	int count, have_point, point_source, bound_source;
	double best_objective, best_bound;
	struct simplex_AutoOutcome outcome;
	struct simplex_AutoAdmission admission[SIMPLEX_AUTO_ROUTE_COUNT];
};

/* Borrowed callback tables, copied at entry. An explicit object requests
 * injected capabilities: NULL tables/functions stay absent, never defaults.
 * cover_enabled separately admits the production restricted-column adapter;
 * injection tests normally leave it zero. No runtime or history lives here. */
struct simplex_AutoBindings {
	const struct simplex_AutoKernels *kernels;
	const struct simplex_AutoRoutes *routes;
	int cover_enabled;
};

/* Sole AUTO controller entry. NULL bindings select production capabilities,
 * including policy-controlled cover. An optional initializer supplies points
 * only. NULL history uses a private journal; it does not disable accounting,
 * audits or candidate retention. All borrowed objects live through the call. */
int simplex_auto_solve(const struct simplex_Problem *problem,
	const struct lp_simplex_Options *options, double *x,
	struct lp_simplex_Result *result, int propagate_bounds,
	const struct simplex_AutoBindings *bindings,
	const struct simplex_AutoInitializer *initializer,
	struct simplex_AutoHistory *history);

#endif
