/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_dual.h
 * @brief Dual revised-simplex iteration and numerical recovery.
 */

#ifndef LP_SIMPLEX_DUAL_INTERNAL_H
#define LP_SIMPLEX_DUAL_INTERNAL_H

#include "simplex_problem.h"
#include "simplex_retry.h"
#include "simplex_problem_identity.h"

#include <lp_simplex/model.h>
#include <lp_simplex/solve.h>

int simplex_dual_solve(const struct lp_Model *model,
	const struct lp_simplex_Options *options, double *x, double *row_dual,
	struct lp_simplex_Result *result, int propagate_bounds);

int simplex_dual_solve_problem(const struct simplex_Problem *problem,
	const struct lp_simplex_Options *options, double *x, double *row_dual,
	struct lp_simplex_Result *result, int propagate_bounds);

/* Caller-owned basis seed for related internal LPs. index has rows entries;
 * status has columns+rows entries using DUAL_STATUS_*. Logical ids are
 * columns+row. A successful solve exports a seed. A resource-refused checked
 * import may retain the input hint, without exporting it; other failures
 * invalidate it.
 * Import rebuilds and certifies the factor/primal state, never trusts a cached
 * inverse or reduced costs. prefer_dual permits moving nonbasics to dual-feasible
 * bounds, retaining a dual-feasible seed even when primal repair is needed.
 * Seeds satisfying neither feasibility condition use a cold start.
 * Bind storage with simplex_dual_seed_init. Reuse on another model requires
 * checked column projection and explicit rebinding; matching dimensions alone
 * never authorize import or cold export into foreign storage. */
struct simplex_DualBasisSeed {
	struct simplex_ProblemIdentity identity;
	int rows, columns;
	int valid, accepted;
	/* valid can describe a retained older hint; exported is this call's
	 * fresh output. Diagnostic source never authorizes import. */
	int exported;
	struct simplex_ExecutionScope export_source;
	unsigned long export_serial;
	int import_attempted, rejection;
	int prefer_dual;
	int *index;
	unsigned char *status;
};

void simplex_dual_seed_init(struct simplex_DualBasisSeed *seed,
	const struct simplex_Problem *problem, int *index,
	unsigned char *status);

int simplex_dual_solve_problem_seeded(const struct simplex_Problem *problem,
	const struct lp_simplex_Options *options, double *x, double *row_dual,
	struct lp_simplex_Result *result, int propagate_bounds,
	struct simplex_DualBasisSeed *seed);

/* A checkpoint belongs to one live, immutable input problem. Unlike the
 * related-LP seed service, it never authorizes implicit index remapping.
 * No factors, bounds, costs or residual estimates are stored here. */
struct simplex_DualCheckpoint {
	struct simplex_DualBasisSeed basis;
};

void simplex_dual_checkpoint_init(struct simplex_DualCheckpoint *checkpoint,
	const struct simplex_Problem *problem, int *index,
	unsigned char *status);
int simplex_dual_checkpoint_rejection(
	const struct simplex_DualCheckpoint *checkpoint,
	const struct simplex_Problem *problem);
int simplex_dual_checkpoint_matches(
	const struct simplex_DualCheckpoint *checkpoint,
	const struct simplex_Problem *problem);
/* Only a numerical import failure or an accepted seed supports retrying that
 * seed with a different factorization backend. Unknown is not permission. */
int simplex_dual_seed_retry_allowed(const struct simplex_DualBasisSeed *seed);
/* Reset invocation diagnostics before any allocation; preserve a valid input
 * hint and its producer. A kernel export replaces that producer explicitly. */
void simplex_dual_seed_begin(struct simplex_DualBasisSeed *seed);
void simplex_dual_seed_export(struct simplex_DualBasisSeed *seed,
	const struct simplex_ExecutionScope *producer);
int simplex_dual_seed_source_event(const struct simplex_DualBasisSeed *seed,
	const struct simplex_ExecutionScope *observer);

/* May export a basis hint on yield/budget exhaustion. Matching is checked on
 * immutable input before internal scaling. Every import refactors/reclassifies.
 * A stale checkpoint is discarded locally and the solve proceeds cold, without
 * exporting into another model's storage. Reinitialize explicitly to rebind.
 * The source and caller-owned arrays must outlive this checkpoint. */
int simplex_dual_solve_problem_checkpoint(const struct simplex_Problem *problem,
	const struct lp_simplex_Options *options, double *x, double *row_dual,
	struct lp_simplex_Result *result, int propagate_bounds,
	struct simplex_DualCheckpoint *checkpoint);

/* Shared start evidence. Basis requests compare the exact input arrays and
 * import policy, never a basis/source event id. Failed requests retain copies
 * owned by the solve context. Callers must not alias seed arrays with outputs. */
enum simplex_DualRetryRepresentation {
	SIMPLEX_DUAL_RETRY_UNSCALED = 0,
	SIMPLEX_DUAL_RETRY_IDENTITY,
	SIMPLEX_DUAL_RETRY_EQUILIBRATED
};

struct simplex_DualRetry {
	struct simplex_SolveContext *context;
	const struct simplex_Execution *history;
	struct simplex_RetryKey key;
};

void simplex_dual_retry_init(struct simplex_DualRetry *retry,
	const struct simplex_Problem *input,
	const struct lp_simplex_Options *options, int propagate,
	int feasibility_only, int seeded);
void simplex_dual_retry_bind_basis(struct simplex_DualRetry *retry,
	const struct simplex_DualBasisSeed *seed);
/* Returns the matching record index plus one when skipped, otherwise zero.
 * Both primary and optional cold backend calls obey this admission rule. */
int simplex_dual_retry_blocked(const struct simplex_DualRetry *retry,
	const struct simplex_Problem *kernel, int backend, int grant);
void simplex_dual_retry_record(const struct simplex_DualRetry *retry,
	int backend, int event, int state,
	const struct lp_simplex_Result *result, int grant, int basis_retry);

/* Internal Phase I service: success means only a verified primal feasible
 * point; result.status remains CondUnsatisfied, never an optimality claim. */
int simplex_dual_find_feasible(const struct simplex_Problem *problem,
	const struct lp_simplex_Options *options, double *x,
	struct lp_simplex_Result *result);

/* A solve-local diagnostic view keeps the stable input as the numerical
 * retry identity. Scope may change history/parent, never context/model. */
int simplex_dual_find_feasible_reported(const struct simplex_Problem *input,
	const struct lp_simplex_Options *options, double *x,
	struct lp_simplex_Result *result,
	const struct simplex_ExecutionScope *execution,
	struct simplex_StageReport *report);

#endif
