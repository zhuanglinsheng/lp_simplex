/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_progress.h
 * @brief Audited input-model progress and producer provenance.
 */

#ifndef LP_SIMPLEX_PROGRESS_INTERNAL_H
#define LP_SIMPLEX_PROGRESS_INTERNAL_H
#include "simplex_problem.h"
#include "simplex_presolve.h"
#include "simplex_scaling.h"
#include "simplex_problem_identity.h"
#include <lp_simplex/solve.h>
#include <time.h>

enum simplex_ProgressTransform {
	SIMPLEX_PROGRESS_SCALING = 1,
	SIMPLEX_PROGRESS_PRESOLVE = 2
};

/* Synchronous borrowed metadata. producer is the emitting algorithm's active
 * scope; coordinates follows reversible bridges. Neither is solver state. */
struct simplex_ProgressSource {
	struct simplex_ExecutionScope producer, coordinates;
	struct simplex_ProblemIdentity input;
	int rows, columns, work;
	unsigned int transforms;
};

/* All events use c + A^T (multiplier_scale * y). Vectors are borrowed.
 * due() gates O(nnz) work; force is reserved for first feasible/final events. */
struct simplex_ProgressSink {
	void *context;
	int (*due)(void *context);
	void (*report)(void *context, const double *x, const double *y,
		double multiplier_scale, const char *stage, int iterations,
		int force, const struct simplex_ProgressSource *source);
};

struct simplex_ProgressTracker {
	struct simplex_ProgressSink sink;
	const struct lp_Model *model;
	const struct lp_simplex_Options *options;
	struct simplex_Problem original;
	double *storage, *point;
	double best_primal, best_dual, primal_error;
	clock_t begun, last;
	int have_primal, point_pending;
	unsigned int polls;
	struct simplex_ExecutionScope execution;
	struct simplex_ExecutionRef parent_actor;
	int point_source;
	unsigned long point_serial;
};

struct simplex_ProgressBridge {
	struct simplex_ProgressSink sink;
	const struct simplex_ProgressSink *parent;
	const struct simplex_Presolve *presolve;
	const struct simplex_Scaling *scaling;
	double *storage, *x, *y;
	int rows, columns;
	struct simplex_ProblemIdentity input;
	struct simplex_ProblemIdentity output;
	struct simplex_ExecutionScope target;
};

enum simplex_ProgressBridgeResult {
	SIMPLEX_PROGRESS_BRIDGE_MEMORY = -1,
	SIMPLEX_PROGRESS_BRIDGE_READY = 0,
	SIMPLEX_PROGRESS_BRIDGE_INVALID = 1
};

int simplex_progress_init(struct simplex_ProgressTracker *t,
	const struct lp_Model *model, const struct lp_simplex_Options *o);
int simplex_progress_init_reported(struct simplex_ProgressTracker *t,
	const struct lp_Model *model, const struct lp_simplex_Options *o,
	const struct simplex_ExecutionScope *execution);
/* A final assembly observation is explicitly attributed to the root scope,
 * rather than guessing the last numerical action. */
void simplex_progress_finish(
	struct simplex_ProgressTracker *t, const double *x, int work);
int simplex_progress_point_source(const struct simplex_ProgressTracker *t,
	const struct simplex_ExecutionScope *execution);
/* Diagnostic only: caller has already committed the archived point/result. */
int simplex_progress_record_restore(const struct simplex_ProgressTracker *t,
	const struct simplex_ExecutionScope *execution,
	const struct lp_simplex_Result *result);
void simplex_progress_destroy(struct simplex_ProgressTracker *t);
int simplex_progress_bridge(struct simplex_ProgressBridge *b,
	const struct simplex_ProgressSink *parent,
	const struct simplex_Presolve *presolve,
	const struct simplex_Scaling *scaling, int rows, int columns,
	const struct simplex_Problem *input,
	const struct simplex_ExecutionScope *target,
	const struct simplex_Problem *output);
void simplex_progress_bridge_destroy(struct simplex_ProgressBridge *b);
void simplex_progress_emit(const struct simplex_Problem *p, const double *x,
	const double *y, double scale, const char *stage, int iterations,
	int force);
/* Preserve an already-audited useful observation without forcing a display
 * notification. Sampling and UI refresh are independent policies. */
void simplex_progress_observe(const struct simplex_Problem *problem,
	const double *x, const double *y, double scale, const char *stage,
	int work);
/* Original-cost box Lagrangian bound. Projects inequality multipliers to
 * their permitted signs. Returns -infinity if any infinite box endpoint
 * prevents a finite, conservatively rounded bound. */
double simplex_progress_dual_bound(
	const struct simplex_Problem *p, const double *y, double scale);
#endif
