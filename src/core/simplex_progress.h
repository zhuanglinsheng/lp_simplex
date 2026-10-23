/* Copyright (C) 2022 Zhuang Linsheng. License: LGPL 3.0. */
#ifndef LP_SIMPLEX_PROGRESS_INTERNAL_H
#define LP_SIMPLEX_PROGRESS_INTERNAL_H
#include "simplex_problem.h"
#include "simplex_presolve.h"
#include "simplex_scaling.h"
#include <lp_simplex/solve.h>
#include <time.h>

/* All events use c + A^T (multiplier_scale * y). Vectors are borrowed.
 * due() gates O(nnz) work; force is reserved for first feasible/final events. */
struct simplex_ProgressSink {
	void *context;
	int (*due)(void *context);
	void (*report)(void *context, const double *x, const double *y,
		double multiplier_scale, const char *stage, int iterations,
		int force);
};
struct simplex_ProgressTracker {
	struct simplex_ProgressSink sink;
	const struct lp_Model *model;
	const struct lp_simplex_Options *options;
	struct simplex_Problem original;
	double *storage, *point;
	double best_primal, best_dual, primal_error;
	clock_t begun, last;
	int have_primal;
};
struct simplex_ProgressBridge {
	struct simplex_ProgressSink sink;
	const struct simplex_ProgressSink *parent;
	const struct simplex_Presolve *presolve;
	const struct simplex_Scaling *scaling;
	double *storage, *x, *y;
	int rows, columns;
};
int simplex_progress_init(struct simplex_ProgressTracker *t,
	const struct lp_Model *model, const struct lp_simplex_Options *o);
void simplex_progress_destroy(struct simplex_ProgressTracker *t);
int simplex_progress_bridge(struct simplex_ProgressBridge *b,
	const struct simplex_ProgressSink *parent,
	const struct simplex_Presolve *presolve,
	const struct simplex_Scaling *scaling, int rows, int columns);
void simplex_progress_bridge_destroy(struct simplex_ProgressBridge *b);
void simplex_progress_emit(const struct simplex_Problem *p, const double *x,
	const double *y, double scale, const char *stage, int iterations,
	int force);
/* Original-cost box Lagrangian bound. Projects inequality multipliers to
 * their permitted signs. Returns -infinity if any infinite box endpoint
 * prevents a finite, conservatively rounded bound. */
double simplex_progress_dual_bound(
	const struct simplex_Problem *p, const double *y, double scale);
#endif
