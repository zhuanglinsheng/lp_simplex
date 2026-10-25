/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_problem.h
 * @brief Canonical sparse problem storage and primal-dual residual checks.
 */

#ifndef LP_SIMPLEX_PROBLEM_INTERNAL_H
#define LP_SIMPLEX_PROBLEM_INTERNAL_H

#include "simplex_csc.h"
#include "simplex_execution_scope.h"

#include <lp_simplex/model.h>

struct simplex_ProgressSink;
struct simplex_StageReport;

/* Canonical immutable input to the optimization kernel.  Public models,
 * presolve output and future decomposition nodes all enter through this
 * boundary; the simplex code never depends on dense-model ownership rules. */
struct simplex_Problem {
	int rows;
	int columns;
	struct simplex_CscMatrix matrix;
	double *objective;
	struct optm_VariableBound *bounds;
	double *rhs;
	unsigned char *row_type;
	void *vector_storage;
	int owns_vectors;
	/* Immutable while a checkpoint borrows this problem. Writers must mark
	 * in-place changes before attempting to reuse any existing checkpoint. */
	unsigned long revision;
	/* Borrowed observer; NULL outside an observed solve. */
	const struct simplex_ProgressSink *progress;
	/* Borrowed per-invocation diagnostics; never owns kernel state. */
	struct simplex_StageReport *report;
	/* Borrowed solve context/model identity, plus optional output journal and
	 * causal parent. Diagnostic storage never owns numerical decision state. */
	struct simplex_ExecutionScope execution;
};

/* Constructors initialize fresh outputs and reject unrepresentable owned
 * storage before allocation. They do not select an algorithm or impose its
 * combined variable-id limit. Borrowed inputs must be well formed/immutable. */
int simplex_problem_from_model(
	struct simplex_Problem *problem, const struct lp_Model *model);

struct simplex_Problem *simplex_problem_create_sparse(
	int rows, int columns, int nonzeros);

/* Allocate vectors and matrix values, borrowing the immutable CSC/CSR indices.
 * The pattern owner must outlive this problem and any views derived from it. */
struct simplex_Problem *simplex_problem_create_with_pattern(
	const struct simplex_CscMatrix *pattern);

/* Absolute primal violation on the unscaled problem; no workspace allocation.
 * Nonfinite coordinates/products/residuals return infinity. */
double simplex_problem_primal_infeasibility(
	const struct simplex_Problem *problem, const double *x);

/* Fresh feasibility quantities on this immutable input. Only violations
 * larger than tolerance contribute to sum/count; maximum remains absolute.
 * Failure leaves all diagnostics unavailable, never partially accumulated. */
struct simplex_PrimalSummary {
	double maximum, sum, count;
};

int simplex_problem_primal_summary(const struct simplex_Problem *problem,
	const double *x, double tolerance,
	struct simplex_PrimalSummary *summary);

/* Absolute box/row KKT sign violation. sign=+1 means c+A^T y,
 * sign=-1 means c-A^T y. Bounds/rows within primal_tolerance are active. */
double simplex_problem_dual_infeasibility(const struct simplex_Problem *problem,
	const double *x, const double *row_dual, double primal_tolerance,
	int sign);

/* Check a separating row combination against immutable rows and bounds.
 * Both orientations are tested; objective perturbations are irrelevant.
 * Failure to verify is a numerical failure, never proof of feasibility. */
int simplex_problem_verify_infeasibility(
	const struct simplex_Problem *problem, const double *row_multiplier);

/* Borrowed A z = rhs, z >= 0 view, using raw CSC coefficients. The same
 * conservative separation verifier applies; no normalized correlations or
 * NNLS stationarity test can substitute for this check. */
int simplex_problem_verify_nonnegative_infeasibility(
	const struct simplex_CscMatrix *matrix, const double *rhs,
	const double *row_multiplier);

/* Reconstruct and verify the basis row in exact rational arithmetic when
 * rounded multipliers cannot decide signs on unbounded basic coordinates. */
int simplex_problem_verify_basis_infeasibility(
	const struct simplex_Problem *problem, const int *basis,
	int leaving_row);

void simplex_problem_mark_changed(struct simplex_Problem *problem);

void simplex_problem_destroy(struct simplex_Problem *problem);

void simplex_problem_free(struct simplex_Problem *problem);

#endif
