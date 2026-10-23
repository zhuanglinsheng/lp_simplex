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

#include <lp_simplex/model.h>


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
	/* Borrowed observer; NULL outside an observed solve. */
	const struct simplex_ProgressSink *progress;
};


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

/* Absolute box/row KKT sign violation. sign=+1 means c+A^T y,
 * sign=-1 means c-A^T y. Bounds/rows within primal_tolerance are active. */
double simplex_problem_dual_infeasibility(const struct simplex_Problem *problem,
	const double *x, const double *row_dual, double primal_tolerance,
	int sign);

void simplex_problem_destroy(struct simplex_Problem *problem);

void simplex_problem_free(struct simplex_Problem *problem);

#endif
