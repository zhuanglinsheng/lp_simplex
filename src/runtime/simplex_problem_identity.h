/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_problem_identity.h
 * @brief Model identity snapshots and immutable-input validation.
 */

#ifndef LP_SIMPLEX_PROBLEM_IDENTITY_INTERNAL_H
#define LP_SIMPLEX_PROBLEM_IDENTITY_INTERNAL_H
struct simplex_Problem;
struct simplex_SolveContext;

/* Borrowed immutable storage. No event ids, factors or residual estimates. */
struct simplex_ProblemIdentity {
	const struct simplex_Problem *input;
	const struct simplex_SolveContext *context;
	unsigned long revision;
	int model, rows, columns, nonzeros;
	const void *values, *starts, *indices, *row_starts, *column_indices,
		*row_values;
	const void *rhs, *cost, *bounds, *row_type;
};
struct simplex_ProblemIdentity simplex_problem_identity(
	const struct simplex_Problem *p);
/* Compare borrowed snapshots without dereferencing their source objects.
 * A caller depth, diagnostic parent or regenerated factor is not identity. */
int simplex_problem_identity_equal(const struct simplex_ProblemIdentity *a,
	const struct simplex_ProblemIdentity *b);
/* Current-invocation audits may observe saturated revisions; reusable state
 * must use matches(), which never accepts a saturated revision. */
int simplex_problem_identity_unchanged(const struct simplex_ProblemIdentity *id,
	const struct simplex_Problem *p);
int simplex_problem_identity_matches(const struct simplex_ProblemIdentity *id,
	const struct simplex_Problem *p);
#endif
