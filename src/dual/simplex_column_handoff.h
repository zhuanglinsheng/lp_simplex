/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_column_handoff.h
 * @brief Validated restricted-column model and candidate transfers.
 */

#ifndef LP_SIMPLEX_COLUMN_HANDOFF_INTERNAL_H
#define LP_SIMPLEX_COLUMN_HANDOFF_INTERNAL_H
#include "simplex_dual.h"
#include "simplex_execution_scope.h"
#include <stddef.h>

enum simplex_ColumnRejection {
	SIMPLEX_COLUMN_OK = 0,
	SIMPLEX_COLUMN_MODEL,
	SIMPLEX_COLUMN_MAP,
	SIMPLEX_COLUMN_RELATION,
	SIMPLEX_COLUMN_LAYOUT,
	SIMPLEX_COLUMN_POINT,
	SIMPLEX_COLUMN_MULTIPLIER,
	SIMPLEX_COLUMN_STORAGE
};

enum simplex_ColumnHandoffKind {
	SIMPLEX_COLUMN_BIND = 0,
	SIMPLEX_COLUMN_BASIS_LIFT,
	SIMPLEX_COLUMN_BASIS_RESTRICT,
	SIMPLEX_COLUMN_POINT_LIFT,
	SIMPLEX_COLUMN_ROWS_LIFT
};

/* Exact column restriction, unchanged rows/costs/bounds, omitted variables
 * fixed at zero. Maps and models are borrowed and immutable until retirement. */
struct simplex_ColumnProjection {
	struct simplex_ProblemIdentity parent, child;
	const int *map, *reverse;
	int valid;
};

/* Basis ids/statuses in parent coordinates; no cached primal/dual/factor.
 * It survives destruction of a restricted model, never of its parent. */
struct simplex_ColumnBasis {
	struct simplex_ProblemIdentity parent;
	int valid, *index;
	/* Output-only provenance; absent/truncated history does not change validity. */
	struct simplex_ExecutionRef producer;
	unsigned char *status;
};

void simplex_column_basis_init(struct simplex_ColumnBasis *basis,
	const struct simplex_Problem *parent, int *index,
	unsigned char *status);
int simplex_column_projection_bind(struct simplex_ColumnProjection *projection,
	const struct simplex_Problem *parent,
	const struct simplex_Problem *child, const int *map, const int *reverse,
	int *row_cursor, int cursor_capacity);
/* Scratch is ceil((parent.columns + parent.rows)/8) bytes, caller-owned.
 * Rejection never changes either seed's arrays or its validity/binding. */
int simplex_column_basis_lift(const struct simplex_ColumnProjection *projection,
	const struct simplex_DualBasisSeed *seed,
	struct simplex_ColumnBasis *basis, unsigned char *scratch,
	size_t bytes);
int simplex_column_basis_restrict(
	const struct simplex_ColumnProjection *projection,
	const struct simplex_ColumnBasis *basis,
	struct simplex_DualBasisSeed *seed, int prefer_dual,
	unsigned char *scratch, size_t bytes);
/* Output is caller-owned staging storage with parent.columns entries. Only a
 * zero return authorizes committing it; input/output arrays must not alias. */
int simplex_column_point_lift(const struct simplex_ColumnProjection *projection,
	const double *point, double *staging, double tolerance);
/* Same row coordinates/sign convention. This only authorizes a multiplier
 * candidate for full-parent pricing/audit, never a parent dual certificate. */
int simplex_column_rows_lift(const struct simplex_ColumnProjection *projection,
	const double *multipliers);
#endif
