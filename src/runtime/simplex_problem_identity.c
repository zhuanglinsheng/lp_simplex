/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_problem_identity.c
 * @brief Model identity snapshots and immutable-input validation.
 */

#include "simplex_problem_identity.h"
#include "simplex_problem.h"
#include <limits.h>
#include <string.h>

struct simplex_ProblemIdentity simplex_problem_identity(
	const struct simplex_Problem *p)
{
	struct simplex_ProblemIdentity id;
	memset(&id, 0, sizeof(id));
	if (!p)
		return id;
	id.input = p;
	id.context = p->execution.context;
	id.revision = p->revision;
	id.model = p->execution.model;
	id.rows = p->rows;
	id.columns = p->columns;
	id.nonzeros = p->matrix.nonzeros;
	id.values = p->matrix.value;
	id.starts = p->matrix.column_start;
	id.indices = p->matrix.row_index;
	id.row_starts = p->matrix.row_start;
	id.column_indices = p->matrix.column_index;
	id.row_values = p->matrix.row_value;
	id.rhs = p->rhs;
	id.cost = p->objective;
	id.bounds = p->bounds;
	id.row_type = p->row_type;
	return id;
}

int simplex_problem_identity_equal(const struct simplex_ProblemIdentity *a,
	const struct simplex_ProblemIdentity *b)
{
	return a && b && a->input == b->input && a->context == b->context &&
	       a->revision == b->revision && a->model == b->model &&
	       a->rows == b->rows && a->columns == b->columns &&
	       a->nonzeros == b->nonzeros && a->values == b->values &&
	       a->starts == b->starts && a->indices == b->indices &&
	       a->row_starts == b->row_starts &&
	       a->column_indices == b->column_indices &&
	       a->row_values == b->row_values && a->rhs == b->rhs &&
	       a->cost == b->cost && a->bounds == b->bounds &&
	       a->row_type == b->row_type;
}

int simplex_problem_identity_unchanged(const struct simplex_ProblemIdentity *id,
	const struct simplex_Problem *p)
{
	return id && p && id->input == p &&
	       id->context == p->execution.context &&
	       id->revision == p->revision && id->model == p->execution.model &&
	       id->rows == p->rows && id->columns == p->columns &&
	       id->nonzeros == p->matrix.nonzeros &&
	       p->matrix.rows == id->rows && p->matrix.columns == id->columns &&
	       id->values == p->matrix.value &&
	       id->starts == p->matrix.column_start &&
	       id->indices == p->matrix.row_index &&
	       id->row_starts == p->matrix.row_start &&
	       id->column_indices == p->matrix.column_index &&
	       id->row_values == p->matrix.row_value && id->rhs == p->rhs &&
	       id->cost == p->objective && id->bounds == p->bounds &&
	       id->row_type == p->row_type && (!id->context || id->model > 0);
}

int simplex_problem_identity_matches(const struct simplex_ProblemIdentity *id,
	const struct simplex_Problem *p)
{
	return simplex_problem_identity_unchanged(id, p) &&
	       id->revision != ULONG_MAX;
}
