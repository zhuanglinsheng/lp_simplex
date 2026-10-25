/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_column_handoff.c
 * @brief Validated restricted-column model and candidate transfers.
 */

#include "simplex_stage.h"
#include "simplex_execution.h"
#include "simplex_column_handoff.h"
#include "simplex_dual_internal.h"
#include <limits.h>
#include <math.h>
#include <string.h>

static int column_trace_link(const struct simplex_ColumnProjection *v, int kind,
	int rejection, int cause, int *out_event)
{
	const struct simplex_ProblemIdentity *to =
		kind == SIMPLEX_COLUMN_BIND ||
				kind == SIMPLEX_COLUMN_BASIS_RESTRICT
			? &v->child
			: &v->parent;
	struct simplex_ExecutionScope scope = to->input->execution;
	int event;
	/* Outputs of a master belong under that invocation even when lifted to
	 * parent coordinates. Diagnostic references do not authorize the handoff. */
	if (kind != SIMPLEX_COLUMN_BIND)
		simplex_execution_inherit_parent(
			&scope, &v->child.input->execution);
	event = simplex_execution_begin(&scope,
		SIMPLEX_EXECUTION_COLUMN_HANDOFF, kind, SIMPLEX_REASON_NONE,
		cause);
	if (event) {
		struct simplex_ExecutionEvent *e =
			scope.history->events + event - 1;
		e->source_model =
			kind == SIMPLEX_COLUMN_BIND ||
					kind == SIMPLEX_COLUMN_BASIS_RESTRICT
				? v->parent.model
				: v->child.model;
		e->handoff_rejection = rejection;
		e->model_rows = to->rows;
		e->model_columns = to->columns;
		e->complete = 1;
		e->stop = rejection ? SIMPLEX_STOP_NUMERICAL
				    : SIMPLEX_STOP_HANDOFF;
		e->failure = rejection ? SIMPLEX_FAILURE_STAGE_CONTRACT
				       : SIMPLEX_FAILURE_UNKNOWN;
	}
	if (out_event)
		*out_event = event;
	return rejection;
}

static int column_trace(
	const struct simplex_ColumnProjection *v, int kind, int rejection)
{
	return column_trace_link(v, kind, rejection, 0, NULL);
}

static int column_maps(const struct simplex_ColumnProjection *v)
{
	int j, k;
	for (j = 0; j < v->child.columns; j++) {
		k = v->map[j];
		if (k < 0 || k >= v->parent.columns || v->reverse[k] != j)
			return 0;
	}
	for (j = 0; j < v->parent.columns; j++) {
		k = v->reverse[j];
		if (k < -1 || k >= v->child.columns ||
			(k >= 0 && v->map[k] != j))
			return 0;
	}
	return 1;
}

static int column_current(const struct simplex_ColumnProjection *v)
{
	return v && v->valid &&
	       simplex_problem_identity_matches(&v->parent, v->parent.input) &&
	       simplex_problem_identity_matches(&v->child, v->child.input);
}

static int column_bounds_equal(
	const struct optm_VariableBound *a, const struct optm_VariableBound *b)
{
	return a->b_type == b->b_type && a->v_type == b->v_type &&
	       ((a->b_type != optm_BOUND_T_LO &&
			a->b_type != optm_BOUND_T_BS) ||
		       a->lb == b->lb) &&
	       ((a->b_type != optm_BOUND_T_UP &&
			a->b_type != optm_BOUND_T_BS) ||
		       a->ub == b->ub);
}

static int column_zero_allowed(const struct optm_VariableBound *b)
{
	return b->b_type >= optm_BOUND_T_FR && b->b_type <= optm_BOUND_T_BS &&
	       ((b->b_type != optm_BOUND_T_LO &&
			b->b_type != optm_BOUND_T_BS) ||
		       b->lb <= 0.) &&
	       ((b->b_type != optm_BOUND_T_UP &&
			b->b_type != optm_BOUND_T_BS) ||
		       b->ub >= 0.);
}

static unsigned char column_endpoint(const struct optm_VariableBound *b)
{
	if (b->b_type == optm_BOUND_T_BS && b->lb == b->ub)
		return DUAL_STATUS_FIXED;
	if ((b->b_type == optm_BOUND_T_LO || b->b_type == optm_BOUND_T_BS) &&
		isfinite(b->lb))
		return DUAL_STATUS_LOWER;
	if ((b->b_type == optm_BOUND_T_UP || b->b_type == optm_BOUND_T_BS) &&
		isfinite(b->ub))
		return DUAL_STATUS_UPPER;
	return DUAL_STATUS_FREE;
}

static int column_csc(const struct simplex_Problem *p)
{
	int j;
	if (p->rows < 0 || p->columns < 0 || p->columns > INT_MAX - p->rows ||
		p->matrix.rows != p->rows || p->matrix.columns != p->columns ||
		p->matrix.nonzeros < 0 || !p->matrix.column_start ||
		!p->matrix.row_start ||
		(p->rows && (!p->rhs || !p->row_type)) ||
		(p->columns && (!p->objective || !p->bounds)) ||
		(p->matrix.nonzeros &&
			(!p->matrix.value || !p->matrix.row_index ||
				!p->matrix.row_value ||
				!p->matrix.column_index)) ||
		p->matrix.column_start[0] != 0 ||
		p->matrix.column_start[p->columns] != p->matrix.nonzeros)
		return 0;
	for (j = 0; j < p->columns; j++)
		if (p->matrix.column_start[j] > p->matrix.column_start[j + 1])
			return 0;
	return 1;
}

void simplex_column_basis_init(struct simplex_ColumnBasis *basis,
	const struct simplex_Problem *parent, int *index, unsigned char *status)
{
	memset(basis, 0, sizeof(*basis));
	basis->parent = simplex_problem_identity(parent);
	basis->index = index;
	basis->status = status;
}

int simplex_column_projection_bind(struct simplex_ColumnProjection *projection,
	const struct simplex_Problem *parent,
	const struct simplex_Problem *child, const int *map, const int *reverse,
	int *cursor, int capacity)
{
	struct simplex_ColumnProjection v;
	int i, j, k, s, t, rejection = SIMPLEX_COLUMN_RELATION;
	memset(&v, 0, sizeof(v));
	if (!parent || !child)
		return SIMPLEX_COLUMN_MODEL;
	v.parent = simplex_problem_identity(parent);
	v.child = simplex_problem_identity(child);
	v.map = map;
	v.reverse = reverse;
	if (parent->execution.context != child->execution.context ||
		!simplex_problem_identity_matches(&v.parent, parent) ||
		!simplex_problem_identity_matches(&v.child, child))
		rejection = SIMPLEX_COLUMN_MODEL;
	else if (!column_csc(parent) || !column_csc(child) ||
		 parent->rows != child->rows ||
		 child->columns > parent->columns)
		rejection = SIMPLEX_COLUMN_RELATION;
	else if ((!map && child->columns) || (!reverse && parent->columns) ||
		 !column_maps(&v))
		rejection = SIMPLEX_COLUMN_MAP;
	else if (capacity < child->rows || (!cursor && child->rows))
		rejection = SIMPLEX_COLUMN_STORAGE;
	else {
		for (i = 0; i < child->rows; i++) {
			if (!isfinite(child->rhs[i]) ||
				child->rhs[i] != parent->rhs[i] ||
				child->row_type[i] != parent->row_type[i] ||
				child->matrix.row_start[i] >
					child->matrix.row_start[i + 1])
				goto done;
			cursor[i] = child->matrix.row_start[i];
		}
		if (child->matrix.row_start[0] != 0 ||
			child->matrix.row_start[child->rows] !=
				child->matrix.nonzeros)
			goto done;
		for (j = 0; j < parent->columns; j++)
			if (reverse[j] < 0 &&
				!column_zero_allowed(parent->bounds + j))
				goto done;
		for (j = 0; j < child->columns; j++) {
			i = map[j];
			s = parent->matrix.column_start[i];
			if (!isfinite(child->objective[j]) ||
				child->objective[j] != parent->objective[i] ||
				!column_bounds_equal(child->bounds + j,
					parent->bounds + i) ||
				child->matrix.column_start[j + 1] -
						child->matrix.column_start[j] !=
					parent->matrix.column_start[i + 1] - s)
				goto done;
			for (k = child->matrix.column_start[j];
				k < child->matrix.column_start[j + 1];
				k++, s++) {
				i = child->matrix.row_index[k];
				if (i < 0 || i >= child->rows ||
					i != parent->matrix.row_index[s] ||
					!isfinite(child->matrix.value[k]) ||
					child->matrix.value[k] !=
						parent->matrix.value[s])
					goto done;
				t = cursor[i];
				if (t >= child->matrix.row_start[i + 1])
					goto done;
				cursor[i]++;
				if (child->matrix.column_index[t] != j ||
					child->matrix.row_value[t] !=
						child->matrix.value[k])
					goto done;
			}
		}
		for (i = 0; i < child->rows; i++)
			if (cursor[i] != child->matrix.row_start[i + 1])
				goto done;
		rejection = SIMPLEX_COLUMN_OK;
		v.valid = 1;
	}
done:
	if (!rejection)
		*projection = v;
	return column_trace(&v, SIMPLEX_COLUMN_BIND, rejection);
}

static int column_layout(int rows, int columns, const int *index,
	const unsigned char *status, unsigned char *seen, size_t bytes)
{
	int i, j, n;
	if (rows < 0 || columns < 0 || columns > INT_MAX - rows || !index ||
		!status)
		return SIMPLEX_COLUMN_LAYOUT;
	n = rows + columns;
	if (!seen || bytes < ((size_t)n + 7) / 8)
		return SIMPLEX_COLUMN_STORAGE;
	memset(seen, 0, ((size_t)n + 7) / 8);
	for (i = 0; i < rows; i++) {
		j = index[i];
		if (j < 0 || j >= n || (seen[j / 8] & (1U << (j % 8))) ||
			status[j] != DUAL_STATUS_BASIC)
			return SIMPLEX_COLUMN_LAYOUT;
		seen[j / 8] |= (unsigned char)(1U << (j % 8));
	}
	for (j = 0; j < n; j++)
		if (status[j] > DUAL_STATUS_FREE ||
			(status[j] == DUAL_STATUS_BASIC) !=
				!!(seen[j / 8] & (1U << (j % 8))))
			return SIMPLEX_COLUMN_LAYOUT;
	return SIMPLEX_COLUMN_OK;
}

int simplex_column_basis_lift(const struct simplex_ColumnProjection *v,
	const struct simplex_DualBasisSeed *seed,
	struct simplex_ColumnBasis *basis, unsigned char *scratch, size_t bytes)
{
	int i, j, rejection = SIMPLEX_COLUMN_MODEL;
	if (!column_current(v))
		return rejection;
	if (!column_maps(v))
		rejection = SIMPLEX_COLUMN_MAP;
	else if (!basis || !seed || !seed->valid ||
		 !simplex_problem_identity_matches(
			 &seed->identity, v->child.input) ||
		 !simplex_problem_identity_matches(
			 &basis->parent, v->parent.input))
		rejection = SIMPLEX_COLUMN_MODEL;
	else if (seed->rows != v->child.rows ||
		 seed->columns != v->child.columns || !basis->index ||
		 !basis->status)
		rejection = SIMPLEX_COLUMN_LAYOUT;
	else if (!(rejection = column_layout(seed->rows, seed->columns,
			   seed->index, seed->status, scratch, bytes))) {
		for (j = 0; j < v->parent.columns; j++)
			basis->status[j] =
				column_endpoint(v->parent.input->bounds + j);
		for (j = 0; j < v->child.columns; j++)
			basis->status[v->map[j]] = seed->status[j];
		for (i = 0; i < v->parent.rows; i++) {
			j = seed->index[i];
			basis->index[i] =
				j < v->child.columns
					? v->map[j]
					: v->parent.columns +
						  (j - v->child.columns);
			basis->status[v->parent.columns + i] =
				seed->status[v->child.columns + i];
		}
		basis->valid = 1;
	}
	{
		int event;
		column_trace_link(
			v, SIMPLEX_COLUMN_BASIS_LIFT, rejection, 0, &event);
		if (!rejection) {
			basis->producer = simplex_execution_ref(
				&v->parent.input->execution, event);
		}
		return rejection;
	}
}

int simplex_column_basis_restrict(const struct simplex_ColumnProjection *v,
	const struct simplex_ColumnBasis *basis,
	struct simplex_DualBasisSeed *seed, int prefer_dual,
	unsigned char *scratch, size_t bytes)
{
	int i, j, rejection = SIMPLEX_COLUMN_MODEL;
	if (!column_current(v))
		return rejection;
	if (!column_maps(v))
		rejection = SIMPLEX_COLUMN_MAP;
	else if (!basis || !basis->valid ||
		 !simplex_problem_identity_matches(
			 &basis->parent, v->parent.input))
		rejection = SIMPLEX_COLUMN_MODEL;
	else if (!seed || !seed->index || !seed->status)
		rejection = SIMPLEX_COLUMN_STORAGE;
	else if (!(rejection = column_layout(v->parent.rows, v->parent.columns,
			   basis->index, basis->status, scratch, bytes))) {
		for (i = 0; i < v->parent.rows; i++) {
			j = basis->index[i];
			if (j < v->parent.columns && v->reverse[j] < 0) {
				rejection = SIMPLEX_COLUMN_LAYOUT;
				break;
			}
		}
		if (!rejection) {
			simplex_dual_seed_init(seed, v->child.input,
				seed->index, seed->status);
			for (i = 0; i < v->parent.rows; i++) {
				j = basis->index[i];
				seed->index[i] =
					j < v->parent.columns
						? v->reverse[j]
						: v->child.columns +
							  (j - v->parent.columns);
				seed->status[v->child.columns + i] =
					basis->status[v->parent.columns + i];
			}
			for (j = 0; j < v->child.columns; j++)
				seed->status[j] = basis->status[v->map[j]];
			seed->valid = 1;
			seed->prefer_dual = prefer_dual;
		}
	}
	return column_trace_link(v, SIMPLEX_COLUMN_BASIS_RESTRICT, rejection,
		basis ? simplex_execution_ref_id(
				&basis->producer, &v->child.input->execution, 1)
		      : 0,
		NULL);
}

int simplex_column_point_lift(const struct simplex_ColumnProjection *v,
	const double *point, double *staging, double tolerance)
{
	int j, rejection = SIMPLEX_COLUMN_POINT;
	if (!column_current(v))
		return SIMPLEX_COLUMN_MODEL;
	if (!column_maps(v))
		return column_trace(
			v, SIMPLEX_COLUMN_POINT_LIFT, SIMPLEX_COLUMN_MAP);
	if (point && staging && point != staging && isfinite(tolerance) &&
		tolerance > 0. &&
		simplex_problem_primal_infeasibility(v->child.input, point) <=
			tolerance) {
		memset(staging, 0, (size_t)v->parent.columns * sizeof(double));
		for (j = 0; j < v->child.columns; j++)
			staging[v->map[j]] = point[j];
		if (simplex_problem_primal_infeasibility(
			    v->parent.input, staging) <= tolerance)
			rejection = SIMPLEX_COLUMN_OK;
	}
	return column_trace(v, SIMPLEX_COLUMN_POINT_LIFT, rejection);
}

int simplex_column_rows_lift(
	const struct simplex_ColumnProjection *v, const double *y)
{
	int i, rejection = SIMPLEX_COLUMN_MULTIPLIER;
	if (!column_current(v))
		return SIMPLEX_COLUMN_MODEL;
	if (!column_maps(v))
		return column_trace(
			v, SIMPLEX_COLUMN_ROWS_LIFT, SIMPLEX_COLUMN_MAP);
	if (y) {
		for (i = 0; i < v->parent.rows; i++)
			if (!isfinite(y[i]))
				break;
		if (i == v->parent.rows)
			rejection = SIMPLEX_COLUMN_OK;
	}
	return column_trace(v, SIMPLEX_COLUMN_ROWS_LIFT, rejection);
}
