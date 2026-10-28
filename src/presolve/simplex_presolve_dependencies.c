/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_presolve_dependencies.c
 * @brief Exact sparse witnesses for overlapping marginal equalities.
 */
#include "simplex_presolve_dependencies.h"
#include "simplex_storage.h"
#include "utils.h"
#include <limits.h>
#include <math.h>

/* (p-1)^2 fits even a 32-bit unsigned long. A nonzero minor modulo this
 * prime proves rational independence; it does not certify the witnesses. */
#define DEPENDENCY_PRIME 65521UL

struct dependency_Term {
	int row;
	unsigned long value;
};

struct dependency_Pivot {
	struct dependency_Term *term;
	int count;
};

struct dependency_Row {
	int anchor, sign, base_sign, color, target, next, witness;
	unsigned long residue;
	struct dependency_Pivot pivot;
};

struct dependency_Column {
	int head, first, second;
	long sum;
};

struct dependency_Workspace {
	const struct lp_Model *model;
	const double *rhs;
	const unsigned char *row_active, *column_active;
	struct dependency_Row *row;
	struct dependency_Column *column;
	int *queue;
	size_t work, bytes, limit;
};

static int dependency_charge(struct dependency_Workspace *w, size_t count)
{
	if (count > w->work) {
		w->work = 0;
		return 0;
	}
	w->work -= count;
	return 1;
}

static void *dependency_allocate(
	struct dependency_Workspace *w, size_t count, size_t size)
{
	size_t bytes;
	void *result;
	if (!simplex_storage_bytes(count, size, w->limit - w->bytes, &bytes))
		return NULL;
	result = lp_simplex_malloc(bytes);
	if (result != NULL) {
		lp_simplex_memset(result, 0, bytes);
		w->bytes += bytes;
	}
	return result;
}

static void dependency_destroy(struct dependency_Workspace *w)
{
	int i;
	if (w->row != NULL)
		for (i = 0; i < w->model->m; i++)
			lp_simplex_free(w->row[i].pivot.term);
	lp_simplex_free(w->row);
	lp_simplex_free(w->column);
	lp_simplex_free(w->queue);
}

/* Only exact signed-unit equations with safely integral adjusted RHS enter
 * this optional rule. A fractional/fixed-column RHS simply loses the search. */
static int dependency_classify(struct dependency_Workspace *w)
{
	const struct lp_Model *m = w->model;
	int i, j, k, positive, negative, pos, neg;
	for (j = 0; j < m->n; j++)
		w->column[j].head = w->column[j].first = w->column[j].second =
			-1;
	for (i = 0; i < m->m; i++) {
		struct dependency_Row *r = w->row + i;
		r->anchor = r->target = r->next = -1;
		if (!w->row_active[i] ||
			m->constraints[i].type != optm_CONS_T_EQ ||
			!isfinite(w->rhs[i]) || fabs(w->rhs[i]) > 1073741823. ||
			floor(w->rhs[i]) != w->rhs[i])
			continue;
		positive = negative = 0;
		pos = neg = -1;
		for (k = m->row_start[i]; k < m->row_start[i + 1]; k++) {
			if (!dependency_charge(w, 1))
				return 0;
			j = m->column_index[k];
			if (!w->column_active[j])
				continue;
			if (m->row_value[k] == 1.) {
				positive++;
				pos = j;
			} else if (m->row_value[k] == -1.) {
				negative++;
				neg = j;
			} else {
				positive = negative = -1;
				break;
			}
		}
		if (positive > 0 && negative == 0)
			r->base_sign = 1;
		else if (negative > 0 && positive == 0)
			r->base_sign = -1;
		else if (positive >= 2 && negative == 1) {
			r->anchor = neg;
			r->sign = 1;
		} else if (negative >= 2 && positive == 1) {
			r->anchor = pos;
			r->sign = -1;
		}
		if (r->anchor >= 0) {
			r->next = w->column[r->anchor].head;
			w->column[r->anchor].head = i;
		}
	}
	return 1;
}

/* An anchor belongs to two homogeneous sum equations. Their incidence graph
 * supplies orientations, without row names, cardinality or RHS templates. */
static int dependency_members(struct dependency_Workspace *w)
{
	const struct lp_Model *m = w->model;
	int j, k, i, count;
	for (j = 0; j < m->n; j++) {
		struct dependency_Column *c = w->column + j;
		if (c->head < 0)
			continue;
		count = 0;
		for (k = m->column_start[j]; k < m->column_start[j + 1]; k++) {
			if (!dependency_charge(w, 1))
				return 0;
			i = m->row_index[k];
			if (w->row[i].base_sign) {
				if (count == 0)
					c->first = i;
				else if (count == 1)
					c->second = i;
				count++;
			}
		}
		if (count != 2)
			c->first = c->second = -1;
	}
	return 1;
}

static int dependency_color(struct dependency_Workspace *w)
{
	const struct lp_Model *m = w->model;
	int start, front, count, i, j, k, other;
	for (start = 0; start < m->m; start++) {
		if (!w->row[start].base_sign || w->row[start].color)
			continue;
		front = 0;
		count = 1;
		w->queue[0] = start;
		w->row[start].color = 1;
		while (front < count) {
			i = w->queue[front++];
			for (k = m->row_start[i]; k < m->row_start[i + 1];
				k++) {
				if (!dependency_charge(w, 1))
					return 0;
				j = m->column_index[k];
				if (!w->column_active[j] ||
					w->column[j].first < 0)
					continue;
				other = w->column[j].first == i
						? w->column[j].second
						: w->column[j].first;
				if (w->row[other].color == 0) {
					w->row[other].color = -w->row[i].color;
					w->queue[count++] = other;
				} else if (w->row[other].color ==
					   w->row[i].color)
					return 0;
			}
		}
	}
	return 1;
}

/* Intersect the two base memberships of every opposite endpoint. A marginal
 * row has one common target; incomplete or mixed neighborhoods are skipped. */
static int dependency_targets(struct dependency_Workspace *w)
{
	const struct lp_Model *m = w->model;
	int i, j, k, t, a, other, first, second, x, valid, found = 0;
	for (i = 0; i < m->m; i++) {
		x = w->row[i].anchor;
		if (x < 0 || w->column[x].first < 0)
			continue;
		first = second = -1;
		valid = 1;
		for (k = m->row_start[i]; k < m->row_start[i + 1]; k++) {
			j = m->column_index[k];
			if (!w->column_active[j] || j == x)
				continue;
			other = -1;
			for (t = m->column_start[j]; t < m->column_start[j + 1];
				t++) {
				if (!dependency_charge(w, 1))
					return 0;
				a = m->row_index[t];
				if (!w->row_active[a])
					continue;
				a = w->row[a].anchor;
				if (a < 0 || w->column[a].first < 0) {
					valid = 0;
					break;
				}
				if (a != x) {
					if (other >= 0 && other != a) {
						valid = 0;
						break;
					}
					other = a;
				}
			}
			if (!valid || other < 0) {
				valid = 0;
				break;
			}
			if (first < 0 && second < 0) {
				first = w->column[other].first;
				second = w->column[other].second;
			} else {
				a = w->column[other].first;
				t = w->column[other].second;
				if (first != a && first != t)
					first = -1;
				if (second != a && second != t)
					second = -1;
				if (first < 0 && second < 0) {
					valid = 0;
					break;
				}
			}
		}
		if (valid && ((first >= 0) != (second >= 0))) {
			w->row[i].target = first >= 0 ? first : second;
			found++;
		}
	}
	return found != 0;
}

static int dependency_add_long(long *sum, long value)
{
	if ((value > 0 && *sum > LONG_MAX - value) ||
		(value < 0 && *sum < LONG_MIN - value))
		return 0;
	*sum += value;
	return 1;
}

/* This exact product is the soundness boundary. A modular zero alone must
 * never establish an LP equation or permit removal. */
static int dependency_verify(struct dependency_Workspace *w)
{
	const struct lp_Model *m = w->model;
	int i, j, k, nonzero = 0, v;
	long rhs = 0;
	for (j = 0; j < m->n; j++)
		w->column[j].sum = 0;
	if (!dependency_charge(w, (size_t)m->n + (size_t)m->m))
		return 0;
	for (i = 0; i < m->m; i++) {
		v = w->row[i].witness;
		if (!v)
			continue;
		nonzero = 1;
		if (!dependency_add_long(&rhs, v * (long)w->rhs[i]))
			return 0;
		for (k = m->row_start[i]; k < m->row_start[i + 1]; k++) {
			if (!dependency_charge(w, 1))
				return 0;
			j = m->column_index[k];
			if (w->column_active[j] &&
				!dependency_add_long(&w->column[j].sum,
					v * (long)m->row_value[k]))
				return 0;
		}
	}
	if (!nonzero || rhs != 0)
		return 0;
	for (j = 0; j < m->n; j++)
		if (w->column[j].sum != 0)
			return 0;
	return 1;
}

static unsigned long dependency_inverse(unsigned long value)
{
	unsigned long power = DEPENDENCY_PRIME - 2, result = 1;
	while (power) {
		if (power & 1UL)
			result = result * value % DEPENDENCY_PRIME;
		value = value * value % DEPENDENCY_PRIME;
		power >>= 1;
	}
	return result;
}

static void dependency_accept(struct dependency_Workspace *w)
{
	const struct lp_Model *m = w->model;
	struct dependency_Term *terms;
	int i, j, k, count, position, tail;
	unsigned long value, inverse, product;
	if (!dependency_verify(w))
		return;
	for (i = 0; i < m->m; i++) {
		int v = w->row[i].witness;
		w->row[i].residue =
			v < 0 ? DEPENDENCY_PRIME - (unsigned long)(-v)
			      : (unsigned long)v;
	}
	for (position = 0; position < m->m; position++) {
		if (!dependency_charge(w, 1))
			return;
		j = w->queue[position];
		value = w->row[j].residue;
		if (!value)
			continue;
		if (w->row[j].pivot.term != NULL) {
			struct dependency_Pivot *pivot = &w->row[j].pivot;
			if (!dependency_charge(w, (size_t)pivot->count))
				return;
			for (k = 0; k < pivot->count; k++) {
				i = pivot->term[k].row;
				product = value * pivot->term[k].value %
					  DEPENDENCY_PRIME;
				w->row[i].residue =
					(w->row[i].residue + DEPENDENCY_PRIME -
						product) %
					DEPENDENCY_PRIME;
			}
			continue;
		}
		count = 0;
		if (!dependency_charge(w, (size_t)(m->m - position)))
			return;
		for (tail = position; tail < m->m; tail++)
			if (w->row[w->queue[tail]].residue)
				count++;
		terms = (struct dependency_Term *)dependency_allocate(
			w, (size_t)count, sizeof(*terms));
		if (terms == NULL) {
			w->work = 0;
			return;
		}
		inverse = dependency_inverse(value);
		count = 0;
		for (tail = position; tail < m->m; tail++) {
			i = w->queue[tail];
			if (w->row[i].residue) {
				terms[count].row = i;
				terms[count++].value = w->row[i].residue *
						       inverse %
						       DEPENDENCY_PRIME;
			}
		}
		w->row[j].pivot.term = terms;
		w->row[j].pivot.count = count;
		return;
	}
}

static void dependency_clear_witness(struct dependency_Workspace *w)
{
	int i;
	if (!dependency_charge(w, (size_t)w->model->m))
		return;
	for (i = 0; i < w->model->m; i++)
		w->row[i].witness = 0;
}

static void dependency_candidates(struct dependency_Workspace *w)
{
	const struct lp_Model *m = w->model;
	int x, i, a, b, first, second;
	for (x = 0; x < m->n && w->work; x++) {
		if (w->column[x].first < 0)
			continue;
		dependency_clear_witness(w);
		for (i = w->column[x].head; i >= 0; i = w->row[i].next) {
			if (!dependency_charge(w, 1))
				return;
			if (w->row[i].target >= 0)
				w->row[i].witness =
					w->row[i].sign *
					w->row[w->row[i].target].color;
		}
		dependency_accept(w);
	}
	for (a = 0; a < m->m && w->work; a++) {
		if (!w->row[a].base_sign)
			continue;
		for (b = a + 1; b < m->m && w->work; b++) {
			if (!dependency_charge(w, 1))
				return;
			if (!w->row[b].base_sign)
				continue;
			dependency_clear_witness(w);
			w->row[a].witness = w->row[a].base_sign;
			w->row[b].witness = -w->row[b].base_sign;
			for (i = 0; i < m->m; i++) {
				if (!dependency_charge(w, 1))
					return;
				x = w->row[i].anchor;
				if (x < 0 || w->row[i].target < 0)
					continue;
				first = w->column[x].first;
				second = w->column[x].second;
				if ((first == a || second == a) &&
					w->row[i].target == b)
					w->row[i].witness = w->row[i].sign;
				else if ((first == b || second == b) &&
					 w->row[i].target == a)
					w->row[i].witness = -w->row[i].sign;
			}
			dependency_accept(w);
		}
	}
	if (!w->work)
		return;
	dependency_clear_witness(w);
	for (i = 0; i < m->m; i++)
		if (w->row[i].base_sign)
			w->row[i].witness =
				w->row[i].base_sign * w->row[i].color;
	dependency_accept(w);
}

int simplex_presolve_dependencies(const struct lp_Model *model,
	const double *rhs, unsigned char *row_active,
	const unsigned char *column_active, int *column_degree,
	size_t work_limit, size_t byte_limit)
{
	struct dependency_Workspace w;
	int i, j, k, count, removed = 0;
	lp_simplex_memset(&w, 0, sizeof(w));
	if (model->m <= 0 || model->n <= 0 || !model->row_start ||
		!model->column_start || work_limit == 0)
		return 0;
	w.model = model;
	w.rhs = rhs;
	w.row_active = row_active;
	w.column_active = column_active;
	w.work = work_limit;
	w.limit = byte_limit;
	w.row = (struct dependency_Row *)dependency_allocate(
		&w, (size_t)model->m, sizeof(*w.row));
	w.column = (struct dependency_Column *)dependency_allocate(
		&w, (size_t)model->n, sizeof(*w.column));
	w.queue = (int *)dependency_allocate(&w, (size_t)model->m, sizeof(int));
	if (w.row && w.column && w.queue && dependency_classify(&w) &&
		dependency_members(&w) && dependency_color(&w) &&
		dependency_targets(&w)) {
		/* Keep homogeneous sum equations when another deletion set exists.
		 * They expose bounds directly and avoid expressing totals only through
		 * chains of signed marginals. This is a row-feature ordering, not a
		 * model identity or a different dependency acceptance criterion. */
		count = 0;
		for (i = 0; i < model->m; i++)
			if (!w.row[i].base_sign)
				w.queue[count++] = i;
		for (i = 0; i < model->m; i++)
			if (w.row[i].base_sign)
				w.queue[count++] = i;
		dependency_candidates(&w);
		/* A pivot minor in the verified witness matrix is invertible over Q.
		 * Delete its original equation columns simultaneously, preserving all
		 * variables and their bounds. No allocation follows this publication. */
		for (i = 0; i < model->m; i++)
			if (w.row[i].pivot.term != NULL) {
				row_active[i] = 0;
				removed++;
				for (k = model->row_start[i];
					k < model->row_start[i + 1]; k++) {
					j = model->column_index[k];
					if (column_active[j])
						column_degree[j]--;
				}
			}
	}
	dependency_destroy(&w);
	return removed;
}
