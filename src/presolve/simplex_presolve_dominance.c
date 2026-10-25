/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_presolve_dominance.c
 * @brief Presolve dominance reductions and recovery metadata.
 */

#include "policy/resources.h"
#include "policy/controls.h"
#include "simplex_presolve_dominance.h"
#include "simplex_presolve_rules.h"
#include "utils.h"
#include "simplex_storage.h"
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define WORD_BITS (sizeof(unsigned long) * CHAR_BIT)

struct support_Fingerprint {
	unsigned long present;
	unsigned long required;
};

struct dominance_ColumnOrder {
	int column, degree;
	double cost;
};

struct dominance_BitStorage {
	size_t words, bits, candidates;
};

struct dominance_PostingStorage {
	size_t order, starts, postings, cursor;
};

/* The allowance bounds all simultaneously owned optional index arrays,
 * including retained cover/fingerprints and temporary ordering/cursors. */
static int dominance_bit_storage(int rows, int columns, size_t retained,
	size_t limit, struct dominance_BitStorage *storage)
{
	struct dominance_BitStorage next = {0};
	size_t remaining;

	*storage = next;
	if (rows <= 0 || columns <= 0 || retained > limit)
		return 0;
	remaining = limit - retained;
	next.words = (size_t)columns / WORD_BITS +
		     ((size_t)columns % WORD_BITS != 0);
	if (!simplex_storage_bytes(next.words, sizeof(unsigned long), remaining,
		    &next.candidates) ||
		!simplex_storage_bytes((size_t)rows, next.candidates,
			remaining - next.candidates, &next.bits))
		return 0;
	*storage = next;
	return 1;
}

static int dominance_posting_storage(int rows, int count, size_t entries,
	size_t retained, size_t limit, struct dominance_PostingStorage *storage)
{
	struct dominance_PostingStorage next = {0};
	size_t remaining;

	*storage = next;
	if (rows <= 0 || count <= 0 || !entries || entries > INT_MAX ||
		retained > limit)
		return 0;
	remaining = limit - retained;
	if (!simplex_storage_bytes((size_t)count,
		    sizeof(struct dominance_ColumnOrder), remaining,
		    &next.order))
		return 0;
	remaining -= next.order;
	if (!simplex_storage_bytes(
		    (size_t)rows + 1, sizeof(int), remaining, &next.starts))
		return 0;
	remaining -= next.starts;
	if (!simplex_storage_bytes(
		    entries, sizeof(int), remaining, &next.postings))
		return 0;
	remaining -= next.postings;
	if (!simplex_storage_bytes(
		    (size_t)rows, sizeof(int), remaining, &next.cursor))
		return 0;
	*storage = next;
	return 1;
}

/* Collisions only admit extra candidates. Every accepted reduction is still
 * checked against the original coefficients, bounds and row directions. */
static unsigned long row_fingerprint(int row)
{
	unsigned long hash = (unsigned long)row + 1;
	hash ^= hash >> 16;
	hash *= 2654435761UL;
	hash ^= hash >> 13;
	return (1UL << (hash % WORD_BITS)) |
	       (1UL << ((hash / WORD_BITS) % WORD_BITS));
}

static int compare_column_order(const void *left, const void *right)
{
	const struct dominance_ColumnOrder *a = left, *b = right;
	if (a->cost != b->cost)
		return a->cost < b->cost ? -1 : 1;
	if (a->degree != b->degree)
		return a->degree > b->degree ? -1 : 1;
	return a->column < b->column ? -1 : a->column != b->column;
}

static int spend(size_t *work, size_t amount)
{
	if (*work < amount) {
		*work = 0;
		return 0;
	}
	*work -= amount;
	return 1;
}

static int nonnegative(const struct optm_VariableBound *b)
{
	return simplex_presolve_has_lower(b) && isfinite(b->lb) && b->lb >= 0.;
}

static int movable(const struct optm_VariableBound *b)
{
	return simplex_presolve_has_lower(b) && isfinite(b->lb) &&
	       (!simplex_presolve_has_upper(b) || b->ub > b->lb);
}

/* Public sparse models may use other layouts. This optional merge-based
 * search must not infer dominance from unsorted or repeated indices. */
static int ordered_storage(const struct lp_Model *m)
{
	int i, k;
	for (i = 0; i < m->m; i++)
		for (k = m->row_start[i]; k < m->row_start[i + 1]; k++)
			if (!isfinite(m->row_value[k]) ||
				(k > m->row_start[i] &&
					m->column_index[k] <=
						m->column_index[k - 1]))
				return 0;
	for (i = 0; i < m->n; i++)
		for (k = m->column_start[i]; k < m->column_start[i + 1]; k++)
			if (!isfinite(m->value[k]) ||
				(k > m->column_start[i] &&
					m->row_index[k] <= m->row_index[k - 1]))
				return 0;
	return 1;
}

static double direction(const struct lp_Model *m, int row)
{
	return m->constraints[row].type == optm_CONS_T_LE ? -1. : 1.;
}

/* Compact movable candidates, ordered by cost and support size. Ordering
 * only accelerates search; it never establishes dominance by itself. */
static void candidate_postings(const struct lp_Model *m,
	const struct optm_VariableBound *bounds, const unsigned char *rows,
	const unsigned char *columns, const int *degree, size_t retained,
	int **starts_out, int **postings_out)
{
	struct dominance_ColumnOrder *order = NULL;
	int *starts = NULL, *postings = NULL, *cursor = NULL;
	int j, k, i, count = 0;
	size_t entries = 0, maximum;
	struct dominance_PostingStorage storage;
	*starts_out = *postings_out = NULL;
	if (m->m <= 0 || m->n <= 0 || retained > DOMINANCE_BYTES)
		return;
	maximum = (DOMINANCE_BYTES - retained) / sizeof(int);
	if (maximum > INT_MAX)
		maximum = INT_MAX;
	for (j = 0; j < m->n; j++)
		if (columns[j] && movable(bounds + j) &&
			isfinite(m->objective[j])) {
			if (degree[j] < 0 ||
				(size_t)degree[j] > maximum - entries)
				return;
			count++;
			entries += (size_t)degree[j];
		}
	if (!dominance_posting_storage(
		    m->m, count, entries, retained, DOMINANCE_BYTES, &storage))
		return;
	order = lp_simplex_malloc(storage.order);
	starts = lp_simplex_malloc(storage.starts);
	postings = lp_simplex_malloc(storage.postings);
	cursor = lp_simplex_malloc(storage.cursor);
	if (!order || !starts || !postings || !cursor)
		goto done;
	count = 0;
	for (j = 0; j < m->n; j++)
		if (columns[j] && movable(bounds + j) &&
			isfinite(m->objective[j])) {
			order[count].column = j;
			order[count].cost = m->objective[j];
			order[count++].degree = degree[j];
		}
	qsort(order, (size_t)count, sizeof(*order), compare_column_order);
	lp_simplex_memset(starts, 0, storage.starts);
	for (j = 0; j < count; j++)
		for (k = m->column_start[order[j].column];
			k < m->column_start[order[j].column + 1]; k++)
			if (rows[m->row_index[k]])
				starts[m->row_index[k] + 1]++;
	for (i = 0; i < m->m; i++)
		starts[i + 1] += starts[i];
	lp_simplex_memcpy(cursor, starts, storage.cursor);
	for (j = 0; j < count; j++)
		for (k = m->column_start[order[j].column];
			k < m->column_start[order[j].column + 1]; k++) {
			i = m->row_index[k];
			if (rows[i])
				postings[cursor[i]++] = order[j].column;
		}
	*starts_out = starts;
	*postings_out = postings;
	starts = postings = NULL;
done:
	lp_simplex_free(cursor);
	lp_simplex_free(starts);
	lp_simplex_free(postings);
	lp_simplex_free(order);
}

/* All variables in a monotone covering row are nonnegative. This is a
 * sufficient condition, not an inference from the instance name or size. */
static int covering_row(const struct lp_Model *m, int row,
	const struct optm_VariableBound *bounds, const unsigned char *active)
{
	int k, j;
	double sign = direction(m, row);
	if (m->constraints[row].type == optm_CONS_T_EQ)
		return 0;
	for (k = m->row_start[row]; k < m->row_start[row + 1]; k++) {
		j = m->column_index[k];
		if (active[j] && (!nonnegative(bounds + j) ||
					 sign * m->row_value[k] < 0.))
			return 0;
	}
	return 1;
}

/* Prove target >= source componentwise on x>=0 and rhs_source>=rhs_target.
 * No feasibility tolerance is used to discard a genuinely stronger row. */
static int row_implies(const struct lp_Model *m, int source, int target,
	const unsigned char *active, size_t *work)
{
	int a = m->row_start[source], b = m->row_start[target];
	double sa = direction(m, source), sb = direction(m, target);
	while (a < m->row_start[source + 1]) {
		int j = m->column_index[a];
		if (!spend(work, 1))
			return 0;
		if (!active[j]) {
			a++;
			continue;
		}
		while (b < m->row_start[target + 1] && m->column_index[b] < j) {
			if (!spend(work, 1))
				return 0;
			b++;
		}
		if (b == m->row_start[target + 1] || m->column_index[b] != j ||
			sa * m->row_value[a] > sb * m->row_value[b])
			return 0;
		a++;
	}
	return 1;
}

/* Certify a*upper>=rhs. fma recovers the sign of product rounding at equality;
 * skip the subnormal range where that residual itself could underflow. */
static int capacity(double a, double upper, double rhs)
{
	double product;
	int ea, eb;
	if (rhs <= 0.)
		return 1;
	if (!(a > 0.) || !isfinite(upper))
		return 0;
	(void)frexp(a, &ea);
	(void)frexp(upper, &eb);
	if (ea + eb < DBL_MIN_EXP + 2 * DBL_MANT_DIG)
		return 0;
	product = a * upper;
	if (!isfinite(product))
		return product > 0.;
	if (product < rhs)
		return 0;
	if (product > rhs)
		return 1;
	return fma(a, upper, -product) >= 0.;
}

/* Transfer delta from source to target. An unbounded target works for
 * arbitrary signed rows (equalities require identical coefficients).
 * A bounded target requires monotone rows and enough capacity to satisfy
 * each affected row alone when the transfer saturates its upper bound. */
static int column_dominates(const struct lp_Model *m, int source, int target,
	const struct optm_VariableBound *bounds, const double *rhs,
	const unsigned char *rows, const unsigned char *cover, size_t *work)
{
	const struct optm_VariableBound *b = bounds + target;
	int bounded = simplex_presolve_has_upper(b) && isfinite(b->ub);
	int a = m->column_start[source], z = m->column_start[target];
	if (!movable(b) || !isfinite(m->objective[target]) ||
		!isfinite(m->objective[source]) ||
		m->objective[target] > m->objective[source])
		return 0;
	if (bounded &&
		(bounds[source].lb != 0. || !nonnegative(b) ||
			m->objective[source] < 0. || m->objective[target] < 0.))
		return 0;
	while (a < m->column_start[source + 1] ||
		z < m->column_start[target + 1]) {
		int ra = a < m->column_start[source + 1] ? m->row_index[a]
							 : m->m;
		int rb = z < m->column_start[target + 1] ? m->row_index[z]
							 : m->m;
		int row = ra < rb ? ra : rb;
		double av = 0., bv = 0., sign;
		if (!spend(work, 1))
			return 0;
		if (ra == row)
			av = m->value[a++];
		if (rb == row)
			bv = m->value[z++];
		if (!rows[row])
			continue;
		sign = direction(m, row);
		if (m->constraints[row].type == optm_CONS_T_EQ) {
			if (av != bv)
				return 0;
		} else if (sign * bv < sign * av)
			return 0;
		if (bounded && (!cover[row] || !capacity(sign * bv, b->ub,
						       sign * rhs[row])))
			return 0;
	}
	return 1;
}

static int posting_lower_bound(
	const struct lp_Model *m, int low, int high, int column, size_t *work)
{
	int end = high;
	while (low < high) {
		int middle = low + (high - low) / 2;
		if (!spend(work, 1))
			return end;
		if (m->column_index[middle] < column)
			low = middle + 1;
		else
			high = middle;
	}
	return low;
}

/* Intersect necessary row supports without a dense bit index. When a row
 * first contains a column beyond target, every intervening column can be
 * skipped in the anchor posting list as well. */
static int next_support_candidate(const struct lp_Model *m, int source,
	int target, const unsigned char *rows, size_t *work)
{
	int k;
	for (k = m->column_start[source]; k < m->column_start[source + 1];
		k++) {
		int row = m->row_index[k], at;
		if (!rows[row] ||
			(m->constraints[row].type != optm_CONS_T_EQ &&
				direction(m, row) * m->value[k] <= 0.))
			continue;
		at = posting_lower_bound(m, m->row_start[row],
			m->row_start[row + 1], target, work);
		if (at == m->row_start[row + 1])
			return m->n;
		if (m->column_index[at] > target)
			return m->column_index[at];
	}
	return target;
}

void simplex_presolve_dominance(const struct lp_Model *m,
	struct optm_VariableBound *bounds, const double *rhs,
	unsigned char *rows, const unsigned char *columns,
	const int *row_degree, int *column_degree, size_t *work,
	int *removed_rows, int *fixed_columns)
{
	unsigned char *cover;
	unsigned long *bits = NULL, *candidates = NULL;
	struct support_Fingerprint *fingerprints = NULL;
	int *posting_starts = NULL, *postings = NULL;
	size_t words = 0, w, cover_bytes, fingerprint_bytes;
	struct dominance_BitStorage bit_storage;
	size_t initial_work = *work;
	int i, j = 0, k, anchor, best;
	*removed_rows = *fixed_columns = 0;
	/* Dense-only models keep the established rules; never build an m*n
	 * matrix here. Sparse CSR/CSC views must have sorted, unique indices. */
	if (!*work || m->m <= 0 || m->n <= 0 || !m->row_start ||
		!m->column_start ||
		!simplex_storage_bytes((size_t)m->m, sizeof(*cover),
			DOMINANCE_BYTES, &cover_bytes) ||
		!ordered_storage(m))
		return;
	cover = (unsigned char *)lp_simplex_malloc(cover_bytes);
	if (!cover)
		return;
	for (i = 0; i < m->m; i++)
		cover[i] = rows[i] && covering_row(m, i, bounds, columns);
	for (i = 0; i < m->m && *work; i++) {
		if (!rows[i] || !cover[i] || row_degree[i] < 2)
			continue;
		anchor = -1;
		best = INT_MAX;
		for (k = m->row_start[i]; k < m->row_start[i + 1]; k++) {
			j = m->column_index[k];
			if (columns[j] && column_degree[j] < best) {
				anchor = j;
				best = column_degree[j];
			}
		}
		if (anchor < 0)
			continue;
		for (k = m->column_start[anchor];
			k < m->column_start[anchor + 1] && *work; k++) {
			int target = m->row_index[k], q;
			if (!spend(work, 1))
				break;
			if (target == i || !rows[target] || !cover[target] ||
				row_degree[target] < row_degree[i] ||
				direction(m, i) * rhs[i] <
					direction(m, target) * rhs[target] ||
				!row_implies(m, i, target, columns, work))
				continue;
			rows[target] = 0;
			(*removed_rows)++;
			for (q = m->row_start[target];
				q < m->row_start[target + 1]; q++)
				if (columns[m->column_index[q]])
					column_degree[m->column_index[q]]--;
		}
	}
	/* Small bit indexes accelerate support intersections. Larger models use
	 * a rarest-row posting list instead, with the same shared work budget. */
	if (!*work) {
		j = 0;
		goto done;
	}
	if (dominance_bit_storage(
		    m->m, m->n, cover_bytes, DOMINANCE_BYTES, &bit_storage)) {
		words = bit_storage.words;
		bits = (unsigned long *)lp_simplex_malloc(bit_storage.bits);
		candidates = (unsigned long *)lp_simplex_malloc(
			bit_storage.candidates);
		if (!bits || !candidates) {
			lp_simplex_free(bits);
			lp_simplex_free(candidates);
			bits = candidates = NULL;
		}
	}
	if (bits) {
		lp_simplex_memset(bits, 0, bit_storage.bits);
		for (i = 0; i < m->m; i++)
			if (rows[i])
				for (k = m->row_start[i];
					k < m->row_start[i + 1]; k++) {
					j = m->column_index[k];
					if (columns[j])
						bits[(size_t)i * words +
							j / WORD_BITS] |=
							1UL << (j % WORD_BITS);
				}
	} else if (simplex_storage_bytes((size_t)m->n, sizeof(*fingerprints),
			   DOMINANCE_BYTES - cover_bytes, &fingerprint_bytes)) {
		fingerprints = (struct support_Fingerprint *)lp_simplex_malloc(
			fingerprint_bytes);
		if (fingerprints) {
			lp_simplex_memset(fingerprints, 0, fingerprint_bytes);
			for (j = 0; j < m->n; j++)
				if (columns[j])
					for (k = m->column_start[j];
						k < m->column_start[j + 1];
						k++) {
						i = m->row_index[k];
						if (!rows[i])
							continue;
						fingerprints[j].present |=
							row_fingerprint(i);
						if (m->constraints[i].type ==
								optm_CONS_T_EQ ||
							direction(m,
								i) * m->value[k] >
								0.)
							fingerprints[j]
								.required |=
								row_fingerprint(
									i);
					}
		}
	}
	if (fingerprints)
		candidate_postings(m, bounds, rows, columns, column_degree,
			cover_bytes + fingerprint_bytes, &posting_starts,
			&postings);
	for (j = 0; j < m->n && *work; j++) {
		if (!columns[j] || !movable(bounds + j))
			continue;
		anchor = -1;
		best = INT_MAX;
		for (k = m->column_start[j]; k < m->column_start[j + 1]; k++) {
			i = m->row_index[k];
			if (rows[i] &&
				(m->constraints[i].type == optm_CONS_T_EQ ||
					direction(m, i) * m->value[k] > 0.) &&
				row_degree[i] < best) {
				anchor = i;
				best = row_degree[i];
			}
		}
		if (anchor < 0)
			continue;
		if (bits) {
			if (!spend(work, (words + 31) / 32))
				break;
			lp_simplex_memcpy(candidates,
				bits + (size_t)anchor * words,
				words * sizeof(*bits));
			for (k = m->column_start[j]; k < m->column_start[j + 1];
				k++) {
				i = m->row_index[k];
				if (!rows[i] || i == anchor ||
					(m->constraints[i].type !=
							optm_CONS_T_EQ &&
						direction(m, i) * m->value[k] <=
							0.))
					continue;
				if (!spend(work, (words + 31) / 32))
					break;
				for (w = 0; w < words; w++)
					candidates[w] &=
						bits[(size_t)i * words + w];
			}
			if (!*work)
				break;
			for (w = 0; w < words && *work; w++) {
				unsigned long mask = candidates[w];
				size_t bit = 0;
				while (mask && *work) {
					int target;
					while (!(mask & 1UL)) {
						mask >>= 1;
						bit++;
					}
					target = (int)(w * WORD_BITS + bit);
					mask >>= 1;
					bit++;
					if (!spend(work, 1))
						break;
					if (target == j || !columns[target] ||
						!column_dominates(m, j, target,
							bounds, rhs, rows,
							cover, work))
						continue;
					simplex_presolve_set_upper(
						bounds + j, bounds[j].lb);
					(*fixed_columns)++;
					break;
				}
				if (!movable(bounds + j))
					break;
			}
		} else {
			int start = postings ? posting_starts[anchor]
					     : m->row_start[anchor];
			int end = postings ? posting_starts[anchor + 1]
					   : m->row_start[anchor + 1];
			for (k = start; k < end && *work; k++) {
				int target = postings ? postings[k]
						      : m->column_index[k];
				int next;
				if (fingerprints) {
					/* Charge cheap word filters in batches, as for the
					 * bit index. Exact candidate checks charge separately. */
					if ((k - start) % 32 == 0 &&
						!spend(work, 1))
						break;
					if (postings && m->objective[target] >
								m->objective[j])
						break;
					if (target == j || !columns[target] ||
						!movable(bounds + target) ||
						m->objective[target] >
							m->objective[j] ||
						(fingerprints[j].required &
							fingerprints[target]
								.present) !=
							fingerprints[j]
								.required)
						continue;
					if (column_dominates(m, j, target,
						    bounds, rhs, rows, cover,
						    work)) {
						simplex_presolve_set_upper(
							bounds + j,
							bounds[j].lb);
						(*fixed_columns)++;
						break;
					}
					continue;
				}
				if (!spend(work, 1))
					break;
				next = next_support_candidate(
					m, j, target, rows, work);
				if (next == m->n)
					break;
				if (next > target) {
					k = posting_lower_bound(m, k + 1,
						    m->row_start[anchor + 1],
						    next, work) -
					    1;
					continue;
				}
				if (target != j && columns[target] &&
					column_dominates(m, j, target, bounds,
						rhs, rows, cover, work)) {
					simplex_presolve_set_upper(
						bounds + j, bounds[j].lb);
					(*fixed_columns)++;
					break;
				}
			}
		}
	}
done:
	lp_simplex_free(candidates);
	lp_simplex_free(bits);
	lp_simplex_free(cover);
	lp_simplex_free(fingerprints);
	lp_simplex_free(postings);
	lp_simplex_free(posting_starts);
	if (getenv(SIMPLEX_ENV_PROFILE))
		fprintf(stderr,
			"presolve dominance: columns-scanned=%d/%d work=%lu remaining=%lu removed-rows=%d fixed-columns=%d\n",
			j, m->n, (unsigned long)(initial_work - *work),
			(unsigned long)*work, *removed_rows, *fixed_columns);
}
