/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_klu_sparse.c
 * @brief Sparse event propagation through KLU block triangular solve graphs.
 */

#include "simplex_klu_sparse.h"

#ifdef LP_SIMPLEX_HAVE_KLU
#include "policy/resources.h"
#include "simplex_storage.h"
#include "simplex_sparse_vector.h"
#include "utils.h"
#include <limits.h>
#include <math.h>

/* BTF gives P*A*Q = L*U + F after row scaling. A normal solve executes
 * each diagonal block's L columns forward, U columns backward, and sends F
 * contributions to earlier blocks. Separate L, U and off-block nodes retain
 * KLU's arithmetic order while forming a DAG. Reversing its edges/order gives the transpose
 * solve, with the same U diagonal divisions and the opposite permutations.
 * Nodes are ranked in the full triangular traversal order. A priority queue
 * visits only activated nodes while retaining that order. No magnitude-based
 * dropping, approximate factorization, or new feasibility decision is used. */
struct simplex_KluSparseImpl {
	int n, nodes, edges;
	int *integer;
	double *real;
	int *lower, *upper, *off, *tag, *column_inverse;
	int *start[2], *index[2];
	double *value[2];
	double *diagonal, *work, *output;
	int *heap, *queued, *visited, *output_index;
	unsigned int *queue_level[7];
	int queue_depth;
	const int32_t *row_permutation, *column_permutation;
	const int32_t *row_inverse;
	const double *scale;
};

void simplex_klu_sparse_destroy(struct simplex_KluSparse *cache)
{
	if (cache->impl != NULL) {
		lp_simplex_free(cache->impl->integer);
		lp_simplex_free(cache->impl->real);
		lp_simplex_free(cache->impl);
		cache->impl = NULL;
	}
}

static int klu_cache_bytes(
	size_t count, size_t width, size_t *remaining, size_t *bytes)
{
	if (!simplex_storage_bytes(count, width, *remaining, bytes))
		return 0;
	*remaining -= *bytes;
	return 1;
}

static int klu_graph_edge(struct simplex_KluSparseImpl *c, int source,
	int target, double value, int fill)
{
	int p;
	if (source < 0 || target >= c->nodes || source >= target ||
		!isfinite(value))
		return 0;
	if (!fill) {
		c->start[0][source + 1]++;
		c->start[1][target + 1]++;
	} else {
		p = c->heap[source]++;
		c->index[0][p] = target;
		c->value[0][p] = value;
		p = c->queued[target]++;
		c->index[1][p] = source;
		c->value[1][p] = value;
	}
	return 1;
}

/* A hierarchical bit queue preserves the exact rank order of a heap while
 * using a bounded number of word operations. Reuse graph-construction scratch;
 * no extra allocation or factor-dependent hash table is needed. */
static int klu_queue_initialize(struct simplex_KluSparseImpl *c)
{
	int size = c->nodes, offset = 0;
	c->queue_depth = 0;
	while (size > 1) {
		int words = size / 32 + (size % 32 != 0);
		if (c->queue_depth >= 7 || offset > c->nodes - words)
			return 0;
		c->queue_level[c->queue_depth++] =
			(unsigned int *)(c->heap + offset);
		offset += words;
		size = words;
	}
	lp_simplex_memset(c->heap, 0, (size_t)offset * sizeof(int));
	return c->queue_depth > 0;
}

int simplex_klu_sparse_build(struct simplex_KluSparse *cache, klu_symbolic *s,
	klu_numeric *f, klu_common *common)
{
	struct simplex_KluSparse next;
	struct simplex_KluSparseImpl *c;
	int n = s->n, nodes, block, i, k, pass, rank = 0, singletons = 0;
	int *temporary = NULL, *lp, *li, *up, *ui, *ip;
	double *temporary_value = NULL, *lx, *ux, *xp;
	size_t edges, integers, reals, remaining = SIMPLEX_KLU_SPARSE_BYTES;
	size_t ib, rb, tib, trb;
	next.impl = NULL;
	simplex_klu_sparse_destroy(cache);
	if (sizeof(unsigned int) * CHAR_BIT != 32 || n < 64 ||
		n > INT_MAX / 3 || f->lnz < n || f->unz < n || f->nzoff < 0 ||
		(size_t)f->lnz > remaining / (sizeof(int) + sizeof(double)) ||
		(size_t)f->unz > remaining / (sizeof(int) + sizeof(double)) ||
		(size_t)f->nzoff > remaining / (sizeof(int) + sizeof(double)))
		return 1;
	for (block = 0; block < s->nblocks; block++)
		singletons += s->R[block + 1] - s->R[block] == 1;
	nodes = 3 * n - 2 * singletons;
	edges = (size_t)f->lnz + f->unz + f->nzoff - 2 * singletons;
	/* Counts are bounded before narrow additions and pointer partitions. */
	if (edges > INT_MAX || remaining < sizeof(*c) ||
		(size_t)n >
			remaining / (23 * sizeof(int) + 7 * sizeof(double)) ||
		edges > remaining / (2 * (sizeof(int) + sizeof(double))))
		return 1;
	remaining -= sizeof(*c);
	integers = (size_t)23 * n + 2 + 2 * edges;
	reals = (size_t)7 * n + 2 * edges;
	if (!klu_cache_bytes(integers, sizeof(int), &remaining, &ib) ||
		!klu_cache_bytes(reals, sizeof(double), &remaining, &rb) ||
		!klu_cache_bytes((size_t)2 * (n + 1) + f->lnz + f->unz,
			sizeof(int), &remaining, &tib) ||
		!klu_cache_bytes((size_t)f->lnz + f->unz, sizeof(double),
			&remaining, &trb))
		return 1;
	c = (struct simplex_KluSparseImpl *)lp_simplex_malloc(sizeof(*c));
	if (c == NULL)
		return 1;
	lp_simplex_memset(c, 0, sizeof(*c));
	next.impl = c;
	c->integer = (int *)lp_simplex_malloc(ib);
	c->real = (double *)lp_simplex_malloc(rb);
	temporary = (int *)lp_simplex_malloc(tib);
	temporary_value = (double *)lp_simplex_malloc(trb);
	if (!c->integer || !c->real || !temporary || !temporary_value)
		goto REFUSED;
	lp_simplex_memset(c->integer, 0, ib);
	lp_simplex_memset(c->real, 0, rb);
	c->n = n;
	c->nodes = nodes;
	c->edges = (int)edges;
	ip = c->integer;
	c->lower = ip;
	ip += n;
	c->upper = ip;
	ip += n;
	c->off = ip;
	ip += n;
	c->column_inverse = ip;
	ip += n;
	c->tag = ip;
	ip += nodes;
	c->start[0] = ip;
	ip += nodes + 1;
	c->start[1] = ip;
	ip += nodes + 1;
	c->index[0] = ip;
	ip += edges;
	c->index[1] = ip;
	ip += edges;
	c->heap = ip;
	ip += nodes;
	c->queued = ip;
	ip += nodes;
	c->visited = ip;
	ip += nodes;
	c->output_index = ip;
	xp = c->real;
	c->value[0] = xp;
	xp += edges;
	c->value[1] = xp;
	xp += edges;
	c->diagonal = xp;
	xp += nodes;
	c->work = xp;
	xp += nodes;
	c->output = xp;
	c->row_permutation = f->Pnum;
	c->row_inverse = f->Pinv;
	c->column_permutation = s->Q;
	for (i = 0; i < n; i++)
		c->column_inverse[s->Q[i]] = i;
	c->scale = f->Rs;
	lp = temporary;
	up = lp + n + 1;
	li = up + n + 1;
	ui = li + f->lnz;
	lx = temporary_value;
	ux = lx + f->lnz;
	if (!klu_extract(f, s, lp, li, lx, up, ui, ux, NULL, NULL, NULL, NULL,
		    NULL, NULL, NULL, common))
		goto REFUSED;
	for (block = s->nblocks - 1; block >= 0; block--) {
		if (s->R[block + 1] - s->R[block] == 1) {
			/* A singleton block has no triangular off-diagonals. Its two
			 * unit transfers can be contracted into the U division without
			 * changing its off-block gather/scatter position. */
			i = s->R[block];
			c->lower[i] = c->upper[i] = c->off[i] = rank;
			c->tag[rank] = 3 * i;
			c->diagonal[rank] = ((double *)f->Udiag)[i];
			if (!isfinite(c->diagonal[rank]) ||
				c->diagonal[rank] == 0.)
				goto REFUSED;
			rank++;
			continue;
		}
		for (i = s->R[block]; i < s->R[block + 1]; i++) {
			c->lower[i] = rank;
			c->tag[rank] = 3 * i;
			c->diagonal[rank++] = 1.;
		}
		for (i = s->R[block + 1] - 1; i >= s->R[block]; i--) {
			c->upper[i] = rank;
			c->tag[rank] = 3 * i + 1;
			c->diagonal[rank] = ((double *)f->Udiag)[i];
			if (!isfinite(c->diagonal[rank]) ||
				c->diagonal[rank] == 0.)
				goto REFUSED;
			rank++;
		}
		/* KLU sends off-block contributions only after both triangular
		 * solves, in ascending column order. Keep that rounding order. */
		for (i = s->R[block]; i < s->R[block + 1]; i++) {
			c->off[i] = rank;
			c->tag[rank] = 3 * i + 2;
			c->diagonal[rank++] = 1.;
		}
	}
	if (rank != nodes)
		goto REFUSED;
	for (pass = 0; pass < 2; pass++) {
		for (i = 0; i < n; i++) {
			if (c->lower[i] != c->upper[i] &&
				!klu_graph_edge(
					c, c->lower[i], c->upper[i], -1., pass))
				goto REFUSED;
			if (c->upper[i] != c->off[i] &&
				!klu_graph_edge(
					c, c->upper[i], c->off[i], -1., pass))
				goto REFUSED;
			for (k = lp[i]; k < lp[i + 1]; k++)
				if (li[k] != i &&
					!klu_graph_edge(c, c->lower[i],
						c->lower[li[k]], lx[k], pass))
					goto REFUSED;
			for (k = up[i]; k < up[i + 1]; k++)
				if (ui[k] != i &&
					!klu_graph_edge(c, c->upper[i],
						c->upper[ui[k]], ux[k], pass))
					goto REFUSED;
			for (k = f->Offp[i]; k < f->Offp[i + 1]; k++)
				if (!klu_graph_edge(c, c->off[i],
					    c->lower[f->Offi[k]],
					    ((double *)f->Offx)[k], pass))
					goto REFUSED;
		}
		if (!pass) {
			for (i = 0; i < nodes; i++) {
				c->start[0][i + 1] += c->start[0][i];
				c->start[1][i + 1] += c->start[1][i];
				c->heap[i] = c->start[0][i];
				c->queued[i] = c->start[1][i];
			}
			if (c->start[0][nodes] != c->edges)
				goto REFUSED;
		}
	}
	if (!klu_queue_initialize(c))
		goto REFUSED;
	lp_simplex_memset(c->queued, 0, (size_t)nodes * sizeof(int));
	lp_simplex_free(temporary);
	lp_simplex_free(temporary_value);
	cache->impl = c;
	return 0;
REFUSED:
	lp_simplex_free(temporary);
	lp_simplex_free(temporary_value);
	simplex_klu_sparse_destroy(&next);
	return 1;
}

static int klu_queue_bit(unsigned int value, int transpose)
{
#if defined(__GNUC__) || defined(__clang__)
	return transpose ? 31 - __builtin_clz(value) : __builtin_ctz(value);
#else
	int bit = 0;
	if (transpose) {
		while (value >>= 1)
			bit++;
	} else {
		while ((value & 1U) == 0) {
			value >>= 1;
			bit++;
		}
	}
	return bit;
#endif
}

static void klu_heap_push(
	struct simplex_KluSparseImpl *c, int node, int *count, int transpose)
{
	int level, coordinate = node;
	(void)transpose;
	(*count)++;
	c->queued[node] = 1;
	for (level = 0; level < c->queue_depth; level++) {
		int word = coordinate / 32;
		unsigned int old = c->queue_level[level][word];
		c->queue_level[level][word] = old | (1U << (coordinate % 32));
		if (old != 0)
			break;
		coordinate = word;
	}
}

static int klu_heap_pop(
	struct simplex_KluSparseImpl *c, int *count, int transpose)
{
	int level, node = 0, coordinate;
	for (level = c->queue_depth - 1; level >= 0; level--)
		node = 32 * node +
		       klu_queue_bit(c->queue_level[level][node], transpose);
	coordinate = node;
	for (level = 0; level < c->queue_depth; level++) {
		int word = coordinate / 32;
		c->queue_level[level][word] &= ~(1U << (coordinate % 32));
		if (c->queue_level[level][word] != 0)
			break;
		coordinate = word;
	}
	(*count)--;
	c->queued[node] = 0;
	return node;
}

static int klu_sparse_solve_impl(struct simplex_KluSparse *cache,
	double *vector, struct simplex_SparseVector *packed, int transpose)
{
	struct simplex_KluSparseImpl *c = cache->impl;
	int count = 0, outputs = 0, visited = 0, i, k, node, variable,
	    result = 1;
	if (c == NULL)
		return 1;
	if (packed && packed->count > c->n / 8)
		return 1;
	for (k = 0; k < (packed ? packed->count : c->n); k++) {
		double value;
		if (packed) {
			variable = packed->index[k];
			i = transpose ? c->column_inverse[variable]
				      : c->row_inverse[variable];
			value = packed->value[k];
		} else {
			i = k;
			value = vector[transpose ? c->column_permutation[i]
						 : c->row_permutation[i]];
		}
		if (value == 0.)
			continue;
		if (count >= c->n / 8)
			goto DECLINED;
		if (!transpose && c->scale)
			value /= c->scale[i];
		node = transpose ? c->off[i] : c->lower[i];
		c->work[node] = value;
		klu_heap_push(c, node, &count, transpose);
	}
	while (count) {
		double value;
		if (visited >= c->nodes / 4)
			goto DECLINED;
		node = klu_heap_pop(c, &count, transpose);
		c->visited[visited++] = node;
		value = c->work[node];
		/* Transposed triangular solves are column dots in KLU. Gather in
		 * that same stored-column order instead of scattering in row order.
		 * Reverse edges below activate dependencies without changing sums. */
		if (transpose)
			for (k = c->start[0][node]; k < c->start[0][node + 1];
				k++)
				value -= c->value[0][k] *
					 c->work[c->index[0][k]];
		value /= c->diagonal[node];
		c->work[node] = value;
		if (!isfinite(value)) {
			result = -1;
			goto DECLINED;
		}
		if (value == 0.)
			continue;
		variable = c->tag[node] / 3;
		if ((!transpose && c->off[variable] == node) ||
			(transpose && c->lower[variable] == node)) {
			double output = value;
			i = variable;
			if (transpose && c->scale)
				output /= c->scale[i];
			if (!isfinite(output)) {
				result = -1;
				goto DECLINED;
			}
			c->output_index[outputs] =
				transpose ? c->row_permutation[i]
					  : c->column_permutation[i];
			c->output[outputs++] = output;
		}
		for (k = c->start[transpose][node];
			k < c->start[transpose][node + 1]; k++) {
			int target = c->index[transpose][k];
			if (!transpose)
				c->work[target] -= c->value[0][k] * value;
			if (!c->queued[target] &&
				(transpose || c->work[target] != 0.))
				klu_heap_push(c, target, &count, transpose);
		}
	}
	if (packed) {
		simplex_sparse_vector_clear(packed);
		for (i = 0; i < outputs; i++)
			(void)simplex_sparse_vector_set(
				packed, c->output_index[i], c->output[i]);
	} else {
		lp_simplex_memset(vector, 0, (size_t)c->n * sizeof(double));
		for (i = 0; i < outputs; i++)
			vector[c->output_index[i]] = c->output[i];
	}
	for (i = 0; i < visited; i++)
		c->work[c->visited[i]] = 0.;
	return 0;
DECLINED:
	for (i = 0; i < visited; i++)
		c->work[c->visited[i]] = 0.;
	while (count) {
		node = klu_heap_pop(c, &count, transpose);
		c->work[node] = 0.;
	}
	return result;
}

int simplex_klu_sparse_solve(
	struct simplex_KluSparse *cache, double *vector, int transpose)
{
	return klu_sparse_solve_impl(cache, vector, NULL, transpose);
}

int simplex_klu_sparse_solve_packed(struct simplex_KluSparse *cache,
	struct simplex_SparseVector *vector, int transpose)
{
	if (cache->impl == NULL || vector->dimension != cache->impl->n)
		return 1;
	return klu_sparse_solve_impl(cache, NULL, vector, transpose);
}
#else
typedef int simplex_klu_sparse_unavailable;
#endif
