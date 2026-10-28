/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file test_klu_sparse.c
 * @brief KLU sparse solve equivalence, original residuals and optional refusal.
 */

#include "simplex_klu_sparse.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

static int requests, fail_at, active;

static void *cache_malloc(size_t bytes)
{
	void *p;
	requests++;
	if (requests == fail_at)
		return NULL;
	p = malloc(bytes);
	active += p != NULL;
	return p;
}

static void cache_free(void *p)
{
	if (p != NULL) {
		active--;
		free(p);
	}
}

#define lp_simplex_malloc cache_malloc
#define lp_simplex_free cache_free
#define simplex_klu_sparse_build cache_build
#define simplex_klu_sparse_destroy cache_destroy
#define simplex_klu_sparse_solve cache_solve
#define simplex_klu_sparse_solve_packed cache_solve_packed
#include "../src/basis/simplex_klu_sparse.c"

#define N 192
#define CHECK(c) \
	do { \
		if (!(c)) { \
			fprintf(stderr, "KLU sparse line %d: %s\n", __LINE__, \
				#c); \
			return 1; \
		} \
	} while (0)

static int check_case(int scaling, int btf, int chain, int singleton)
{
	double dense[N * N], values[12 * N], expected[N], actual[N], rhs[N];
	int starts[N + 1], indices[12 * N], i, j, r, c, nz = 0;
	int transpose, pass, code, accepted = 0, allocations, failure;
	struct simplex_KluSparse cache = {NULL};
	struct simplex_SparseVector packed;
	klu_symbolic *symbolic;
	klu_numeric *numeric;
	klu_common common;
	static const double block[9] = {0., 2., .25, 3., 1., -.5, .125, 4., 5.};
	memset(dense, 0, sizeof(dense));
	for (i = 0; i < N; i++)
		for (j = 0; j < N; j++) {
			double value = i / 3 == j / 3
					       ? block[3 * (i % 3) + j % 3]
					       : 0.;
			if (singleton && i / 3 == j / 3 && (i / 3) % 5 == 0)
				value = i == j ? 2. + i % 3 : 0.;
			if (i / 3 + 1 == j / 3 && (chain || (i / 3) % 4 != 3) &&
				i % 3 == j % 3)
				value = .0625;
			r = (5 * i + 7) % N;
			c = (7 * j + 11) % N;
			dense[r * N + c] = ldexp(value, (i % 13) - (j % 11));
		}
	for (j = 0; j < N; j++) {
		starts[j] = nz;
		for (i = 0; i < N; i++)
			if (dense[i * N + j] != 0.) {
				CHECK(nz < 12 * N);
				indices[nz] = i;
				values[nz++] = dense[i * N + j];
			}
	}
	starts[N] = nz;
	klu_defaults(&common);
	common.scale = scaling;
	common.btf = btf;
	symbolic = klu_analyze(N, starts, indices, &common);
	CHECK(symbolic);
	numeric = klu_factor(starts, indices, values, symbolic, &common);
	CHECK(numeric);
	requests = fail_at = 0;
	CHECK(!cache_build(&cache, symbolic, numeric, &common));
	allocations = requests;
	CHECK(!simplex_sparse_vector_create(&packed, N));
	memset(rhs, 0, sizeof(rhs));
	rhs[5] = NAN;
	memcpy(actual, rhs, sizeof(rhs));
	CHECK(cache_solve(&cache, actual, 0) == -1);
	CHECK(!memcmp(actual, rhs, sizeof(rhs)));
	for (transpose = 0; transpose < 2; transpose++)
		for (pass = 0; pass < 48; pass++) {
			double scale = pass % 3 == 2 ? 1e-200 : 1.;
			memset(rhs, 0, sizeof(rhs));
			if (pass % 4 == 0)
				for (i = 0; i < N; i++)
					rhs[i] = scale * ((i % 7) - 3);
			else
				rhs[(pass * 17 + 13) % N] =
					pass % 2 ? scale : -scale;
			if (pass % 11 == 10)
				memset(rhs, 0, sizeof(rhs));
			memcpy(expected, rhs, sizeof(rhs));
			memcpy(actual, rhs, sizeof(rhs));
			CHECK(transpose ? klu_tsolve(symbolic, numeric, N, 1,
						  expected, &common)
					: klu_solve(symbolic, numeric, N, 1,
						  expected, &common));
			code = cache_solve(&cache, actual, transpose);
			CHECK(code >= 0);
			CHECK(!simplex_sparse_vector_pack(&packed, rhs, 0.));
			CHECK(cache_solve_packed(&cache, &packed, transpose) ==
				code);
			for (i = 0; i < N; i++)
				CHECK(packed.dense[i] == actual[i]);
			if (code == 1) {
				CHECK(!memcmp(actual, rhs, sizeof(rhs)));
				continue;
			}
			accepted++;
			for (i = 0; i < N; i++) {
				long double residual = rhs[i],
					    magnitude =
						    fabsl((long double)rhs[i]);
				CHECK(isfinite(actual[i]));
				CHECK(fabs((actual[i] - expected[i]) / scale) <=
					1e-10 * (1. + fabs(expected[i] /
							      scale)));
				for (j = 0; j < N; j++) {
					double a = transpose ? dense[j * N + i]
							     : dense[i * N + j];
					long double term =
						(long double)a * actual[j];
					residual -= term;
					/* Match row-scaled auditing, with the test's tiny RHS
					 * scale instead of a unit floor that would mask underflow. */
					magnitude +=
						fabs(a) *
						fmax(scale, fabs(actual[j]));
				}
				CHECK(fabsl(residual) <= 1e-10L * magnitude);
			}
		}
	CHECK(accepted > 0);
	simplex_sparse_vector_destroy(&packed);
	cache_destroy(&cache);
	CHECK(!active);
	for (failure = 1; failure <= allocations; failure++) {
		requests = 0;
		fail_at = failure;
		CHECK(cache_build(&cache, symbolic, numeric, &common) == 1);
		CHECK(!cache.impl && !active);
		memset(actual, 0, sizeof(actual));
		actual[7] = 1.;
		CHECK(klu_solve(symbolic, numeric, N, 1, actual, &common));
	}
	fail_at = 0;
	klu_free_numeric(&numeric, &common);
	klu_free_symbolic(&symbolic, &common);
	printf("KLU sparse scale=%d BTF=%d chain=%d accepted=%d allocations=%d\n",
		scaling, btf, chain, accepted, allocations);
	return 0;
}

static int check_rank_queue(void)
{
	const int sizes[] = {64, 65, 1024, 1025, 32768, 32769, 1048577};
	int test, transpose, i;
	for (test = 0; test < 7; test++) {
		struct simplex_KluSparseImpl c = {0};
		c.nodes = sizes[test];
		c.heap = malloc((size_t)c.nodes * sizeof(int));
		c.queued = calloc((size_t)c.nodes, sizeof(int));
		CHECK(c.heap && c.queued && klu_queue_initialize(&c));
		for (transpose = 0; transpose < 2; transpose++) {
			int count = 0, expected = transpose ? c.nodes - 1 : 0;
			for (i = 0; i < 4096; i++) {
				int node = (i * 7919 + 31) % c.nodes;
				if (!c.queued[node])
					klu_heap_push(
						&c, node, &count, transpose);
			}
			while (count) {
				while (!c.queued[expected])
					expected += transpose ? -1 : 1;
				CHECK(klu_heap_pop(&c, &count, transpose) ==
					expected);
			}
			for (i = 0; i < c.nodes; i++)
				CHECK(!c.queued[i]);
		}
		free(c.heap);
		free(c.queued);
	}
	return 0;
}

int main(void)
{
	int scaling, btf, chain, singleton;
	if (check_rank_queue())
		return 1;
	for (scaling = 0; scaling <= 2; scaling++)
		for (btf = 0; btf <= 1; btf++)
			for (chain = 0; chain <= 1; chain++)
				for (singleton = 0; singleton < 2; singleton++)
					if (check_case(scaling, btf, chain,
						    singleton))
						return 1;
	return 0;
}
