/**
 * @file test_pan_qr_storage.c
 * @brief Production dynamic R implementation, bounded allocation fault
 * injection.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>

static void *live[64];
static int requests, fail_at, active, bad_free;

void *qr_storage_malloc(size_t bytes)
{
	int i;
	void *p;
	requests++;
	if (requests == fail_at || bytes > 1024 * 1024)
		return NULL;
	p = malloc(bytes ? bytes : 1);
	if (!p)
		return NULL;
	for (i = 0; i < 64; i++)
		if (!live[i]) {
			live[i] = p;
			active++;
			return p;
		}
	free(p);
	return NULL;
}

void qr_storage_free(void *p)
{
	int i;
	if (!p)
		return;
	for (i = 0; i < 64; i++)
		if (live[i] == p) {
			live[i] = NULL;
			active--;
			free(p);
			return;
		}
	bad_free++;
}

#define lp_simplex_malloc qr_storage_malloc
#define lp_simplex_free qr_storage_free
#define pan_dynamic_sparse_qr_initialize qr_storage_initialize
#define pan_dynamic_sparse_qr_clear qr_storage_clear
#define pan_dynamic_sparse_qr_destroy qr_storage_destroy
#define pan_dynamic_sparse_qr_append qr_storage_append
#define pan_dynamic_sparse_qr_load_csc qr_storage_load_csc
#define pan_dynamic_sparse_qr_delete qr_storage_delete
#define pan_dynamic_sparse_qr_move_last qr_storage_move_last
#define pan_dynamic_sparse_qr_solve_upper_transpose \
	qr_storage_solve_upper_transpose
#define pan_dynamic_sparse_qr_solve_upper qr_storage_solve_upper
#define pan_dynamic_sparse_qr_solve qr_storage_solve
#define pan_dynamic_sparse_qr_validate qr_storage_validate
#ifndef QR_STORAGE_SOURCE
#define QR_STORAGE_SOURCE "../src/bda/linalg/pan_dynamic_sparse_qr.c"
#endif
#include QR_STORAGE_SOURCE

#undef lp_simplex_malloc
#undef lp_simplex_free
#define CHECK(c) \
	do { \
		if (!(c)) { \
			fprintf(stderr, "QR storage line %d: %s\n", __LINE__, \
				#c); \
			return 1; \
		} \
	} while (0)

static int seed(struct pan_DynamicSparseQR *f, int count)
{
	double border[16];
	int i;
	requests = fail_at = 0;
	CHECK(!qr_storage_initialize(f));
	for (i = 0; i < 16; i++)
		border[i] = .125;
	for (i = 0; i < count; i++)
		CHECK(!qr_storage_append(f, i, border, 1., 1e-10));
	return 0;
}

static int solves(struct pan_DynamicSparseQR *f, int count)
{
	double right[16];
	int i, j;
	CHECK(f->count == count);
	CHECK(!qr_storage_validate(f));
	for (i = 0; i < count; i++) {
		CHECK(f->tag[i] == i);
		for (j = 0; j < count; j++)
			CHECK(pan_get(f, i, j) ==
				(j == i ? 1. : (j > i ? .125 : 0.)));
		right[i] = 1. + .125 * (count - 1 - i);
	}
	CHECK(!qr_storage_solve_upper(f, right));
	for (i = 0; i < count; i++)
		CHECK(fabs(right[i] - 1.) < 1e-12);
	return 0;
}

static int allocation_failures(void)
{
	struct pan_DynamicSparseQR f, before;
	double border[16];
	int failure, i, existing;
	for (i = 0; i < 16; i++)
		border[i] = .125;
	/* The 14 buffers are staged: failure preserves an existing R and cache. */
	for (failure = 1; failure <= 14; failure++) {
		CHECK(!seed(&f, 8) && !solves(&f, 8));
		before = f;
		existing = active;
		requests = 0;
		fail_at = failure;
		CHECK(qr_storage_append(&f, 8, border, 1., 1e-10));
		CHECK(active == existing && !memcmp(&before, &f, sizeof(f)));
		fail_at = 0;
		CHECK(!solves(&f, 8));
		CHECK(!qr_storage_append(&f, 8, border, 1., 1e-10));
		CHECK(!solves(&f, 9));
		qr_storage_destroy(&f);
		CHECK(!active && !bad_free);
	}
	/* After four inserts, fail either the next border node or the diagonal. */
	for (failure = 0; failure < 2; failure++) {
		CHECK(!seed(&f, 7) && !solves(&f, 7));
		CHECK(f.node_count == 28 && f.node_capacity == 32);
		for (i = 0; i < 7; i++)
			border[i] = failure && i >= 4 ? 0. : .125;
		requests = 0;
		fail_at = 1;
		CHECK(qr_storage_append(&f, 7, border, 1., 1e-10));
		CHECK(requests == 1);
		fail_at = 0;
		CHECK(!solves(&f, 7));
		CHECK(f.column_head[7] == -1 && f.column_tail[7] == -1 &&
			!f.column_count[7]);
		for (i = 0; i < 7; i++)
			border[i] = .125;
		CHECK(!qr_storage_append(&f, 7, border, 1., 1e-10));
		CHECK(!solves(&f, 8));
		qr_storage_destroy(&f);
		CHECK(!active && !bad_free);
	}
	/* Both packed-buffer allocations must preserve a usable old view. */
	for (failure = 1; failure <= 2; failure++) {
		CHECK(!seed(&f, 8) && !solves(&f, 8));
		before = f;
		existing = active;
		requests = 0;
		fail_at = failure;
		CHECK(pan_packed_reserve(&f, f.packed_capacity + 1));
		CHECK(active == existing && !memcmp(&before, &f, sizeof(f)));
		fail_at = 0;
		CHECK(!solves(&f, 8));
		qr_storage_destroy(&f);
		CHECK(!active && !bad_free);
	}
	return 0;
}

static int extreme_reserve(void)
{
	struct pan_DynamicSparseQR f;
	CHECK(!qr_storage_initialize(&f));
	requests = fail_at = 0;
	CHECK(pan_factor_reserve(&f, INT_MAX));
	qr_storage_destroy(&f);
	CHECK(!active && !bad_free);
	return 0;
}
#ifndef QR_STORAGE_LEGACY
static int arithmetic(void)
{
	struct pan_DynamicSparseQR f;
	struct pan_DynamicSparseRow rows[2];
	double border = 1.;
	int next;
	CHECK(pan_reserve_capacity(
		      8, 9, 8, sizeof(double), (size_t)-1, &next) &&
		next == 16);
	CHECK(pan_reserve_capacity(INT_MAX - 1, INT_MAX, 8, sizeof(double),
		      (size_t)-1, &next) == (sizeof(size_t) >= 8));
	CHECK(!pan_reserve_capacity(
		      0, INT_MAX, 8, 16, (size_t)4294967295UL, &next) &&
		!next);
	CHECK(!pan_reserve_capacity(0, -1, 8, 16, (size_t)-1, &next));
	CHECK(!qr_storage_initialize(&f));
	f.count = INT_MAX;
	requests = 0;
	CHECK(qr_storage_append(&f, 0, &border, 1., 1e-10) && !requests);
	f.count = 0;
	f.node_count = INT_MAX;
	CHECK(pan_node_allocate(&f) < 0 && !requests);
	f.node_count = 0;
	f.count = 2;
	f.row = rows;
	rows[0].count = INT_MAX;
	rows[1].count = 1;
	CHECK(pan_pack_numeric(&f) && !requests);
	return 0;
}
#endif
int main(int argc, char **argv)
{
	if (argc > 1 && !strcmp(argv[1], "overflow"))
		return extreme_reserve();
#ifndef QR_STORAGE_LEGACY
	if (arithmetic())
		return 1;
#endif
	return allocation_failures() || extreme_reserve();
}
