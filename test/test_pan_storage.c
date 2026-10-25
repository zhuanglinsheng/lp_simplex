/**
 * @file test_pan_storage.c
 * @brief Production equality-form body with bounded allocation fault injection.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include "simplex_pan_storage.h"

static void *live[16];
static int requests, fail_at, active, bad_free;

void *pan_test_malloc(size_t bytes)
{
	int i;
	void *p;
	requests++;
	if (requests == fail_at || bytes > 1024 * 1024)
		return NULL;
	p = malloc(bytes ? bytes : 1);
	if (!p)
		return NULL;
	for (i = 0; i < 16; i++)
		if (!live[i]) {
			live[i] = p;
			active++;
			return p;
		}
	free(p);
	return NULL;
}

void pan_test_free(void *p)
{
	int i;
	if (!p)
		return;
	for (i = 0; i < 16; i++)
		if (live[i] == p) {
			live[i] = NULL;
			active--;
			free(p);
			return;
		}
	bad_free++;
}

#define lp_simplex_malloc pan_test_malloc
#define lp_simplex_free pan_test_free
#define simplex_pan_standard_create pan_test_standard_create
#define simplex_pan_standard_destroy pan_test_standard_destroy
#define simplex_pan_standard_recover pan_test_standard_recover
#define simplex_pan_standard_encode pan_test_standard_encode
#define simplex_pan_column_dot pan_test_column_dot
#define simplex_pan_residual pan_test_residual
#define simplex_pan_original_form_smaller pan_test_original_form_smaller
#ifndef PAN_STANDARD_SOURCE
#define PAN_STANDARD_SOURCE "../src/bda/simplex_pan_standard.c"
#endif
#include PAN_STANDARD_SOURCE

#undef lp_simplex_malloc
#undef lp_simplex_free
#include "simplex_pan_bounded.h"

#define CHECK(c) \
	do { \
		if (!(c)) { \
			fprintf(stderr, "Pan storage line %d: %s\n", __LINE__, \
				#c); \
			return 1; \
		} \
	} while (0)

static int arithmetic(void)
{
	int count = INT_MAX - 1;
	size_t narrow = (size_t)4294967295UL;
	CHECK(pan_storage_add(&count, 1) && count == INT_MAX);
	CHECK(!pan_storage_add(&count, 1) && count == INT_MAX);
	CHECK(!pan_storage_add(&count, -1) && count == INT_MAX);
	CHECK(pan_storage_fits(2, 4, 3, 5, 0, narrow));
	CHECK(!pan_storage_fits(1, INT_MAX, 1, 0, 0, narrow));
	CHECK(!pan_storage_fits(1, 1, 1, INT_MAX, 1, narrow));
	CHECK(!pan_storage_fits(INT_MAX, 1, 1, 0, 1, narrow));
	CHECK(!pan_storage_fits(1, 1, -1, 0, 0, narrow));
	return 0;
}

static int conversion(void)
{
	struct simplex_Problem p;
	struct simplex_PanStandard a, zero;
	struct simplex_PanBounded b;
	struct optm_VariableBound bounds[3];
	int starts[] = {0, 1, 2, 3}, index[] = {0, 1, 0};
	double values[] = {1., 2., 3.}, rhs[] = {4., 2.};
	double objective[] = {2., 1., 3.};
	double x[] = {0., 1., 2.}, z[5], recovered[3], residual[3];
	unsigned char types[] = {optm_CONS_T_LE, optm_CONS_T_EQ};
	int failure;
	memset(&p, 0, sizeof(p));
	memset(bounds, 0, sizeof(bounds));
	memset(&zero, 0, sizeof(zero));
	p.rows = 2;
	p.columns = 3;
	p.matrix.column_start = starts;
	p.matrix.row_index = index;
	p.matrix.value = values;
	p.matrix.nonzeros = 3;
	p.bounds = bounds;
	p.rhs = rhs;
	p.objective = objective;
	p.row_type = types;
	bounds[0].b_type = optm_BOUND_T_FR;
	bounds[1].b_type = optm_BOUND_T_BS;
	bounds[1].lb = 0.;
	bounds[1].ub = 2.;
	bounds[2].b_type = optm_BOUND_T_BS;
	bounds[2].lb = bounds[2].ub = 2.;
	/* Ten allocations, including normalization and original-space maps. */
	for (failure = 1; failure <= 11; failure++) {
		requests = bad_free = 0;
		fail_at = failure;
		if (failure <= 10) {
			CHECK(pan_test_standard_create(&a, &p));
			CHECK(!memcmp(&a, &zero, sizeof(a)));
		} else {
			CHECK(!pan_test_standard_create(&a, &p));
			CHECK(a.rows == 3 && a.columns == 5 && a.nonzeros == 6);
			fail_at = 0;
			/* Fixed x[2] shifts row zero by 6; choose a feasible free x. */
			x[0] = -3.;
			CHECK(!pan_test_standard_encode(
				&a, x, z, residual, 1e-10));
			pan_test_standard_recover(&a, z, recovered);
			CHECK(recovered[0] == x[0] && recovered[1] == 1. &&
				recovered[2] == 2.);
		}
		pan_test_standard_destroy(&a);
		CHECK(!active && !bad_free);
	}
	/* One column with a huge support is enough to trigger free-column
	 * duplication overflow, without allocating or reading that support. */
	p.rows = p.columns = 1;
	p.matrix.nonzeros = starts[1] = INT_MAX;
	types[0] = optm_CONS_T_EQ;
	requests = 0;
	fail_at = 0;
	CHECK(pan_test_standard_create(&a, &p));
	CHECK(!requests && !active && !memcmp(&a, &zero, sizeof(a)));
	/* Direct allocation also checks the column-start sentinel. On a wide
	 * host the bounded allocator rejects the request; on a narrow host the
	 * byte preflight must reject it before any allocation. */
	CHECK(pan_standard_allocate(&a, 1, INT_MAX, 1, 0));
	pan_test_standard_destroy(&a);
	CHECK(!active && !bad_free);
	/* An inequality slack alone would overflow the native nnz count. */
	types[0] = optm_CONS_T_LE;
	CHECK(simplex_pan_bounded_create(&b, &p));
	simplex_pan_bounded_destroy(&b);
	CHECK(!active && !bad_free);
	return 0;
}

int main(void)
{
	return arithmetic() || conversion();
}
