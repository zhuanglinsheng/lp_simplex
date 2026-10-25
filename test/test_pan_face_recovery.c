/**
 * @file test_pan_face_recovery.c
 * @brief Exact native-face requests, optional storage and no-op suppression.
 */

#include <lp_simplex/lp_simplex.h>
#include "utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

static void *live[8];
static int active, requests, fail_at, bad_free;

void *face_test_malloc(size_t bytes)
{
	void *p;
	int i;
	requests++;
	if (requests == fail_at || bytes > 1024 * 1024)
		return NULL;
	p = malloc(bytes ? bytes : 1);
	if (!p)
		return NULL;
	memset(p, 0xa5, bytes);
	for (i = 0; i < 8; i++)
		if (!live[i]) {
			live[i] = p;
			active++;
			return p;
		}
	free(p);
	return NULL;
}

void face_test_free(void *p)
{
	int i;
	if (!p)
		return;
	for (i = 0; i < 8; i++)
		if (live[i] == p) {
			live[i] = NULL;
			active--;
			free(p);
			return;
		}
	bad_free++;
}

#define lp_simplex_malloc face_test_malloc
#define lp_simplex_free face_test_free
#define simplex_pan_face_recover tested_face_recover
#define simplex_pan_face_recover_reported tested_face_recover_reported
#define simplex_pan_face_recovery_destroy tested_face_destroy
#include "../src/bda/simplex_pan_face_recovery.c"
#undef lp_simplex_malloc
#undef lp_simplex_free

#define CHECK(c) \
	do { \
		if (!(c)) { \
			fprintf(stderr, "face recovery line %d: %s\n", \
				__LINE__, #c); \
			return 1; \
		} \
	} while (0)

static int cases(void)
{
	struct lp_Model *m = lp_model_create(1, 2);
	struct simplex_Problem p;
	struct simplex_PanBounded b;
	struct simplex_PanFaceRecovery c;
	struct simplex_PanFaceResult face_result;
	struct lp_simplex_Options o;
	double x[3] = {.5, .5, 0.}, y[1] = {0.}, saved[3], r[1], work[1];
	int f, before;
	size_t xb, yb;
	CHECK(m);
	m->constraints[0].type = optm_CONS_T_GE;
	m->constraints[0].rhs = 1.;
	m->constraints[0].coef[0] = m->constraints[0].coef[1] = 1.;
	m->objective[1] = 1.;
	CHECK(!simplex_problem_from_model(&p, m));
	CHECK(!simplex_pan_bounded_create(&b, &p));
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_PAN_BDA);
	CHECK(face_snapshot_bytes(1, 3, 32, &xb, &yb) && xb == 24 && yb == 8);
	CHECK(!face_snapshot_bytes(1, 3, 31, &xb, &yb) && !xb && !yb);
	CHECK(!face_snapshot_bytes(-1, 3, 32, &xb, &yb));
	CHECK(!face_snapshot_bytes(1, INT_MAX, (size_t)4294967295UL, &xb, &yb));
	for (f = 1; f <= 2; f++) {
		memset(&c, 0, sizeof(c));
		requests = 0;
		fail_at = f;
		CHECK(tested_face_recover_reported(&c, &b, x, y, &o, saved, r,
			      work, &face_result) == (f == 1 ? -1 : 0));
		CHECK(face_result.action ==
			(f == 1 ? SIMPLEX_PAN_FACE_RESOURCE
				: SIMPLEX_PAN_FACE_NO_CHANGE));
		CHECK(!face_result.changed && !face_result.remembered);
		CHECK(!c.known && x[0] == .5 && x[1] == .5 && !y[0]);
		CHECK(!active && !bad_free);
		requests = fail_at = 0;
		CHECK(!tested_face_recover_reported(
			&c, &b, x, y, &o, saved, r, work, &face_result));
		CHECK(face_result.action == SIMPLEX_PAN_FACE_NO_CHANGE &&
			face_result.remembered);
		CHECK(c.known && requests == 2);
		before = requests;
		saved[0] = r[0] = work[0] = 123.;
		CHECK(!tested_face_recover_reported(
			&c, &b, x, y, &o, saved, r, work, &face_result));
		CHECK(face_result.action == SIMPLEX_PAN_FACE_EXACT_SKIP &&
			face_result.remembered && !face_result.changed);
		CHECK(c.skipped == 1 && requests == before);
		CHECK(saved[0] == 123. && r[0] == 123. && work[0] == 123.);
		tested_face_destroy(&c);
		tested_face_destroy(&c);
		CHECK(!active && !bad_free);
	}
	memset(&c, 0, sizeof(c));
	requests = fail_at = 0;
	CHECK(!tested_face_recover(&c, &b, x, y, &o, saved, r, work));
	CHECK(c.attempts == 1);
	/* Diagnostic counters do not authorize or invalidate the retained request. */
	c.skipped = LONG_MAX;
	before = requests;
	CHECK(!tested_face_recover(&c, &b, x, y, &o, saved, r, work));
	CHECK(c.skipped == LONG_MAX && requests == before);
	c.skipped = 0;
	x[0] = nextafter(.5, 1.);
	CHECK(!tested_face_recover(&c, &b, x, y, &o, saved, r, work));
	CHECK(c.attempts == 2 && !c.skipped);
	y[0] = -0.;
	CHECK(!tested_face_recover(&c, &b, x, y, &o, saved, r, work));
	CHECK(c.attempts == 3);
	o.dual_tolerance *= 2.;
	CHECK(!tested_face_recover(&c, &b, x, y, &o, saved, r, work));
	CHECK(c.attempts == 4);
	o.primal_tolerance *= 2.;
	CHECK(!tested_face_recover(&c, &b, x, y, &o, saved, r, work));
	CHECK(c.attempts == 5);
	simplex_problem_mark_changed(&p);
	CHECK(!tested_face_recover(&c, &b, x, y, &o, saved, r, work));
	CHECK(c.attempts == 6);
	b.a.service_revision++;
	CHECK(!tested_face_recover(&c, &b, x, y, &o, saved, r, work));
	CHECK(c.attempts == 7);
	x[0] = 1.;
	x[1] = 1e-9;
	CHECK(tested_face_recover(&c, &b, x, y, &o, saved, r, work) == 1);
	CHECK(!c.known && x[1] == 0. && x[0] == 1.);
	CHECK(simplex_pan_bounded_primal(&b, x, r, work) == 0.);
	CHECK(!y[0] && signbit(y[0]));
	tested_face_destroy(&c);
	CHECK(!active && !bad_free);
	/* Without an owner the same audited action remains available. */
	x[0] = x[1] = .5;
	requests = 0;
	CHECK(!tested_face_recover(NULL, &b, x, y, &o, saved, r, work));
	CHECK(requests == 1 && !active);
	simplex_pan_bounded_destroy(&b);
	simplex_problem_destroy(&p);
	lp_model_free(m);
	return 0;
}

int main(void)
{
	return cases();
}
