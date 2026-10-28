/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_pan_face_recovery.c
 * @brief Audited native-face correction and exact no-change request memory.
 */

#include "simplex_pan_face_recovery.h"
#include "simplex_pan_cone.h"
#include "simplex_storage.h"
#include "policy/pan.h"
#include "policy/resources.h"
#include "policy/controls.h"
#include "utils.h"
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

struct face_Column {
	int column;
	double price;
};

static int face_order(const void *first, const void *second)
{
	const struct face_Column *a = (const struct face_Column *)first;
	const struct face_Column *b = (const struct face_Column *)second;

	if (a->price > b->price)
		return -1;
	if (a->price < b->price)
		return 1;
	return a->column < b->column ? -1 : (a->column > b->column);
}

/* Use a dual prediction to identify a numerical active face. This permits a
 * bounded objective correction, but never relaxes the terminal KKT/gap test.
 * Aggregate absolute row changes protect small essential coordinates. */
static int face_recover(const struct simplex_PanBounded *b, double *x,
	const double *y, const struct lp_simplex_Options *o, double *saved,
	double *r, double *work)
{
	struct face_Column *order;
	long double before, cost = 0.;
	double cost_budget;
	size_t bytes;
	int i, j, k, proposal, count = 0, changed = 0;

	if (b->a.columns < 0 || !simplex_storage_bytes((size_t)b->a.columns + 1,
					sizeof(*order), (size_t)-1, &bytes))
		return -1;
	order = (struct face_Column *)lp_simplex_malloc(bytes);
	if (!order)
		return -1;
	before = simplex_pan_cone_cost(&b->a, x);
	cost_budget = o->dual_tolerance * SIMPLEX_PAN_FACE_COST_FRACTION *
		      (1. + fabsl(before));
	for (j = 0; j < b->a.columns; j++)
		if ((x[j] > b->lower[j] && x[j] < b->upper[j])) {
			double g = b->a.objective[j] +
				   simplex_pan_column_dot(&b->a, j, y);

			if (!isfinite(g) ||
				fabs(g) <=
					o->dual_tolerance *
						(1. + fabs(b->a.objective[j])))
				continue;
			order[count].column = j;
			order[count++].price =
				fabs(g) / fmax(b->a.column_norm[j],
						  SIMPLEX_PAN_FACE_NORM_FLOOR);
		}
	qsort(order, (size_t)count, sizeof(*order), face_order);
	lp_simplex_memcpy(saved, x, (size_t)b->a.columns * sizeof(double));
	/* A relative-activity batch can cancel row changes that a strictly
	 * coordinatewise absolute budget must refuse. It has no acceptance
	 * authority: the complete compensated absolute audit below is unchanged.
	 * Refusal retries the original conservative proposal from saved. */
	for (proposal = 0; proposal < 2; proposal++) {
		changed = 0;
		cost = 0.;
		(void)simplex_pan_bounded_primal(b, x, r, work);
		if (proposal == 0) {
			for (i = 0; i < b->a.rows; i++)
				r[i] = fabs(b->a.rhs[i]);
			for (j = 0; j < b->a.columns; j++)
				for (k = b->a.column_start[j];
					k < b->a.column_start[j + 1]; k++)
					r[b->a.row_index[k]] +=
						fabs(b->a.value[k]) *
						fmax(1., fabs(x[j]));
		}
		for (i = 0; i < b->a.rows; i++)
			work[i] = 0.;
		for (i = 0; i < count; i++) {
			double g, target, delta;
			j = order[i].column;
			g = b->a.objective[j] +
			    simplex_pan_column_dot(&b->a, j, y);
			target = g > 0. ? b->lower[j] : b->upper[j];
			if (!isfinite(target))
				continue;
			delta = target - x[j];
			if (!isfinite(delta) ||
				fabsl(cost + (long double)b->a.objective[j] *
						     delta) > cost_budget)
				continue;
			for (k = b->a.column_start[j];
				k < b->a.column_start[j + 1]; k++) {
				int row = b->a.row_index[k];
				double allowance =
					proposal == 0
						? o->primal_tolerance *
							  SIMPLEX_PAN_FACE_COST_FRACTION *
							  r[row]
						: o->primal_tolerance *
								  SIMPLEX_PAN_FEASIBLE_MARGIN -
							  fabs(r[row]);
				if (!isfinite(allowance) ||
					work[row] + fabs(b->a.value[k] *
							    delta) >
						allowance)
					break;
			}
			if (k < b->a.column_start[j + 1])
				continue;
			for (k = b->a.column_start[j];
				k < b->a.column_start[j + 1]; k++)
				work[b->a.row_index[k]] +=
					fabs(b->a.value[k] * delta);
			cost += (long double)b->a.objective[j] * delta;
			x[j] = target;
			changed++;
		}
		if (!changed)
			continue;
		if (simplex_pan_bounded_primal(b, x, r, work) <=
				o->primal_tolerance *
					SIMPLEX_PAN_FEASIBLE_MARGIN &&
			simplex_pan_cone_cost(&b->a, x) <= before + cost_budget)
			break;
		lp_simplex_memcpy(
			x, saved, (size_t)b->a.columns * sizeof(double));
		changed = 0;
	}

	if (changed && getenv(SIMPLEX_ENV_PROFILE))
		fprintf(stderr,
			"pan-bda native face recovery: bounds=%d cost-change=%.6g\n",
			changed, (double)cost);
	lp_simplex_free(order);
	return changed;
}

static int face_same_model(const struct simplex_PanFaceRecovery *c,
	const struct simplex_PanBounded *b)
{
	const struct simplex_PanStandard *a = &b->a, *p = &c->representation;
	return c->model == b && b->problem &&
	       simplex_problem_identity_matches(&c->input, b->problem) &&
	       c->lower == b->lower && c->upper == b->upper &&
	       p->rows == a->rows && p->columns == a->columns &&
	       p->nonzeros == a->nonzeros &&
	       p->service_revision == a->service_revision &&
	       p->service_model == a->service_model &&
	       p->services == a->services && p->value == a->value &&
	       p->column_start == a->column_start &&
	       p->row_index == a->row_index && p->objective == a->objective &&
	       p->rhs == a->rhs && p->column_norm == a->column_norm;
}

static int face_snapshot_bytes(int rows, int columns, size_t limit,
	size_t *columns_bytes, size_t *rows_bytes)
{
	size_t xb, yb;
	*columns_bytes = *rows_bytes = 0;
	if (rows < 0 || columns < 0 ||
		!simplex_storage_bytes(
			(size_t)columns, sizeof(double), limit, &xb) ||
		!simplex_storage_bytes(
			(size_t)rows, sizeof(double), limit - xb, &yb) ||
		(!xb && !yb))
		return 0;
	*columns_bytes = xb;
	*rows_bytes = yb;
	return 1;
}

static void face_remember(struct simplex_PanFaceRecovery *c,
	const struct simplex_PanBounded *b, const double *x, const double *y,
	const struct lp_simplex_Options *o)
{
	size_t xb, yb, total;
	int j;
	if (!b->problem || b->a.columns < 0 || b->a.rows < 0 ||
		!isfinite(o->primal_tolerance) || o->primal_tolerance <= 0. ||
		!isfinite(o->dual_tolerance) || o->dual_tolerance <= 0.)
		goto unknown;
	if (!face_snapshot_bytes(b->a.rows, b->a.columns,
		    SIMPLEX_PAN_FACE_RETRY_BYTES, &xb, &yb))
		goto unknown;
	c->input = simplex_problem_identity(b->problem);
	if (!simplex_problem_identity_matches(&c->input, b->problem))
		goto unknown;
	for (j = 0; j < b->a.columns; j++)
		if (!isfinite(x[j]))
			goto unknown;
	for (j = 0; j < b->a.rows; j++)
		if (!isfinite(y[j]))
			goto unknown;
	total = xb + yb;
	if (!total)
		goto unknown;
	if (c->columns_bytes != xb || c->rows_bytes != yb) {
		lp_simplex_free(c->storage);
		c->storage = NULL;
		c->columns_bytes = xb;
		c->rows_bytes = yb;
	}
	if (!c->storage)
		c->storage = (double *)lp_simplex_malloc(total);
	if (!c->storage)
		goto unknown;
	if (xb)
		lp_simplex_memcpy(c->storage, x, xb);
	if (yb)
		lp_simplex_memcpy(c->storage + b->a.columns, y, yb);
	c->model = b;
	c->representation = b->a;
	c->lower = b->lower;
	c->upper = b->upper;
	c->primal_tolerance = o->primal_tolerance;
	c->dual_tolerance = o->dual_tolerance;
	c->known = 1;
	return;
unknown:
	if (c->unremembered < LONG_MAX)
		c->unremembered++;
}

int simplex_pan_face_recover_reported(struct simplex_PanFaceRecovery *c,
	const struct simplex_PanBounded *b, double *x, const double *y,
	const struct lp_simplex_Options *o, double *saved, double *r,
	double *work, struct simplex_PanFaceResult *report)
{
	int result;
	if (report)
		lp_simplex_memset(report, 0, sizeof(*report));
	if (c && c->known && face_same_model(c, b) &&
		c->primal_tolerance == o->primal_tolerance &&
		c->dual_tolerance == o->dual_tolerance &&
		(!c->columns_bytes ||
			!lp_simplex_memcmp(c->storage, x, c->columns_bytes)) &&
		(!c->rows_bytes || !lp_simplex_memcmp(c->storage + b->a.columns,
					   y, c->rows_bytes))) {
		if (c->skipped < LONG_MAX)
			c->skipped++;
		if (report) {
			report->action = SIMPLEX_PAN_FACE_EXACT_SKIP;
			report->remembered = 1;
		}
		return 0;
	}
	if (c) {
		c->known = 0;
		if (c->attempts < LONG_MAX)
			c->attempts++;
	}
	result = face_recover(b, x, y, o, saved, r, work);
	if (c && result == 0)
		face_remember(c, b, x, y, o);
	if (report) {
		report->action = result < 0   ? SIMPLEX_PAN_FACE_RESOURCE
				 : result > 0 ? SIMPLEX_PAN_FACE_ADJUSTED
					      : SIMPLEX_PAN_FACE_NO_CHANGE;
		report->changed = result > 0 ? result : 0;
		report->remembered = c && c->known;
	}
	return result;
}

int simplex_pan_face_recover(struct simplex_PanFaceRecovery *c,
	const struct simplex_PanBounded *b, double *x, const double *y,
	const struct lp_simplex_Options *o, double *saved, double *r,
	double *work)
{
	return simplex_pan_face_recover_reported(
		c, b, x, y, o, saved, r, work, NULL);
}

void simplex_pan_face_recovery_destroy(struct simplex_PanFaceRecovery *c)
{
	if (!c)
		return;
	if (getenv(SIMPLEX_ENV_PROFILE))
		fprintf(stderr,
			"pan-bda face recovery: attempts=%ld unchanged-skips=%ld unremembered=%ld\n",
			c->attempts, c->skipped, c->unremembered);
	lp_simplex_free(c->storage);
	lp_simplex_memset(c, 0, sizeof(*c));
}
