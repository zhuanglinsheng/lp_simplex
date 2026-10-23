/* Copyright (C) 2022 Zhuang Linsheng. License: LGPL 3.0. */
#include "simplex_progress.h"
#include "utils.h"
#include <math.h>
#include <float.h>
#include <string.h>

/* Outward long-double intervals also work where long double == double.
 * Exact zero operations are kept exact, notably for unbounded/free columns. */
static long double down_add(long double a, long double b)
{
	long double sum, v, error;
	if (b == 0.)
		return a;
	if (a == 0.)
		return b;
	sum = a + b;
	if (!isfinite(sum))
		return -HUGE_VALL;
	v = sum - a;
	error = (a - (sum - v)) + (b - v);
	return error >= 0. ? sum : nextafterl(sum, -HUGE_VALL);
}
static long double down_product(long double a, long double b)
{
	long double product;
	if (a == 0. || b == 0.)
		return 0.;
	product = a * b;
	if (!isfinite(product))
		return -HUGE_VALL;
	/* Close to underflow, even the FMA residual can round to zero. */
	if (fabsl(product) >= 4.L * (LDBL_MIN / LDBL_EPSILON) &&
		fmal(a, b, -product) >= 0.)
		return product;
	return nextafterl(product, -HUGE_VALL);
}
static long double multiplier(
	const struct simplex_Problem *p, const double *y, double scale, int i)
{
	long double v = y ? (long double)y[i] * scale : 0.;
	if (p->row_type[i] == optm_CONS_T_LE && v < 0.)
		v = 0.;
	if (p->row_type[i] == optm_CONS_T_GE && v > 0.)
		v = 0.;
	return v;
}

double simplex_progress_dual_bound(
	const struct simplex_Problem *p, const double *y, double scale)
{
	long double total = 0.;
	int i, j, k;
	if (!isfinite(scale))
		return -HUGE_VAL;
	for (i = 0; i < p->rows; i++) {
		long double v = multiplier(p, y, scale, i);
		if (!isfinite(v))
			return -HUGE_VAL;
		total = down_add(total, down_product(-p->rhs[i], v));
	}
	for (j = 0; j < p->columns; j++) {
		const struct optm_VariableBound *b = p->bounds + j;
		long double lo = p->objective[j], hi = lo, term;
		double lower = (b->b_type == optm_BOUND_T_LO ||
				       b->b_type == optm_BOUND_T_BS)
				       ? b->lb
				       : -HUGE_VAL;
		double upper = (b->b_type == optm_BOUND_T_UP ||
				       b->b_type == optm_BOUND_T_BS)
				       ? b->ub
				       : HUGE_VAL;
		for (k = p->matrix.column_start[j];
			k < p->matrix.column_start[j + 1]; k++) {
			long double v =
				multiplier(p, y, scale, p->matrix.row_index[k]);
			long double a = p->matrix.value[k];
			lo = down_add(lo, down_product(a, v));
			hi = -down_add(-hi, down_product(-a, v));
		}
		if (!isfinite(lo) || !isfinite(hi))
			return -HUGE_VAL;
		/* Minimise g*x over both the box and the enclosing g interval.
		 * Crossing zero with an infinite endpoint cannot yield a bound. */
		if ((!isfinite(lower) && hi > 0.) ||
			(!isfinite(upper) && lo < 0.))
			return -HUGE_VAL;
		term = 0.;
		if (isfinite(lower)) {
			term = down_product(lower, lower >= 0. ? lo : hi);
			if (!isfinite(upper) && lo <= 0. && hi >= 0.)
				term = fminl(term, 0.);
		}
		if (isfinite(upper)) {
			long double other =
				down_product(upper, upper >= 0. ? lo : hi);
			term = isfinite(lower) ? fminl(term, other) : other;
			if (!isfinite(lower) && lo <= 0. && hi >= 0.)
				term = fminl(term, 0.);
		}
		total = down_add(total, term);
	}
	if (!isfinite(total) || !isfinite((double)total))
		return -HUGE_VAL;
	{
		double rounded = (double)total;
		return (long double)rounded > total
			       ? nextafter(rounded, -HUGE_VAL)
			       : rounded;
	}
}

static int tracker_due(void *context)
{
	struct simplex_ProgressTracker *t = context;
	return (double)(clock() - t->last) / CLOCKS_PER_SEC >=
	       t->options->progress_interval;
}
static void tracker_report(void *context, const double *x, const double *y,
	double scale, const char *stage, int iterations, int force)
{
	struct simplex_ProgressTracker *t = context;
	struct lp_simplex_Progress event;
	long double objective = t->model->objective_offset;
	double primal = HUGE_VAL, bound, candidate;
	int j, improved = 0;
	(void)force;
	if (x) {
		primal = simplex_problem_primal_infeasibility(&t->original, x);
		if (isfinite(primal) &&
			primal <= t->options->primal_tolerance) {
			for (j = 0; j < t->original.columns; j++)
				objective = -down_add(-objective,
					down_product(-t->original.objective[j],
						x[j]));
			candidate = (double)objective;
			if ((long double)candidate < objective)
				candidate = nextafter(candidate, HUGE_VAL);
			if (isfinite(objective) && isfinite(candidate) &&
				(!t->have_primal ||
					candidate < t->best_primal)) {
				t->best_primal = candidate;
				t->primal_error = primal;
				t->have_primal = improved = 1;
				memcpy(t->point, x,
					(size_t)t->model->n * sizeof(double));
			}
		}
	}
	bound = simplex_progress_dual_bound(&t->original, y, scale);
	if (isfinite(bound)) {
		long double v = down_add(bound, t->model->objective_offset);
		bound = (double)v;
		if ((long double)bound > v)
			bound = nextafter(bound, -HUGE_VAL);
		if (isfinite(bound) && bound > t->best_dual)
			t->best_dual = bound;
	}
	event.stage = stage;
	event.stage_iterations = iterations;
	event.cpu_seconds = (double)(clock() - t->begun) / CLOCKS_PER_SEC;
	if (t->model->objective_sense == 1.) {
		event.lower_bound = t->best_dual;
		event.upper_bound = t->best_primal;
	} else {
		event.lower_bound = -t->best_primal;
		event.upper_bound = -t->best_dual;
	}
	event.absolute_gap = event.relative_gap = HUGE_VAL;
	/* A crossed interval is unavailable rather than a false zero gap.
	 * Finite-tolerance primal feasibility can cause a tiny crossing. */
	if (isfinite(event.lower_bound) && isfinite(event.upper_bound) &&
		event.lower_bound <= event.upper_bound) {
		event.absolute_gap = event.upper_bound - event.lower_bound;
		event.relative_gap = event.absolute_gap /
				     fmax(1., fmax(fabs(event.lower_bound),
						      fabs(event.upper_bound)));
	}
	event.point = improved ? t->point : NULL;
	event.primal_infeasibility = t->primal_error;
	t->last = clock();
	t->options->progress_callback(t->options->progress_context, &event);
}
int simplex_progress_init(struct simplex_ProgressTracker *t,
	const struct lp_Model *model, const struct lp_simplex_Options *o)
{
	memset(t, 0, sizeof(*t));
	if (!o->progress_callback)
		return 0;
	t->model = model;
	t->options = o;
	t->best_primal = t->primal_error = HUGE_VAL;
	t->best_dual = -HUGE_VAL;
	t->begun = t->last = clock();
	if (simplex_problem_from_model(&t->original, model))
		return -1;
	t->storage = lp_simplex_malloc((size_t)model->n * sizeof(double));
	if (!t->storage) {
		simplex_problem_destroy(&t->original);
		return -1;
	}
	t->point = t->storage;
	t->sink.context = t;
	t->sink.due = tracker_due;
	t->sink.report = tracker_report;
	tracker_report(t, NULL, NULL, 1., "start", 0, 1);
	return 0;
}
void simplex_progress_destroy(struct simplex_ProgressTracker *t)
{
	if (!t->sink.report)
		return;
	simplex_problem_destroy(&t->original);
	lp_simplex_free(t->storage);
}
static int bridge_due(void *context)
{
	struct simplex_ProgressBridge *b = context;
	return b->parent->due(b->parent->context);
}
static void bridge_report(void *context, const double *x, const double *y,
	double scale, const char *stage, int iterations, int force)
{
	struct simplex_ProgressBridge *b = context;
	int i;
	if (b->presolve) {
		if (x)
			simplex_presolve_postsolve(b->presolve, x, b->x);
		memset(b->y, 0, (size_t)b->rows * sizeof(double));
		if (y)
			for (i = 0; i < b->presolve->reduced->rows; i++)
				b->y[b->presolve->row_map[i]] = y[i];
	} else {
		if (x)
			for (i = 0; i < b->columns; i++)
				b->x[i] = x[i] * b->scaling->column[i];
		if (y)
			for (i = 0; i < b->rows; i++)
				b->y[i] = y[i] * b->scaling->row[i];
	}
	b->parent->report(b->parent->context, x ? b->x : NULL, y ? b->y : NULL,
		scale, stage, iterations, force);
}
int simplex_progress_bridge(struct simplex_ProgressBridge *b,
	const struct simplex_ProgressSink *parent,
	const struct simplex_Presolve *presolve,
	const struct simplex_Scaling *scaling, int rows, int columns)
{
	memset(b, 0, sizeof(*b));
	if (!parent)
		return 0;
	b->storage =
		lp_simplex_malloc(((size_t)rows + columns) * sizeof(double));
	if (!b->storage)
		return -1;
	b->parent = parent;
	b->presolve = presolve;
	b->scaling = scaling;
	b->rows = rows;
	b->columns = columns;
	b->x = b->storage;
	b->y = b->x + columns;
	b->sink.context = b;
	b->sink.due = bridge_due;
	b->sink.report = bridge_report;
	return 0;
}
void simplex_progress_bridge_destroy(struct simplex_ProgressBridge *b)
{
	lp_simplex_free(b->storage);
}
void simplex_progress_emit(const struct simplex_Problem *p, const double *x,
	const double *y, double scale, const char *stage, int iterations,
	int force)
{
	const struct simplex_ProgressSink *s = p->progress;
	if (s && (force || s->due(s->context)))
		s->report(s->context, x, y, scale, stage, iterations, force);
}
