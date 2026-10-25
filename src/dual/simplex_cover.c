/**
 * @file simplex_cover.c
 * @brief Restricted-column solving for nonnegative covering LPs.
 */

#include "simplex_stage.h"
#include "simplex_solve_context.h"
#include "simplex_execution.h"
#include "policy/resources.h"
#include "policy/factorization.h"
#include "policy/controls.h"
#include "simplex_cover.h"
#include "simplex_dual.h"
#include "simplex_column_handoff.h"
#include "simplex_dual_internal.h"
#include "linalg.h"
#include "simplex_progress.h"
#include "utils.h"
#include <limits.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

struct cover_Candidate {
	int column, material;
	double score;
};

struct cover_Work {
	int *indices, *map, *reverse, *global_basis;
	unsigned char *flags, *selected, *global_status, *scratch;
	size_t scratch_bytes;
	int contract_failed;
	double *numbers, *small_x, *candidate, *dual, *certificate_dual;
	struct cover_Candidate *heap;
	struct simplex_DualBasisSeed seed;
	struct simplex_ColumnProjection projection;
	struct simplex_ColumnBasis parent_basis;
};

static int cover_allowed(const struct simplex_Problem *p)
{
	int i, j, k;
	long double n = p->columns, m = p->rows;
	if (p->rows <= 0 || p->columns <= 0 || p->columns > INT_MAX - p->rows ||
		(2 * n + 2 * m) * sizeof(int) + 3 * n + 2 * m +
				(2 * n + 2 * m) * sizeof(double) +
				(n + m + 7) / 8 +
				m * sizeof(struct cover_Candidate) >
			COVER_WORKSPACE_BYTES)
		return 0;
	for (i = 0; i < p->rows; i++)
		if (p->row_type[i] != optm_CONS_T_GE || !isfinite(p->rhs[i]) ||
			p->rhs[i] < 0.)
			return 0;
	for (j = 0; j < p->columns; j++) {
		const struct optm_VariableBound *b = p->bounds + j;
		if (!isfinite(p->objective[j]) || p->objective[j] < 0. ||
			b->lb != 0. ||
			(b->b_type != optm_BOUND_T_LO &&
				b->b_type != optm_BOUND_T_BS) ||
			(b->b_type == optm_BOUND_T_BS &&
				(!isfinite(b->ub) || b->ub < 0.)))
			return 0;
	}
	for (k = 0; k < p->matrix.nonzeros; k++)
		if (!isfinite(p->matrix.value[k]) || p->matrix.value[k] < 0.)
			return 0;
	return 1;
}

static void cover_destroy(struct cover_Work *w)
{
	lp_simplex_free(w->indices);
	lp_simplex_free(w->flags);
	lp_simplex_free(w->numbers);
	lp_simplex_free(w->heap);
}

static int cover_create(struct cover_Work *w, int n, int m)
{
	lp_simplex_memset(w, 0, sizeof(*w));
	w->indices =
		(int *)lp_simplex_malloc(((size_t)2 * n + 2 * m) * sizeof(int));
	w->scratch_bytes = ((size_t)n + m + 7) / 8;
	w->flags = (unsigned char *)lp_simplex_malloc(
		(size_t)3 * n + 2 * m + w->scratch_bytes);
	w->numbers = (double *)lp_simplex_malloc(
		((size_t)2 * n + 2 * (size_t)m) * sizeof(double));
	w->heap = (struct cover_Candidate *)lp_simplex_malloc(
		(size_t)m * sizeof(*w->heap));
	if (!w->indices || !w->flags || !w->numbers || !w->heap) {
		cover_destroy(w);
		return -1;
	}
	w->map = w->indices;
	w->reverse = w->map + n;
	w->global_basis = w->reverse + n;
	w->seed.index = w->global_basis + m;
	w->selected = w->flags;
	w->global_status = w->selected + n;
	w->seed.status = w->global_status + n + m;
	w->scratch = w->seed.status + n + m;
	w->small_x = w->numbers;
	w->candidate = w->small_x + n;
	w->dual = w->candidate + n;
	w->certificate_dual = w->dual + m;
	w->seed.rows = m;
	w->seed.prefer_dual = 1;
	lp_simplex_memset(w->selected, 0, (size_t)n);
	lp_simplex_memset(w->global_status, DUAL_STATUS_LOWER, (size_t)n + m);
	return 0;
}

/* One finite-capacity column per positive row is a sufficient feasibility
 * construction for nonnegative covering rows. Refuse this optional route when
 * no such column exists; multiple smaller capacities may still be feasible. */
static int cover_initial_columns(
	const struct simplex_Problem *p, struct cover_Work *w)
{
	int i, j, k, best;
	for (i = 0; i < p->rows; i++) {
		double score = HUGE_VAL;
		if (p->rhs[i] == 0.)
			continue;
		best = -1;
		for (k = p->matrix.row_start[i]; k < p->matrix.row_start[i + 1];
			k++) {
			double a = p->matrix.row_value[k], ratio;
			j = p->matrix.column_index[k];
			if (a <= 0. ||
				(p->bounds[j].b_type == optm_BOUND_T_BS &&
					a * p->bounds[j].ub < p->rhs[i]))
				continue;
			ratio = p->objective[j] / a;
			if (ratio < score) {
				score = ratio;
				best = j;
			}
		}
		if (best < 0)
			return -1;
		w->selected[best] = 1;
	}
	return 0;
}

static struct simplex_Problem *cover_master(
	const struct simplex_Problem *p, struct cover_Work *w)
{
	struct simplex_Problem *q;
	int i, j, k, t, n = 0, nz = 0, *next, events;
	for (j = 0; j < p->columns; j++)
		w->reverse[j] = -1;
	for (j = 0; j < p->columns; j++)
		if (w->selected[j]) {
			w->map[n] = j;
			w->reverse[j] = n++;
			nz += p->matrix.column_start[j + 1] -
			      p->matrix.column_start[j];
		}
	q = simplex_problem_create_sparse(p->rows, n, nz);
	if (!q)
		return NULL;
	next = (int *)lp_simplex_malloc(((size_t)p->rows + 1) * sizeof(int));
	if (!next) {
		simplex_problem_free(q);
		return NULL;
	}
	lp_simplex_memcpy(q->rhs, p->rhs, (size_t)p->rows * sizeof(double));
	lp_simplex_memcpy(q->row_type, p->row_type, (size_t)p->rows);
	lp_simplex_memset(
		q->matrix.row_start, 0, ((size_t)p->rows + 1) * sizeof(int));
	t = 0;
	for (j = 0; j < n; j++) {
		i = w->map[j];
		q->objective[j] = p->objective[i];
		q->bounds[j] = p->bounds[i];
		q->matrix.column_start[j] = t;
		for (k = p->matrix.column_start[i];
			k < p->matrix.column_start[i + 1]; k++) {
			q->matrix.row_index[t] = p->matrix.row_index[k];
			q->matrix.value[t++] = p->matrix.value[k];
			q->matrix.row_start[p->matrix.row_index[k] + 1]++;
		}
	}
	q->matrix.column_start[n] = t;
	for (i = 0; i < p->rows; i++) {
		q->matrix.row_start[i + 1] += q->matrix.row_start[i];
		next[i] = q->matrix.row_start[i];
	}
	for (j = 0; j < n; j++)
		for (k = q->matrix.column_start[j];
			k < q->matrix.column_start[j + 1]; k++) {
			t = next[q->matrix.row_index[k]]++;
			q->matrix.column_index[t] = j;
			q->matrix.row_value[t] = q->matrix.value[k];
		}
	events = p->execution.history ? p->execution.history->count : 0;
	q->execution = simplex_solve_context_derive(
		&p->execution, SIMPLEX_EXECUTION_COLUMN_MODEL);
	if (q->execution.history && q->execution.history->count > events) {
		q->execution.history->events[events].model_rows = q->rows;
		q->execution.history->events[events].model_columns = q->columns;
	}
	if (simplex_column_projection_bind(
		    &w->projection, p, q, w->map, w->reverse, next, p->rows)) {
		w->contract_failed = 1;
		lp_simplex_free(next);
		simplex_problem_free(q);
		return NULL;
	}
	lp_simplex_free(next);
	return q;
}

/* Bounded max-heap retains the strongest rows-sized batch in O(n log m)
 * time and O(m) memory; the root is the worst retained score/id pair. */
static int cover_worse(struct cover_Candidate a, struct cover_Candidate b)
{
	if (a.material != b.material)
		return a.material < b.material;
	return a.score > b.score || (a.score == b.score && a.column > b.column);
}

static void cover_offer(struct cover_Candidate *heap, int *count, int capacity,
	struct cover_Candidate candidate)
{
	int i, parent, child;
	if (*count < capacity) {
		i = (*count)++;
		while (i > 0) {
			parent = (i - 1) / 2;
			if (!cover_worse(candidate, heap[parent]))
				break;
			heap[i] = heap[parent];
			i = parent;
		}
		heap[i] = candidate;
		return;
	}
	if (!cover_worse(heap[0], candidate))
		return;
	i = 0;
	while ((child = 2 * i + 1) < *count) {
		if (child + 1 < *count &&
			cover_worse(heap[child + 1], heap[child]))
			child++;
		if (!cover_worse(heap[child], candidate))
			break;
		heap[i] = heap[child];
		i = child;
	}
	heap[i] = candidate;
}

/* Every feasible covering point can be clipped to this capacity without
 * losing a row: beyond max_i(rhs_i/a_ij), column j alone covers each of its
 * incident rows. Nonnegative costs make clipping objective-nonincreasing.
 * Thus the capped box has the same optimum, even for original LO columns. */
static double cover_certificate_capacity(const struct simplex_Problem *p, int j)
{
	double cap = 0.;
	int k;
	for (k = p->matrix.column_start[j]; k < p->matrix.column_start[j + 1];
		k++)
		if (p->matrix.value[k] > 0.)
			cap = fmax(cap,
				p->rhs[p->matrix.row_index[k]] == 0.
					? 0.
					: nextafter(p->rhs[p->matrix.row_index
								    [k]] /
							    p->matrix.value[k],
						  HUGE_VAL));
	if (p->bounds[j].b_type == optm_BOUND_T_BS)
		cap = fmin(cap, p->bounds[j].ub);
	return cap;
}

int simplex_cover_solve(const struct simplex_Problem *p,
	const struct lp_simplex_Options *o, double *x, double *row_dual,
	struct lp_simplex_Result *result)
{
	struct cover_Work w;
	struct simplex_Problem owned_input = *p;
	struct simplex_SolveContext local_context;
	struct simplex_StageReport master_report;
	int master_event = 0, previous_columns, prefer_dual;
	int owns_context = !p->execution.context;
	struct simplex_Problem *q = NULL;
	struct lp_simplex_Options stage = *o;
	struct lp_simplex_Result current;
	int i, j, code, count, used = 0, rounds = 0, have_point = 0,
			       refinements = 0, material_count = 0,
			       projected = 0, polish_attempted = 0;
	long double objective, bound, gap, accuracy, capacity;
	long double magnitude, uncertainty, gamma, roundoff, candidate_floor;
	double rounding_floor = nextafter(0., 1.);

	lp_simplex_memset(result, 0, sizeof(*result));
	result->status = lp_simplex_CondUnsatisfied;
	result->objective = NAN;
	result->primal_infeasibility = result->dual_infeasibility = HUGE_VAL;
	result->dual_residual_space = lp_simplex_RESIDUAL_ORIGINAL;
	if (p->report) {
		simplex_stage_reset(p->report);
		p->report->stop = SIMPLEX_STOP_NOT_APPLICABLE;
	}
	if (!cover_allowed(p))
		return lp_simplex_EXIT_FAILURE;
	if (!p->execution.context) {
		simplex_solve_context_init(
			&local_context, NULL, &owned_input.execution);
		p = &owned_input;
	}
	lp_simplex_memset(x, 0, (size_t)p->columns * sizeof(double));
	if (cover_create(&w, p->columns, p->rows)) {
		result->status = lp_simplex_MemoryAllocError;
		if (p->report)
			p->report->stop = SIMPLEX_STOP_RESOURCE;
		if (owns_context)
			simplex_solve_context_destroy(&local_context);
		return lp_simplex_EXIT_FAILURE;
	}
	simplex_column_basis_init(
		&w.parent_basis, p, w.global_basis, w.global_status);
	if (cover_initial_columns(p, &w))
		goto DONE;
	stage.algorithm = lp_simplex_ALGORITHM_DUAL_REVISED;
	while (used < o->iteration_limit) {
		if (!simplex_problem_identity_matches(
			    &w.parent_basis.parent, p)) {
			w.contract_failed = 1;
			result->status = lp_simplex_PrecisionError;
			goto DONE;
		}
		q = cover_master(p, &w);
		if (!q) {
			result->status = w.contract_failed
						 ? lp_simplex_PrecisionError
						 : lp_simplex_MemoryAllocError;
			goto DONE;
		}
		stage.iteration_limit = o->iteration_limit - used;
		simplex_stage_reset(&master_report);
		q->report = &master_report;
		master_event = simplex_execution_begin(&q->execution,
			SIMPLEX_EXECUTION_COVER_MASTER, 0,
			polish_attempted ? SIMPLEX_REASON_DUAL_REFINEMENT
					 : SIMPLEX_REASON_COLUMN_ADMISSION,
			0);
		simplex_execution_set_parent(&q->execution, master_event);
		if (master_event)
			q->execution.history->events[master_event - 1].grant =
				stage.iteration_limit;
		previous_columns = w.seed.columns;
		prefer_dual = q->columns - previous_columns >=
			      p->rows / 4 + (p->rows % 4 != 0);
		simplex_dual_seed_init(&w.seed, q, w.seed.index, w.seed.status);
		if (w.parent_basis.valid &&
			simplex_column_basis_restrict(&w.projection,
				&w.parent_basis, &w.seed, prefer_dual,
				w.scratch, w.scratch_bytes)) {
			w.contract_failed = 1;
			result->status = lp_simplex_PrecisionError;
			goto DONE;
		}
		lp_simplex_memset(&current, 0, sizeof(current));
		code = simplex_dual_solve_problem_seeded(
			q, &stage, w.small_x, w.dual, &current, 0, &w.seed);
		master_event = simplex_execution_parent(&q->execution);
		simplex_execution_end(q->execution.history, master_event,
			&current, &master_report);
		master_event = 0;
		if (current.iterations < 0 ||
			current.iterations > stage.iteration_limit) {
			used = o->iteration_limit;
			result->status = lp_simplex_PrecisionError;
			goto DONE;
		}
		used += current.iterations;
		if (code || current.status != lp_simplex_Success) {
			result->status =
				current.status == lp_simplex_ExceedIterLimit ||
						current.status ==
							lp_simplex_MemoryAllocError
					? current.status
					: lp_simplex_PrecisionError;
			goto DONE;
		}
		if (simplex_column_point_lift(&w.projection, w.small_x,
			    w.candidate, o->primal_tolerance) ||
			simplex_column_rows_lift(&w.projection, w.dual) ||
			simplex_column_basis_lift(&w.projection, &w.seed,
				&w.parent_basis, w.scratch, w.scratch_bytes)) {
			w.contract_failed = 1;
			result->status = lp_simplex_PrecisionError;
			goto DONE;
		}
		lp_simplex_memcpy(
			x, w.candidate, (size_t)p->columns * sizeof(double));
		have_point = 1;
		if (getenv(SIMPLEX_ENV_PROFILE))
			fprintf(stderr,
				"cover columns: round=%d selected=%d/%d work=%d seed=%d dual=%d objective=%.12g\n",
				rounds, q->columns, p->columns, used,
				w.seed.accepted, w.seed.prefer_dual,
				current.objective);
		w.projection.valid = 0;
		simplex_problem_free(q);
		q = NULL;
		simplex_progress_emit(
			p, x, w.dual, -1., "cover-columns", used, 0);
		if (used == o->iteration_limit)
			break;
		used++;
		rounds++;
		/* Price every column and form a box-dual lower bound for an
		 * equivalent capped covering domain. Local KKT tolerance alone
		 * does not control the objective over a wide variable range. */
		/* Candidate admission may ignore aggregate-negligible noise, while
		 * the objective certificate below still accounts for every column. */
		candidate_floor = o->dual_tolerance *
				  (1. + fabsl((long double)current.objective)) /
				  (16. * (1. + p->columns));
		objective = bound = magnitude = uncertainty = 0.;
		projected = 0;
		for (i = 0; i < p->rows; i++) {
			long double term;
			if (!isfinite(w.dual[i])) {
				result->status = lp_simplex_PrecisionError;
				goto DONE;
			}
			projected |= w.dual[i] < 0.;
			w.certificate_dual[i] = fmax(0., w.dual[i]);
			term = (long double)p->rhs[i] * w.certificate_dual[i];
			bound += term;
			magnitude += fabsl(term);
		}
		count = material_count = 0;
		for (j = 0; j < p->columns; j++) {
			if (!isfinite(x[j])) {
				result->status = lp_simplex_PrecisionError;
				goto DONE;
			}
			objective += (long double)p->objective[j] * x[j];
			{
				double price_dot = lp_simplex_linalg_sparse_dot(
					p->matrix.column_start[j],
					p->matrix.column_start[j + 1],
					p->matrix.row_index, p->matrix.value,
					w.dual);
				double dot =
					projected
						? lp_simplex_linalg_sparse_dot(
							  p->matrix.column_start
								  [j],
							  p->matrix.column_start
								  [j + 1],
							  p->matrix.row_index,
							  p->matrix.value,
							  w.certificate_dual)
						: price_dot;
				double pricing_reduced =
					p->objective[j] - price_dot;
				double reduced = p->objective[j] - dot;
				int material =
					pricing_reduced < -o->dual_tolerance;
				double cap = 0.;
				double error =
					((double)(p->matrix
							  .column_start[j + 1] -
						  p->matrix.column_start[j]) +
						2.) *
					(8. * DBL_EPSILON *
							(fabs(p->objective[j]) +
								fabs(dot)) +
						8. * rounding_floor);
				if (!isfinite(reduced) || !isfinite(error) ||
					!isfinite(pricing_reduced)) {
					result->status =
						lp_simplex_PrecisionError;
					goto DONE;
				}
				if (reduced < error) {
					cap = cover_certificate_capacity(p, j);
					uncertainty += (long double)cap * error;
					if (reduced < 0.) {
						long double term =
							(long double)reduced *
							cap;
						bound += term;
						magnitude += fabsl(term);
					}
				}
				if (!w.selected[j] &&
					(material ||
						(reduced < 0. &&
							-(long double)reduced *
									cap >
								candidate_floor)) &&
					!(p->bounds[j].b_type ==
							optm_BOUND_T_BS &&
						p->bounds[j].ub == 0.)) {
					struct cover_Candidate candidate;
					candidate.column = j;
					candidate.material = material;
					material_count += material;
					candidate.score =
						(material ? pricing_reduced
							  : reduced) /
						(1. + p->objective[j]);
					cover_offer(w.heap, &count, p->rows,
						candidate);
				}
			}
		}
		gap = fabsl(objective - bound);
		accuracy = o->dual_tolerance * (1. + fabsl(objective));
		gamma = ((long double)p->columns + p->rows + 8.) * DBL_EPSILON;
		roundoff =
			gamma < SIMPLEX_COVER_ROUNDING_LIMIT
				? 2. * uncertainty +
					  8. *
						  ((long double)p->columns +
							  p->rows + 8.) *
						  rounding_floor +
					  4. * gamma / (1. - 2. * gamma) *
						  (magnitude + fabsl(objective))
				: HUGE_VAL;
		if (isfinite(gap) && isfinite(roundoff) &&
			gap + roundoff <= accuracy) {
			result->primal_infeasibility =
				simplex_problem_primal_infeasibility(p, x);
			result->dual_infeasibility =
				simplex_problem_dual_infeasibility(p, x,
					w.certificate_dual, o->primal_tolerance,
					-1);
			result->status =
				result->primal_infeasibility <=
							o->primal_tolerance &&
						result->dual_infeasibility <=
							o->dual_tolerance
					? lp_simplex_Success
					: lp_simplex_PrecisionError;
			if (result->status == lp_simplex_Success)
				goto DONE;
		}
		if (getenv(SIMPLEX_ENV_PROFILE))
			fprintf(stderr,
				"cover certificate: work=%d gap=%.6Lg rounding=%.6Lg budget=%.6Lg material=%d pool=%d polish=%d\n",
				used, gap, roundoff, accuracy, material_count,
				count, polish_attempted);
		if (!material_count && (!polish_attempted || !count)) {
			/* The restricted solve can satisfy its local tolerance
			 * while failing the full objective budget, particularly
			 * for LO columns. Reuse its basis at a tighter tolerance. */
			double tighter;
			if (++refinements > 3) {
				result->status = lp_simplex_PrecisionError;
				goto DONE;
			}
			capacity = 0.;
			for (j = 0; j < p->columns; j++)
				if (w.selected[j])
					capacity += cover_certificate_capacity(
						p, j);
			tighter = (double)(accuracy / (4. * (1. + capacity)));
			tighter = fmin(tighter,
				stage.dual_tolerance *
					SIMPLEX_COVER_RETRY_TOLERANCE_FACTOR);
			if (!isfinite(tighter) || tighter <= 0. ||
				tighter >= stage.dual_tolerance) {
				result->status = lp_simplex_PrecisionError;
				goto DONE;
			}
			stage.dual_tolerance = tighter;
			polish_attempted = 1;
			continue;
		}
		if (!count) {
			result->status = lp_simplex_PrecisionError;
			goto DONE;
		}
		refinements = polish_attempted = 0;
		/* Fill the master with material pricing violations first. Small
		 * objective-budget candidates get their own batch after a dual polish,
		 * rather than changing the size and starting basis of a material batch. */
		for (i = 0; i < count; i++)
			if (!material_count || w.heap[i].material)
				w.selected[w.heap[i].column] = 1;
	}
	result->status = lp_simplex_ExceedIterLimit;
DONE:
	if (master_event && q) {
		master_event = simplex_execution_parent(&q->execution);
		simplex_stage_reset(&master_report);
		master_report.stop = SIMPLEX_STOP_NUMERICAL;
		master_report.failure = SIMPLEX_FAILURE_STAGE_CONTRACT;
		simplex_execution_end(q->execution.history, master_event, NULL,
			&master_report);
	}
	if (!simplex_problem_identity_matches(&w.parent_basis.parent, p)) {
		result->status = lp_simplex_PrecisionError;
		have_point = 0;
		w.contract_failed = 1;
	}
	if (have_point) {
		objective = 0.;
		for (j = 0; j < p->columns; j++)
			objective += (long double)p->objective[j] * x[j];
		result->objective = (double)objective;
		if (!isfinite(result->objective))
			result->status = lp_simplex_PrecisionError;
	}
	if (row_dual && result->status == lp_simplex_Success)
		lp_simplex_memcpy(row_dual, w.certificate_dual,
			(size_t)p->rows * sizeof(double));
	result->iterations = used;
	if (have_point)
		result->primal_infeasibility =
			simplex_problem_primal_infeasibility(p, x);
	if (p->report) {
		simplex_stage_reset(p->report);
		p->report->stop =
			simplex_stage_stop_from_status(result->status);
		p->report->point_available = have_point;
		if (w.contract_failed)
			p->report->failure = SIMPLEX_FAILURE_STAGE_CONTRACT;
		simplex_stage_observe(
			p, o, have_point ? x : NULL, NULL, -1, used, p->report);
	}
	w.projection.valid = 0;
	simplex_problem_free(q);
	cover_destroy(&w);
	if (owns_context)
		simplex_solve_context_destroy(&local_context);
	return result->status == lp_simplex_Success ? lp_simplex_EXIT_SUCCESS
						    : lp_simplex_EXIT_FAILURE;
}
