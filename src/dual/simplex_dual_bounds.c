/* Bound initialization policy used before the dual crash. */
#include "simplex_dual_bounds.h"
#include "utils.h"


struct dual_BoundQueue {
	int *row;
	int *queued;
	int capacity;
	int head;
	int tail;
	int count;
};


static int dual_bound_has_lower(const double value)
{
	return value > __lp_simplex_NINF__;
}


static int dual_bound_has_upper(const double value)
{
	return value < __lp_simplex_INF__;
}


static void dual_bound_queue_push(
		struct dual_BoundQueue *queue, const int row)
{
	if (queue->queued[row] || queue->count >= queue->capacity)
		return;
	queue->row[queue->tail] = row;
	queue->tail = (queue->tail + 1) % queue->capacity;
	queue->count++;
	queue->queued[row] = 1;
}


static int dual_bound_queue_pop(struct dual_BoundQueue *queue)
{
	int row = queue->row[queue->head];
	queue->head = (queue->head + 1) % queue->capacity;
	queue->count--;
	queue->queued[row] = 0;
	return row;
}


static void dual_bound_queue_column(
		const struct simplex_DualState *state,
		struct dual_BoundQueue *queue, const int column)
{
	int k;
	for (k = state->matrix.column_start[column];
	     k < state->matrix.column_start[column + 1]; k++)
		dual_bound_queue_push(queue, state->matrix.row_index[k]);
}


static void dual_initialize_column_bounds(
		struct simplex_DualState *state, const struct simplex_Problem *problem)
{
	int i, j;
	for (j = 0; j < problem->columns; j++) {
		const struct optm_VariableBound *bound = problem->bounds + j;
		state->lower[j] = __lp_simplex_NINF__;
		state->upper[j] = __lp_simplex_INF__;
		if (bound->b_type == optm_BOUND_T_LO || bound->b_type == optm_BOUND_T_BS)
			state->lower[j] = bound->lb;
		if (bound->b_type == optm_BOUND_T_UP || bound->b_type == optm_BOUND_T_BS)
			state->upper[j] = bound->ub;
		state->cost[j] = problem->objective[j];
	}
	for (i = 0; i < problem->rows; i++) {
		int variable = problem->columns + i;
		int type = problem->row_type[i];
		state->lower[variable] = __lp_simplex_NINF__;
		state->upper[variable] = __lp_simplex_INF__;
		if (type == optm_CONS_T_EQ || type == optm_CONS_T_GE)
			state->lower[variable] = problem->rhs[i];
		if (type == optm_CONS_T_EQ || type == optm_CONS_T_LE)
			state->upper[variable] = problem->rhs[i];
		state->cost[variable] = 0.;
	}
}


static void dual_propagate_row_bounds(
		struct simplex_DualState *state, const struct simplex_Problem *problem,
		struct dual_BoundQueue *queue, const int row)
{
	int k;
	int type = problem->row_type[row];
	double rhs = problem->rhs[row];
	double finite_minimum = 0., finite_maximum = 0.;
	int minimum_infinite = 0, maximum_infinite = 0;
	for (k = state->matrix.row_start[row];
	     k < state->matrix.row_start[row + 1]; k++) {
		int j = state->matrix.column_index[k];
		double a = state->matrix.row_value[k];
		if ((a > 0. && !dual_bound_has_lower(state->lower[j])) ||
		    (a < 0. && !dual_bound_has_upper(state->upper[j])))
			minimum_infinite++;
		else
			finite_minimum += a * (a > 0.
				? state->lower[j] : state->upper[j]);
		if ((a > 0. && !dual_bound_has_upper(state->upper[j])) ||
		    (a < 0. && !dual_bound_has_lower(state->lower[j])))
			maximum_infinite++;
		else
			finite_maximum += a * (a > 0.
				? state->upper[j] : state->lower[j]);
	}
	for (k = state->matrix.row_start[row];
	     k < state->matrix.row_start[row + 1]; k++) {
		int changed = 0;
		int j = state->matrix.column_index[k];
		double a = state->matrix.row_value[k];
		int own_minimum_infinite, own_maximum_infinite;
		double own_minimum = 0., own_maximum = 0.;
		own_minimum_infinite =
			(a > 0. && !dual_bound_has_lower(state->lower[j])) ||
			(a < 0. && !dual_bound_has_upper(state->upper[j]));
		own_maximum_infinite =
			(a > 0. && !dual_bound_has_upper(state->upper[j])) ||
			(a < 0. && !dual_bound_has_lower(state->lower[j]));
		if (!own_minimum_infinite)
			own_minimum = a * (a > 0.
				? state->lower[j] : state->upper[j]);
		if (!own_maximum_infinite)
			own_maximum = a * (a > 0.
				? state->upper[j] : state->lower[j]);
		if (type != optm_CONS_T_GE &&
		    minimum_infinite - own_minimum_infinite == 0) {
			double bound = (rhs - (finite_minimum - own_minimum)) / a;
			if (a > 0. && bound < state->upper[j]) {
				state->upper[j] = bound;
				changed = 1;
			} else if (a < 0. && bound > state->lower[j]) {
				state->lower[j] = bound;
				changed = 1;
			}
		}
		if (type != optm_CONS_T_LE &&
		    maximum_infinite - own_maximum_infinite == 0) {
			double bound = (rhs - (finite_maximum - own_maximum)) / a;
			if (a > 0. && bound > state->lower[j]) {
				state->lower[j] = bound;
				changed = 1;
			} else if (a < 0. && bound < state->upper[j]) {
				state->upper[j] = bound;
				changed = 1;
			}
		}
		if (changed)
			dual_bound_queue_column(state, queue, j);
	}
}


static void dual_propagate_bounds(
		struct simplex_DualState *state, const struct simplex_Problem *problem)
{
	struct dual_BoundQueue queue;
	long processed = 0;
	long budget = 8L * problem->rows;
	int row;
	/* Pricing is not initialized yet, so its integer workspace is a safe
	 * same-lifetime queue owner and avoids another solve-time allocation. */
	queue.row = state->candidate_index;
	queue.queued = state->candidate_sign;
	queue.capacity = problem->rows;
	queue.head = 0;
	queue.tail = 0;
	queue.count = 0;
	lp_simplex_memset(queue.queued, 0,
		(size_t)problem->rows * sizeof(int));
	for (row = 0; row < problem->rows; row++)
		dual_bound_queue_push(&queue, row);
	while (queue.count > 0 && processed++ < budget) {
		row = dual_bound_queue_pop(&queue);
		dual_propagate_row_bounds(state, problem, &queue, row);
	}
}


void simplex_dual_initialize_bounds(
		struct simplex_DualState *state,
		const struct simplex_Problem *problem, const int propagate)
{
	dual_initialize_column_bounds(state, problem);
	if (propagate)
		dual_propagate_bounds(state, problem);
}
