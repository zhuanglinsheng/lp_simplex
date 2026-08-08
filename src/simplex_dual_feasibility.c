/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
/* Hypersparse maintenance of primal-infeasible basic rows. */
#include "simplex_dual_feasibility.h"
#include "simplex_dual_internal.h"
#include "utils.h"

#include <lp_simplex/status.h>


static int feasibility_better(
		const struct simplex_DualFeasibility *feasibility,
		const int left, const int right)
{
	return feasibility->score[left] > feasibility->score[right] ||
		(feasibility->score[left] == feasibility->score[right] &&
		 left < right);
}


static void feasibility_swap(int *heap, int *slot, const int a, const int b)
{
	int row = heap[a];
	heap[a] = heap[b];
	heap[b] = row;
	slot[heap[a]] = a;
	slot[heap[b]] = b;
}


static void feasibility_sift_up(
		const struct simplex_DualFeasibility *feasibility,
		int *heap, int *slot, int position)
{
	while (position > 0) {
		int parent = (position - 1) / 2;
		if (!feasibility_better(feasibility, heap[position], heap[parent]))
			break;
		feasibility_swap(heap, slot, position, parent);
		position = parent;
	}
}


static void feasibility_sift_down(
		const struct simplex_DualFeasibility *feasibility,
		int *heap, int *slot, const int count, int position)
{
	for (;;) {
		int child = 2 * position + 1;
		if (child >= count)
			break;
		if (child + 1 < count && feasibility_better(feasibility,
				heap[child + 1], heap[child]))
			child++;
		if (!feasibility_better(feasibility, heap[child], heap[position]))
			break;
		feasibility_swap(heap, slot, position, child);
		position = child;
	}
}


static void feasibility_set(
		struct simplex_DualFeasibility *feasibility,
		int *heap, int *slot, int *count, const int row, const int active)
{
	int position = slot[row];
	if (!active) {
		int moved;
		if (position < 0)
			return;
		(*count)--;
		slot[row] = -1;
		if (position != *count) {
			moved = heap[*count];
			heap[position] = moved;
			slot[moved] = position;
			feasibility_sift_up(feasibility, heap, slot, position);
			position = slot[moved];
			feasibility_sift_down(feasibility, heap, slot, *count, position);
		}
		return;
	}
	if (position < 0) {
		position = (*count)++;
		heap[position] = row;
		slot[row] = position;
		feasibility_sift_up(feasibility, heap, slot, position);
		return;
	}
	feasibility_sift_up(feasibility, heap, slot, position);
	position = slot[row];
	feasibility_sift_down(feasibility, heap, slot, *count, position);
}


static double feasibility_row_score(
		const struct simplex_DualState *state, const int row)
{
	double violation = 0.;
	if (state->basic_value[row] < state->basic_lower[row] -
	    state->options->primal_tolerance)
		violation = state->basic_lower[row] - state->basic_value[row];
	else if (state->basic_value[row] > state->basic_upper[row] +
		 state->options->primal_tolerance)
		violation = state->basic_value[row] - state->basic_upper[row];
	return violation > state->options->primal_tolerance
		? violation * violation / __lp_simplex_MAX__(state->edge_weight[row],
			1e-12) : 0.;
}


int simplex_dual_feasibility_create(
		struct simplex_DualFeasibility *feasibility, const int rows)
{
	lp_simplex_memset(feasibility, 0, sizeof(*feasibility));
	feasibility->rows = rows;
	feasibility->heap = (int *)lp_simplex_malloc((size_t)rows * sizeof(int));
	feasibility->slot = (int *)lp_simplex_malloc((size_t)rows * sizeof(int));
	feasibility->structural_heap = (int *)lp_simplex_malloc(
		(size_t)rows * sizeof(int));
	feasibility->structural_slot = (int *)lp_simplex_malloc(
		(size_t)rows * sizeof(int));
	feasibility->score = (double *)lp_simplex_malloc(
		(size_t)rows * sizeof(double));
	if (feasibility->heap == NULL || feasibility->slot == NULL ||
	    feasibility->structural_heap == NULL ||
	    feasibility->structural_slot == NULL || feasibility->score == NULL) {
		simplex_dual_feasibility_destroy(feasibility);
		return lp_simplex_EXIT_FAILURE;
	}
	return lp_simplex_EXIT_SUCCESS;
}


void simplex_dual_feasibility_destroy(
		struct simplex_DualFeasibility *feasibility)
{
	if (feasibility == NULL)
		return;
	lp_simplex_free(feasibility->heap);
	lp_simplex_free(feasibility->slot);
	lp_simplex_free(feasibility->structural_heap);
	lp_simplex_free(feasibility->structural_slot);
	lp_simplex_free(feasibility->score);
	lp_simplex_memset(feasibility, 0, sizeof(*feasibility));
}


void simplex_dual_feasibility_update(
		struct simplex_DualState *state, const int row)
{
	struct simplex_DualFeasibility *feasibility = &state->feasibility;
	feasibility->score[row] = feasibility_row_score(state, row);
	feasibility_set(feasibility, feasibility->heap, feasibility->slot,
		&feasibility->count, row, feasibility->score[row] > 0.);
	feasibility_set(feasibility, feasibility->structural_heap,
		feasibility->structural_slot, &feasibility->structural_count, row,
		feasibility->score[row] > 0. &&
		state->basis[row] < state->structural);
	feasibility->incremental_updates++;
}


void simplex_dual_feasibility_rebuild(struct simplex_DualState *state)
{
	int row;
	struct simplex_DualFeasibility *feasibility = &state->feasibility;
	feasibility->count = 0;
	feasibility->structural_count = 0;
	for (row = 0; row < feasibility->rows; row++) {
		feasibility->slot[row] = -1;
		feasibility->structural_slot[row] = -1;
	}
	for (row = 0; row < feasibility->rows; row++) {
		feasibility->score[row] = feasibility_row_score(state, row);
		if (feasibility->score[row] <= 0.)
			continue;
		feasibility->slot[row] = feasibility->count;
		feasibility->heap[feasibility->count++] = row;
		if (state->basis[row] < state->structural) {
			feasibility->structural_slot[row] =
				feasibility->structural_count;
			feasibility->structural_heap[
				feasibility->structural_count++] = row;
		}
	}
	for (row = feasibility->count / 2; row > 0; row--)
		feasibility_sift_down(feasibility, feasibility->heap,
			feasibility->slot, feasibility->count, row - 1);
	for (row = feasibility->structural_count / 2; row > 0; row--)
		feasibility_sift_down(feasibility, feasibility->structural_heap,
			feasibility->structural_slot,
			feasibility->structural_count, row - 1);
	feasibility->rebuilds++;
}


void simplex_dual_feasibility_update_packed(
		struct simplex_DualState *state,
		const struct simplex_SparseVector *changed_rows)
{
	int k;
	/* Heap updates win for truly hypersparse directions.  Once the direction
	 * touches a material fraction of the basis, Floyd heap construction is
	 * cheaper and much more cache friendly than O(nnz log m) repairs. */
	if (changed_rows->count * 4 >= state->rows) {
		simplex_dual_feasibility_rebuild(state);
		return;
	}
	for (k = 0; k < changed_rows->count; k++)
		simplex_dual_feasibility_update(state, changed_rows->index[k]);
}


int simplex_dual_feasibility_choose(
		struct simplex_DualState *state, double *target, int *kappa,
		double *maximum, const int prefer_structural)
{
	struct simplex_DualFeasibility *feasibility = &state->feasibility;
	int row;
	if (feasibility->count == 0)
		return -1;
	row = feasibility->heap[0];
	if (prefer_structural && feasibility->structural_count > 0) {
		int structural = feasibility->structural_heap[0];
		if (feasibility->score[structural] >= .8 * feasibility->score[row])
			row = structural;
	}
	if (state->basic_value[row] < state->basic_lower[row]) {
		*target = state->basic_lower[row];
		*kappa = 1;
		*maximum = state->basic_lower[row] - state->basic_value[row];
	} else {
		*target = state->basic_upper[row];
		*kappa = -1;
		*maximum = state->basic_value[row] - state->basic_upper[row];
	}
	return row;
}
