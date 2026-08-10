/* Dynamic sparse R with stable node handles and simultaneous row/column views. */
#include "pan_dynamic_sparse_qr.h"
#include "utils.h"

#include <lp_simplex/status.h>

#include <math.h>


static int pan_factor_reserve(
		struct pan_DynamicSparseQR *factor, const int capacity)
{
	struct pan_DynamicSparseRow *row;
	double *rotation_first, *rotation_second, *event_cosine, *event_sine;
	double *packed_diagonal;
	int *tag, *union_column, *scratch_stamp, *event_row;
	int *packed_start;
	int *column_head, *column_tail, *column_count;
	int next;
	if (capacity <= factor->capacity)
		return lp_simplex_EXIT_SUCCESS;
	next = factor->capacity > 0 ? 2 * factor->capacity : 8;
	if (next < capacity)
		next = capacity;
	row = (struct pan_DynamicSparseRow *)lp_simplex_malloc(
		(size_t)next * sizeof(*row));
	tag = (int *)lp_simplex_malloc((size_t)next * sizeof(int));
	column_head = (int *)lp_simplex_malloc((size_t)next * sizeof(int));
	column_tail = (int *)lp_simplex_malloc((size_t)next * sizeof(int));
	column_count = (int *)lp_simplex_malloc((size_t)next * sizeof(int));
	union_column = (int *)lp_simplex_malloc((size_t)(3 * next) * sizeof(int));
	scratch_stamp = (int *)lp_simplex_malloc((size_t)next * sizeof(int));
	rotation_first = (double *)lp_simplex_malloc(
		(size_t)(2 * next) * sizeof(double));
	rotation_second = (double *)lp_simplex_malloc(
		(size_t)(2 * next) * sizeof(double));
	event_row = (int *)lp_simplex_malloc((size_t)next * sizeof(int));
	event_cosine = (double *)lp_simplex_malloc((size_t)next * sizeof(double));
	event_sine = (double *)lp_simplex_malloc((size_t)next * sizeof(double));
	packed_start = (int *)lp_simplex_malloc(
		(size_t)(next + 1) * sizeof(int));
	packed_diagonal = (double *)lp_simplex_malloc(
		(size_t)next * sizeof(double));
	if (row == NULL || tag == NULL || column_head == NULL ||
	    column_tail == NULL || column_count == NULL || union_column == NULL ||
	    scratch_stamp == NULL ||
	    rotation_first == NULL || rotation_second == NULL ||
	    event_row == NULL || event_cosine == NULL || event_sine == NULL ||
	    packed_start == NULL || packed_diagonal == NULL) {
		lp_simplex_free(row); lp_simplex_free(tag);
		lp_simplex_free(column_head); lp_simplex_free(column_tail);
		lp_simplex_free(column_count); lp_simplex_free(union_column);
		lp_simplex_free(scratch_stamp);
		lp_simplex_free(rotation_first); lp_simplex_free(rotation_second);
		lp_simplex_free(event_row); lp_simplex_free(event_cosine);
		lp_simplex_free(event_sine);
		lp_simplex_free(packed_start); lp_simplex_free(packed_diagonal);
		return lp_simplex_EXIT_FAILURE;
	}
	lp_simplex_memset(row, 0xff, (size_t)next * sizeof(*row));
	lp_simplex_memset(column_head, 0xff, (size_t)next * sizeof(int));
	lp_simplex_memset(column_tail, 0xff, (size_t)next * sizeof(int));
	lp_simplex_memset(column_count, 0, (size_t)next * sizeof(int));
	lp_simplex_memset(scratch_stamp, 0, (size_t)next * sizeof(int));
	if (factor->count > 0) {
		lp_simplex_memcpy(row, factor->row,
			(size_t)factor->count * sizeof(*row));
		lp_simplex_memcpy(tag, factor->tag,
			(size_t)factor->count * sizeof(int));
		lp_simplex_memcpy(column_head, factor->column_head,
			(size_t)factor->count * sizeof(int));
		lp_simplex_memcpy(column_tail, factor->column_tail,
			(size_t)factor->count * sizeof(int));
		lp_simplex_memcpy(column_count, factor->column_count,
			(size_t)factor->count * sizeof(int));
	}
	lp_simplex_free(factor->row); lp_simplex_free(factor->tag);
	lp_simplex_free(factor->column_head); lp_simplex_free(factor->column_tail);
	lp_simplex_free(factor->column_count); lp_simplex_free(factor->union_column);
	lp_simplex_free(factor->scratch_stamp);
	lp_simplex_free(factor->rotation_first); lp_simplex_free(factor->rotation_second);
	lp_simplex_free(factor->event_row); lp_simplex_free(factor->event_cosine);
	lp_simplex_free(factor->event_sine);
	lp_simplex_free(factor->packed_start);
	lp_simplex_free(factor->packed_diagonal);
	factor->row = row; factor->tag = tag;
	factor->column_head = column_head; factor->column_tail = column_tail;
	factor->column_count = column_count; factor->union_column = union_column;
	factor->scratch_stamp = scratch_stamp;
	factor->rotation_first = rotation_first;
	factor->rotation_second = rotation_second; factor->event_row = event_row;
	factor->event_cosine = event_cosine; factor->event_sine = event_sine;
	factor->packed_start = packed_start;
	factor->packed_diagonal = packed_diagonal;
	factor->capacity = next;
	return lp_simplex_EXIT_SUCCESS;
}


static int pan_packed_reserve(
		struct pan_DynamicSparseQR *factor, const int capacity)
{
	int *column;
	double *value;
	int next;
	if (capacity <= factor->packed_capacity)
		return lp_simplex_EXIT_SUCCESS;
	next = factor->packed_capacity > 0 ? 2 * factor->packed_capacity : 64;
	while (next < capacity)
		next *= 2;
	column = (int *)lp_simplex_malloc((size_t)next * sizeof(int));
	value = (double *)lp_simplex_malloc((size_t)next * sizeof(double));
	if (column == NULL || value == NULL) {
		lp_simplex_free(column);
		lp_simplex_free(value);
		return lp_simplex_EXIT_FAILURE;
	}
	lp_simplex_free(factor->packed_column);
	lp_simplex_free(factor->packed_value);
	factor->packed_column = column;
	factor->packed_value = value;
	factor->packed_capacity = next;
	return lp_simplex_EXIT_SUCCESS;
}


static int pan_pack_numeric(struct pan_DynamicSparseQR *factor)
{
	int nonzeros = 0;
	int row;
	for (row = 0; row < factor->count; row++)
		nonzeros += factor->row[row].count;
	if (pan_packed_reserve(factor, nonzeros) == lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	nonzeros = 0;
	for (row = 0; row < factor->count; row++) {
		int node = factor->row[row].head;
		double diagonal = 0.;
		factor->packed_start[row] = nonzeros;
		while (node >= 0) {
			int column = factor->node[node].column;
			factor->packed_column[nonzeros] = column;
			factor->packed_value[nonzeros++] = factor->node[node].value;
			if (column == row)
				diagonal = factor->node[node].value;
			node = factor->node[node].row_next;
		}
		factor->packed_diagonal[row] = diagonal;
	}
	factor->packed_start[factor->count] = nonzeros;
	factor->packed_valid = 1;
	return lp_simplex_EXIT_SUCCESS;
}


static int pan_node_reserve(
		struct pan_DynamicSparseQR *factor, const int capacity)
{
	struct pan_DynamicSparseNode *node;
	int next;
	if (capacity <= factor->node_capacity)
		return lp_simplex_EXIT_SUCCESS;
	next = factor->node_capacity > 0 ? 2 * factor->node_capacity : 32;
	while (next < capacity)
		next *= 2;
	node = (struct pan_DynamicSparseNode *)lp_simplex_malloc(
		(size_t)next * sizeof(*node));
	if (node == NULL)
		return lp_simplex_EXIT_FAILURE;
	if (factor->node_count > 0)
		lp_simplex_memcpy(node, factor->node,
			(size_t)factor->node_count * sizeof(*node));
	lp_simplex_free(factor->node);
	factor->node = node;
	factor->node_capacity = next;
	return lp_simplex_EXIT_SUCCESS;
}


static int pan_node_allocate(struct pan_DynamicSparseQR *factor)
{
	int index;
	if (factor->free_node >= 0) {
		index = factor->free_node;
		factor->free_node = factor->node[index].owner.free_next;
		return index;
	}
	if (pan_node_reserve(factor, factor->node_count + 1) ==
	    lp_simplex_EXIT_FAILURE)
		return -1;
	return factor->node_count++;
}


static void pan_node_unlink(struct pan_DynamicSparseQR *factor, const int index)
{
	struct pan_DynamicSparseNode *node = factor->node + index;
	struct pan_DynamicSparseRow *row = factor->row + node->owner.row;
	if (node->row_previous >= 0)
		factor->node[node->row_previous].row_next = node->row_next;
	else
		row->head = node->row_next;
	if (node->row_next >= 0)
		factor->node[node->row_next].row_previous = node->row_previous;
	else
		row->tail = node->row_previous;
	if (node->column_previous >= 0)
		factor->node[node->column_previous].column_next = node->column_next;
	else
		factor->column_head[node->column] = node->column_next;
	if (node->column_next >= 0)
		factor->node[node->column_next].column_previous = node->column_previous;
	else
		factor->column_tail[node->column] = node->column_previous;
	row->count--;
	factor->column_count[node->column]--;
	node->owner.free_next = factor->free_node;
	factor->free_node = index;
}


static int pan_row_find(
		const struct pan_DynamicSparseQR *factor,
		const int row, const int column)
{
	int node = factor->row[row].head;
	while (node >= 0 && factor->node[node].column < column)
		node = factor->node[node].row_next;
	return node >= 0 && factor->node[node].column == column ? node : -1;
}


static double pan_get(
		const struct pan_DynamicSparseQR *factor,
		const int row, const int column)
{
	int node = pan_row_find(factor, row, column);
	return node >= 0 ? factor->node[node].value : 0.;
}


static int pan_append_node(struct pan_DynamicSparseQR *factor,
		const int row_index, const int column, const double value)
{
	struct pan_DynamicSparseRow *row = factor->row + row_index;
	int index;
	if (value == 0.)
		return lp_simplex_EXIT_SUCCESS;
	index = pan_node_allocate(factor);
	if (index < 0)
		return lp_simplex_EXIT_FAILURE;
	factor->node[index].owner.row = row_index;
	factor->node[index].column = column;
	factor->node[index].value = value;
	factor->node[index].row_previous = row->tail;
	factor->node[index].row_next = -1;
	if (row->tail >= 0)
		factor->node[row->tail].row_next = index;
	else
		row->head = index;
	row->tail = index;
	factor->node[index].column_previous = factor->column_tail[column];
	factor->node[index].column_next = -1;
	if (factor->column_tail[column] >= 0)
		factor->node[factor->column_tail[column]].column_next = index;
	else
		factor->column_head[column] = index;
	factor->column_tail[column] = index;
	row->count++;
	factor->column_count[column]++;
	return lp_simplex_EXIT_SUCCESS;
}


static int pan_set(struct pan_DynamicSparseQR *factor,
		const int row_index, const int column, const double value)
{
	struct pan_DynamicSparseRow *row = factor->row + row_index;
	int next = row->head;
	int previous = -1;
	int index;
	if (value != 0. && (row->tail < 0 ||
	    factor->node[row->tail].column < column))
		return pan_append_node(factor, row_index, column, value);
	while (next >= 0 && factor->node[next].column < column) {
		previous = next;
		next = factor->node[next].row_next;
	}
	if (next >= 0 && factor->node[next].column == column) {
		if (value == 0.)
			pan_node_unlink(factor, next);
		else
			factor->node[next].value = value;
		return lp_simplex_EXIT_SUCCESS;
	}
	if (value == 0.)
		return lp_simplex_EXIT_SUCCESS;
	index = pan_node_allocate(factor);
	if (index < 0)
		return lp_simplex_EXIT_FAILURE;
	factor->node[index].owner.row = row_index;
	factor->node[index].column = column;
	factor->node[index].value = value;
	factor->node[index].row_previous = previous;
	factor->node[index].row_next = next;
	if (previous >= 0)
		factor->node[previous].row_next = index;
	else
		row->head = index;
	if (next >= 0)
		factor->node[next].row_previous = index;
	else
		row->tail = index;
	factor->node[index].column_previous = factor->column_tail[column];
	factor->node[index].column_next = -1;
	if (factor->column_tail[column] >= 0)
		factor->node[factor->column_tail[column]].column_next = index;
	else
		factor->column_head[column] = index;
	factor->column_tail[column] = index;
	row->count++;
	factor->column_count[column]++;
	return lp_simplex_EXIT_SUCCESS;
}


static void pan_clear_row(struct pan_DynamicSparseQR *factor, const int row)
{
	int node = factor->row[row].head;
	while (node >= 0) {
		int next = factor->node[node].row_next;
		pan_node_unlink(factor, node);
		node = next;
	}
}


static void pan_release_detached_node(
		struct pan_DynamicSparseQR *factor, const int index)
{
	struct pan_DynamicSparseNode *node = factor->node + index;
	if (node->column_previous >= 0)
		factor->node[node->column_previous].column_next = node->column_next;
	else
		factor->column_head[node->column] = node->column_next;
	if (node->column_next >= 0)
		factor->node[node->column_next].column_previous = node->column_previous;
	else
		factor->column_tail[node->column] = node->column_previous;
	factor->column_count[node->column]--;
	node->owner.free_next = factor->free_node;
	factor->free_node = index;
}


static void pan_move_node_column(struct pan_DynamicSparseQR *factor,
		const int index, const int column)
{
	struct pan_DynamicSparseNode *node = factor->node + index;
	int old_column = node->column;
	if (node->column_previous >= 0)
		factor->node[node->column_previous].column_next = node->column_next;
	else
		factor->column_head[old_column] = node->column_next;
	if (node->column_next >= 0)
		factor->node[node->column_next].column_previous =
			node->column_previous;
	else
		factor->column_tail[old_column] = node->column_previous;
	factor->column_count[old_column]--;
	node->column = column;
	node->column_previous = factor->column_tail[column];
	node->column_next = -1;
	if (factor->column_tail[column] >= 0)
		factor->node[factor->column_tail[column]].column_next = index;
	else
		factor->column_head[column] = index;
	factor->column_tail[column] = index;
	factor->column_count[column]++;
}


static int pan_new_detached_node(
		struct pan_DynamicSparseQR *factor, const int column)
{
	int index = pan_node_allocate(factor);
	if (index < 0)
		return -1;
	factor->node[index].column = column;
	factor->node[index].column_previous = factor->column_tail[column];
	factor->node[index].column_next = -1;
	if (factor->column_tail[column] >= 0)
		factor->node[factor->column_tail[column]].column_next = index;
	else
		factor->column_head[column] = index;
	factor->column_tail[column] = index;
	factor->column_count[column]++;
	return index;
}


static double pan_stable_hypot(const double first, const double second)
{
	double a = fabs(first), b = fabs(second), swap;
	if (a < b) { swap = a; a = b; b = swap; }
	if (a == 0.) return 0.;
	b /= a;
	return a * sqrt(1. + b * b);
}


static int pan_rotate_rows(struct pan_DynamicSparseQR *factor,
		const int first, const int second, const double c, const double s)
{
	int a = factor->row[first].head, b = factor->row[second].head;
	int *first_node = factor->union_column + factor->capacity;
	int *second_node = factor->union_column + 2 * factor->capacity;
	int total = 0, i;
	int old_first = factor->row[first].count;
	int old_second = factor->row[second].count;
	while (a >= 0 || b >= 0) {
		int column;
		double av = 0., bv = 0.;
		if (b < 0 || (a >= 0 && factor->node[a].column < factor->node[b].column)) {
			column = factor->node[a].column; av = factor->node[a].value;
			first_node[total] = a; second_node[total] = -1;
			a = factor->node[a].row_next;
		} else if (a < 0 || factor->node[b].column < factor->node[a].column) {
			column = factor->node[b].column; bv = factor->node[b].value;
			first_node[total] = -1; second_node[total] = b;
			b = factor->node[b].row_next;
		} else {
			column = factor->node[a].column; av = factor->node[a].value;
			bv = factor->node[b].value;
			first_node[total] = a; second_node[total] = b;
			a = factor->node[a].row_next; b = factor->node[b].row_next;
		}
		factor->union_column[total] = column;
		factor->rotation_first[total] = c * av + s * bv;
		factor->rotation_second[total++] = -s * av + c * bv;
	}
	factor->row[first].head = factor->row[first].tail = -1;
	factor->row[second].head = factor->row[second].tail = -1;
	factor->row[first].count = factor->row[second].count = 0;
	for (i = 0; i < total; i++) {
		int top = -1, bottom = -1;
		int one = first_node[i], two = second_node[i];
		if (factor->rotation_first[i] != 0.) {
			top = one >= 0 ? one : two;
			if (top < 0)
				top = pan_new_detached_node(factor,
					factor->union_column[i]);
		}
		if (factor->rotation_second[i] != 0.) {
			bottom = two >= 0 && two != top ? two :
				(one >= 0 && one != top ? one : -1);
			if (bottom < 0)
				bottom = pan_new_detached_node(factor,
					factor->union_column[i]);
		}
		if ((factor->rotation_first[i] != 0. && top < 0) ||
		    (factor->rotation_second[i] != 0. && bottom < 0))
			return lp_simplex_EXIT_FAILURE;
		if (one >= 0 && one != top && one != bottom)
			pan_release_detached_node(factor, one);
		if (two >= 0 && two != top && two != bottom)
			pan_release_detached_node(factor, two);
		if (top >= 0) {
			struct pan_DynamicSparseRow *row = factor->row + first;
			factor->node[top].owner.row = first;
			factor->node[top].value = factor->rotation_first[i];
			factor->node[top].row_previous = row->tail;
			factor->node[top].row_next = -1;
			if (row->tail >= 0) factor->node[row->tail].row_next = top;
			else row->head = top;
			row->tail = top; row->count++;
		}
		if (bottom >= 0) {
			struct pan_DynamicSparseRow *row = factor->row + second;
			factor->node[bottom].owner.row = second;
			factor->node[bottom].value = factor->rotation_second[i];
			factor->node[bottom].row_previous = row->tail;
			factor->node[bottom].row_next = -1;
			if (row->tail >= 0) factor->node[row->tail].row_next = bottom;
			else row->head = bottom;
			row->tail = bottom; row->count++;
		}
	}
	if (factor->row[first].count > old_first ||
	    factor->row[second].count > old_second)
		factor->symbolic_updates++;
	factor->rotations++;
	factor->event_row[factor->event_count] = first;
	factor->event_cosine[factor->event_count] = c;
	factor->event_sine[factor->event_count++] = s;
	return lp_simplex_EXIT_SUCCESS;
}


int pan_dynamic_sparse_qr_initialize(struct pan_DynamicSparseQR *factor)
{
	lp_simplex_memset(factor, 0, sizeof(*factor));
	factor->free_node = -1;
	factor->packed_valid = 0;
	return lp_simplex_EXIT_SUCCESS;
}


void pan_dynamic_sparse_qr_clear(struct pan_DynamicSparseQR *factor)
{
	int i;
	for (i = 0; i < factor->count; i++) {
		factor->row[i].head = factor->row[i].tail = -1;
		factor->row[i].count = 0;
		factor->column_head[i] = factor->column_tail[i] = -1;
		factor->column_count[i] = 0;
	}
	factor->count = 0;
	factor->node_count = 0;
	factor->free_node = -1;
	factor->packed_valid = 0;
}


void pan_dynamic_sparse_qr_destroy(struct pan_DynamicSparseQR *factor)
{
	lp_simplex_free(factor->row); lp_simplex_free(factor->node);
	lp_simplex_free(factor->column_head); lp_simplex_free(factor->column_tail);
	lp_simplex_free(factor->column_count); lp_simplex_free(factor->tag);
	lp_simplex_free(factor->union_column); lp_simplex_free(factor->scratch_stamp);
	lp_simplex_free(factor->rotation_first);
	lp_simplex_free(factor->rotation_second); lp_simplex_free(factor->event_row);
	lp_simplex_free(factor->event_cosine); lp_simplex_free(factor->event_sine);
	lp_simplex_free(factor->packed_start); lp_simplex_free(factor->packed_column);
	lp_simplex_free(factor->packed_value); lp_simplex_free(factor->packed_diagonal);
	lp_simplex_memset(factor, 0, sizeof(*factor));
}


int pan_dynamic_sparse_qr_append(struct pan_DynamicSparseQR *factor,
		const int tag, const double *border, const double diagonal,
		const double rank_tolerance)
{
	int i, count = factor->count;
	factor->packed_valid = 0;
	if ((count > 0 && border == NULL) || fabs(diagonal) <= rank_tolerance ||
	    pan_factor_reserve(factor, count + 1) == lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	factor->row[count].head = factor->row[count].tail = -1;
	factor->row[count].count = 0;
	factor->column_head[count] = factor->column_tail[count] = -1;
	factor->column_count[count] = 0;
	for (i = 0; i < count; i++)
		if (pan_set(factor, i, count, border[i]) == lp_simplex_EXIT_FAILURE)
			return lp_simplex_EXIT_FAILURE;
	if (pan_set(factor, count, count, diagonal) == lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	factor->tag[count] = tag; factor->count++; factor->extensions++;
	return lp_simplex_EXIT_SUCCESS;
}


int pan_dynamic_sparse_qr_load_csc(struct pan_DynamicSparseQR *factor,
		const int *tag, const int count, const int *start, const int *row,
		const double *value, const double rank_tolerance)
{
	int column, k;
	pan_dynamic_sparse_qr_clear(factor);
	if (count < 0 || pan_factor_reserve(factor, count) == lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	for (column = 0; column < count; column++) {
		factor->row[column].head = factor->row[column].tail = -1;
		factor->row[column].count = 0;
		factor->column_head[column] = factor->column_tail[column] = -1;
		factor->column_count[column] = 0;
	}
	factor->count = count;
	for (column = 0; column < count; column++)
		for (k = start[column]; k < start[column + 1]; k++)
			if (row[k] < 0 || row[k] > column ||
			    pan_set(factor, row[k], column, value[k]) ==
			    lp_simplex_EXIT_FAILURE)
				goto failure;
	for (column = 0; column < count; column++) {
		if (fabs(pan_get(factor, column, column)) <= rank_tolerance)
			goto failure;
		factor->tag[column] = tag[column];
	}
	return lp_simplex_EXIT_SUCCESS;
failure:
	pan_dynamic_sparse_qr_clear(factor);
	return lp_simplex_EXIT_FAILURE;
}


static void pan_delete_column(struct pan_DynamicSparseQR *factor,
		const int column)
{
	int node = factor->column_head[column];
	int c;
	while (node >= 0) {
		int next = factor->node[node].column_next;
		pan_node_unlink(factor, node); node = next;
	}
	for (c = column + 1; c < factor->count; c++) {
		node = factor->column_head[c];
		while (node >= 0) {
			factor->node[node].column--;
			node = factor->node[node].column_next;
		}
		factor->column_head[c - 1] = factor->column_head[c];
		factor->column_tail[c - 1] = factor->column_tail[c];
		factor->column_count[c - 1] = factor->column_count[c];
	}
	factor->column_head[factor->count - 1] = -1;
	factor->column_tail[factor->count - 1] = -1;
	factor->column_count[factor->count - 1] = 0;
}


int pan_dynamic_sparse_qr_delete(struct pan_DynamicSparseQR *factor,
		const int position, const double rank_tolerance)
{
	int i, count = factor->count;
	factor->event_count = 0;
	factor->packed_valid = 0;
	if (position < 0 || position >= count) return lp_simplex_EXIT_FAILURE;
	pan_delete_column(factor, position);
	for (i = position; i < count - 1; i++) {
		double a = pan_get(factor, i, i), b = pan_get(factor, i + 1, i);
		double radius;
		if (b == 0.) continue;
		radius = pan_stable_hypot(a, b);
		if (radius == 0. || pan_rotate_rows(factor, i, i + 1,
			a / radius, b / radius) == lp_simplex_EXIT_FAILURE ||
		    pan_set(factor, i + 1, i, 0.) == lp_simplex_EXIT_FAILURE)
			return lp_simplex_EXIT_FAILURE;
	}
	pan_clear_row(factor, count - 1);
	for (i = position; i < count - 1; i++) factor->tag[i] = factor->tag[i + 1];
	factor->count--;
	for (i = 0; i < factor->count; i++)
		if (fabs(pan_get(factor, i, i)) <= rank_tolerance)
			return lp_simplex_EXIT_FAILURE;
	factor->downdates++;
	return lp_simplex_EXIT_SUCCESS;
}


static int pan_swap_adjacent(struct pan_DynamicSparseQR *factor,
		const int position, const double rank_tolerance)
{
	int row, tag, node, touched_count = 0, i;
	int *first_node = factor->union_column;
	int *second_node = factor->union_column + factor->capacity;
	int *touched_row = factor->union_column + 2 * factor->capacity;
	double a, b, radius;
	if (factor->scratch_epoch == 2147483647) {
		lp_simplex_memset(factor->scratch_stamp, 0,
			(size_t)factor->capacity * sizeof(int));
		factor->scratch_epoch = 0;
	}
	factor->scratch_epoch++;
	node = factor->column_head[position];
	while (node >= 0) {
		row = factor->node[node].owner.row;
		factor->scratch_stamp[row] = factor->scratch_epoch;
		first_node[row] = node;
		second_node[row] = -1;
		touched_row[touched_count++] = row;
		node = factor->node[node].column_next;
	}
	node = factor->column_head[position + 1];
	while (node >= 0) {
		row = factor->node[node].owner.row;
		if (factor->scratch_stamp[row] != factor->scratch_epoch) {
			factor->scratch_stamp[row] = factor->scratch_epoch;
			first_node[row] = -1;
			touched_row[touched_count++] = row;
		}
		second_node[row] = node;
		node = factor->node[node].column_next;
	}
	for (i = 0; i < touched_count; i++) {
		int first;
		int second;
		row = touched_row[i];
		first = first_node[row];
		second = second_node[row];
		if (first >= 0 && second >= 0) {
			double value = factor->node[first].value;
			factor->node[first].value = factor->node[second].value;
			factor->node[second].value = value;
		} else if (first >= 0) {
			pan_move_node_column(factor, first, position + 1);
		} else if (second >= 0) {
			pan_move_node_column(factor, second, position);
		}
	}
	a = pan_get(factor, position, position);
	b = pan_get(factor, position + 1, position);
	radius = pan_stable_hypot(a, b);
	if (radius <= rank_tolerance || pan_rotate_rows(factor, position,
		position + 1, a / radius, b / radius) == lp_simplex_EXIT_FAILURE ||
	    pan_set(factor, position + 1, position, 0.) == lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	tag = factor->tag[position]; factor->tag[position] = factor->tag[position + 1];
	factor->tag[position + 1] = tag;
	return lp_simplex_EXIT_SUCCESS;
}


int pan_dynamic_sparse_qr_move_last(struct pan_DynamicSparseQR *factor,
		const int position, const double rank_tolerance)
{
	int i;
	factor->event_count = 0;
	factor->packed_valid = 0;
	if (position < 0 || position >= factor->count) return lp_simplex_EXIT_FAILURE;
	for (i = factor->count - 2; i >= position; i--)
		if (pan_swap_adjacent(factor, i, rank_tolerance) ==
		    lp_simplex_EXIT_FAILURE) return lp_simplex_EXIT_FAILURE;
	return lp_simplex_EXIT_SUCCESS;
}


int pan_dynamic_sparse_qr_solve_upper_transpose(
		struct pan_DynamicSparseQR *factor, double *right)
{
	int i;
	if (!factor->packed_valid &&
	    pan_pack_numeric(factor) == lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	for (i = 0; i < factor->count; i++) {
		int k;
		double diagonal = factor->packed_diagonal[i], value;
		if (diagonal == 0.) return lp_simplex_EXIT_FAILURE;
		value = right[i] / diagonal; right[i] = value;
		for (k = factor->packed_start[i]; k < factor->packed_start[i + 1]; k++)
			if (factor->packed_column[k] > i)
				right[factor->packed_column[k]] -=
					factor->packed_value[k] * value;
	}
	return lp_simplex_EXIT_SUCCESS;
}


int pan_dynamic_sparse_qr_solve_upper(
		struct pan_DynamicSparseQR *factor, double *right)
{
	int i;
	if (!factor->packed_valid &&
	    pan_pack_numeric(factor) == lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	for (i = factor->count - 1; i >= 0; i--) {
		int k;
		double diagonal = factor->packed_diagonal[i], value = right[i];
		for (k = factor->packed_start[i]; k < factor->packed_start[i + 1]; k++) {
			int column = factor->packed_column[k];
			if (column > i)
				value -= factor->packed_value[k] * right[column];
		}
		if (diagonal == 0.) return lp_simplex_EXIT_FAILURE;
		right[i] = value / diagonal;
	}
	return lp_simplex_EXIT_SUCCESS;
}


int pan_dynamic_sparse_qr_solve(
		struct pan_DynamicSparseQR *factor, double *right)
{
	if (pan_dynamic_sparse_qr_solve_upper_transpose(factor, right) ==
	    lp_simplex_EXIT_FAILURE) return lp_simplex_EXIT_FAILURE;
	return pan_dynamic_sparse_qr_solve_upper(factor, right);
}


int pan_dynamic_sparse_qr_validate(
		const struct pan_DynamicSparseQR *factor)
{
	unsigned char *seen;
	int row, column, active = 0;
	seen = (unsigned char *)lp_simplex_malloc(
		(size_t)(factor->node_count > 0 ? factor->node_count : 1));
	if (seen == NULL)
		return lp_simplex_EXIT_FAILURE;
	lp_simplex_memset(seen, 0,
		(size_t)(factor->node_count > 0 ? factor->node_count : 1));
	for (row = 0; row < factor->count; row++) {
		int node = factor->row[row].head;
		int previous = -1;
		int previous_column = -1;
		int count = 0;
		while (node >= 0) {
			if (node >= factor->node_count || seen[node] != 0 ||
			    factor->node[node].owner.row != row ||
			    factor->node[node].row_previous != previous ||
			    factor->node[node].column <= previous_column)
				goto failure;
			seen[node] = 1;
			previous_column = factor->node[node].column;
			previous = node;
			node = factor->node[node].row_next;
			count++;
		}
		if (previous != factor->row[row].tail ||
		    count != factor->row[row].count)
			goto failure;
		active += count;
	}
	for (column = 0; column < factor->count; column++) {
		int node = factor->column_head[column];
		int previous = -1;
		int count = 0;
		while (node >= 0) {
			if (node >= factor->node_count || seen[node] != 1 ||
			    factor->node[node].column != column ||
			    factor->node[node].column_previous != previous)
				goto failure;
			seen[node] = 2;
			previous = node;
			node = factor->node[node].column_next;
			count++;
		}
		if (previous != factor->column_tail[column] ||
		    count != factor->column_count[column])
			goto failure;
	}
	for (row = 0; row < factor->node_count; row++)
		if (seen[row] == 1)
			goto failure;
	lp_simplex_free(seen);
	return active >= factor->count ? lp_simplex_EXIT_SUCCESS :
		lp_simplex_EXIT_FAILURE;
failure:
	lp_simplex_free(seen);
	return lp_simplex_EXIT_FAILURE;
}
