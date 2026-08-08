/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#include "simplex_presolve_queue.h"
#include "utils.h"

#include <lp_simplex/status.h>


int simplex_presolve_queue_init(
		struct simplex_PresolveQueue *queue, const int capacity)
{
	queue->item = NULL;
	queue->present = NULL;
	queue->capacity = capacity;
	queue->head = 0;
	queue->count = 0;
	if (capacity <= 0)
		return lp_simplex_EXIT_SUCCESS;
	queue->item = (int *)lp_simplex_malloc(
		(size_t)capacity * sizeof(int));
	queue->present = (unsigned char *)lp_simplex_malloc(
		(size_t)capacity * sizeof(unsigned char));
	if (queue->item == NULL || queue->present == NULL) {
		simplex_presolve_queue_destroy(queue);
		return lp_simplex_EXIT_FAILURE;
	}
	lp_simplex_memset(queue->present, 0,
		(size_t)capacity * sizeof(unsigned char));
	return lp_simplex_EXIT_SUCCESS;
}


void simplex_presolve_queue_clear(struct simplex_PresolveQueue *queue)
{
	queue->head = 0;
	queue->count = 0;
	if (queue->present != NULL)
		lp_simplex_memset(queue->present, 0,
			(size_t)queue->capacity * sizeof(unsigned char));
}


int simplex_presolve_queue_push(
		struct simplex_PresolveQueue *queue, const int item)
{
	int tail;
	if (item < 0 || item >= queue->capacity)
		return lp_simplex_EXIT_FAILURE;
	if (queue->present[item])
		return lp_simplex_EXIT_SUCCESS;
	if (queue->count == queue->capacity)
		return lp_simplex_EXIT_FAILURE;
	tail = (queue->head + queue->count) % queue->capacity;
	queue->item[tail] = item;
	queue->present[item] = 1;
	queue->count++;
	return lp_simplex_EXIT_SUCCESS;
}


int simplex_presolve_queue_pop(
		struct simplex_PresolveQueue *queue, int *item)
{
	if (queue->count == 0)
		return 0;
	*item = queue->item[queue->head];
	queue->head = (queue->head + 1) % queue->capacity;
	queue->count--;
	queue->present[*item] = 0;
	return 1;
}


void simplex_presolve_queue_destroy(struct simplex_PresolveQueue *queue)
{
	if (queue == NULL)
		return;
	lp_simplex_free(queue->item);
	lp_simplex_free(queue->present);
	queue->item = NULL;
	queue->present = NULL;
	queue->capacity = 0;
	queue->head = 0;
	queue->count = 0;
}
