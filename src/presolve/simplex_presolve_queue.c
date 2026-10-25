/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_presolve_queue.c
 * @brief Bounded work queues with duplicate suppression for presolve.
 */

#include "simplex_presolve_queue.h"
#include "utils.h"
#include "simplex_storage.h"

#include <lp_simplex/status.h>

int simplex_presolve_queue_init(
	struct simplex_PresolveQueue *queue, const int capacity)
{
	size_t item_bytes, present_bytes;

	lp_simplex_memset(queue, 0, sizeof(*queue));
	if (capacity < 0 ||
		!simplex_storage_bytes((size_t)capacity, sizeof(int),
			(size_t)-1, &item_bytes) ||
		!simplex_storage_bytes((size_t)capacity, sizeof(unsigned char),
			(size_t)-1, &present_bytes))
		return lp_simplex_EXIT_FAILURE;
	if (!capacity)
		return lp_simplex_EXIT_SUCCESS;
	queue->capacity = capacity;
	queue->item = (int *)lp_simplex_malloc(item_bytes);
	queue->present = (unsigned char *)lp_simplex_malloc(present_bytes);
	if (queue->item == NULL || queue->present == NULL) {
		simplex_presolve_queue_destroy(queue);
		return lp_simplex_EXIT_FAILURE;
	}
	lp_simplex_memset(queue->present, 0, present_bytes);
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
	/* Both summands can be valid indices while their sum exceeds INT_MAX. */
	if (queue->count >= queue->capacity - queue->head)
		tail = queue->count - (queue->capacity - queue->head);
	else
		tail = queue->head + queue->count;
	queue->item[tail] = item;
	queue->present[item] = 1;
	queue->count++;
	return lp_simplex_EXIT_SUCCESS;
}

int simplex_presolve_queue_pop(struct simplex_PresolveQueue *queue, int *item)
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
