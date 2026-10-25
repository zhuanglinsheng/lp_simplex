/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_presolve_queue.h
 * @brief Bounded work queues with duplicate suppression for presolve.
 */

#ifndef LP_SIMPLEX_PRESOLVE_QUEUE_INTERNAL_H
#define LP_SIMPLEX_PRESOLVE_QUEUE_INTERNAL_H

struct simplex_PresolveQueue {
	int *item;
	unsigned char *present;
	int capacity;
	int head;
	int count;
};

/* Fresh output; zero capacity is empty, negative or unrepresentable capacity
 * is refused. Failed allocation leaves a cleared, safely destructible queue. */
int simplex_presolve_queue_init(
	struct simplex_PresolveQueue *queue, int capacity);

void simplex_presolve_queue_clear(struct simplex_PresolveQueue *queue);

int simplex_presolve_queue_push(struct simplex_PresolveQueue *queue, int item);

int simplex_presolve_queue_pop(struct simplex_PresolveQueue *queue, int *item);

void simplex_presolve_queue_destroy(struct simplex_PresolveQueue *queue);

#endif
