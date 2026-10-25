/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_retry.c
 * @brief Exact failed-start records independent of diagnostic history.
 */

#include "simplex_retry.h"
#include "simplex_execution.h"
#include "policy/resources.h"
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

static const struct simplex_ExecutionEvent *retry_actor(
	const struct simplex_Execution *history, int event)
{
	return history && history->count >= 0 &&
			       history->count <= SIMPLEX_EXECUTION_CAPACITY &&
			       event > 0 && event <= history->count
		       ? history->events + event - 1
		       : NULL;
}

/* Check against the fixed cap before multiplying or adding dimensions. */
static size_t retry_basis_bytes(const struct simplex_RetryKey *k)
{
	size_t statuses, rows;
	if (k->input.rows < 0 || k->input.columns < 0)
		return 0;
	rows = (size_t)k->input.rows;
	if (rows > SIMPLEX_RETRY_BASIS_BYTES ||
		(size_t)k->input.columns > SIMPLEX_RETRY_BASIS_BYTES - rows)
		return 0;
	statuses = rows + (size_t)k->input.columns;
	if (rows > (SIMPLEX_RETRY_BASIS_BYTES - statuses) / sizeof(int))
		return 0;
	return rows * sizeof(int) + statuses;
}

static int retry_key_known(const struct simplex_RetryKey *k)
{
	return k && k->known && k->model > 0 && k->action > 0 &&
	       k->backend > 0 && k->input.input && k->input.context &&
	       k->input.model == k->model && k->input.revision == k->revision &&
	       (!k->seeded ||
		       (k->basis && k->status && retry_basis_bytes(k))) &&
	       k->representation >= 0 && k->revision != ULONG_MAX &&
	       isfinite(k->primal_tolerance) && k->primal_tolerance > 0. &&
	       isfinite(k->dual_tolerance) && k->dual_tolerance > 0. &&
	       isfinite(k->pivot_tolerance) && k->pivot_tolerance > 0.;
}

static int retry_key_equal(
	const struct simplex_RetryKey *a, const struct simplex_RetryKey *b)
{
	return simplex_problem_identity_equal(&a->input, &b->input) &&
	       a->model == b->model && a->revision == b->revision &&
	       a->action == b->action && a->backend == b->backend &&
	       a->representation == b->representation &&
	       a->controls == b->controls && a->algorithm == b->algorithm &&
	       a->pricing == b->pricing && a->propagate == b->propagate &&
	       a->feasibility_only == b->feasibility_only &&
	       a->primal_tolerance == b->primal_tolerance &&
	       a->dual_tolerance == b->dual_tolerance &&
	       a->pivot_tolerance == b->pivot_tolerance &&
	       a->seeded == b->seeded &&
	       (!a->seeded ||
		       (a->prefer_dual == b->prefer_dual &&
			       !memcmp(a->basis, b->basis,
				       (size_t)a->input.rows * sizeof(int)) &&
			       !memcmp(a->status, b->status,
				       (size_t)a->input.rows +
					       a->input.columns)));
}

int simplex_retry_failed(const struct simplex_RetryLedger *ledger,
	const struct simplex_RetryKey *key)
{
	int i;
	if (!ledger || !retry_key_known(key))
		return 0;
	for (i = 0; i < ledger->count; i++)
		if (retry_key_equal(&ledger->records[i].key, key))
			return i + 1;
	return 0;
}

void simplex_retry_record(struct simplex_RetryLedger *ledger,
	const struct simplex_RetryKey *key,
	const struct simplex_Execution *history, int event, int state,
	const struct lp_simplex_Result *result, int grant, int basis_retry)
{
	struct simplex_RetryRecord *record;
	const struct simplex_ExecutionEvent *actor;
	size_t bytes;
	void *storage = NULL;
	if (!ledger || !retry_key_known(key) || !result ||
		state != lp_simplex_EXIT_FAILURE || result->iterations < 0 ||
		result->iterations >= grant ||
		(result->status != lp_simplex_PrecisionError &&
			result->status != lp_simplex_Singularity) ||
		simplex_retry_failed(ledger, key))
		return;
	bytes = key->seeded ? retry_basis_bytes(key) : 0;
	if (ledger->count == SIMPLEX_RETRY_CAPACITY ||
		bytes > SIMPLEX_RETRY_BASIS_BYTES - ledger->basis_bytes ||
		(bytes && !(storage = malloc(bytes)))) {
		/* Missing history is not permission to invent a previous failure.
		 * Existing kernel/route limits still bound untracked recovery. */
		if (ledger->unrecorded < INT_MAX)
			ledger->unrecorded++;
		return;
	}
	record = ledger->records + ledger->count++;
	record->key = *key;
	record->basis_storage = storage;
	if (bytes) {
		memcpy(storage, key->basis,
			(size_t)key->input.rows * sizeof(int));
		record->key.basis = (const int *)storage;
		record->key.status = (const unsigned char *)storage +
				     (size_t)key->input.rows * sizeof(int);
		memcpy((void *)record->key.status, key->status,
			(size_t)key->input.rows + key->input.columns);
		ledger->basis_bytes += bytes;
	}
	record->history = history;
	record->event = event;
	record->status = result->status;
	record->work = result->iterations;
	record->basis_retry = basis_retry != 0;
	/* Trace validation cannot grant or erase a completed numerical failure. */
	record->event_serial = 0;
	record->event_model = record->event_detail = 0;
	actor = retry_actor(history, event);
	if (actor && actor->serial && actor->context == key->input.context &&
		actor->complete && actor->result_available &&
		actor->kind == key->action && actor->status == result->status &&
		actor->work == result->iterations) {
		record->event_serial = actor->serial;
		record->event_model = actor->model;
		record->event_detail = actor->detail;
	}
}

void simplex_retry_destroy(struct simplex_RetryLedger *ledger)
{
	int i;
	for (i = 0; i < ledger->count; i++)
		free(ledger->records[i].basis_storage);
	memset(ledger, 0, sizeof(*ledger));
}

int simplex_retry_cause(const struct simplex_RetryRecord *record,
	const struct simplex_Execution *history)
{
	const struct simplex_ExecutionEvent *actor;
	if (!record || !history || record->history != history ||
		!record->event_serial)
		return 0;
	actor = retry_actor(history, record->event);
	return actor && actor->context == record->key.input.context &&
			       actor->serial == record->event_serial &&
			       actor->complete && actor->result_available &&
			       actor->kind == record->key.action &&
			       actor->model == record->event_model &&
			       actor->detail == record->event_detail &&
			       actor->status == record->status &&
			       actor->work == record->work
		       ? record->event
		       : 0;
}
