/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_retry.h
 * @brief Exact failed-start records independent of diagnostic history.
 */

#ifndef LP_SIMPLEX_RETRY_INTERNAL_H
#define LP_SIMPLEX_RETRY_INTERNAL_H

#include <lp_simplex/solve.h>
#include "simplex_problem_identity.h"
#include <stddef.h>

/* This ledger is numerical decision state, independent of the bounded trace.
 * Only exact, known starts are comparable. No vector fingerprints, CPU
 * times or weighted progress scores participate in admission. */
#define SIMPLEX_RETRY_CAPACITY 32
struct simplex_Execution;

struct simplex_RetryKey {
	struct simplex_ProblemIdentity input;
	int known, model, action, backend, representation;
	int algorithm, pricing, propagate, feasibility_only;
	unsigned long revision, controls;
	/* Borrowed during admission; records own exact copies. Cold starts have
	 * seeded=0, no arrays. Only the consumed truth value of prefer_dual matters. */
	int seeded, prefer_dual;
	const int *basis;
	const unsigned char *status;
	double primal_tolerance, dual_tolerance, pivot_tolerance;
};

struct simplex_RetryRecord {
	struct simplex_RetryKey key;
	/* Diagnostic ownership only; never part of numerical equality. */
	const struct simplex_Execution *history;
	int event, status, work, basis_retry;
	unsigned long event_serial;
	int event_model, event_detail;
	void *basis_storage;
};

struct simplex_RetryLedger {
	struct simplex_RetryRecord records[SIMPLEX_RETRY_CAPACITY];
	int count, unrecorded;
	size_t basis_bytes;
};

/* Returns the matching record plus one; zero means no known prior failure.
 * Event ids may be unavailable after trace truncation and are never keys. */
int simplex_retry_failed(const struct simplex_RetryLedger *ledger,
	const struct simplex_RetryKey *key);
void simplex_retry_record(struct simplex_RetryLedger *ledger,
	const struct simplex_RetryKey *key,
	const struct simplex_Execution *history, int event, int state,
	const struct lp_simplex_Result *result, int grant, int basis_retry);
void simplex_retry_destroy(struct simplex_RetryLedger *ledger);
/* References match a completed actor's context, token and failure payload.
 * Missing/reset/replaced journals lose provenance, never numerical evidence. */
int simplex_retry_cause(const struct simplex_RetryRecord *record,
	const struct simplex_Execution *history);
#endif
