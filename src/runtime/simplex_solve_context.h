/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_solve_context.h
 * @brief Solve-local model registry, identity and shared execution context.
 */

#ifndef LP_SIMPLEX_SOLVE_CONTEXT_INTERNAL_H
#define LP_SIMPLEX_SOLVE_CONTEXT_INTERNAL_H

#include "simplex_retry.h"
#include "simplex_run_config.h"

/* Numerical decision state shared by one live solve. Its lifetime and model
 * identities do not depend on whether an output journal exists or has room. */
struct simplex_SolveContext {
	int models;
	/* Diagnostic actor tokens survive journal resets; saturation loses only provenance. */
	unsigned long execution_serial;
	struct simplex_RetryLedger retries;
	int configured;
	struct simplex_RunConfig controls;
};

void simplex_solve_context_configure(
	struct simplex_SolveContext *context, int experimental);
/* Release owned failure snapshots before reinitializing or ending the owner.
 * Borrowed model/basis/journal storage is never freed here. */
void simplex_solve_context_destroy(struct simplex_SolveContext *context);
struct simplex_Execution;
struct simplex_ExecutionScope;
void simplex_solve_context_init(struct simplex_SolveContext *context,
	struct simplex_Execution *history,
	struct simplex_ExecutionScope *scope);
struct simplex_ExecutionScope simplex_solve_context_derive(
	const struct simplex_ExecutionScope *source, int kind);
#endif
