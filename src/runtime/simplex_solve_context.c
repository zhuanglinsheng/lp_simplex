/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_solve_context.c
 * @brief Solve-local model registry, identity and shared execution context.
 */

#include "simplex_solve_context.h"
#include "simplex_execution.h"
#include <limits.h>
#include <string.h>

void simplex_solve_context_destroy(struct simplex_SolveContext *context)
{
	simplex_retry_destroy(&context->retries);
}

void simplex_solve_context_configure(
	struct simplex_SolveContext *context, int experimental)
{
	simplex_run_config_init(&context->controls, experimental);
	context->configured = 1;
}

void simplex_solve_context_init(struct simplex_SolveContext *context,
	struct simplex_Execution *history, struct simplex_ExecutionScope *scope)
{
	memset(context, 0, sizeof(*context));
	context->models = 1;
	if (history)
		memset(history, 0, sizeof(*history));
	scope->context = context;
	scope->history = history;
	scope->model = 1;
	simplex_execution_set_parent(scope, 0);
}

struct simplex_ExecutionScope simplex_solve_context_derive(
	const struct simplex_ExecutionScope *source, int kind)
{
	struct simplex_ExecutionScope target = *source;
	struct simplex_SolveContext *context = source->context;
	int event;
	/* An unbound private call has no shared model namespace. Do not invent
	 * comparable identity from its dimensions or diagnostic storage. */
	target.model =
		context && context->models < INT_MAX ? ++context->models : 0;
	event = simplex_execution_begin(
		&target, kind, 0, SIMPLEX_REASON_NONE, 0);
	if (event) {
		source->history->events[event - 1].source_model = source->model;
		source->history->events[event - 1].complete = 1;
	}
	return target;
}
