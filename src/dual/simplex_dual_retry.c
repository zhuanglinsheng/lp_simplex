/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_dual_retry.c
 * @brief Exact Dual backend request comparison and failure admission.
 */

#include "simplex_stage.h"
#include "simplex_solve_context.h"
#include "simplex_execution.h"
#include "simplex_dual.h"
#include "simplex_basis.h"
#include "policy/controls.h"
#include <limits.h>
#include <stdlib.h>
#include <string.h>

void simplex_dual_retry_init(struct simplex_DualRetry *retry,
	const struct simplex_Problem *input,
	const struct lp_simplex_Options *options, int propagate,
	int feasibility_only, int seeded)
{
	struct simplex_RetryKey *key = &retry->key;
	struct simplex_RunConfig local;
	const struct simplex_RunConfig *controls =
		simplex_run_config_get(input->execution.context, &local);
	memset(retry, 0, sizeof(*retry));
	retry->context = input->execution.context;
	retry->history = input->execution.history;
	key->known = !seeded;
	key->input = simplex_problem_identity(input);
	key->controls = simplex_run_dual_signature(&controls->dual);
	key->model = input->execution.model;
	key->representation = SIMPLEX_DUAL_RETRY_UNSCALED;
	key->revision = input->revision;
	key->action = SIMPLEX_EXECUTION_BACKEND;
	key->algorithm = options->algorithm;
	key->pricing = options->pricing;
	key->primal_tolerance = options->primal_tolerance;
	key->dual_tolerance = options->dual_tolerance;
	key->pivot_tolerance = options->pivot_tolerance;
	key->propagate = propagate;
	key->feasibility_only = feasibility_only;
}

void simplex_dual_retry_bind_basis(struct simplex_DualRetry *retry,
	const struct simplex_DualBasisSeed *seed)
{
	struct simplex_RetryKey *key = &retry->key;
	key->known = 1;
	key->seeded = key->prefer_dual = 0;
	key->basis = NULL;
	key->status = NULL;
	if (!seed || !seed->valid)
		return;
	key->known =
		seed->rows == key->input.rows &&
		seed->columns == key->input.columns && seed->index &&
		seed->status &&
		simplex_problem_identity_equal(&seed->identity, &key->input);
	key->seeded = 1;
	key->prefer_dual = seed->prefer_dual != 0;
	key->basis = seed->index;
	key->status = seed->status;
}

static struct simplex_RetryKey dual_retry_backend(
	const struct simplex_DualRetry *retry, int backend)
{
	struct simplex_RetryKey key = retry->key;
#ifdef LP_SIMPLEX_HAVE_KLU
	key.backend = backend == SIMPLEX_BASIS_BACKEND_DEFAULT
			      ? SIMPLEX_BASIS_BACKEND_KLU
			      : backend;
#else
	key.backend = backend == SIMPLEX_BASIS_BACKEND_DEFAULT
			      ? SIMPLEX_BASIS_BACKEND_INTERNAL
			      : backend;
#endif
	return key;
}

int simplex_dual_retry_blocked(const struct simplex_DualRetry *retry,
	const struct simplex_Problem *kernel, int backend, int grant)
{
	struct simplex_RetryKey key = dual_retry_backend(retry, backend);
	struct simplex_SolveContext *context = retry->context;
	struct simplex_StageReport skipped;
	int prior, event;
	if (!context || context != kernel->execution.context ||
		!simplex_problem_identity_matches(&key.input, key.input.input))
		return 0;
	prior = simplex_retry_failed(&context->retries, &key);
	if (!prior || grant <= context->retries.records[prior - 1].work)
		/* A shorter prefix may hand off a basis before the known breakdown.
		 * Prior failure does not classify that budget-limited action. */
		return 0;
	event = simplex_execution_begin(&kernel->execution,
		SIMPLEX_EXECUTION_RECOVERY_SKIP, backend,
		SIMPLEX_REASON_REPEATED_FAILURE,
		simplex_retry_cause(context->retries.records + prior - 1,
			kernel->execution.history));
	simplex_stage_reset(&skipped);
	skipped.stop = SIMPLEX_STOP_HANDOFF;
	skipped.recovery_skipped = 1;
	/* No result or work is fabricated for an action that was not executed. */
	simplex_execution_end(kernel->execution.history, event, NULL, &skipped);
	if (kernel->report && kernel->report->recovery_skipped < INT_MAX)
		kernel->report->recovery_skipped++;
	return prior;
}

void simplex_dual_retry_record(const struct simplex_DualRetry *retry,
	int backend, int event, int state,
	const struct lp_simplex_Result *result, int grant, int basis_retry)
{
	struct simplex_RetryKey key = dual_retry_backend(retry, backend);
	if (retry->context &&
		simplex_problem_identity_matches(&key.input, key.input.input))
		simplex_retry_record(&retry->context->retries, &key,
			retry->history, event, state, result, grant,
			basis_retry);
}
