/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_dual_backend.c
 * @brief Dual backend attempts, bounded numerical recovery and point selection.
 */

#include "simplex_execution.h"
#include "simplex_dual_backend.h"
#include "simplex_dual_kernel.h"
#include "simplex_dual.h"
#include "simplex_basis.h"
#include "simplex_solve_context.h"
#include "simplex_stage.h"
#include "policy/controls.h"
#include "utils.h"

#include <lp_simplex/status.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* Retry a default-backend numerical breakdown once with internal sparse LU,
 * using the same model and the remaining grant. Different pivoting may help;
 * no fallback is authorized by a budget exit or a mathematical classification.
 * Point selection assembles both attempts without performing another solve. */
int simplex_dual_backend_solve(const struct simplex_Problem *problem,
	const struct lp_simplex_Options *options, double *x, double *dual,
	struct lp_simplex_Result *result, const int propagate_bounds,
	const int columns, const int rows, const int feasibility_only,
	struct simplex_DualBasisSeed *seed,
	const struct simplex_DualRetry *retry)
{
	struct simplex_Problem attempt = *problem;
	struct simplex_DualRetry invocation = *retry;
	int event, state, j, prior, basis_retry = 0;

#ifdef LP_SIMPLEX_HAVE_KLU
	struct lp_simplex_Result first;
	struct simplex_StageReport first_report;
	struct lp_simplex_Options retry_options;
	struct simplex_DualBasisSeed next_seed;
	struct simplex_ExecutionRef first_actor;
	int retry_event, selection, chosen_source, retry_iterations,
		had_seed = seed && seed->valid;
	double *saved;
#endif

	simplex_dual_retry_bind_basis(&invocation, seed);
	prior = simplex_dual_retry_blocked(&invocation, problem,
		SIMPLEX_BASIS_BACKEND_DEFAULT, options->iteration_limit);
	if (prior) {
		const struct simplex_RetryRecord *failure =
			retry->context->retries.records + prior - 1;
		/* No backend ran and no point, multiplier or basis was produced.
		 * Carry only the known numerical failure to the enclosing adapter;
		 * a different backend remains eligible below. The skip event has no
		 * fabricated kernel result. */
		event = simplex_retry_cause(
			failure, problem->execution.history);
		basis_retry = failure->basis_retry;
		result->status = failure->status;
		result->iterations = 0;
		result->objective = NAN;
		result->primal_infeasibility = HUGE_VAL;
		result->dual_infeasibility = NAN;
		result->dual_residual_space = lp_simplex_RESIDUAL_UNAVAILABLE;
		if (seed) {
			seed->valid = seed->exported = seed->accepted =
				seed->import_attempted = 0;
			seed->rejection = SIMPLEX_BASIS_NOT_REJECTED;
		}
		for (j = 0; j < columns; j++)
			x[j] = NAN;
		if (dual)
			for (j = 0; j < rows; j++)
				dual[j] = NAN;
		if (problem->report)
			problem->report->stop = SIMPLEX_STOP_HANDOFF;
		state = lp_simplex_EXIT_FAILURE;
	} else {
		event = simplex_execution_begin(&problem->execution,
			SIMPLEX_EXECUTION_BACKEND,
			SIMPLEX_BASIS_BACKEND_DEFAULT, SIMPLEX_REASON_NONE, 0);
		if (event)
			problem->execution.history->events[event - 1]
				.basis_source = simplex_dual_seed_source_event(
				seed, &problem->execution);
		simplex_execution_set_parent(&attempt.execution, event);
		state = simplex_dual_kernel_solve(&attempt, options, x, dual,
			result, propagate_bounds, SIMPLEX_BASIS_BACKEND_DEFAULT,
			feasibility_only, seed);
		event = simplex_execution_parent(&attempt.execution);
		simplex_execution_end(problem->execution.history, event, result,
			problem->report);
		if (event)
			problem->execution.history->events[event - 1]
				.basis_exported =
				seed && seed->valid && seed->exported;
		basis_retry = simplex_dual_seed_retry_allowed(seed);
		simplex_dual_retry_record(&invocation,
			SIMPLEX_BASIS_BACKEND_DEFAULT, event, state, result,
			options->iteration_limit, basis_retry);
	}
#ifdef LP_SIMPLEX_HAVE_KLU
	first_actor = simplex_execution_ref(&problem->execution, event);
	if (state != lp_simplex_EXIT_FAILURE ||
		(result->status != lp_simplex_PrecisionError &&
			result->status != lp_simplex_Singularity))
		return state;
	first = *result;
	if (problem->report)
		first_report = *problem->report;
	/* The caller's pivot budget covers both attempts. */
	if (first.iterations >= options->iteration_limit)
		return state;
	/* The next backend receives either this exact hint or a cold request.
	 * Do not compare a cold fallback under the failed hint's key. */
	if (seed) {
		next_seed = *seed;
		next_seed.valid = had_seed && basis_retry;
	}
	invocation = *retry;
	simplex_dual_retry_bind_basis(&invocation, seed ? &next_seed : NULL);
	if (simplex_dual_retry_blocked(&invocation, problem,
		    SIMPLEX_BASIS_BACKEND_INTERNAL,
		    options->iteration_limit - first.iterations))
		return state;
	saved = (double *)lp_simplex_malloc(
		((size_t)columns + (dual != NULL ? (size_t)rows : 0)) *
		sizeof(double));
	if (saved == NULL)
		return state;
	lp_simplex_memcpy(saved, x, (size_t)columns * sizeof(double));
	lp_simplex_memset(x, 0, (size_t)columns * sizeof(double));
	if (dual != NULL) {
		lp_simplex_memcpy(
			saved + columns, dual, (size_t)rows * sizeof(double));
		lp_simplex_memset(dual, 0, (size_t)rows * sizeof(double));
	}
	retry_options = *options;
	retry_options.iteration_limit -= first.iterations;
	lp_simplex_memset(result, 0, sizeof(*result));
	/* Retain the original hint only when an accepted import or a numerical
	 * import failure supports another factorization. A duplicate index, bad
	 * status or classified feasibility failure supplies no such evidence. */
	if (seed)
		seed->valid = next_seed.valid;
	if (problem->report)
		simplex_stage_reset(problem->report);
	retry_event = simplex_execution_begin(&problem->execution,
		SIMPLEX_EXECUTION_BACKEND, SIMPLEX_BASIS_BACKEND_INTERNAL,
		SIMPLEX_REASON_NUMERICAL_FAILURE, event);
	if (retry_event)
		problem->execution.history->events[retry_event - 1]
			.basis_source = simplex_dual_seed_source_event(
			seed, &problem->execution);
	simplex_execution_set_parent(&attempt.execution, retry_event);
	state = simplex_dual_kernel_solve(&attempt, &retry_options, x, dual,
		result, propagate_bounds, SIMPLEX_BASIS_BACKEND_INTERNAL,
		feasibility_only, seed);
	retry_event = simplex_execution_parent(&attempt.execution);
	simplex_execution_end(problem->execution.history, retry_event, result,
		problem->report);
	if (retry_event)
		problem->execution.history->events[retry_event - 1]
			.basis_exported = seed && seed->valid && seed->exported;
	simplex_dual_retry_record(&invocation, SIMPLEX_BASIS_BACKEND_INTERNAL,
		retry_event, state, result, retry_options.iteration_limit,
		simplex_dual_seed_retry_allowed(seed));
	if (problem->report) {
		problem->report->recovery_skipped +=
			first_report.recovery_skipped;
		problem->report->backend_retries =
			first_report.backend_retries + 1;
		problem->report->work.factorizations +=
			first_report.work.factorizations;
		problem->report->work.ftran += first_report.work.ftran;
		problem->report->work.btran += first_report.work.btran;
		problem->report->work.cpu_seconds +=
			first_report.work.cpu_seconds;
	}
	if (state != lp_simplex_EXIT_SUCCESS &&
		result->status != lp_simplex_Infeasibility &&
		result->status != lp_simplex_Unboundedness) {
		retry_iterations = result->iterations;
		selection =
			simplex_stage_prefer_point(problem, options, x, saved);
		event = simplex_execution_ref_id(
			&first_actor, &problem->execution, 1);
		retry_event = simplex_execution_parent(&attempt.execution);
		chosen_source = selection ? retry_event : event;
		if (!selection) {
			struct simplex_StageWork total_work;
			if (problem->report) {
				total_work = problem->report->work;
				*problem->report = first_report;
				problem->report->work = total_work;
				problem->report->backend_retries++;
			}
			*result = first;
			lp_simplex_memcpy(
				x, saved, (size_t)columns * sizeof(double));
			if (dual != NULL)
				lp_simplex_memcpy(dual, saved + columns,
					(size_t)rows * sizeof(double));
		}
		/* This assembles the whole Dual call. It is a sibling of the two
		 * local attempts, with separate causal and point-producing sources. */
		result->iterations = first.iterations + retry_iterations;
		selection = simplex_execution_begin(&problem->execution,
			SIMPLEX_EXECUTION_POINT_SELECTION, 0,
			selection ? SIMPLEX_REASON_KEEP_RETRY_POINT
				  : SIMPLEX_REASON_RESTORE_SAVED_POINT,
			retry_event);
		if (selection)
			problem->execution.history->events[selection - 1]
				.point_source = chosen_source;
		simplex_execution_end(problem->execution.history, selection,
			result, problem->report);
		lp_simplex_free(saved);
		return lp_simplex_EXIT_FAILURE;
	}
	lp_simplex_free(saved);
	result->iterations += first.iterations;
	if (getenv(SIMPLEX_ENV_PROFILE) != NULL)
		fprintf(stderr,
			"dual profile: recovered on the internal "
			"factorization after %d iterations\n",
			first.iterations);
#else
	(void)columns;
	(void)rows;
#endif
	return state;
}
