/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_dual_checkpoint.c
 * @brief Dual basis checkpoint export, import and producer validation.
 */

#include "simplex_execution.h"
#include "simplex_dual.h"
#include <limits.h>
#include <string.h>

void simplex_dual_seed_begin(struct simplex_DualBasisSeed *seed)
{
	if (!seed)
		return;
	seed->exported = seed->accepted = seed->import_attempted = 0;
	seed->rejection = SIMPLEX_BASIS_NOT_REJECTED;
	if (!seed->valid) {
		memset(&seed->export_source, 0, sizeof(seed->export_source));
		seed->export_serial = 0;
	}
}

void simplex_dual_seed_export(struct simplex_DualBasisSeed *seed,
	const struct simplex_ExecutionScope *producer)
{
	int actor = simplex_execution_parent(producer);
	seed->valid = seed->exported = 1;
	seed->export_source = *producer;
	seed->export_serial =
		actor > 0 ? producer->history->events[actor - 1].serial : 0;
}

int simplex_dual_seed_source_event(const struct simplex_DualBasisSeed *seed,
	const struct simplex_ExecutionScope *observer)
{
	const struct simplex_ExecutionScope *source;
	const struct simplex_ExecutionEvent *event;
	if (!seed || !seed->valid || !seed->export_serial || !observer)
		return 0;
	source = &seed->export_source;
	if (!source->context || source->context != observer->context ||
		!source->history || source->history != observer->history ||
		!simplex_execution_parent(source))
		return 0;
	event = source->history->events + source->parent - 1;
	return event->context == source->context &&
			       event->serial == seed->export_serial &&
			       event->complete && event->result_available &&
			       event->basis_exported &&
			       event->kind == SIMPLEX_EXECUTION_BACKEND &&
			       event->model == source->model &&
			       (event->status == lp_simplex_Success ||
				       event->status ==
					       lp_simplex_ExceedIterLimit)
		       ? source->parent
		       : 0;
}

void simplex_dual_seed_init(struct simplex_DualBasisSeed *seed,
	const struct simplex_Problem *problem, int *index,
	unsigned char *status)
{
	memset(seed, 0, sizeof(*seed));
	seed->identity = simplex_problem_identity(problem);
	seed->rows = problem->rows;
	seed->columns = problem->columns;
	seed->index = index;
	seed->status = status;
	seed->prefer_dual = 1;
}

void simplex_dual_checkpoint_init(struct simplex_DualCheckpoint *checkpoint,
	const struct simplex_Problem *problem, int *index,
	unsigned char *status)
{
	simplex_dual_seed_init(&checkpoint->basis, problem, index, status);
}

int simplex_dual_checkpoint_rejection(
	const struct simplex_DualCheckpoint *checkpoint,
	const struct simplex_Problem *problem)
{
	if (!checkpoint || !simplex_problem_identity_matches(
				   &checkpoint->basis.identity, problem))
		return SIMPLEX_BASIS_MODEL_REJECTED;
	if (checkpoint->basis.rows != problem->rows ||
		checkpoint->basis.columns != problem->columns ||
		!checkpoint->basis.index || !checkpoint->basis.status)
		return SIMPLEX_BASIS_LAYOUT_REJECTED;
	return SIMPLEX_BASIS_NOT_REJECTED;
}

int simplex_dual_checkpoint_matches(
	const struct simplex_DualCheckpoint *checkpoint,
	const struct simplex_Problem *problem)
{
	return simplex_dual_checkpoint_rejection(checkpoint, problem) ==
	       SIMPLEX_BASIS_NOT_REJECTED;
}

int simplex_dual_seed_retry_allowed(const struct simplex_DualBasisSeed *seed)
{
	return seed && seed->import_attempted &&
	       (seed->accepted ||
		       seed->rejection ==
			       SIMPLEX_BASIS_FACTORIZATION_REJECTED ||
		       seed->rejection == SIMPLEX_BASIS_RESIDUAL_REJECTED);
}
