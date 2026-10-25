/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_execution.h
 * @brief Bounded execution journal and validated causal references.
 */

#ifndef LP_SIMPLEX_EXECUTION_INTERNAL_H
#define LP_SIMPLEX_EXECUTION_INTERNAL_H

#include "simplex_execution_scope.h"
#include "simplex_stage.h"
#include "simplex_solve_context.h"

/* Bounded, solve-local causal history. No allocation, callbacks or wall-clock
 * decisions. Exhaustion loses diagnostics explicitly, never numerical state.
 * Work is inclusive within an event; nested event counts must not be summed. */
#define SIMPLEX_EXECUTION_CAPACITY 128

enum simplex_ExecutionKind {
	SIMPLEX_EXECUTION_SOLVE = 1,
	SIMPLEX_EXECUTION_PRESOLVE,
	SIMPLEX_EXECUTION_REDUCED_SOLVE,
	SIMPLEX_EXECUTION_SCALING,
	SIMPLEX_EXECUTION_AUTO,
	SIMPLEX_EXECUTION_BACKEND,
	SIMPLEX_EXECUTION_STABLE_INITIALIZE,
	SIMPLEX_EXECUTION_POSTSOLVE,
	SIMPLEX_EXECUTION_ORIGINAL_RETRY,
	SIMPLEX_EXECUTION_POINT_SELECTION,
	SIMPLEX_EXECUTION_RECOVERY_SKIP,
	SIMPLEX_EXECUTION_PAN_SERVICE,
	SIMPLEX_EXECUTION_PAN_MODEL,
	SIMPLEX_EXECUTION_COLUMN_MODEL,
	SIMPLEX_EXECUTION_COVER_MASTER,
	SIMPLEX_EXECUTION_COLUMN_HANDOFF,
	SIMPLEX_EXECUTION_AUTO_FINISH,
	SIMPLEX_EXECUTION_POINT_ARCHIVE,
	SIMPLEX_EXECUTION_RESULT_AUDIT,
	SIMPLEX_EXECUTION_TABLEAU,
	SIMPLEX_EXECUTION_TABLEAU_MODEL,
	SIMPLEX_EXECUTION_TABLEAU_STANDARD,
	SIMPLEX_EXECUTION_SINGLETON,
	SIMPLEX_EXECUTION_SINGLETON_MODEL,
	SIMPLEX_EXECUTION_PAN_TRIGGER,
	SIMPLEX_EXECUTION_PAN_FACE_OPPORTUNITY
};

enum simplex_ExecutionReason {
	SIMPLEX_REASON_NONE = 0,
	SIMPLEX_REASON_NUMERICAL_FAILURE,
	SIMPLEX_REASON_INVALID_POSTSOLVE,
	SIMPLEX_REASON_STARTUP_FAILURE,
	SIMPLEX_REASON_FEASIBLE_CANDIDATE,
	SIMPLEX_REASON_KEEP_RETRY_POINT,
	SIMPLEX_REASON_RESTORE_SAVED_POINT,
	SIMPLEX_REASON_RESTORE_ARCHIVED_POINT,
	SIMPLEX_REASON_REPEATED_FAILURE,
	SIMPLEX_REASON_DUAL_HANDOFF,
	SIMPLEX_REASON_IMPLICIT_ORACLE_FAILURE,
	SIMPLEX_REASON_BASIS_REBUILD,
	SIMPLEX_REASON_ALTERNATIVE_REPRESENTATION,
	SIMPLEX_REASON_RECOVERY_EXHAUSTED,
	SIMPLEX_REASON_COVER_UNCERTIFIED,
	SIMPLEX_REASON_ALTERNATIVE_ALGORITHM,
	SIMPLEX_REASON_STAGE_CONTRACT,
	SIMPLEX_REASON_COLUMN_ADMISSION,
	SIMPLEX_REASON_DUAL_REFINEMENT,
	SIMPLEX_REASON_TERMINAL_CERTIFICATE,
	SIMPLEX_REASON_TOTAL_BUDGET,
	SIMPLEX_REASON_ROUTES_EXHAUSTED,
	SIMPLEX_REASON_INPUT_CHANGED,
	SIMPLEX_REASON_REQUIRED_WORKSPACE,
	SIMPLEX_REASON_FINAL_AUDIT,
	SIMPLEX_REASON_INFEASIBILITY_CONFLICT,
	SIMPLEX_REASON_REQUIRED_CAPABILITY,
	SIMPLEX_REASON_ALTERNATIVE_DIRECTION,
	SIMPLEX_REASON_PDHG_STAGNATION,
	SIMPLEX_REASON_FACE_INITIAL,
	SIMPLEX_REASON_FACE_LOW_PROGRESS,
	SIMPLEX_REASON_FACE_MULTIPLIER_SEARCH,
	SIMPLEX_REASON_PRIMAL_MARGIN
};

struct simplex_ExecutionEvent {
	struct simplex_PanFaceReport face;
	struct simplex_SingletonReport singleton;
	struct simplex_TableauReport tableau;
	/* Token namespaces belong to live contexts, not reusable journal buffers. */
	const struct simplex_SolveContext *context;
	int kind, detail, reason;
	int model, source_model, parent, cause;
	int point_source, basis_source, grant;
	int prerequisite, purpose, on_failure;
	int complete, result_available, status, work, stop, failure;
	int basis_import_attempted, basis_rejection, basis_exported;
	int recovery_skipped;
	int service_output, service_code, service_valid, service_raw_stop;
	int service_reported_end, service_reported_work, service_relative_work;
	int service_implementation, model_rows, model_columns;
	int handoff_rejection;
	int experimental_environment;
	int observation_work;
	unsigned long serial;
	char observation_stage[48];
	long bounds_epoch, cost_epoch;
	struct simplex_Observation primal, primal_audit;
};

struct simplex_Execution {
	struct simplex_ExecutionEvent events[SIMPLEX_EXECUTION_CAPACITY];
	int count, dropped;
};

void simplex_execution_set_parent(struct simplex_ExecutionScope *scope, int id);
int simplex_execution_parent(const struct simplex_ExecutionScope *scope);
void simplex_execution_inherit_parent(struct simplex_ExecutionScope *target,
	const struct simplex_ExecutionScope *source);
struct simplex_ExecutionRef simplex_execution_ref(
	const struct simplex_ExecutionScope *scope, int id);
int simplex_execution_ref_id(const struct simplex_ExecutionRef *reference,
	const struct simplex_ExecutionScope *scope, int completed);
int simplex_execution_begin(const struct simplex_ExecutionScope *scope,
	int kind, int detail, int reason, int cause);
void simplex_execution_end(struct simplex_Execution *history, int id,
	const struct lp_simplex_Result *result,
	const struct simplex_StageReport *report);
void simplex_execution_print(const struct simplex_Execution *history);
/* Internal inspection entry point; public solve API remains unchanged. */
struct lp_Model;
int simplex_solve_reported(const struct lp_Model *model,
	const struct lp_simplex_Options *options, double *x,
	struct lp_simplex_Result *result, struct simplex_Execution *history);
int simplex_solve_reported_experimental(const struct lp_Model *model,
	const struct lp_simplex_Options *options, double *x,
	struct lp_simplex_Result *result, struct simplex_Execution *history);
#endif
