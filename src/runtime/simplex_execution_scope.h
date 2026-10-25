/**
 * @file simplex_execution_scope.h
 * @brief Borrowed execution context and causal references for model scopes.
 */

#ifndef LP_SIMPLEX_EXECUTION_SCOPE_INTERNAL_H
#define LP_SIMPLEX_EXECUTION_SCOPE_INTERNAL_H

struct simplex_Execution;
struct simplex_SolveContext;

/* A diagnostic actor reference, never a model/state identity. The token
 * survives index reuse in a live context; its journal must remain alive. */
struct simplex_ExecutionRef {
	const struct simplex_Execution *history;
	const struct simplex_SolveContext *context;
	unsigned long serial;
	int id, kind, detail, model;
};

/* Parent provenance is bound when assigned, not rediscovered from an index
 * after a callback. Neither field grants model or numerical permissions. */
struct simplex_ExecutionScope {
	struct simplex_Execution *history;
	struct simplex_SolveContext *context;
	int model, parent;
	struct simplex_ExecutionRef parent_actor;
};

#endif
