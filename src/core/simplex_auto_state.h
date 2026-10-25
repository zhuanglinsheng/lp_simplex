/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_auto_state.h
 * @brief Controller-owned AUTO runtime, candidates and attempt outputs.
 */

#ifndef LP_SIMPLEX_AUTO_STATE_INTERNAL_H
#define LP_SIMPLEX_AUTO_STATE_INTERNAL_H

#include "simplex_auto_policy.h"
#include "simplex_dual.h"

struct simplex_AutoCandidate {
	double *point;
	int valid;
	double objective;
	struct lp_simplex_Result result;
	struct simplex_ProblemIdentity identity;
	struct simplex_AutoOrigin origin;
};

/* One invocation's raw output, accounting/audit result and independent audit.
 * The raw kernel claim is preserved even when the controller rejects it. */
struct simplex_AutoAttempt {
	struct lp_simplex_Result result, raw_result;
	struct simplex_StageReport report, audit;
	struct simplex_ProblemIdentity identity;
	int state, feasible, terminal, work_valid;
	double objective;
};

struct simplex_AutoRuntime {
	struct simplex_ExecutionRef parent_actor;
	int used, cursor, route, recovery_count;
	int allowed[SIMPLEX_AUTO_ROUTE_COUNT];
	struct simplex_AutoCandidate best;
	struct simplex_DualCheckpoint checkpoint;
	struct simplex_AutoPolicyEvidence evidence;
	struct simplex_ProblemIdentity input;
};

#endif
