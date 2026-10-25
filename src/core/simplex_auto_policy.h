/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_auto_policy.h
 * @brief Read-only AUTO triggers, route selection and work grants.
 */

#ifndef LP_SIMPLEX_AUTO_POLICY_INTERNAL_H
#define LP_SIMPLEX_AUTO_POLICY_INTERNAL_H

#include "simplex_auto.h"

/* Policy quantities, not correctness tolerances. One observation quantum
 * gives a recovery kernel a chance to measure progress. Reserving at most
 * half the remainder keeps the live primary path the largest allocation. */
int simplex_auto_recovery_reserve(int remaining, int recovery_paths);
/* Feasibility is an optional search, bounded by eight observations and
 * its caller's grant, rather than another fraction of the total budget. */
int simplex_auto_initializer_limit(int available);

int simplex_auto_policy_initializer_grant(int remaining, int recovery_paths);
/* Return the actual grant; zero means no optional cover attempt. Dual reserves
 * one observation per recovery route plus its final recovery, if any. */
int simplex_auto_policy_cover_grant(int rows, int columns, int remaining);
int simplex_auto_policy_dual_grant(int remaining, int recovery_routes);

/* Controller-produced evidence, separate from the display journal. Presence
 * flags describe audited state, not whether a diagnostic origin was retained. */
struct simplex_AutoPolicyEvidence {
	int have_point, have_basis, dual_yielded, cover_yielded;
	int dual_contract_failure;
	int standard_basis_uncertified, standard_basis_contract_failure;
	int native_oracle_failure, standard_oracle_failure;
	struct simplex_AutoOrigin point, basis, dual, cover, retryable;
	struct simplex_AutoOrigin native_oracle, standard_oracle;
	struct simplex_AutoOrigin standard_basis;
};

struct simplex_AutoTrigger {
	int reason;
	struct simplex_AutoOrigin origin;
};
struct simplex_AutoTrigger simplex_auto_policy_trigger(
	int route, const struct simplex_AutoPolicyEvidence *evidence);

/* Proposal only: no pointers into mutable AUTO runtime, candidates or history.
 * index is a priority cursor, not an enum value. A feasible initializer may
 * supersede singleton search through a policy cursor proposal. */
struct simplex_AutoDecision {
	int index, route, grant, recovery_paths, reserve;
};
struct simplex_AutoDecision simplex_auto_policy_next(
	const int allowed[SIMPLEX_AUTO_ROUTE_COUNT], int previous,
	int remaining);
/* A committed feasible initializer may bypass singleton feasibility search.
 * Cursor positions belong to policy, never hard-coded by the controller. */
int simplex_auto_policy_seed_cursor(int previous, int accepted);
int simplex_auto_policy_oracle_retry(
	int route, const struct simplex_StageReport *report);
int simplex_auto_policy_purpose(
	int route, const struct simplex_AutoPolicyEvidence *evidence);
int simplex_auto_policy_prerequisite(int route, int have_point, int have_basis);
int simplex_auto_policy_uses_point(int route);

#endif
