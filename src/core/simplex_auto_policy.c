/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_auto_policy.c
 * @brief Read-only AUTO decisions. No kernel calls, allocation or state
 * commits.
 */

#include "simplex_stage.h"
#include "simplex_execution.h"
#include "simplex_auto_policy.h"
#include <limits.h>

static const int auto_order[SIMPLEX_AUTO_ROUTE_COUNT] = {SIMPLEX_AUTO_SINGLETON,
	SIMPLEX_AUTO_HYBRID, SIMPLEX_AUTO_EXPLICIT, SIMPLEX_AUTO_CLASSIC,
	SIMPLEX_AUTO_STANDARD_IMPLICIT, SIMPLEX_AUTO_STANDARD_EXPLICIT,
	SIMPLEX_AUTO_TABLEAU};

/* An observation is useful only if some work remains for another method.
 * Saturate before multiplication; this policy also covers tiny user budgets. */
int simplex_auto_recovery_reserve(int remaining, int recovery_paths)
{
	int cap;
	if (remaining <= 1 || recovery_paths <= 0)
		return 0;
	cap = remaining / 2;
	if (recovery_paths > cap / SIMPLEX_AUTO_PROBE_CAP)
		return cap;
	return recovery_paths * SIMPLEX_AUTO_PROBE_CAP;
}

int simplex_auto_initializer_limit(int available)
{
	int cap =
		SIMPLEX_AUTO_INITIALIZER_OBSERVATIONS * SIMPLEX_AUTO_PROBE_CAP;
	if (available <= 0)
		return 0;
	return available < cap ? available : cap;
}

int simplex_auto_policy_initializer_grant(int remaining, int recovery_paths)
{
	return simplex_auto_initializer_limit(
		remaining -
		simplex_auto_recovery_reserve(remaining, recovery_paths));
}

/* One action description is used for both the journal and trace. The fixed
 * priority order remains a policy; these fields state its actual prerequisite
 * and intended effect without inventing a numerical failure diagnosis. */
int simplex_auto_policy_purpose(
	int route, const struct simplex_AutoPolicyEvidence *e)
{
	if (route == SIMPLEX_AUTO_COVER)
		return SIMPLEX_AUTO_RESTRICT_COLUMNS;
	if (route == SIMPLEX_AUTO_INITIALIZER)
		return SIMPLEX_AUTO_SEEK_FEASIBLE;
	if (route == SIMPLEX_AUTO_DUAL_FALLBACK && e->have_basis)
		return SIMPLEX_AUTO_REBUILD_BASIS;
	if (route == SIMPLEX_AUTO_DUAL_EXPLORE ||
		route == SIMPLEX_AUTO_DUAL_FALLBACK)
		return SIMPLEX_AUTO_PRIMARY_SEARCH;
	if (route == SIMPLEX_AUTO_EXPLICIT ||
		route == SIMPLEX_AUTO_STANDARD_EXPLICIT)
		return SIMPLEX_AUTO_REPLACE_ORACLE;
	if (route == SIMPLEX_AUTO_TABLEAU)
		return SIMPLEX_AUTO_ALTERNATIVE_ALGORITHM;
	if (route == SIMPLEX_AUTO_STANDARD_IMPLICIT &&
		e->standard_basis_uncertified)
		return SIMPLEX_AUTO_REPLACE_DIRECTION;
	if (route == SIMPLEX_AUTO_HYBRID)
		return e->have_point ? SIMPLEX_AUTO_OPTIMIZE_POINT
				     : SIMPLEX_AUTO_SEEK_FEASIBLE;
	return SIMPLEX_AUTO_REFORMULATE;
}

int simplex_auto_policy_uses_point(int route)
{
	return route == SIMPLEX_AUTO_HYBRID || route == SIMPLEX_AUTO_EXPLICIT ||
	       route == SIMPLEX_AUTO_CLASSIC ||
	       route == SIMPLEX_AUTO_STANDARD_IMPLICIT ||
	       route == SIMPLEX_AUTO_STANDARD_EXPLICIT;
}

int simplex_auto_policy_prerequisite(int route, int have_point, int have_basis)
{
	if (route == SIMPLEX_AUTO_EXPLICIT ||
		route == SIMPLEX_AUTO_STANDARD_EXPLICIT)
		return SIMPLEX_AUTO_IMPLICIT_ORACLE_FAILURE;
	if (route == SIMPLEX_AUTO_DUAL_FALLBACK && have_basis)
		return SIMPLEX_AUTO_BASIS_HINT;
	if (route == SIMPLEX_AUTO_DUAL_EXPLORE ||
		route == SIMPLEX_AUTO_DUAL_FALLBACK)
		return SIMPLEX_AUTO_UNCERTIFIED;
	return have_point && simplex_auto_policy_uses_point(route)
		       ? SIMPLEX_AUTO_AUDITED_POINT
		       : SIMPLEX_AUTO_UNCERTIFIED;
}

int simplex_auto_policy_cover_grant(int rows, int columns, int remaining)
{
	return rows > 0 &&
			       columns / (2 * SIMPLEX_AUTO_COVER_COLUMN_MARGIN) >=
				       rows &&
			       remaining / (2 * SIMPLEX_AUTO_COVER_MIN_SWEEPS) >=
				       rows
		       ? remaining / 2
		       : 0;
}

int simplex_auto_policy_dual_grant(int remaining, int recovery_routes)
{
	int paths = recovery_routes <= 0	 ? 0
		    : recovery_routes == INT_MAX ? INT_MAX
						 : recovery_routes + 1;
	if (remaining <= 0)
		return 0;
	return remaining - simplex_auto_recovery_reserve(remaining, paths);
}

int simplex_auto_policy_seed_cursor(int previous, int accepted)
{
	int i;
	if (accepted)
		for (i = 0; i < SIMPLEX_AUTO_ROUTE_COUNT; i++)
			if (auto_order[i] == SIMPLEX_AUTO_SINGLETON &&
				previous < i)
				return i;
	return previous;
}

int simplex_auto_policy_oracle_retry(
	int route, const struct simplex_StageReport *report)
{
	return (route == SIMPLEX_AUTO_HYBRID ||
		       route == SIMPLEX_AUTO_STANDARD_IMPLICIT) &&
	       simplex_stage_oracle_retry(report);
}

/* A read-only scheduling query. The caller commits index/route/grant. */
struct simplex_AutoDecision simplex_auto_policy_next(
	const int allowed[SIMPLEX_AUTO_ROUTE_COUNT], int previous,
	int remaining)
{
	struct simplex_AutoDecision d;
	int j;

	d.index = previous < -1 ? -1 : previous;
	if (d.index < SIMPLEX_AUTO_ROUTE_COUNT)
		d.index++;
	while (d.index < SIMPLEX_AUTO_ROUTE_COUNT &&
		!allowed[auto_order[d.index]])
		d.index++;
	d.route = SIMPLEX_AUTO_DUAL_FALLBACK;
	d.grant = remaining > 0 ? remaining : 0;
	d.recovery_paths = d.reserve = 0;
	if (d.index >= SIMPLEX_AUTO_ROUTE_COUNT || remaining <= 0)
		return d;
	d.route = auto_order[d.index];
	d.recovery_paths = 1;
	for (j = d.index + 1; j < SIMPLEX_AUTO_ROUTE_COUNT; j++)
		d.recovery_paths += allowed[auto_order[j]] != 0;
	if (d.route == SIMPLEX_AUTO_HYBRID) {
		d.reserve = simplex_auto_recovery_reserve(
			remaining, d.recovery_paths);
		d.grant = remaining - d.reserve;
	} else
		d.grant = remaining / (d.recovery_paths + 1);
	if (d.grant < 1)
		d.grant = 1;
	return d;
}

/* Causality follows the selected mechanism's enabling evidence. A point may
 * come from an initializer while the failed dual still motivates the handoff. */
struct simplex_AutoTrigger simplex_auto_policy_trigger(
	int route, const struct simplex_AutoPolicyEvidence *e)
{
	struct simplex_AutoTrigger t = {0};
	t.reason = SIMPLEX_REASON_NONE;
	t.origin.record = -1;
	t.origin.event = 0;
	if (route == SIMPLEX_AUTO_COVER)
		return t;
	if (route == SIMPLEX_AUTO_DUAL_EXPLORE) {
		if (e->cover_yielded) {
			t.reason = SIMPLEX_REASON_COVER_UNCERTIFIED;
			t.origin = e->cover;
		}
		return t;
	}
	if (route == SIMPLEX_AUTO_EXPLICIT && e->native_oracle_failure) {
		t.reason = SIMPLEX_REASON_IMPLICIT_ORACLE_FAILURE;
		t.origin = e->native_oracle;
	} else if (route == SIMPLEX_AUTO_STANDARD_EXPLICIT &&
		   e->standard_oracle_failure) {
		t.reason = SIMPLEX_REASON_IMPLICIT_ORACLE_FAILURE;
		t.origin = e->standard_oracle;
	} else if (route == SIMPLEX_AUTO_DUAL_FALLBACK && e->have_basis) {
		t.reason = SIMPLEX_REASON_BASIS_REBUILD;
		t.origin = e->basis;
	} else if (route == SIMPLEX_AUTO_DUAL_FALLBACK) {
		t.reason = SIMPLEX_REASON_RECOVERY_EXHAUSTED;
		t.origin = e->retryable;
	} else if ((route == SIMPLEX_AUTO_INITIALIZER ||
			   route == SIMPLEX_AUTO_HYBRID ||
			   route == SIMPLEX_AUTO_SINGLETON) &&
		   e->dual_yielded) {
		t.reason = e->dual_contract_failure
				   ? SIMPLEX_REASON_STAGE_CONTRACT
				   : SIMPLEX_REASON_DUAL_HANDOFF;
		t.origin = e->dual;
	} else if (route == SIMPLEX_AUTO_TABLEAU) {
		t.reason = SIMPLEX_REASON_ALTERNATIVE_ALGORITHM;
		t.origin = e->retryable;
	} else if (route == SIMPLEX_AUTO_STANDARD_IMPLICIT &&
		   e->standard_basis_uncertified) {
		t.reason = e->standard_basis_contract_failure
				   ? SIMPLEX_REASON_STAGE_CONTRACT
				   : SIMPLEX_REASON_ALTERNATIVE_DIRECTION;
		t.origin = e->standard_basis;
	} else if (route == SIMPLEX_AUTO_CLASSIC ||
		   route == SIMPLEX_AUTO_STANDARD_IMPLICIT) {
		t.reason = SIMPLEX_REASON_ALTERNATIVE_REPRESENTATION;
		t.origin = e->retryable;
	}
	return t;
}
