/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_stage.c
 * @brief Stage work accounting, evidence and handoff validation.
 */

#include "simplex_stage.h"
#include "simplex_problem.h"
#include "simplex_progress.h"
#include <limits.h>
#include "utils.h"
#include <math.h>

void simplex_stage_reset(struct simplex_StageReport *r)
{
	lp_simplex_memset(r, 0, sizeof(*r));
	r->phase = -1;
	r->primal.value = r->dual.value = r->working_primal.value =
		r->working_dual.value = r->first_order_merit.value =
			r->lower_bound.value = r->feasibility_merit.value =
				r->violation_count.value =
					r->feasible_objective.value = NAN;
}

int simplex_stage_stop_from_status(int status)
{
	switch (status) {
	case lp_simplex_Success:
		return SIMPLEX_STOP_OPTIMAL;
	case lp_simplex_Infeasibility:
		return SIMPLEX_STOP_INFEASIBLE;
	case lp_simplex_Unboundedness:
		return SIMPLEX_STOP_UNBOUNDED;
	case lp_simplex_ExceedIterLimit:
		return SIMPLEX_STOP_BUDGET;
	case lp_simplex_Degeneracy:
		return SIMPLEX_STOP_STAGNATION;
	case lp_simplex_MemoryAllocError:
		return SIMPLEX_STOP_RESOURCE;
	case lp_simplex_PrecisionError:
	case lp_simplex_Singularity:
		return SIMPLEX_STOP_NUMERICAL;
	default:
		return SIMPLEX_STOP_UNKNOWN;
	}
}

const char *simplex_stage_stop_name(int stop)
{
	switch (stop) {
	case SIMPLEX_STOP_OPTIMAL:
		return "optimal";
	case SIMPLEX_STOP_INFEASIBLE:
		return "certified-infeasible";
	case SIMPLEX_STOP_UNBOUNDED:
		return "certified-unbounded";
	case SIMPLEX_STOP_BUDGET:
		return "work-budget";
	case SIMPLEX_STOP_STAGNATION:
		return "observed-stagnation";
	case SIMPLEX_STOP_CYCLE:
		return "repeated-state";
	case SIMPLEX_STOP_NUMERICAL:
		return "numerical-failure";
	case SIMPLEX_STOP_REASSESS:
		return "measurement-changed";
	case SIMPLEX_STOP_HANDOFF:
		return "candidate-handoff";
	case SIMPLEX_STOP_RESOURCE:
		return "local-resource-failure";
	case SIMPLEX_STOP_NOT_APPLICABLE:
		return "not-applicable";
	default:
		return "unknown";
	}
}

int simplex_stage_observation_valid(const struct simplex_Observation *v,
	int space, long epoch, int minimum_quality)
{
	return v->quality != SIMPLEX_EVIDENCE_UNKNOWN &&
	       v->quality <= SIMPLEX_EVIDENCE_VERIFIED &&
	       v->space != SIMPLEX_SPACE_UNKNOWN &&
	       v->quality >= minimum_quality && v->space == space &&
	       v->epoch == epoch && isfinite(v->value) && v->value >= 0.;
}

int simplex_stage_oracle_retry(const struct simplex_StageReport *r)
{
	return r->stop == SIMPLEX_STOP_NUMERICAL &&
	       r->failure == SIMPLEX_FAILURE_IMPLICIT_ORACLE;
}

int simplex_stage_prefer_point(const struct simplex_Problem *p,
	const struct lp_simplex_Options *o, const double *candidate,
	const double *incumbent)
{
	double primal = simplex_problem_primal_infeasibility(p, candidate);
	long double cost = 0., saved_cost = 0.;
	int j;

	if (!isfinite(primal) || primal > o->primal_tolerance)
		return 0;
	for (j = 0; j < p->columns; j++)
		cost += (long double)p->objective[j] * candidate[j];
	if (!isfinite(cost))
		return 0;
	primal = simplex_problem_primal_infeasibility(p, incumbent);
	if (!isfinite(primal) || primal > o->primal_tolerance)
		return 1;
	for (j = 0; j < p->columns; j++)
		saved_cost += (long double)p->objective[j] * incumbent[j];
	return !isfinite(saved_cost) || cost < saved_cost;
}

int simplex_stage_observation_current(const struct simplex_Observation *v,
	int space, long epoch, int work, int minimum_quality)
{
	int signed_value = v->kind == SIMPLEX_OBSERVATION_LOWER_BOUND ||
			   v->kind == SIMPLEX_OBSERVATION_FEASIBLE_OBJECTIVE;
	return v->quality >= minimum_quality &&
	       v->quality > SIMPLEX_EVIDENCE_UNKNOWN &&
	       v->quality <= SIMPLEX_EVIDENCE_VERIFIED && v->space == space &&
	       space != SIMPLEX_SPACE_UNKNOWN && v->epoch == epoch &&
	       v->work == work && isfinite(v->value) &&
	       (signed_value || v->value >= 0.);
}

int simplex_stage_progress_observe(struct simplex_ObservationHistory *h,
	const struct simplex_Observation *v, int work, double relative,
	double noise)
{
	double gain;
	int increasing = v->kind == SIMPLEX_OBSERVATION_LOWER_BOUND;

	if (work < 0 || v->work < 0)
		return SIMPLEX_OBSERVATION_UNAVAILABLE;
	if (v->work != work)
		return SIMPLEX_OBSERVATION_STALE;
	if (v->quality < SIMPLEX_EVIDENCE_RECOMPUTED ||
		v->quality > SIMPLEX_EVIDENCE_VERIFIED ||
		(v->space < SIMPLEX_SPACE_INPUT ||
			v->space > SIMPLEX_SPACE_SCALED_INPUT) ||
		v->kind <= SIMPLEX_OBSERVATION_UNKNOWN ||
		v->kind > SIMPLEX_OBSERVATION_LOWER_BOUND || v->work < 0 ||
		!isfinite(v->value) || !isfinite(relative) || relative < 0. ||
		!isfinite(noise) || noise < 0.)
		return SIMPLEX_OBSERVATION_UNAVAILABLE;
	if (v->kind != SIMPLEX_OBSERVATION_FEASIBLE_OBJECTIVE && !increasing &&
		v->value < 0.)
		return SIMPLEX_OBSERVATION_UNAVAILABLE;
	if (h->initialized && v->work < h->last_work)
		return SIMPLEX_OBSERVATION_STALE;
	if (!h->initialized || h->best.space != v->space ||
		h->best.epoch != v->epoch || h->best.kind != v->kind) {
		h->initialized = 1;
		h->last_work = v->work;
		h->best = *v;
		return SIMPLEX_OBSERVATION_BASELINE;
	}
	if (v->work == h->last_work)
		return SIMPLEX_OBSERVATION_STALE;
	h->last_work = v->work;
	gain = increasing ? v->value - h->best.value : h->best.value - v->value;
	if (gain > noise && gain > relative * fabs(h->best.value)) {
		h->best = *v;
		return SIMPLEX_OBSERVATION_IMPROVED;
	}
	return SIMPLEX_OBSERVATION_UNCHANGED;
}

static void stage_observation(struct simplex_Observation *v, double value,
	int kind, int work, int quality)
{
	v->value = isfinite(value) ? value : NAN;
	v->kind = kind;
	v->work = work;
	v->space = SIMPLEX_SPACE_INPUT;
	v->epoch = 0;
	v->quality = isfinite(value) ? quality : SIMPLEX_EVIDENCE_UNKNOWN;
}

void simplex_stage_observe(const struct simplex_Problem *p,
	const struct lp_simplex_Options *o, const double *x, const double *y,
	int sign, int work, struct simplex_StageReport *r)
{
	struct simplex_PrimalSummary primal;
	long double cost = 0.;
	double maximum = NAN, sum = NAN, count = NAN, objective = NAN;
	double dual = NAN, bound = NAN;
	int j;

	if (x && !simplex_problem_primal_summary(
			 p, x, o->primal_tolerance, &primal)) {
		maximum = primal.maximum;
		sum = primal.sum;
		count = primal.count;
		if (maximum <= o->primal_tolerance) {
			for (j = 0; j < p->columns; j++)
				cost += (long double)p->objective[j] * x[j];
			if (isfinite(cost))
				objective = (double)cost;
		}
	}
	if (y && (sign == 1 || sign == -1)) {
		if (x)
			dual = simplex_problem_dual_infeasibility(
				p, x, y, o->primal_tolerance, sign);
		bound = simplex_progress_dual_bound(p, y, sign);
	}
	stage_observation(&r->primal, maximum, SIMPLEX_OBSERVATION_PRIMAL_MAX,
		work, SIMPLEX_EVIDENCE_RECOMPUTED);
	stage_observation(&r->feasibility_merit, sum,
		SIMPLEX_OBSERVATION_PRIMAL_SUM, work,
		SIMPLEX_EVIDENCE_RECOMPUTED);
	stage_observation(&r->violation_count, count,
		SIMPLEX_OBSERVATION_VIOLATION_COUNT, work,
		SIMPLEX_EVIDENCE_RECOMPUTED);
	stage_observation(&r->feasible_objective, objective,
		SIMPLEX_OBSERVATION_FEASIBLE_OBJECTIVE, work,
		SIMPLEX_EVIDENCE_RECOMPUTED);
	stage_observation(&r->dual, dual, SIMPLEX_OBSERVATION_DUAL_MAX, work,
		SIMPLEX_EVIDENCE_RECOMPUTED);
	stage_observation(&r->lower_bound, bound,
		SIMPLEX_OBSERVATION_LOWER_BOUND, work,
		SIMPLEX_EVIDENCE_VERIFIED);
	if (r->input_audits < INT_MAX)
		r->input_audits++;
	if (!r->input_progress_space)
		r->input_progress_space = SIMPLEX_SPACE_INPUT;
}
