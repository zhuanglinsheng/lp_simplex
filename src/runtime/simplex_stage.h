/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_stage.h
 * @brief Stage work accounting, evidence and handoff validation.
 */

#ifndef LP_SIMPLEX_STAGE_INTERNAL_H
#define LP_SIMPLEX_STAGE_INTERNAL_H

#include <lp_simplex/solve.h>
#include "simplex_pan_service_types.h"
#include "simplex_tableau_types.h"
#include "simplex_singleton_types.h"

/* All zero values mean unknown, never a zero residual or a proof. Input is
 * the immutable problem passed to this stage, possibly a presolved model.
 * Working observations cannot be compared across changes of epoch. */
enum simplex_EvidenceQuality {
	SIMPLEX_EVIDENCE_UNKNOWN = 0,
	SIMPLEX_EVIDENCE_ESTIMATED,
	SIMPLEX_EVIDENCE_RECOMPUTED,
	SIMPLEX_EVIDENCE_VERIFIED
};

enum simplex_EvidenceSpace {
	SIMPLEX_SPACE_UNKNOWN = 0,
	SIMPLEX_SPACE_INPUT,
	SIMPLEX_SPACE_WORKING,
	SIMPLEX_SPACE_SCALED_INPUT
};

enum simplex_StageStop {
	SIMPLEX_STOP_UNKNOWN = 0,
	SIMPLEX_STOP_OPTIMAL,
	SIMPLEX_STOP_INFEASIBLE,
	SIMPLEX_STOP_UNBOUNDED,
	SIMPLEX_STOP_BUDGET,
	SIMPLEX_STOP_STAGNATION,
	SIMPLEX_STOP_CYCLE,
	SIMPLEX_STOP_NUMERICAL,
	SIMPLEX_STOP_RESOURCE,
	SIMPLEX_STOP_HANDOFF,
	SIMPLEX_STOP_REASSESS,
	SIMPLEX_STOP_NOT_APPLICABLE
};

enum simplex_StageFailure {
	SIMPLEX_FAILURE_UNKNOWN = 0,
	SIMPLEX_FAILURE_FACTORIZATION,
	/* Repeated/conflicting implicit cone oracle work, not a model form. */
	SIMPLEX_FAILURE_IMPLICIT_ORACLE,
	SIMPLEX_FAILURE_DIRECTION,
	SIMPLEX_FAILURE_FEASIBILITY,
	SIMPLEX_FAILURE_CERTIFICATION,
	SIMPLEX_FAILURE_STAGE_CONTRACT
};

/* A rejected basis hint is not a mathematical infeasibility certificate.
 * Layout/state failures cannot be repaired by changing factor backends. */
enum simplex_BasisRejection {
	SIMPLEX_BASIS_NOT_REJECTED = 0,
	SIMPLEX_BASIS_LAYOUT_REJECTED,
	SIMPLEX_BASIS_FACTORIZATION_REJECTED,
	SIMPLEX_BASIS_RESIDUAL_REJECTED,
	SIMPLEX_BASIS_FEASIBILITY_REJECTED,
	SIMPLEX_BASIS_MODEL_REJECTED,
	SIMPLEX_BASIS_RESOURCE_REJECTED
};

enum simplex_ObservationKind {
	SIMPLEX_OBSERVATION_UNKNOWN = 0,
	SIMPLEX_OBSERVATION_PRIMAL_MAX,
	SIMPLEX_OBSERVATION_PRIMAL_SUM,
	SIMPLEX_OBSERVATION_VIOLATION_COUNT,
	SIMPLEX_OBSERVATION_DUAL_MAX,
	SIMPLEX_OBSERVATION_FEASIBLE_OBJECTIVE,
	SIMPLEX_OBSERVATION_LOWER_BOUND
};

enum simplex_ObservationChange {
	SIMPLEX_OBSERVATION_UNAVAILABLE = 0,
	SIMPLEX_OBSERVATION_BASELINE,
	SIMPLEX_OBSERVATION_IMPROVED,
	SIMPLEX_OBSERVATION_UNCHANGED,
	SIMPLEX_OBSERVATION_STALE
};

struct simplex_Observation {
	double value;
	int quality, space, kind;
	/* Work position of this observation, not of the latest estimate. */
	int work;
	long epoch;
};

/* A solve-local baseline for one quantity. Metric/model changes establish a
 * baseline, never progress. last_work also advances on rebounds, so an older
 * sample cannot masquerade as a new observation of the historical best. */
struct simplex_ObservationHistory {
	int initialized, last_work;
	struct simplex_Observation best;
};

int simplex_stage_progress_observe(struct simplex_ObservationHistory *history,
	const struct simplex_Observation *sample, int work, double relative,
	double noise);
int simplex_stage_observation_current(const struct simplex_Observation *sample,
	int space, long epoch, int work, int minimum_quality);

#define SIMPLEX_INPUT_PRIMAL_PROGRESS 1U
#define SIMPLEX_INPUT_COUNT_PROGRESS 2U
#define SIMPLEX_INPUT_OBJECTIVE_PROGRESS 4U
#define SIMPLEX_INPUT_BOUND_PROGRESS 8U
#define SIMPLEX_INPUT_DUAL_PROGRESS 16U
#define SIMPLEX_INPUT_FEASIBLE_AVAILABLE 32U
#define SIMPLEX_INPUT_BOUND_AVAILABLE 64U

struct simplex_StageWork {
	/* Counts are independent quantities, not equivalent pivot charges.
	 * Unknown counts have a clear availability flag. */
	int basis_counts_available;
	long factorizations, ftran, btran;
	int first_order_steps, repair_attempts, repair_successes;
	int minimum_norm_qr_available;
	long minimum_norm_qr_calls;
	double cpu_seconds;
};

struct simplex_StageReport {
	struct simplex_SingletonReport singleton;
	struct simplex_TableauReport tableau;
	struct simplex_PanServiceSummary services;
	int stop, failure;
	int phase, perturbed_bounds, perturbed_costs;
	long bounds_epoch, cost_epoch;
	struct simplex_Observation primal, dual, working_primal, working_dual;
	/* Conservative box Lagrangian bound on this stage's immutable input.
	 * A finite bound can be negative; it is not a residual. */
	struct simplex_Observation lower_bound;
	struct simplex_Observation feasibility_merit, violation_count;
	struct simplex_Observation feasible_objective;
	unsigned int input_improvements;
	int input_audits, input_progress_space;
	/* Exported basis indices may be a hint even when no reliable factor
	 * survives. Import must reconstruct and audit before using them. */
	int basis_available, factor_available, point_available;
	int basis_import_attempted, basis_rejection;
	/* Backend invocation skipped due to an exact prior cold-start failure. */
	int recovery_skipped;
	long repeated_states, rebuilds, backend_retries;
	struct simplex_StageWork work;
	/* PDHG's internal handoff is separate from the enclosing stage stop. */
	int first_order_stop;
	struct simplex_Observation first_order_merit;
};
struct simplex_Problem;
/* Refresh input observations only. Numerical work state, stop causes and
 * progress history remain owned by the caller. Missing vectors stay unknown. */
void simplex_stage_observe(const struct simplex_Problem *problem,
	const struct lp_simplex_Options *options, const double *x,
	const double *y, int multiplier_sign, int work,
	struct simplex_StageReport *report);
int simplex_stage_prefer_point(const struct simplex_Problem *problem,
	const struct lp_simplex_Options *options, const double *candidate,
	const double *incumbent);
void simplex_stage_reset(struct simplex_StageReport *report);
int simplex_stage_stop_from_status(int status);
const char *simplex_stage_stop_name(int stop);
int simplex_stage_observation_valid(const struct simplex_Observation *v,
	int space, long epoch, int minimum_quality);
/* Numerical ambiguity alone does not establish implicit-oracle failure. */
int simplex_stage_oracle_retry(const struct simplex_StageReport *r);

#endif
