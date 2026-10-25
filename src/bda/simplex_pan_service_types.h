/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_pan_service_types.h
 * @brief BDA service roles, requests and outcome contracts.
 */

#ifndef LP_SIMPLEX_PAN_SERVICE_TYPES_H
#define LP_SIMPLEX_PAN_SERVICE_TYPES_H

#include "simplex_execution_scope.h"
#include "simplex_pan_face_types.h"

/* Services are roles within BDA, not additional public LP algorithms. */
enum simplex_PanServiceKind {
	SIMPLEX_PAN_SERVICE_NNLS = 0,
	SIMPLEX_PAN_SERVICE_DUAL_FEASIBILITY,
	SIMPLEX_PAN_SERVICE_PDHG,
	SIMPLEX_PAN_SERVICE_CONE,
	SIMPLEX_PAN_SERVICE_RESTORE,
	SIMPLEX_PAN_SERVICE_DUAL_REFINE,
	SIMPLEX_PAN_SERVICE_PREDICT,
	SIMPLEX_PAN_SERVICE_FACE_RECOVERY,
	SIMPLEX_PAN_SERVICE_COUNT
};

enum simplex_PanServiceOutput {
	SIMPLEX_PAN_OUTPUT_NONE = 0,
	SIMPLEX_PAN_OUTPUT_NONNEGATIVE_EQUALITY_POINT,
	SIMPLEX_PAN_OUTPUT_INPUT_FEASIBLE_POINT,
	SIMPLEX_PAN_OUTPUT_NATIVE_FEASIBLE_POINT,
	SIMPLEX_PAN_OUTPUT_DIRECTION_PROPOSAL,
	SIMPLEX_PAN_OUTPUT_STATIONARY_RESIDUAL,
	SIMPLEX_PAN_OUTPUT_CHECKED_DIRECTION,
	SIMPLEX_PAN_OUTPUT_CHECKED_DUAL,
	SIMPLEX_PAN_OUTPUT_NATIVE_KKT,
	SIMPLEX_PAN_OUTPUT_NATIVE_FACE_ADJUSTMENT
};
struct simplex_PanStandard;
struct simplex_Problem;

enum simplex_PanServiceTriggerKind {
	SIMPLEX_PAN_TRIGGER_UNKNOWN = 0,
	SIMPLEX_PAN_TRIGGER_PDHG_STAGNATION,
	SIMPLEX_PAN_TRIGGER_PRIMAL_MARGIN
};

/* Producer observations only: neither a permission nor a reusable failure key.
 * Previous residual may be infinity when no repair has yet been attempted. */
struct simplex_PanServiceTrigger {
	int kind, work;
	long epoch;
	double primal, previous_primal;
};

/* Invocation-local output, not a reusable checkpoint. Raw work/status/output
 * remain visible even when contract validation rejects the handoff. */
struct simplex_PanServiceOutcome {
	struct simplex_ExecutionRef actor;
	struct simplex_PanFaceReport face;
	struct simplex_PanServiceTrigger trigger;
	int service, valid, stop, raw_stop, output, raw_output;
	int raw_code, raw_status, before, grant, reported_end, end, event;
	const struct simplex_PanStandard *representation;
	const struct simplex_Problem *input;
	unsigned long revision;
	int rows, columns, input_model, working_model;
	unsigned long representation_revision;
	int reported_work, relative_work, implementation;
};

/* Output-only aggregates. Nested work is inclusive and must not be summed
 * across services. Zero-initialization means not observed, not a free call. */
struct simplex_PanServiceSummary {
	long calls[SIMPLEX_PAN_SERVICE_COUNT];
	long work[SIMPLEX_PAN_SERVICE_COUNT];
	int stop[SIMPLEX_PAN_SERVICE_COUNT], output[SIMPLEX_PAN_SERVICE_COUNT];
	int implementation[SIMPLEX_PAN_SERVICE_COUNT];
	long invalid_contracts;
};

const char *simplex_pan_service_name(int service);
const char *simplex_pan_service_output_name(int output);
#endif
