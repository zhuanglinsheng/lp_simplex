/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#ifndef LP_SIMPLEX_DUAL_STATE_INTERNAL_H
#define LP_SIMPLEX_DUAL_STATE_INTERNAL_H

#include "simplex_basis.h"
#include "simplex_csc.h"
#include "simplex_degeneracy.h"
#include "simplex_dual_feasibility.h"
#include "simplex_sparse_vector.h"

#include <lp_simplex/solve.h>


#define DUAL_STATUS_BASIC 0
#define DUAL_STATUS_LOWER 1
#define DUAL_STATUS_UPPER 2
#define DUAL_STATUS_FIXED 3
#define DUAL_STATUS_FREE  4


struct simplex_DualProfile {
	int enabled;
	double ratio_seconds;
	double total_seconds;
	long rejected_relative;
	long rejected_ftran;
	long bound_flips;
	long flip_batches;
};


/* Owners for the dense solver workspaces.  The state keeps typed aliases for
 * hot-path readability; destruction only releases these four roots. */
struct simplex_DualWorkspace {
	int *index;
	unsigned char *status;
	double *column;
	double *row;
};


/* Shared only by the dual orchestration modules; not part of the public API. */
struct simplex_DualState {
	int rows;
	int structural;
	int variables;
	int iterations;
	int reinversions;
	int structural_basic;
	const struct lp_simplex_Options *options;
	struct simplex_CscMatrix matrix;
	struct simplex_Basis factor;
	struct simplex_DualWorkspace workspace;
	int *basis;
	int *position;
	unsigned char *status;
	double *lower;
	double *upper;
	double *cost;
	double *value;
	double *basic_value;
	double *basic_lower;
	double *basic_upper;
	double *feasibility_tolerance;
	double *reduced;
	double *pi;
	double *rho;
	struct simplex_SparseVector packed_rho;
	struct simplex_SparseVector packed_alpha;
	int packed_alpha_valid;
	struct simplex_SparseVector packed_direction;
	struct simplex_SparseVector flip_rhs;
	double *alpha;
	double *breakpoint;
	double *direction;
	double *work;
	double *edge_weight;
	struct simplex_DualFeasibility feasibility;
	int *candidate_sign;
	int *candidate_index;
	int candidate_count;
	int *nonbasic_index;
	int *nonbasic_slot;
	int nonbasic_count;
	struct simplex_DegeneracyControl degeneracy;
	int pan_enabled;
	int stable_candidate_order;
	int regular_columns;
	int numerically_stressed;
	int modal_column_degree;
	struct simplex_DualProfile profile;
	double ratio_minimum;
	int ratio_minimum_valid;
	int pricing_validation_countdown;
	/* A leaving row whose complete Harris set had no stable actual FTRAN
	 * pivot is skipped once so Pan can move along another face direction. */
	int pan_deferred_row;
};

#endif
