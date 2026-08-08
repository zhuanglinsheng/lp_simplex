#ifndef LP_SIMPLEX_DUAL_STATE_INTERNAL_H
#define LP_SIMPLEX_DUAL_STATE_INTERNAL_H

#include "simplex_basis.h"
#include "simplex_csc.h"
#include "simplex_degeneracy.h"
#include "simplex_sparse_vector.h"
#include <lp_simplex/solve.h>

#define DUAL_STATUS_BASIC 0
#define DUAL_STATUS_LOWER 1
#define DUAL_STATUS_UPPER 2
#define DUAL_STATUS_FIXED 3
#define DUAL_STATUS_FREE  4

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
	double *reduced;
	double *pi;
	double *rho;
	struct simplex_SparseVector packed_rho;
	struct simplex_SparseVector packed_alpha;
	int packed_alpha_valid;
	struct simplex_SparseVector packed_direction;
	double *alpha;
	double *breakpoint;
	double *direction;
	double *work;
	double *edge_weight;
	int *candidate_sign;
	int *candidate_index;
	int candidate_count;
	int *nonbasic_index;
	int *nonbasic_slot;
	int nonbasic_count;
	struct simplex_DegeneracyControl degeneracy;
	int pan_enabled;
	int profile_enabled;
	double profile_ratio_seconds;
	double profile_total_seconds;
	double ratio_minimum;
	int ratio_minimum_valid;
};

#endif
