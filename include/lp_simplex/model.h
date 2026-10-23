/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file model.h
 * @brief Dense and sparse linear-program models and their storage lifecycle.
 */

#ifndef LP_SIMPLEX_MODEL_H
#define LP_SIMPLEX_MODEL_H


#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/** Continuous variable; the only variable type accepted by the solver. */
#define optm_VAR_T_REAL 0
/** Integer declaration retained by the MPS reader; solving is unsupported. */
#define optm_VAR_T_INT 1
/** Binary declaration retained by the MPS reader; solving is unsupported. */
#define optm_VAR_T_BIN 2

/** No lower or upper bound; lb and ub are ignored. */
#define optm_BOUND_T_FR 0
/** x <= ub, with no lower bound; lb is ignored. */
#define optm_BOUND_T_UP 1
/** x >= lb, with no upper bound; ub is ignored. */
#define optm_BOUND_T_LO 2
/** lb <= x <= ub; equal bounds fix the variable. */
#define optm_BOUND_T_BS 3

/** Row activity equals rhs. */
#define optm_CONS_T_EQ 0
/** Row activity is greater than or equal to rhs. */
#define optm_CONS_T_GE 1
/** Row activity is less than or equal to rhs. */
#define optm_CONS_T_LE 2


/**
 * @brief Variable metadata and bounds, indexed in the original variable order.
 */
struct optm_VariableBound {
	/**< Optional name: at most 8 characters plus a terminating NUL. */
	char name[9];

	/**< Lower bound, used only by LO and BS. Default: 0. */
	double lb;

	/**< Upper bound, used only by UP and BS. Default: positive infinity. */
	double ub;

	/**< One of optm_BOUND_T_*; default: optm_BOUND_T_LO. */
	int b_type;

	/**< One of optm_VAR_T_*; default: optm_VAR_T_REAL. */
	int v_type;
};

/**
 * @brief One row: sum_j coef[j] * x[j] {=, >=, <=} rhs.
 */
struct optm_LinearConstraint {
	/**< Optional name: at most 8 characters plus a terminating NUL. */
	char name[9];

	/**< Borrowed dense row of n entries; NULL in a sparse-only model. */
	double *coef;

	/**< Finite right-hand side; default: 0. */
	double rhs;

	/**< One of optm_CONS_T_*; default: optm_CONS_T_EQ. */
	int type;
};

/**
 * @brief Continuous LP input, with dense and/or sparse matrix storage.
 *
 * Row and column indices are zero-based. Constructors own every allocated
 * array; release the whole model with lp_model_free(), never individual rows.
 * Do not shallow-copy a model or replace its owned pointers with borrowed data.
 * Keep dimensions and allocation sizes unchanged after construction.
 *
 * The solver minimises objective_offset + sum_j objective[j] * x[j]. To
 * maximise c*x + c0, store -c and -c0 and set objective_sense to -1. The result
 * restores the caller's sense by multiplying the stored objective by this
 * field. Set it to +1 for minimisation. MPS input handles this automatically.
 *
 * Dense models own an m-by-n row-major coefficients array; constraints[i].coef
 * aliases its i-th row. After editing matrix entries, call
 * lp_model_build_sparse() and check its return value before solving.
 *
 * Sparse-only models have coefficients == NULL and all row coef pointers NULL.
 * Both CSC and CSR must describe the same matrix. CSC column j occupies
 * [column_start[j], column_start[j+1]); entry k is (row_index[k], j, value[k]).
 * CSR row i occupies [row_start[i], row_start[i+1]); entry k is
 * (i, column_index[k], row_value[k]). Both start arrays begin at 0, are
 * nondecreasing, and end at nnz; all indices must be in range. Keep both views
 * consistent when editing; the solver does not reconcile them with each other
 * or with a dense matrix. Use one entry per nonzero matrix coordinate.
 */
struct lp_Model {
	/**< Number of constraints; strictly positive. */
	int m;

	/**< Number of variables; strictly positive. */
	int n;

	/**< Used sparse entries in each view, not the allocated capacity. */
	int nnz;

	/**< Constant in the stored minimisation objective; default: 0. */
	double objective_offset;

	/**< +1 for minimise, -1 for maximise; default: +1. */
	double objective_sense;

	/**< n stored minimisation coefficients; initially zero. */
	double *objective;

	/**< Owned m*n dense matrix; NULL for sparse-only storage. */
	double *coefficients;

	/** m row descriptors; dense coef pointers alias coefficients. */
	struct optm_LinearConstraint *constraints;

	/**< n real-variable descriptors; defaults: x[j] >= 0. */
	struct optm_VariableBound *bounds;

	/**< CSC column offsets, length n+1. */
	int *column_start;

	/**< CSC zero-based row indices, length at least nnz. */
	int *row_index;

	/**< CSC matrix values, length at least nnz. */
	double *value;

	/**< CSR row offsets, length m+1. */
	int *row_start;

	/**< CSR zero-based column indices, length at least nnz. */
	int *column_index;

	/**< CSR matrix values, length at least nnz. */
	double *row_value;
};


/**
 * @brief Allocate a dense model with zero costs, matrix, RHS and empty names.
 *
 * Rows initially have equality type; variables are real with lower bound 0
 * and no upper bound. The objective is minimisation with offset 0. Sparse
 * arrays are initially NULL and nnz is 0.
 *
 * @param m Positive number of constraints.
 * @param n Positive number of variables.
 * @return Caller-owned model, or NULL for invalid dimensions/allocation failure.
 * @see lp_model_build_sparse(), lp_model_free()
 */
struct lp_Model *lp_model_create(int m, int n);

/**
 * @brief Allocate metadata and capacity for both sparse views, with no dense A.
 *
 * Defaults match lp_model_create(). Both start arrays are zeroed, nnz is 0,
 * and coefficient/index buffers are uninitialised. Fill CSC and CSR offsets,
 * indices and values, then set nnz to the used count, which must not exceed
 * nonzero_capacity. Fill row types, RHS, costs and variable bounds as needed.
 * Zero capacity is allowed: value/index buffers are then NULL.
 * lp_model_build_sparse() cannot be used on this sparse-only model.
 *
 * @param m Positive number of constraints.
 * @param n Positive number of variables.
 * @param nonzero_capacity Nonnegative entry capacity of each sparse view.
 * @return Caller-owned model, or NULL for invalid arguments/allocation failure.
 * @see lp_model_free()
 */
struct lp_Model *lp_model_create_sparse(int m, int n, int nonzero_capacity);

/**
 * @brief Release the model and all its owned arrays; NULL is a no-op.
 * @param model Model returned by a constructor or lp_read_mps(). All pointers
 * into it become invalid after this call.
 */
void lp_model_free(struct lp_Model *model);

/**
 * @brief Replace CSC and CSR arrays using the current dense constraint rows.
 *
 * Every coefficient unequal to zero is included, with no tolerance filtering.
 * On success nnz and both sparse views are updated; dense storage is retained.
 * Previous sparse pointers become invalid, including on allocation failure.
 * After failure, free the model or successfully rebuild before solving it.
 *
 * @param model Non-NULL dense model with its original row aliases intact.
 * @return 0 on success; -1 for NULL/sparse-only input or allocation failure.
 */
int lp_model_build_sparse(struct lp_Model *model);


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif
