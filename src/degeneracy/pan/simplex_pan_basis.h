/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_pan_basis.h
 * @brief Dynamic deficient-basis storage and numerical operations for Pan/BDA.
 */

#ifndef LP_SIMPLEX_PAN_BASIS_INTERNAL_H
#define LP_SIMPLEX_PAN_BASIS_INTERNAL_H

#include "simplex_pan_standard.h"


/* Full-column-rank, possibly deficient basis.  SuiteSparseQR is the preferred
 * backend; a corrected semi-normal-equation implementation remains available
 * when SPQR is not present at build time. */
struct simplex_PanBasis {
	const struct simplex_PanStandard *standard;
	int *column;
	int *position;
	int count;
	int capacity;
	int nnls_block_limit, block_crash, singleton_crash, full_refinement;
	int implicit_bulk;
	int prefer_updates; /* Per-solve implicit cone policy. */
	double *factor;
	double *column_norm;
	double *basis_work;
	double *row_work;
	double *right_work;
	void *sparse_backend;
	double rank_tolerance;
	long factorizations;
	long updates;
	long refinements;
	long rank_repairs;
	long stability_drops;
	double reciprocal_condition;
	double projection_solve_error;
	long rank_disagreements;
	long extensions;
	long downdates;
	long rotations;
	long symbolic_updates;
	long anchor_rebuilds;
	long local_retriangularizations;
	long block_householders;
	long maximum_front;
	long factor_nonzeros;
	long numeric_packs;
	long numeric_packed_nodes;
	long border_solves;
	double orthogonal_error;
	double backward_error;
	double factor_flops;
	double factor_seconds;
	double least_squares_seconds;
	double minimum_norm_seconds;
	double projection_seconds;
};


int simplex_pan_basis_create(struct simplex_PanBasis *basis,
	const struct simplex_PanStandard *standard, double rank_tolerance);

void simplex_pan_basis_destroy(struct simplex_PanBasis *basis);

int simplex_pan_basis_add(struct simplex_PanBasis *basis, int column);

void simplex_pan_basis_remove_marked(
	struct simplex_PanBasis *basis, const unsigned char *remove);

int simplex_pan_basis_factorize(struct simplex_PanBasis *basis);

/* Unconditionally build a fresh rank-revealing factor from the active CSC
 * columns. Also used to bootstrap a batch of columns without serial updates. */
int simplex_pan_basis_rebuild(struct simplex_PanBasis *basis);

/* Rebuild with rank revelation, retaining all positive-valued columns.
 * Only exactly-zero active columns may be discarded. */
int simplex_pan_basis_repair_rank(
	struct simplex_PanBasis *basis, const double *value);

/* Optional one-shot independent completion from the whole sparse matrix,
 * retaining the positive support of the current feasible point. */
int simplex_pan_basis_complete(
	struct simplex_PanBasis *basis, const double *value);
/* Rank-reveal a supplied candidate set on an empty basis. Only this call
 * may see an overcomplete list; no triangular solve precedes rank selection.
 * required protects caller-designated columns, independently of primal mass. */
int simplex_pan_basis_select_columns(struct simplex_PanBasis *basis,
	const int *column, int count, const double *required);

/* Drop only zero-valued working columns until a sparse triangular condition
 * estimate is usable. Keep the newest entering column as well as all positive
 * columns. This is a working-set stability heuristic, not an LP certificate. */
int simplex_pan_basis_stabilize(struct simplex_PanBasis *basis,
	const double *value, int protected_column);

/* Cheap pricing screens: a unit vector orthogonal to B, and a normalized
 * multiplier with positive B^T components only on zero-valued active columns. */
int simplex_pan_basis_pricing_probes(struct simplex_PanBasis *basis,
	const double *value, double *span_probe, double *edge_probe);

/* Replace one active column without changing the basis dimension. */
int simplex_pan_basis_replace_factorized(
	struct simplex_PanBasis *basis, int position, int column);

int simplex_pan_basis_least_squares(struct simplex_PanBasis *basis,
	const double *right, double *solution, double *residual_norm);

/* Rank-revealing QR fallback used for certification and recovery when the
 * corrected semi-normal equations cannot meet the requested tolerance. */
int simplex_pan_basis_least_squares_qr(struct simplex_PanBasis *basis,
	const double *right, double *solution, double *residual_norm);

int simplex_pan_basis_minimum_norm(struct simplex_PanBasis *basis,
	const double *right, double *solution, double *residual_norm);

int simplex_pan_basis_column_projection(struct simplex_PanBasis *basis,
	int column, double *coefficient, double *residual_norm);

double simplex_pan_column_norm(
	const struct simplex_PanStandard *standard, int column);

/* Rectangular least-squares Newton correction for the dual hinge objective.
 * Factors the transpose of all free/violated normalized columns, without an
 * objective row or normal equations. Returns failure when SPQR is absent. */
int simplex_pan_dual_newton_qr(const struct simplex_PanStandard *standard,
	const double *point, const unsigned char *free_mask, const double *map,
	const double *multiplier, double scale, double *correction);

#endif
