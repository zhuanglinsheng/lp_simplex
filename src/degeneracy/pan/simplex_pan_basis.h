/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
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
	long extensions;
	long downdates;
	long rotations;
	long symbolic_updates;
	long anchor_rebuilds;
	long local_retriangularizations;
	long block_householders;
	long maximum_front;
	long factor_nonzeros;
	double orthogonal_error;
	double backward_error;
	double factor_flops;
};


int simplex_pan_basis_create(
		struct simplex_PanBasis *basis,
		const struct simplex_PanStandard *standard,
		double rank_tolerance);

void simplex_pan_basis_destroy(struct simplex_PanBasis *basis);

int simplex_pan_basis_add(struct simplex_PanBasis *basis, int column);

void simplex_pan_basis_remove_marked(
		struct simplex_PanBasis *basis, const unsigned char *remove);

int simplex_pan_basis_factorize(struct simplex_PanBasis *basis);

/* Replace one active column without changing the basis dimension. */
int simplex_pan_basis_replace_factorized(
		struct simplex_PanBasis *basis, int position, int column);

int simplex_pan_basis_least_squares(
		struct simplex_PanBasis *basis, const double *right,
		double *solution, double *residual_norm);

/* Rank-revealing QR fallback used for certification and recovery when the
 * corrected semi-normal equations cannot meet the requested tolerance. */
int simplex_pan_basis_least_squares_qr(
		struct simplex_PanBasis *basis, const double *right,
		double *solution, double *residual_norm);

int simplex_pan_basis_minimum_norm(
		struct simplex_PanBasis *basis, const double *right,
		double *solution, double *residual_norm);

int simplex_pan_basis_column_projection(
		struct simplex_PanBasis *basis, int column,
		double *coefficient, double *residual_norm);

double simplex_pan_column_norm(
		const struct simplex_PanStandard *standard, int column);

#endif
