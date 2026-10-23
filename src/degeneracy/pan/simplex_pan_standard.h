/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_pan_standard.h
 * @brief Sparse equality-form conversion and original-variable recovery for
 * Pan/BDA.
 */

#ifndef LP_SIMPLEX_PAN_STANDARD_INTERNAL_H
#define LP_SIMPLEX_PAN_STANDARD_INTERNAL_H

#include "simplex_problem.h"


/* Equality-form problem used by Pan's deficient-basis algorithm:
 *
 *     minimize   objective' z + objective_offset
 *     subject to matrix z = rhs, z >= 0.
 *
 * Structural transformations are recorded so z can be recovered without
 * materializing a dense copy of the public model. */
struct simplex_PanStandard {
	int rows;
	int columns;
	int original_columns;
	int nonzeros;
	int *column_start;
	int *row_index;
	double *value;
	double *normalized_value;
	double *column_norm;
	double *rhs;
	double *objective;
	int *original_column;
	double *original_scale;
	double *original_shift;
	double objective_offset;
	/* Optional original-objective audit for an equivalent dual-centered view.
	 * The pointed-to coefficients are borrowed from the enclosing LP. */
	const double *audit_objective;
	double audit_objective_offset;
};


int simplex_pan_standard_create(struct simplex_PanStandard *standard,
	const struct simplex_Problem *problem);

void simplex_pan_standard_destroy(struct simplex_PanStandard *standard);

void simplex_pan_standard_recover(
	const struct simplex_PanStandard *standard, const double *z, double *x);

/* Encode an original-space point, construct inequality/upper slacks, and
 * verify nonnegativity and every equality against the actual sparse matrix. */
int simplex_pan_standard_encode(const struct simplex_PanStandard *standard,
	const double *x, double *z, double *residual, double tolerance);

/* Compensated sparse dot, including FMA product residuals. Does not rely on
 * long double having more precision than double on the host platform. */
double simplex_pan_column_dot(const struct simplex_PanStandard *standard,
	int column, const double *vector);

/* Componentwise residual with compensated sums and product errors. */
double simplex_pan_residual(const struct simplex_PanStandard *standard,
	const double *value, double *residual, double *correction);

#endif
