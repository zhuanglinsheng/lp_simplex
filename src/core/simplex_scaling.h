/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_scaling.h
 * @brief Reversible power-of-two equilibration for sparse LP kernels.
 */

#ifndef LP_SIMPLEX_SCALING_INTERNAL_H
#define LP_SIMPLEX_SCALING_INTERNAL_H

#include "simplex_problem.h"


struct simplex_Scaling {
	int rows;
	int columns;
	double *storage;
	double *row;
	double *column;
	int active;
};


/* The returned immutable problem borrows the input's sparsity pattern, and
 * all input storage when scaling is inactive.  Free it before the input. */
struct simplex_Problem *simplex_scaling_create_problem(
	const struct simplex_Problem *problem, struct simplex_Scaling *scaling);

void simplex_scaling_recover_primal(
	const struct simplex_Scaling *scaling, double *x);

void simplex_scaling_recover_dual(
	const struct simplex_Scaling *scaling, double *dual);

void simplex_scaling_destroy(struct simplex_Scaling *scaling);

#endif
