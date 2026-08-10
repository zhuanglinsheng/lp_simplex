/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#ifndef LP_SIMPLEX_PROBLEM_INTERNAL_H
#define LP_SIMPLEX_PROBLEM_INTERNAL_H

#include "simplex_csc.h"

#include <lp_simplex/model.h>


/* Canonical immutable input to the optimization kernel.  Public models,
 * presolve output and future decomposition nodes all enter through this
 * boundary; the simplex code never depends on dense-model ownership rules. */
struct simplex_Problem {
	int rows;
	int columns;
	struct simplex_CscMatrix matrix;
	double *objective;
	struct optm_VariableBound *bounds;
	double *rhs;
	unsigned char *row_type;
	void *vector_storage;
	int owns_vectors;
};


int simplex_problem_from_model(
		struct simplex_Problem *problem, const struct lp_Model *model);

struct simplex_Problem *simplex_problem_create_sparse(
		int rows, int columns, int nonzeros);

void simplex_problem_destroy(struct simplex_Problem *problem);

void simplex_problem_free(struct simplex_Problem *problem);

#endif
