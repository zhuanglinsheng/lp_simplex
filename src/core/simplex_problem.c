/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
/* Canonical sparse problem storage shared by presolve and simplex. */
#include "simplex_problem.h"
#include "utils.h"

#include <lp_simplex/status.h>


static int simplex_problem_allocate_vectors(
		struct simplex_Problem *problem, const int rows, const int columns)
{
	problem->objective = (double *)lp_simplex_malloc(
		(size_t)columns * sizeof(double));
	problem->bounds = (struct optm_VariableBound *)lp_simplex_malloc(
		(size_t)columns * sizeof(*problem->bounds));
	problem->rhs = (double *)lp_simplex_malloc((size_t)rows * sizeof(double));
	problem->row_type = (unsigned char *)lp_simplex_malloc(
		(size_t)rows * sizeof(unsigned char));
	return problem->objective != NULL && problem->bounds != NULL &&
		problem->rhs != NULL && problem->row_type != NULL
		? lp_simplex_EXIT_SUCCESS : lp_simplex_EXIT_FAILURE;
}


int simplex_problem_from_model(
		struct simplex_Problem *problem, const struct lp_Model *model)
{
	int i;
	lp_simplex_memset(problem, 0, sizeof(*problem));
	problem->rows = model->m;
	problem->columns = model->n;
	problem->owns_vectors = 1;
	if (simplex_problem_allocate_vectors(problem, model->m, model->n) ==
	    lp_simplex_EXIT_FAILURE ||
	    simplex_csc_from_model(model, &problem->matrix) ==
	    lp_simplex_EXIT_FAILURE) {
		simplex_problem_destroy(problem);
		return lp_simplex_EXIT_FAILURE;
	}
	lp_simplex_memcpy(problem->objective, model->objective,
		(size_t)model->n * sizeof(double));
	lp_simplex_memcpy(problem->bounds, model->bounds,
		(size_t)model->n * sizeof(*problem->bounds));
	for (i = 0; i < model->m; i++) {
		problem->rhs[i] = model->constraints[i].rhs;
		problem->row_type[i] =
			(unsigned char)model->constraints[i].type;
	}
	return lp_simplex_EXIT_SUCCESS;
}


struct simplex_Problem *simplex_problem_create_sparse(
		const int rows, const int columns, const int nonzeros)
{
	struct simplex_Problem *problem = (struct simplex_Problem *)
		lp_simplex_malloc(sizeof(*problem));
	if (problem == NULL)
		return NULL;
	lp_simplex_memset(problem, 0, sizeof(*problem));
	problem->rows = rows;
	problem->columns = columns;
	problem->owns_vectors = 1;
	problem->matrix.rows = rows;
	problem->matrix.columns = columns;
	problem->matrix.nonzeros = nonzeros;
	problem->matrix.owns_storage = 1;
	if (simplex_problem_allocate_vectors(problem, rows, columns) ==
	    lp_simplex_EXIT_FAILURE)
		goto failure;
	problem->matrix.column_start = (int *)lp_simplex_malloc(
		(size_t)(columns + 1) * sizeof(int));
	problem->matrix.row_start = (int *)lp_simplex_malloc(
		(size_t)(rows + 1) * sizeof(int));
	if (nonzeros > 0) {
		problem->matrix.row_index = (int *)lp_simplex_malloc(
			(size_t)nonzeros * sizeof(int));
		problem->matrix.value = (double *)lp_simplex_malloc(
			(size_t)nonzeros * sizeof(double));
		problem->matrix.column_index = (int *)lp_simplex_malloc(
			(size_t)nonzeros * sizeof(int));
		problem->matrix.row_value = (double *)lp_simplex_malloc(
			(size_t)nonzeros * sizeof(double));
	}
	if (problem->matrix.column_start == NULL ||
	    problem->matrix.row_start == NULL ||
	    (nonzeros > 0 && (problem->matrix.row_index == NULL ||
	     problem->matrix.value == NULL ||
	     problem->matrix.column_index == NULL ||
	     problem->matrix.row_value == NULL)))
		goto failure;
	return problem;
failure:
	simplex_problem_free(problem);
	return NULL;
}


void simplex_problem_destroy(struct simplex_Problem *problem)
{
	if (problem == NULL)
		return;
	simplex_csc_destroy(&problem->matrix);
	if (problem->owns_vectors) {
		lp_simplex_free(problem->objective);
		lp_simplex_free(problem->bounds);
		lp_simplex_free(problem->rhs);
		lp_simplex_free(problem->row_type);
	}
	lp_simplex_memset(problem, 0, sizeof(*problem));
}


void simplex_problem_free(struct simplex_Problem *problem)
{
	if (problem == NULL)
		return;
	simplex_problem_destroy(problem);
	lp_simplex_free(problem);
}
