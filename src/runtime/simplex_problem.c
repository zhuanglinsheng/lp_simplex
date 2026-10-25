/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_problem.c
 * @brief Canonical sparse problem storage shared by presolve and simplex.
 */

#include "simplex_problem.h"
#include "simplex_csc_storage.h"
#include "utils.h"

#include <lp_simplex/status.h>

#include <math.h>
#include <limits.h>

void simplex_problem_mark_changed(struct simplex_Problem *problem)
{
	/* Saturation disables checkpoint matching rather than wrapping to a
	 * revision that an old checkpoint might still carry. */
	if (problem->revision < ULONG_MAX)
		problem->revision++;
}

struct problem_VectorStorage {
	size_t numeric, bounds, total;
};

static int problem_vector_storage(int rows, int columns, size_t limit,
	struct problem_VectorStorage *layout)
{
	struct problem_VectorStorage next = {0};
	size_t row, column;
	*layout = next;
	if (rows < 0 || columns < 0 ||
		!simplex_storage_bytes(
			(size_t)rows, sizeof(double), limit, &row) ||
		!simplex_storage_bytes(
			(size_t)columns, sizeof(double), limit, &column) ||
		row > limit - column ||
		!simplex_storage_bytes((size_t)columns,
			sizeof(struct optm_VariableBound), limit, &next.bounds))
		return 0;
	next.numeric = row + column;
	if (next.numeric > limit - next.bounds)
		return 0;
	next.total = next.numeric + next.bounds;
	if ((size_t)rows > limit - next.total)
		return 0;
	next.total += (size_t)rows;
	*layout = next;
	return 1;
}

static int simplex_problem_allocate_vectors(struct simplex_Problem *problem,
	const struct problem_VectorStorage *layout)
{
	unsigned char *storage =
		(unsigned char *)lp_simplex_malloc(layout->total);

	if (storage == NULL)
		return lp_simplex_EXIT_FAILURE;
	problem->vector_storage = storage;
	problem->objective = (double *)storage;
	problem->rhs = problem->objective + problem->columns;
	problem->bounds =
		(struct optm_VariableBound *)(storage + layout->numeric);
	problem->row_type = storage + layout->numeric + layout->bounds;
	return lp_simplex_EXIT_SUCCESS;
}

int simplex_problem_from_model(
	struct simplex_Problem *problem, const struct lp_Model *model)
{
	int i;
	struct problem_VectorStorage layout;

	lp_simplex_memset(problem, 0, sizeof(*problem));
	if (!problem_vector_storage(model->m, model->n, (size_t)-1, &layout))
		return lp_simplex_EXIT_FAILURE;
	problem->rows = model->m;
	problem->columns = model->n;
	problem->owns_vectors = 1;
	if (simplex_problem_allocate_vectors(problem, &layout) ==
			lp_simplex_EXIT_FAILURE ||
		simplex_csc_from_model(model, &problem->matrix) ==
			lp_simplex_EXIT_FAILURE) {
		simplex_problem_destroy(problem);
		return lp_simplex_EXIT_FAILURE;
	}
	lp_simplex_memcpy(problem->objective, model->objective,
		(size_t)model->n * sizeof(double));
	if (model->bounds != NULL)
		lp_simplex_memcpy(problem->bounds, model->bounds,
			(size_t)model->n * sizeof(*problem->bounds));
	else
		for (i = 0; i < model->n; i++) {
			problem->bounds[i].lb = 0.;
			problem->bounds[i].ub = __lp_simplex_INF__;
			problem->bounds[i].b_type = optm_BOUND_T_LO;
			problem->bounds[i].v_type = optm_VAR_T_REAL;
		}
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
	struct problem_VectorStorage layout;
	struct simplex_CscStorage arrays;
	struct simplex_Problem *problem;
	if (!problem_vector_storage(rows, columns, (size_t)-1, &layout) ||
		!simplex_csc_storage(
			rows, columns, nonzeros, (size_t)-1, &arrays))
		return NULL;
	problem = (struct simplex_Problem *)lp_simplex_malloc(sizeof(*problem));

	if (problem == NULL)
		return NULL;
	lp_simplex_memset(problem, 0, sizeof(*problem));
	problem->rows = rows;
	problem->columns = columns;
	problem->owns_vectors = 1;
	problem->matrix.rows = rows;
	problem->matrix.columns = columns;
	problem->matrix.nonzeros = nonzeros;
	problem->matrix.owns_storage = SIMPLEX_CSC_OWNS_ALL;
	if (simplex_problem_allocate_vectors(problem, &layout) ==
		lp_simplex_EXIT_FAILURE)
		goto failure;
	problem->matrix.column_start = (int *)lp_simplex_malloc(arrays.columns);
	problem->matrix.row_start = (int *)lp_simplex_malloc(arrays.rows);
	if (nonzeros > 0) {
		problem->matrix.row_index =
			(int *)lp_simplex_malloc(arrays.indices);
		problem->matrix.value =
			(double *)lp_simplex_malloc(arrays.values);
		problem->matrix.column_index =
			(int *)lp_simplex_malloc(arrays.indices);
		problem->matrix.row_value =
			(double *)lp_simplex_malloc(arrays.values);
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

struct simplex_Problem *simplex_problem_create_with_pattern(
	const struct simplex_CscMatrix *pattern)
{
	struct problem_VectorStorage layout;
	struct simplex_Problem *problem;
	size_t values;
	if (!problem_vector_storage(
		    pattern->rows, pattern->columns, (size_t)-1, &layout) ||
		pattern->nonzeros < 0 ||
		!simplex_storage_bytes((size_t)pattern->nonzeros,
			2 * sizeof(double), (size_t)-1, &values))
		return NULL;
	problem = (struct simplex_Problem *)lp_simplex_malloc(sizeof(*problem));

	if (problem == NULL)
		return NULL;
	lp_simplex_memset(problem, 0, sizeof(*problem));
	problem->rows = pattern->rows;
	problem->columns = pattern->columns;
	problem->owns_vectors = 1;
	problem->matrix = *pattern;
	problem->matrix.owns_storage = SIMPLEX_CSC_OWNS_VALUES;
	problem->matrix.value = NULL;
	problem->matrix.row_value = NULL;
	if (simplex_problem_allocate_vectors(problem, &layout) ==
		lp_simplex_EXIT_FAILURE)
		goto failure;
	if (pattern->nonzeros > 0) {
		/* CSC and CSR numeric arrays have the same lifetime.  A single
		 * allocation also avoids retaining two separately rounded heap blocks. */
		problem->matrix.value = (double *)lp_simplex_malloc(values);
		if (problem->matrix.value == NULL)
			goto failure;
		problem->matrix.row_value =
			problem->matrix.value + pattern->nonzeros;
	}
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
		lp_simplex_free(problem->vector_storage);
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

/* Check CSR rows directly instead of building a Pan representation merely to
 * measure feasibility. Compensated products protect cancellation-heavy rows. */
static double problem_primal_audit(const struct simplex_Problem *p,
	const double *x, double tolerance,
	struct simplex_PrimalSummary *summary)
{
	double maximum = 0.;
	int i, j, k;

	if (summary)
		summary->sum = summary->count = 0.;

	for (j = 0; j < p->columns; j++) {
		int type = p->bounds[j].b_type;
		double violation = 0.;

		if (!isfinite(x[j]))
			return HUGE_VAL;
		if ((type == optm_BOUND_T_LO || type == optm_BOUND_T_BS) &&
			x[j] < p->bounds[j].lb)
			violation = p->bounds[j].lb - x[j];
		if ((type == optm_BOUND_T_UP || type == optm_BOUND_T_BS) &&
			x[j] > p->bounds[j].ub)
			violation = fmax(violation, x[j] - p->bounds[j].ub);
		if (!isfinite(violation))
			return HUGE_VAL;
		maximum = fmax(maximum, violation);
		if (summary && violation > tolerance) {
			summary->sum += violation;
			summary->count += 1.;
		}
	}
	for (i = 0; i < p->rows; i++) {
		double sum = 0., correction = 0., violation;

		for (k = p->matrix.row_start[i]; k < p->matrix.row_start[i + 1];
			k++) {
			double a = p->matrix.row_value[k],
			       v = x[p->matrix.column_index[k]];
			double product = a * v, next = sum + product;

			if (!isfinite(product) || !isfinite(next))
				return HUGE_VAL;
			correction +=
				(fabs(sum) >= fabs(product)
						? (sum - next) + product
						: (product - next) + sum) +
				fma(a, v, -product);
			sum = next;
		}
		violation = (sum + correction) - p->rhs[i];
		if (!isfinite(violation))
			return HUGE_VAL;
		if (p->row_type[i] == optm_CONS_T_EQ)
			violation = fabs(violation);
		else if (p->row_type[i] == optm_CONS_T_GE)
			violation = fmax(0., -violation);
		else
			violation = fmax(0., violation);
		maximum = fmax(maximum, violation);
		if (summary && violation > tolerance) {
			summary->sum += violation;
			summary->count += 1.;
		}
	}
	return maximum;
}

double simplex_problem_primal_infeasibility(
	const struct simplex_Problem *p, const double *x)
{
	return problem_primal_audit(p, x, 0., NULL);
}

int simplex_problem_primal_summary(const struct simplex_Problem *p,
	const double *x, double tolerance,
	struct simplex_PrimalSummary *summary)
{
	if (!isfinite(tolerance) || tolerance < 0.)
		goto unavailable;
	summary->maximum = problem_primal_audit(p, x, tolerance, summary);
	if (isfinite(summary->maximum) && isfinite(summary->sum))
		return lp_simplex_EXIT_SUCCESS;
unavailable:
	summary->maximum = summary->sum = summary->count = NAN;
	return lp_simplex_EXIT_FAILURE;
}

static double problem_product_sum(
	const double *a, const int *index, int begin, int end, const double *x)
{
	double sum = 0., correction = 0.;
	int k;

	for (k = begin; k < end; k++) {
		double v = x[index[k]], product = a[k] * v,
		       next = sum + product;

		if (!isfinite(product) || !isfinite(next))
			return HUGE_VAL;
		correction +=
			(fabs(sum) >= fabs(product) ? (sum - next) + product
						    : (product - next) + sum) +
			fma(a[k], v, -product);
		sum = next;
	}
	return sum + correction;
}

double simplex_problem_dual_infeasibility(const struct simplex_Problem *p,
	const double *x, const double *y, double tolerance, int sign)
{
	double maximum = 0.;
	int i, j;

	for (i = 0; i < p->rows; i++)
		if (!isfinite(y[i]))
			return HUGE_VAL;
	for (j = 0; j < p->columns; j++) {
		const struct optm_VariableBound *b = p->bounds + j;
		double g = p->objective[j] +
			   sign * problem_product_sum(p->matrix.value,
					  p->matrix.row_index,
					  p->matrix.column_start[j],
					  p->matrix.column_start[j + 1], y);
		double v;
		int lower = (b->b_type == optm_BOUND_T_LO ||
			     b->b_type == optm_BOUND_T_BS);
		int upper = (b->b_type == optm_BOUND_T_UP ||
			     b->b_type == optm_BOUND_T_BS);

		if (!isfinite(x[j]) || !isfinite(g))
			return HUGE_VAL;
		if (lower && upper && b->lb == b->ub)
			v = 0.;
		else if (lower && x[j] <= b->lb + tolerance)
			v = fmax(0., -g);
		else if (upper && x[j] >= b->ub - tolerance)
			v = fmax(0., g);
		else
			v = fabs(g);
		maximum = fmax(maximum, v);
	}
	for (i = 0; i < p->rows; i++) {
		double lambda = sign * y[i], v = 0.;
		double activity = problem_product_sum(p->matrix.row_value,
			p->matrix.column_index, p->matrix.row_start[i],
			p->matrix.row_start[i + 1], x);

		if (!isfinite(activity))
			return HUGE_VAL;
		if (p->row_type[i] == optm_CONS_T_LE)
			v = p->rhs[i] - activity > tolerance
				    ? fabs(lambda)
				    : fmax(0., -lambda);
		else if (p->row_type[i] == optm_CONS_T_GE)
			v = activity - p->rhs[i] > tolerance ? fabs(lambda)
							     : fmax(0., lambda);
		maximum = fmax(maximum, v);
	}
	return maximum;
}
