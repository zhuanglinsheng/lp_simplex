/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex.c
 * @brief Public solve dispatch, input validation and original-model result
 * checks.
 */

#include "simplex_dual.h"
#include "simplex_progress.h"
#include "simplex_pan.h"
#include "simplex_auto.h"
#include "simplex_pan_phase1.h"
#include "simplex_pan_basis.h"
#include "simplex_presolve.h"
#include "simplex_singleton_dual.h"
#include "simplex_tableau_solver.h"
#include "utils.h"

#include <lp_simplex/solve.h>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>


#define SIMPLEX_DEFAULT_ITERATION_LIMIT 100000
#define SIMPLEX_DEFAULT_PRIMAL_TOLERANCE 1e-7
#define SIMPLEX_DEFAULT_DUAL_TOLERANCE 1e-7
#define SIMPLEX_DEFAULT_PIVOT_TOLERANCE 1e-9


void lp_simplex_default_options(
	struct lp_simplex_Options *options, const int algorithm)
{
	if (options == NULL)
		return;

	options->algorithm = algorithm;
	options->pricing = (algorithm == lp_simplex_ALGORITHM_DUAL_REVISED ||
				   algorithm == lp_simplex_ALGORITHM_AUTO)
				   ? lp_simplex_PRICING_DUAL_STEEPEST_EDGE
			   : algorithm == lp_simplex_ALGORITHM_PAN_BDA
				   ? lp_simplex_PRICING_DANTZIG
				   : lp_simplex_PRICING_BLAND;
	options->iteration_limit =
		algorithm == lp_simplex_ALGORITHM_AUTO
			? LP_SIMPLEX_DEFAULT_AUTO_ITERATION_LIMIT
			: SIMPLEX_DEFAULT_ITERATION_LIMIT;
	options->presolve = 1;
	options->progress_callback = NULL;
	options->progress_context = NULL;
	options->progress_interval = 5.;
	options->primal_tolerance = SIMPLEX_DEFAULT_PRIMAL_TOLERANCE;
	options->dual_tolerance = SIMPLEX_DEFAULT_DUAL_TOLERANCE;
	options->pivot_tolerance = SIMPLEX_DEFAULT_PIVOT_TOLERANCE;
}


static int simplex_validate_sparse_storage(const struct lp_Model *model)
{
	int i, k;

	if (model->nnz < 0 || model->column_start == NULL ||
		model->row_start == NULL ||
		(model->nnz > 0 &&
			(model->row_index == NULL || model->value == NULL ||
				model->column_index == NULL ||
				model->row_value == NULL)) ||
		model->column_start[0] != 0 ||
		model->column_start[model->n] != model->nnz ||
		model->row_start[0] != 0 ||
		model->row_start[model->m] != model->nnz)
		return 0;
	for (i = 0; i < model->n; i++)
		if (model->column_start[i] > model->column_start[i + 1])
			return 0;
	for (i = 0; i < model->m; i++)
		if (model->row_start[i] > model->row_start[i + 1])
			return 0;
	for (k = 0; k < model->nnz; k++)
		if (model->row_index[k] < 0 ||
			model->row_index[k] >= model->m ||
			model->column_index[k] < 0 ||
			model->column_index[k] >= model->n)
			return 0;
	return 1;
}


static int simplex_validate_model(const struct lp_Model *model)
{
	int i, j;
	int has_sparse, sparse = 0;

	if (model == NULL || model->objective == NULL ||
		model->constraints == NULL || model->m <= 0 || model->n <= 0)
		return 0;
	has_sparse = model->column_start != NULL || model->row_start != NULL ||
		     model->row_index != NULL || model->value != NULL ||
		     model->column_index != NULL || model->row_value != NULL ||
		     model->nnz != 0;
	if (has_sparse) {
		if (!simplex_validate_sparse_storage(model))
			return 0;
		sparse = 1;
	}
	for (i = 0; i < model->m; i++) {
		if ((model->constraints[i].coef == NULL && !sparse) ||
			model->constraints[i].type < optm_CONS_T_EQ ||
			model->constraints[i].type > optm_CONS_T_LE)
			return 0;
	}
	if (model->bounds == NULL)
		return 1;
	for (j = 0; j < model->n; j++) {
		if (model->bounds[j].b_type < optm_BOUND_T_FR ||
			model->bounds[j].b_type > optm_BOUND_T_BS ||
			model->bounds[j].v_type != optm_VAR_T_REAL ||
			(model->bounds[j].b_type == optm_BOUND_T_BS &&
				model->bounds[j].lb > model->bounds[j].ub))
			return 0;
	}
	return 1;
}


static struct lp_Model *simplex_materialize_dense(const struct lp_Model *model)
{
	struct lp_Model *dense = lp_model_create(model->m, model->n);
	int i, j, k;

	if (dense == NULL)
		return NULL;
	lp_simplex_memcpy(dense->objective, model->objective,
		(size_t)model->n * sizeof(double));
	dense->objective_offset = model->objective_offset;
	dense->objective_sense = model->objective_sense;
	if (model->bounds != NULL)
		lp_simplex_memcpy(dense->bounds, model->bounds,
			(size_t)model->n * sizeof(*dense->bounds));
	for (i = 0; i < model->m; i++) {
		lp_simplex_memcpy(dense->constraints[i].name,
			model->constraints[i].name,
			sizeof(dense->constraints[i].name));
		dense->constraints[i].rhs = model->constraints[i].rhs;
		dense->constraints[i].type = model->constraints[i].type;
	}
	if (model->column_start != NULL)
		for (j = 0; j < model->n; j++)
			for (k = model->column_start[j];
				k < model->column_start[j + 1]; k++)
				dense->constraints[model->row_index[k]]
					.coef[j] = model->value[k];
	else
		for (i = 0; i < model->m; i++)
			lp_simplex_memcpy(dense->constraints[i].coef,
				model->constraints[i].coef,
				(size_t)model->n * sizeof(double));
	return dense;
}


static int simplex_validate_options(const struct lp_simplex_Options *options)
{
	if (options == NULL || options->iteration_limit <= 0 ||
		(options->presolve != 0 && options->presolve != 1) ||
		!isfinite(options->primal_tolerance) ||
		!isfinite(options->dual_tolerance) ||
		!isfinite(options->pivot_tolerance) ||
		options->primal_tolerance <= 0. ||
		options->dual_tolerance <= 0. ||
		options->pivot_tolerance <= 0. ||
		(options->progress_callback &&
			(!isfinite(options->progress_interval) ||
				options->progress_interval <= 0.)))
		return 0;
	if (options->algorithm == lp_simplex_ALGORITHM_TABLEAU)
		return options->pricing == lp_simplex_PRICING_BLAND ||
		       options->pricing == lp_simplex_PRICING_DANTZIG;
	if (options->algorithm == lp_simplex_ALGORITHM_DUAL_REVISED ||
		options->algorithm == lp_simplex_ALGORITHM_AUTO)
		return options->pricing ==
		       lp_simplex_PRICING_DUAL_STEEPEST_EDGE;
	if (options->algorithm == lp_simplex_ALGORITHM_PAN_BDA)
		return options->pricing == lp_simplex_PRICING_DANTZIG ||
		       options->pricing == lp_simplex_PRICING_PAN_NORMALIZED;
	return 0;
}


static double simplex_primal_infeasibility(
	const struct lp_Model *model, const double *x)
{
	double maximum = 0.;
	int i, j, k;

	for (j = 0; j < model->n; j++) {
		double violation = 0.;
		/* A non-finite primal entry cannot be feasible.  Testing it here
		 * also keeps NaN out of the comparisons below, where every
		 * "violation > maximum" test would silently be false. */
		if (!__lp_simplex_ISFINITE__(x[j]))
			return __lp_simplex_INF__;
		if (model->bounds == NULL) {
			if (x[j] < 0.)
				violation = -x[j];
		} else {
			int type = model->bounds[j].b_type;

			if ((type == optm_BOUND_T_LO ||
				    type == optm_BOUND_T_BS) &&
				x[j] < model->bounds[j].lb)
				violation = model->bounds[j].lb - x[j];
			if ((type == optm_BOUND_T_UP ||
				    type == optm_BOUND_T_BS) &&
				x[j] > model->bounds[j].ub &&
				x[j] - model->bounds[j].ub > violation)
				violation = x[j] - model->bounds[j].ub;
		}
		if (!__lp_simplex_ISFINITE__(violation))
			return __lp_simplex_INF__;
		if (violation > maximum)
			maximum = violation;
	}
	for (i = 0; i < model->m; i++) {
		double activity = 0., correction = 0.;
		double violation;

		if (!__lp_simplex_ISFINITE__(model->constraints[i].rhs))
			return __lp_simplex_INF__;
		if (model->row_start != NULL) {
			for (k = model->row_start[i];
				k < model->row_start[i + 1]; k++) {
				double a = model->row_value[k],
				       v = x[model->column_index[k]];
				double product = a * v,
				       next = activity + product;

				if (!isfinite(product) || !isfinite(next))
					return HUGE_VAL;
				correction +=
					(fabs(activity) >= fabs(product)
							? (activity - next) +
								  product
							: (product - next) +
								  activity) +
					fma(a, v, -product);
				activity = next;
			}
		} else {
			for (j = 0; j < model->n; j++) {
				double a = model->constraints[i].coef[j],
				       product = a * x[j],
				       next = activity + product;

				if (!isfinite(product) || !isfinite(next))
					return HUGE_VAL;
				correction +=
					(fabs(activity) >= fabs(product)
							? (activity - next) +
								  product
							: (product - next) +
								  activity) +
					fma(a, x[j], -product);
				activity = next;
			}
		}
		activity += correction;
		if (!__lp_simplex_ISFINITE__(activity))
			return __lp_simplex_INF__;
		violation = activity - model->constraints[i].rhs;
		if (model->constraints[i].type == optm_CONS_T_EQ)
			violation = __lp_simplex_ABS__(violation);
		else if (model->constraints[i].type == optm_CONS_T_GE)
			violation = violation < 0. ? -violation : 0.;
		else
			violation = violation > 0. ? violation : 0.;
		if (!__lp_simplex_ISFINITE__(violation))
			return __lp_simplex_INF__;
		if (violation > maximum)
			maximum = violation;
	}
	return maximum;
}


static double simplex_original_objective(
	const struct lp_Model *model, const double *x)
{
	double sum = 0., correction = 0.;
	int j;

	for (j = 0; j < model->n; j++) {
		double a = model->objective[j], product = a * x[j],
		       next = sum + product;

		if (!isfinite(x[j]) || !isfinite(product) || !isfinite(next))
			return NAN;
		correction +=
			(fabs(sum) >= fabs(product) ? (sum - next) + product
						    : (product - next) + sum) +
			fma(a, x[j], -product);
		sum = next;
	}
	return ((sum + correction) + model->objective_offset) *
	       model->objective_sense;
}


/**
 * Single gate in front of every "solved" return.
 *
 * A solver may only claim `lp_simplex_Success` when the vector it hands back
 * is finite and feasible for the caller's own model.  Without this check the
 * presolve path can report success while its reduced problem overflowed
 * (objective inf/NaN) and the tableau path reports success while its
 * recovered vector violates the constraints outright.
 */
static int simplex_finalize_result(const struct lp_Model *model,
	const struct lp_simplex_Options *options, const double *x,
	struct lp_simplex_Result *result, int state)
{
	result->primal_infeasibility = simplex_primal_infeasibility(model, x);

	result->objective = simplex_original_objective(model, x);
	if (result->dual_residual_space == lp_simplex_RESIDUAL_UNAVAILABLE)
		result->dual_infeasibility = NAN;
	if (result->status == lp_simplex_Success &&
		(state != lp_simplex_EXIT_SUCCESS ||
			!isfinite(result->objective) ||
			!isfinite(result->primal_infeasibility) ||
			(result->dual_residual_space !=
					lp_simplex_RESIDUAL_UNAVAILABLE &&
				!isfinite(result->dual_infeasibility)) ||
			result->primal_infeasibility >
				options->primal_tolerance))
		result->status = lp_simplex_PrecisionError;
	return result->status == lp_simplex_Success ? lp_simplex_EXIT_SUCCESS
						    : lp_simplex_EXIT_FAILURE;
}


static int simplex_solve_dual_raw(const struct lp_Model *model,
	const struct lp_simplex_Options *options, double *x,
	struct lp_simplex_Result *result,
	const struct simplex_ProgressSink *progress)
{
	int applicable = 0, state;

	if (model->coefficients != NULL) {
		state = simplex_singleton_dual_solve(
			model, options, x, result, &applicable);
		if (applicable)
			return state;
	}
	/* Raw solves are also the numerical fallback for a rejected presolve
	 * result.  Do not derive additional working bounds here: on unscaled models
	 * with extreme coefficient ranges (ETAMACRO/PILOT) repeated activity
	 * cancellation can make those optional crash bounds less reliable than the
	 * original model.  Presolved problems retain propagation below. */
	if (progress) {
		struct simplex_Problem problem;
		if (simplex_problem_from_model(&problem, model)) {
			result->status = lp_simplex_MemoryAllocError;
			return lp_simplex_EXIT_FAILURE;
		}
		problem.progress = progress;
		state = simplex_dual_solve_problem(
			&problem, options, x, NULL, result, 0);
		simplex_problem_destroy(&problem);
		return state;
	}
	return simplex_dual_solve(model, options, x, NULL, result, 0);
}


static int simplex_solve_pan_raw(const struct lp_Model *model,
	const struct lp_simplex_Options *options, double *x,
	struct lp_simplex_Result *result,
	const struct simplex_ProgressSink *progress)
{
	struct simplex_Problem problem;
	int state;

	if (simplex_problem_from_model(&problem, model) ==
		lp_simplex_EXIT_FAILURE) {
		result->status = lp_simplex_MemoryAllocError;
		return lp_simplex_EXIT_FAILURE;
	}
	problem.progress = progress;
	state = options->algorithm == lp_simplex_ALGORITHM_AUTO
			? simplex_auto_solve_problem(
				  &problem, options, x, result, 0)
			: simplex_pan_solve_problem(
				  &problem, options, x, result);
	simplex_problem_destroy(&problem);
	return state;
}


static int simplex_solve_sparse_raw(const struct lp_Model *model,
	const struct lp_simplex_Options *options, double *x,
	struct lp_simplex_Result *result,
	const struct simplex_ProgressSink *progress)
{
	return (options->algorithm == lp_simplex_ALGORITHM_PAN_BDA ||
		       options->algorithm == lp_simplex_ALGORITHM_AUTO)
		       ? simplex_solve_pan_raw(
				 model, options, x, result, progress)
		       : simplex_solve_dual_raw(
				 model, options, x, result, progress);
}


/* Presolve can remove structural rows yet add many implied finite bounds.
 * Pan currently expands every width into an equality and a slack. Compare
 * the actual equality-form dimensions before choosing its input model. */
static int simplex_pan_original_form_smaller(
	const struct lp_Model *model, const struct simplex_Problem *reduced)
{
	size_t original_rows = (size_t)model->m,
	       original_columns = (size_t)model->n;
	size_t reduced_rows = (size_t)reduced->rows,
	       reduced_columns = (size_t)reduced->columns;
	int i, j;

	for (i = 0; i < model->m; i++)
		if (model->constraints[i].type != optm_CONS_T_EQ)
			original_columns++;
	for (i = 0; i < reduced->rows; i++)
		if (reduced->row_type[i] != optm_CONS_T_EQ)
			reduced_columns++;
	for (j = 0; j < model->n; j++) {
		const struct optm_VariableBound *bound = model->bounds + j;

		if (bound->b_type == optm_BOUND_T_FR)
			original_columns++;
		if (bound->b_type == optm_BOUND_T_BS) {
			if (bound->ub == bound->lb)
				original_columns--;
			else if (bound->ub > bound->lb) {
				original_rows++;
				original_columns++;
			}
		}
	}
	for (j = 0; j < reduced->columns; j++) {
		const struct optm_VariableBound *bound = reduced->bounds + j;

		if (bound->b_type == optm_BOUND_T_FR)
			reduced_columns++;
		if (bound->b_type == optm_BOUND_T_BS) {
			if (bound->ub == bound->lb)
				reduced_columns--;
			else if (bound->ub > bound->lb) {
				reduced_rows++;
				reduced_columns++;
			}
		}
	}
	if (reduced_rows > original_rows &&
		reduced_columns >= original_columns) {
		if (getenv("LP_SIMPLEX_PROFILE"))
			fprintf(stderr,
				"presolve: Pan equality form original=%lu/%lu reduced=%lu/%lu selecting original\n",
				(unsigned long)original_rows,
				(unsigned long)original_columns,
				(unsigned long)reduced_rows,
				(unsigned long)reduced_columns);
		return 1;
	}
	return 0;
}


/* Presolve may replace a cheap equality feasibility problem by many finite
 * width rows. Compare that expansion before choosing an initializer. NNLS
 * returns only a candidate: validate both the original and reduced models. */
static int simplex_auto_original_seed(const struct lp_Model *model,
	const struct simplex_Presolve *presolve,
	const struct lp_simplex_Options *o, double *original, double *reduced,
	int *iterations)
{
	struct simplex_Problem raw;
	struct simplex_PanStandard a;
	struct simplex_PanBasis basis;
	const struct simplex_Problem *p = presolve->reduced;
	size_t expanded_rows = (size_t)p->rows;
	double *z;
	int j, limit, status, state, valid = 0;

	*iterations = 0;
	for (j = 0; j < p->columns; j++)
		if (p->bounds[j].b_type == optm_BOUND_T_BS &&
			p->bounds[j].ub > p->bounds[j].lb)
			expanded_rows++;
	if (expanded_rows <= 2 * (size_t)model->m ||
		o->iteration_limit < SIMPLEX_AUTO_PROBE_CAP ||
		!simplex_pan_original_form_smaller(model, p) ||
		simplex_problem_from_model(&raw, model))
		return 0;
	if (!simplex_auto_route_allowed(&raw, SIMPLEX_AUTO_CLASSIC)) {
		simplex_problem_destroy(&raw);
		return 0;
	}
	if (simplex_pan_standard_create(&a, &raw)) {
		simplex_problem_destroy(&raw);
		return 0;
	}
	if (simplex_pan_basis_create(&basis, &a, o->pivot_tolerance)) {
		simplex_pan_standard_destroy(&a);
		simplex_problem_destroy(&raw);
		return 0;
	}
	z = (double *)lp_simplex_malloc(
		((size_t)a.columns + 1) * sizeof(double));
	if (z) {
		lp_simplex_memset(
			z, 0, ((size_t)a.columns + 1) * sizeof(double));
		basis.singleton_crash = 1;
		limit = o->iteration_limit / 16;
		if (limit > 8 * SIMPLEX_AUTO_PROBE_CAP)
			limit = 8 * SIMPLEX_AUTO_PROBE_CAP;
		state = simplex_pan_phase1_run(&basis, z, o->primal_tolerance,
			o->dual_tolerance, limit, iterations, &status);
		simplex_pan_standard_recover(&a, z, original);
		if (!state && simplex_problem_primal_infeasibility(
				      &raw, original) <= o->primal_tolerance) {
			for (j = 0; j < p->columns; j++)
				reduced[j] = original[presolve->column_map[j]];
			valid = simplex_problem_primal_infeasibility(
					p, reduced) <= o->primal_tolerance;
		}
		lp_simplex_free(z);
	}
	if (getenv("LP_SIMPLEX_PROFILE"))
		fprintf(stderr,
			"auto: original-feasibility expanded-rows=%lu original-rows=%d work=%d accepted=%d\n",
			(unsigned long)expanded_rows, raw.rows, *iterations,
			valid);
	simplex_pan_basis_destroy(&basis);
	simplex_pan_standard_destroy(&a);
	simplex_problem_destroy(&raw);
	return valid;
}


struct simplex_AutoSeedContext {
	const struct lp_Model *model;
	const struct simplex_Presolve *presolve;
	double *original;
};
static int simplex_auto_initialize(void *context,
	const struct lp_simplex_Options *o, double *seed, int *work)
{
	const struct simplex_AutoSeedContext *c =
		(const struct simplex_AutoSeedContext *)context;
	return simplex_auto_original_seed(
		c->model, c->presolve, o, c->original, seed, work);
}


static int simplex_solve_sparse_presolved(const struct lp_Model *model,
	const struct lp_simplex_Options *options, double *x,
	struct lp_simplex_Result *result,
	const struct simplex_ProgressSink *progress)
{
	struct simplex_Presolve presolve;
	struct simplex_ProgressBridge bridge;
	double *reduced_x = NULL;
	clock_t presolve_started;
	int j, state;

	if (!options->presolve ||
		getenv("LP_SIMPLEX_DISABLE_PRESOLVE") != NULL ||
		model->bounds == NULL) {
		state = simplex_solve_sparse_raw(
			model, options, x, result, progress);
		return state;
	}
	lp_simplex_memset(&bridge, 0, sizeof(bridge));
	presolve_started = clock();
	if (simplex_presolve_run(&presolve, model, options->primal_tolerance) ==
		lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	if (getenv("LP_SIMPLEX_PROFILE") != NULL)
		simplex_presolve_print_profile(
			&presolve, (double)(clock() - presolve_started) /
					   (double)CLOCKS_PER_SEC);
	if (presolve.status != lp_simplex_CondUnsatisfied) {
		simplex_presolve_postsolve(&presolve, NULL, x);
		result->status = presolve.status;
		state = presolve.status == lp_simplex_Success
				? lp_simplex_EXIT_SUCCESS
				: lp_simplex_EXIT_FAILURE;
	} else if (options->algorithm == lp_simplex_ALGORITHM_PAN_BDA &&
		   getenv("LP_SIMPLEX_PAN_NATIVE_BOUNDS") == NULL &&
		   ((presolve.reduced->rows == model->m &&
			    presolve.reduced->columns == model->n) ||
			   simplex_pan_original_form_smaller(
				   model, presolve.reduced))) {
		/* With no structural reductions, implied bound tightening only adds
		 * redundant width rows to Pan's equality form. Keep the original LP. */
		state = simplex_solve_pan_raw(
			model, options, x, result, progress);
	} else {
		reduced_x = (double *)lp_simplex_malloc(
			(size_t)presolve.reduced->columns * sizeof(double));
		if (reduced_x == NULL) {
			simplex_presolve_destroy(&presolve);
			return lp_simplex_EXIT_FAILURE;
		}
		for (j = 0; j < presolve.reduced->columns; j++)
			reduced_x[j] = NAN;
		if (simplex_progress_bridge(&bridge, progress, &presolve, NULL,
			    model->m, model->n)) {
			lp_simplex_free(reduced_x);
			simplex_presolve_destroy(&presolve);
			result->status = lp_simplex_MemoryAllocError;
			return lp_simplex_EXIT_FAILURE;
		}
		presolve.reduced->progress = progress ? &bridge.sink : NULL;
		/* Presolve uses outward-rounded bounds for safe elimination.  The crash
		 * phase then derives exact working bounds on the reduced matrix; these
		 * bounds guide initialization but do not authorize more elimination. */
		if (options->algorithm == lp_simplex_ALGORITHM_AUTO) {
			struct simplex_AutoSeedContext context;
			struct simplex_AutoInitializer initializer;

			context.model = model;
			context.presolve = &presolve;
			context.original = x;
			initializer.context = &context;
			initializer.run = simplex_auto_initialize;
			state = simplex_auto_solve_initialized(presolve.reduced,
				options, reduced_x, result, 1, &initializer);
		} else if (options->algorithm == lp_simplex_ALGORITHM_PAN_BDA)
			state = simplex_pan_solve_problem(
				presolve.reduced, options, reduced_x, result);
		else
			state = simplex_dual_solve_problem(presolve.reduced,
				options, reduced_x, NULL, result, 1);
		if (result->dual_residual_space == lp_simplex_RESIDUAL_ORIGINAL)
			result->dual_residual_space =
				lp_simplex_RESIDUAL_PRESOLVED;
		simplex_presolve_postsolve(&presolve, reduced_x, x);
	}
	/* A reduced model is only an optimization, never a weaker correctness
	 * contract.  Outward-rounded substitutions can still magnify cancellation
	 * during postsolve on badly scaled models.  Retry the original sparse model
	 * when postsolve violates the requested primal tolerance, or when a crash on
	 * the reduced model failed before making a simplex iteration.  The latter
	 * is a cheap and useful discriminator: continuing failures are not hidden,
	 * while models such as SHARE1B recover without disabling presolve globally. */
	{
		double postsolve_infeasibility =
			result->status == lp_simplex_Success
				? simplex_primal_infeasibility(model, x)
				: 0.;
		int invalid_postsolve =
			result->status == lp_simplex_Success &&
			postsolve_infeasibility > options->primal_tolerance;

		if (invalid_postsolve && getenv("LP_SIMPLEX_PROFILE") != NULL)
			fprintf(stderr,
				"presolve: rejected postsolve primal=%.6g\n",
				postsolve_infeasibility);
		if (invalid_postsolve ||
			(result->status == lp_simplex_PrecisionError &&
				(result->iterations == 0 ||
					result->primal_infeasibility <=
						options->primal_tolerance))) {
			double *saved_x = (double *)lp_simplex_malloc(
				(size_t)model->n * sizeof(double));
			struct lp_simplex_Result saved_result = *result;
			struct lp_simplex_Result retry_result;
			int retry_state;
			struct lp_simplex_Options retry_options = *options;

			retry_options.iteration_limit -= result->iterations;
			if (saved_x != NULL &&
				retry_options.iteration_limit > 0) {
				lp_simplex_memcpy(saved_x, x,
					(size_t)model->n * sizeof(double));
				lp_simplex_memset(
					&retry_result, 0, sizeof(retry_result));
				retry_state = simplex_solve_sparse_raw(model,
					&retry_options, x, &retry_result,
					progress);
				retry_result.iterations +=
					saved_result.iterations;
				if (retry_result.status == lp_simplex_Success) {
					*result = retry_result;
					state = retry_state;
				} else {
					lp_simplex_memcpy(x, saved_x,
						(size_t)model->n *
							sizeof(double));
					*result = saved_result;
					result->iterations =
						retry_result.iterations;
				}
			}
			lp_simplex_free(saved_x);
		}
	}

	simplex_progress_bridge_destroy(&bridge);
	lp_simplex_free(reduced_x);
	simplex_presolve_destroy(&presolve);
	return state;
}


static int simplex_solve_impl(const struct lp_Model *model,
	const struct lp_simplex_Options *options, double *x,
	struct lp_simplex_Result *result,
	const struct simplex_ProgressSink *progress)
{
	int state, j;
	int status = lp_simplex_CondUnsatisfied;
	double objective = 0.;

	if (result == NULL)
		return lp_simplex_EXIT_FAILURE;
	result->status = lp_simplex_CondUnsatisfied;
	result->iterations = 0;
	result->objective = NAN;
	result->primal_infeasibility = HUGE_VAL;
	result->dual_infeasibility = NAN;
	result->dual_residual_space = lp_simplex_RESIDUAL_UNAVAILABLE;
	if (x == NULL || !simplex_validate_model(model) ||
		!simplex_validate_options(options)) {
		/*
		 * Integer variables are parsed (MPS INTORG/MARKER, BV/LI/UI) but not
		 * supported.  Say so, instead of failing with a bare status code that
		 * looks like a malformed model.
		 */
		if (x != NULL && model != NULL && model->bounds != NULL) {
			int j;

			for (j = 0; j < model->n; j++)
				if (model->bounds[j].v_type !=
					optm_VAR_T_REAL) {
					fprintf(stderr,
						"lp_simplex: model declares "
						"integer or binary variables; "
						"only continuous LPs are supported\n");
					break;
				}
		}
		return lp_simplex_EXIT_FAILURE;
	}

	for (j = 0; j < model->n; j++)
		x[j] = NAN;
	if (options->algorithm == lp_simplex_ALGORITHM_DUAL_REVISED ||
		options->algorithm == lp_simplex_ALGORITHM_PAN_BDA ||
		options->algorithm == lp_simplex_ALGORITHM_AUTO) {
		state = simplex_solve_sparse_presolved(
			model, options, x, result, progress);
		return simplex_finalize_result(
			model, options, x, result, state);
	}
	if (model->coefficients == NULL) {
		struct lp_Model *dense = simplex_materialize_dense(model);

		if (dense == NULL) {
			result->status = lp_simplex_MemoryAllocError;
			return lp_simplex_EXIT_FAILURE;
		}
		state = simplex_solve_impl(dense, options, x, result, progress);
		lp_model_free(dense);
		return state;
	}

	state = simplex_tableau_solve_model(
		model, options, x, &objective, &status, &result->iterations);
	result->status = status;
	result->objective = objective + model->objective_offset;
	return simplex_finalize_result(model, options, x, result, state);
}


int lp_simplex_solve(const struct lp_Model *model,
	const struct lp_simplex_Options *options, double *x,
	struct lp_simplex_Result *result)
{
	struct simplex_ProgressTracker tracker;
	int state, j;
	if (!options || !options->progress_callback)
		return simplex_solve_impl(model, options, x, result, NULL);
	/* Preserve the existing invalid-input contract without observing it. */
	if (!result || !x || !simplex_validate_model(model) ||
		!simplex_validate_options(options))
		return simplex_solve_impl(model, options, x, result, NULL);
	if (simplex_progress_init(&tracker, model, options)) {
		for (j = 0; j < model->n; j++)
			x[j] = NAN;
		lp_simplex_memset(result, 0, sizeof(*result));
		result->status = lp_simplex_MemoryAllocError;
		result->objective = result->dual_infeasibility = NAN;
		result->primal_infeasibility = HUGE_VAL;
		return lp_simplex_EXIT_FAILURE;
	}
	state = simplex_solve_impl(model, options, x, result, &tracker.sink);
	tracker.sink.report(tracker.sink.context, x, NULL, 1., "finish",
		result->iterations, 1);
	simplex_progress_destroy(&tracker);
	return state;
}
