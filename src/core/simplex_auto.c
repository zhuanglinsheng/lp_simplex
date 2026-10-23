/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_auto.c
 * @brief AUTO owns stage selection and the cumulative work budget.
 *
 * Kernels own their numerical safeguards. Iterations count pivots/first-order
 * steps/NNLS work, not seconds; factorization cost is deliberately not called
 * an iteration.
 */

#include "simplex_auto.h"
#include "simplex_dual.h"
#include "simplex_pan_bounded.h"
#include "simplex_pan.h"
#include "simplex_pan_config.h"
#include "simplex_singleton_dual.h"
#include "simplex_tableau_solver.h"
#include "utils.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>


static void auto_reset(struct lp_simplex_Result *r)
{
	lp_simplex_memset(r, 0, sizeof(*r));
	r->status = lp_simplex_CondUnsatisfied;
	r->objective = r->primal_infeasibility = r->dual_infeasibility =
		HUGE_VAL;
}


/* Working bounds and scaled residuals never authorize passing a seed. */
static int auto_feasible(const struct simplex_Problem *p, const double *x,
	const struct lp_simplex_Options *o, double *objective)
{
	long double cost = 0.;
	int j;

	if (simplex_problem_primal_infeasibility(p, x) > o->primal_tolerance)
		return 0;
	for (j = 0; j < p->columns; j++)
		cost += (long double)p->objective[j] * x[j];
	*objective = (double)cost;
	return isfinite(*objective);
}


static int auto_terminal(const struct lp_simplex_Result *r)
{
	return r->status == lp_simplex_Infeasibility ||
	       r->status == lp_simplex_Unboundedness ||
	       r->status == lp_simplex_MemoryAllocError;
}


static const char *auto_reason(int status)
{
	switch (status) {
	case lp_simplex_Success:
		return "optimal";
	case lp_simplex_Infeasibility:
		return "infeasible";
	case lp_simplex_Unboundedness:
		return "unbounded";
	case lp_simplex_MemoryAllocError:
		return "memory-failure";
	case lp_simplex_ExceedIterLimit:
		return "work-budget";
	case lp_simplex_Degeneracy:
		return "stagnation";
	case lp_simplex_Singularity:
		return "singularity";
	default:
		return "numerical-or-audit-failure";
	}
}


static void auto_trace(
	const char *stage, const char *reason, int budget, int used, int seed)
{
	if (getenv("LP_SIMPLEX_PROFILE"))
		fprintf(stderr,
			"auto: stage=%s reason=%s budget=%d used=%d feasible-seed=%d\n",
			stage, reason, budget, used, seed);
}


/* A workspace admission limit, not a claim about total RSS or sparse fill.
 * All arithmetic is promoted before multiplication and before integer casts. */
int simplex_auto_route_allowed(const struct simplex_Problem *p, int route)
{
	long double rows = p->rows, cols = p->columns, cells, nonzeros;
	int i, j, positive, negative;

	if (route < 0 || route >= SIMPLEX_AUTO_ROUTE_COUNT || p->rows <= 0 ||
		p->columns <= 0)
		return 0;
	if (route == SIMPLEX_AUTO_CLASSIC ||
		route == SIMPLEX_AUTO_STANDARD_IMPLICIT ||
		route == SIMPLEX_AUTO_STANDARD_EXPLICIT ||
		route == SIMPLEX_AUTO_TABLEAU ||
		route == SIMPLEX_AUTO_SINGLETON) {
		for (j = 0; j < p->columns; j++) {
			if (p->bounds[j].b_type == optm_BOUND_T_FR)
				cols++;
			if (p->bounds[j].b_type == optm_BOUND_T_BS)
				rows++;
		}
	}
	if (rows >= INT_MAX / 4 || cols >= INT_MAX / 4)
		return 0;
	if (route == SIMPLEX_AUTO_CLASSIC ||
		route == SIMPLEX_AUTO_STANDARD_IMPLICIT ||
		route == SIMPLEX_AUTO_STANDARD_EXPLICIT) {
		/* Free columns duplicate nonzeros; widths add a row, a coefficient
		 * and a slack. Budget materialized storage before performing any
		 * conversion. Sparse QR fill is separate and not predictable here. */
		nonzeros =
			2.L * p->matrix.nonzeros + 2.L * p->columns + p->rows;
		if (nonzeros >= INT_MAX || cols + rows >= INT_MAX / 2)
			return 0;
		cells = 8 * (rows + cols + 1) + 4 * nonzeros;
		if (cells * sizeof(double) > SIMPLEX_AUTO_CONVERSION_BYTES)
			return 0;
	}
	if (route == SIMPLEX_AUTO_TABLEAU || route == SIMPLEX_AUTO_SINGLETON) {
		/* Original dense model + transformed matrix + worst Phase-I table,
		 * plus conservative vector/index storage. Upper bounds add rows,
		 * free variables add columns; neither may be ignored. */
		cells = (long double)p->rows * p->columns + rows * cols +
			(rows + 1) * (cols + 2 * rows + 2) +
			16 * (rows + cols + 1);
		if (cells * sizeof(double) > SIMPLEX_AUTO_DENSE_BYTES)
			return 0;
	}
#ifndef LP_SIMPLEX_HAVE_SPQR
	if (route == SIMPLEX_AUTO_HYBRID || route == SIMPLEX_AUTO_EXPLICIT ||
		route == SIMPLEX_AUTO_CLASSIC ||
		route == SIMPLEX_AUTO_STANDARD_IMPLICIT ||
		route == SIMPLEX_AUTO_STANDARD_EXPLICIT) {
		/* Dense QR stores Q, R, factors and reconstruction scratch. */
		if (8 * (rows + 1) * (rows + 1) * sizeof(double) >
			SIMPLEX_AUTO_DENSE_BYTES)
			return 0;
	}
#endif
	if (route != SIMPLEX_AUTO_SINGLETON)
		return 1;
	/* The transformation is admissible only with a paired positive-cost
	 * singleton for each equality and nonnegative, lower-only variables. */
	for (j = 0; j < p->columns; j++)
		if (p->bounds[j].b_type != optm_BOUND_T_LO ||
			p->bounds[j].lb != 0.)
			return 0;
	for (i = 0; i < p->rows; i++) {
		if (p->row_type[i] != optm_CONS_T_EQ)
			return 0;
		positive = negative = 0;
		for (j = 0; j < p->columns; j++) {
			int k = p->matrix.column_start[j];

			if (p->matrix.column_start[j + 1] - k != 1 ||
				p->matrix.row_index[k] != i ||
				p->objective[j] <= 0.)
				continue;
			positive += p->matrix.value[k] > 0.;
			negative += p->matrix.value[k] < 0.;
		}
		if (positive != 1 || negative != 1)
			return 0;
	}
	return 1;
}


static const char *auto_route_name(int route)
{
	static const char *names[] = {"native-hybrid", "native-explicit",
		"classic-pan", "standard-implicit", "standard-explicit",
		"singleton-dual", "tableau"};
	return names[route];
}


static struct lp_Model *auto_dense_model(const struct simplex_Problem *p)
{
	struct lp_Model *m = lp_model_create(p->rows, p->columns);
	int i, j, k;

	if (!m)
		return NULL;
	lp_simplex_memcpy(
		m->bounds, p->bounds, (size_t)p->columns * sizeof(*m->bounds));
	lp_simplex_memcpy(m->objective, p->objective,
		(size_t)p->columns * sizeof(double));
	for (i = 0; i < p->rows; i++) {
		m->constraints[i].type = p->row_type[i];
		m->constraints[i].rhs = p->rhs[i];
	}
	for (j = 0; j < p->columns; j++)
		for (k = p->matrix.column_start[j];
			k < p->matrix.column_start[j + 1]; k++)
			m->constraints[p->matrix.row_index[k]].coef[j] +=
				p->matrix.value[k];
	return m;
}


int simplex_auto_route_solve(int route, const struct simplex_Problem *p,
	const struct lp_simplex_Options *o, const double *seed, double *x,
	struct lp_simplex_Result *r)
{
	struct lp_simplex_Options stage = *o;
	struct lp_Model *m;
	int state, applicable = 0;

	if (!simplex_auto_route_allowed(p, route)) {
		r->status = lp_simplex_CondUnsatisfied;
		r->iterations = 0;
		return lp_simplex_EXIT_FAILURE;
	}
	stage.algorithm = lp_simplex_ALGORITHM_AUTO;
	if (route == SIMPLEX_AUTO_HYBRID)
		return simplex_pan_bounded_solve(p, &stage, seed, x, r);
	if (route == SIMPLEX_AUTO_EXPLICIT) {
		stage.algorithm = SIMPLEX_PAN_AUTO_EXPLICIT;
		return simplex_pan_bounded_solve(p, &stage, seed, x, r);
	}
	if (route == SIMPLEX_AUTO_CLASSIC ||
		route == SIMPLEX_AUTO_STANDARD_IMPLICIT ||
		route == SIMPLEX_AUTO_STANDARD_EXPLICIT) {
		/* AUTO in this entry explicitly selects classic standard-form Pan
		 * and isolates it from manual environment-controlled experiments. */
		if (route == SIMPLEX_AUTO_STANDARD_IMPLICIT)
			stage.algorithm = SIMPLEX_PAN_AUTO_STANDARD_IMPLICIT;
		if (route == SIMPLEX_AUTO_STANDARD_EXPLICIT)
			stage.algorithm = SIMPLEX_PAN_AUTO_STANDARD_EXPLICIT;
		stage.pricing = lp_simplex_PRICING_PAN_NORMALIZED;
		return simplex_pan_solve_from_feasible(p, &stage, seed, x, r);
	}
	m = auto_dense_model(p);
	if (!m) {
		r->status = lp_simplex_MemoryAllocError;
		return lp_simplex_EXIT_FAILURE;
	}
	if (route == SIMPLEX_AUTO_SINGLETON) {
		stage.algorithm = lp_simplex_ALGORITHM_DUAL_REVISED;
		state = simplex_singleton_dual_solve(
			m, &stage, x, r, &applicable);
		if (!applicable && r->status != lp_simplex_MemoryAllocError) {
			r->status = lp_simplex_CondUnsatisfied;
			state = lp_simplex_EXIT_FAILURE;
		}
	} else {
		stage.algorithm = lp_simplex_ALGORITHM_TABLEAU;
		stage.pricing = lp_simplex_PRICING_BLAND;
		state = simplex_tableau_solve_model(m, &stage, x, &r->objective,
			&r->status, &r->iterations);
		/* Tableau owns its Phase-II optimality check but does not export
		 * multipliers. Never invent a zero dual diagnostic. */
		r->dual_infeasibility = NAN;
		r->dual_residual_space = lp_simplex_RESIDUAL_UNAVAILABLE;
	}
	lp_model_free(m);
	return state;
}


static int auto_run_routes(const struct simplex_Problem *p,
	const struct lp_simplex_Options *o, double *x,
	struct lp_simplex_Result *result, int propagate_bounds,
	const struct simplex_AutoKernels *kernels,
	const struct simplex_AutoRoutes *routes,
	const struct simplex_AutoInitializer *initializer)
{
	static const int order[SIMPLEX_AUTO_ROUTE_COUNT] = {
		SIMPLEX_AUTO_SINGLETON, SIMPLEX_AUTO_HYBRID,
		SIMPLEX_AUTO_EXPLICIT, SIMPLEX_AUTO_CLASSIC,
		SIMPLEX_AUTO_STANDARD_IMPLICIT, SIMPLEX_AUTO_STANDARD_EXPLICIT,
		SIMPLEX_AUTO_TABLEAU};
	struct lp_simplex_Options stage = *o;
	struct lp_simplex_Result current, best_result;
	double *best, cost, best_cost = HUGE_VAL;
	int j, state = lp_simplex_EXIT_FAILURE, used = 0, remaining, probe,
	       reserve, have_best = 0, feasible = 0, index, route = -1,
	       allowed[SIMPLEX_AUTO_ROUTE_COUNT], count = 0, slots;
	const char *reason;

	auto_reset(result);
	auto_reset(&best_result);
	for (j = 0; j < SIMPLEX_AUTO_ROUTE_COUNT; j++) {
		allowed[j] = (j == SIMPLEX_AUTO_HYBRID || routes != NULL) &&
			     simplex_auto_route_allowed(p, j);
		count += allowed[j];
	}
	best = (double *)lp_simplex_malloc(
		((size_t)p->columns + 1) * sizeof(double));
	if (!best) {
		result->status = lp_simplex_MemoryAllocError;
		return lp_simplex_EXIT_FAILURE;
	}
	/* Progress renews the live basis. Keep a recovery reserve, rather than
	 * imposing an absolute pivot cutoff on a progressing solve. */
	probe = o->iteration_limit - o->iteration_limit / 8;
	if (!count) {
		probe = o->iteration_limit;
		stage.algorithm = lp_simplex_ALGORITHM_DUAL_REVISED;
	}
	stage.iteration_limit = probe;
	auto_reset(&current);
	index = -1;
	auto_trace("dual-explore", "start", probe, used, 0);
	for (j = 0; j < p->columns; j++)
		x[j] = HUGE_VAL;
	state = kernels->dual(p, &stage, x, NULL, &current, propagate_bounds);
	for (;;) {
		/* A malformed kernel result must not overflow/reset the shared budget
		 * or authorize success. Charge its entire granted lease. */
		if (current.iterations < 0 ||
			current.iterations > stage.iteration_limit) {
			current.iterations = stage.iteration_limit;
			current.status = lp_simplex_PrecisionError;
			state = lp_simplex_EXIT_FAILURE;
		}
		used += current.iterations;
		feasible = auto_feasible(p, x, o, &cost);
		if (state == 0 && current.status == lp_simplex_Success &&
			feasible &&
			(isfinite(current.dual_infeasibility) ||
				(route == SIMPLEX_AUTO_TABLEAU &&
					current.dual_residual_space ==
						lp_simplex_RESIDUAL_UNAVAILABLE &&
					isnan(current.dual_infeasibility))))
			goto done;
		if (current.status == lp_simplex_Success) {
			current.status = lp_simplex_PrecisionError;
			state = lp_simplex_EXIT_FAILURE;
		}
		if (current.status == lp_simplex_Infeasibility &&
			(feasible || have_best))
			current.status = lp_simplex_PrecisionError;
		if (current.status == lp_simplex_MemoryAllocError ||
			auto_terminal(&current))
			goto done;
		if (feasible && (!have_best || cost < best_cost)) {
			lp_simplex_memcpy(
				best, x, (size_t)p->columns * sizeof(double));
			have_best = 1;
			best_cost = cost;
			best_result = current;
		}
		remaining = o->iteration_limit - used;
		if (remaining <= 0 || route == -2 || !count)
			break;
		if (route == -1 && !have_best && initializer) {
			struct lp_simplex_Options initialization = *o;
			int work = 0, accepted;

			initialization.iteration_limit = remaining;
			accepted = initializer->run(initializer->context,
				&initialization, best, &work);
			if (work < 0 || work > remaining) {
				work = remaining;
				accepted = 0;
			}
			used += work;
			remaining -= work;
			if (accepted && auto_feasible(p, best, o, &cost)) {
				have_best = 1;
				best_cost = cost;
				best_result.objective = cost;
				best_result.dual_infeasibility = NAN;
				best_result.dual_residual_space =
					lp_simplex_RESIDUAL_UNAVAILABLE;
				best_result.primal_infeasibility =
					simplex_problem_primal_infeasibility(
						p, best);
				/* A validated point supersedes singleton feasibility work. */
				index = 0;
			}
			auto_trace("feasibility-initializer",
				have_best ? "validated-seed"
					  : "no-validated-seed",
				work, used, have_best);
			if (remaining <= 0)
				break;
		}
		reason = route == -1 &&
					 current.status ==
						 lp_simplex_ExceedIterLimit &&
					 current.iterations <
						 stage.iteration_limit
				 ? "dual-progress-yield"
				 : auto_reason(current.status);
		/* Explicit QR is a different numerical representation, so only
		 * retry it after failure of the implicit representation, not merely
		 * because a work quota expired. Each route is attempted at most once. */
		if (route == SIMPLEX_AUTO_HYBRID &&
			current.status != lp_simplex_PrecisionError &&
			current.status != lp_simplex_Singularity)
			allowed[SIMPLEX_AUTO_EXPLICIT] = 0;
		if (route == SIMPLEX_AUTO_STANDARD_IMPLICIT &&
			current.status != lp_simplex_PrecisionError &&
			current.status != lp_simplex_Singularity)
			allowed[SIMPLEX_AUTO_STANDARD_EXPLICIT] = 0;
		while (++index < SIMPLEX_AUTO_ROUTE_COUNT) {
			route = order[index];
			if (allowed[route])
				break;
			auto_trace(auto_route_name(route),
				"ineligible-or-not-triggered", 0, used,
				have_best);
		}
		stage = *o;
		if (index >= SIMPLEX_AUTO_ROUTE_COUNT) {
			route = -2;
			stage.algorithm = lp_simplex_ALGORITHM_DUAL_REVISED;
			stage.iteration_limit = remaining;
			auto_trace("dual-fallback", reason, remaining, used, 0);
		} else {
			/* Reserve one share for every untried admissible path and one
			 * final cold dual solve. Hybrid gets the main allocation; rescue
			 * paths share the reserve rather than restarting the budget. */
			slots = 1;
			for (j = index + 1; j < SIMPLEX_AUTO_ROUTE_COUNT; j++)
				slots += allowed[order[j]];
			if (route == SIMPLEX_AUTO_HYBRID) {
				reserve = remaining / 8;
				if (reserve < slots * SIMPLEX_AUTO_PROBE_CAP)
					reserve =
						slots * SIMPLEX_AUTO_PROBE_CAP;
				if (reserve > remaining / 2)
					reserve = remaining / 2;
				stage.iteration_limit = remaining - reserve;
			} else
				stage.iteration_limit = remaining / (slots + 1);
			if (stage.iteration_limit < 1)
				stage.iteration_limit = 1;
			auto_trace(auto_route_name(route), reason,
				stage.iteration_limit, used, have_best);
		}
		auto_reset(&current);
		for (j = 0; j < p->columns; j++)
			x[j] = HUGE_VAL;
		if (route == -2)
			state = kernels->dual(
				p, &stage, x, NULL, &current, propagate_bounds);
		else if (route == SIMPLEX_AUTO_HYBRID)
			state = kernels->hybrid(p, &stage,
				have_best ? best : NULL, x, &current);
		else
			state = routes->solve(route, p, &stage,
				have_best ? best : NULL, x, &current);
	}
	if (have_best) {
		lp_simplex_memcpy(x, best, (size_t)p->columns * sizeof(double));
		current = best_result;
	}
	current.status = used >= o->iteration_limit ? lp_simplex_ExceedIterLimit
						    : lp_simplex_PrecisionError;
	state = lp_simplex_EXIT_FAILURE;
done:
	if (current.status != lp_simplex_Success && have_best &&
		(!feasible || current.status == lp_simplex_MemoryAllocError)) {
		lp_simplex_memcpy(x, best, (size_t)p->columns * sizeof(double));
		current.objective = best_cost;
		current.primal_infeasibility = best_result.primal_infeasibility;
		current.dual_infeasibility = best_result.dual_infeasibility;
		current.dual_residual_space = best_result.dual_residual_space;
	}
	*result = current;
	result->iterations = used;
	auto_trace("finish", auto_reason(result->status), 0, used, 0);
	lp_simplex_free(best);
	return state;
}


int simplex_auto_run_routes(const struct simplex_Problem *p,
	const struct lp_simplex_Options *o, double *x,
	struct lp_simplex_Result *result, int propagate_bounds,
	const struct simplex_AutoKernels *kernels,
	const struct simplex_AutoRoutes *routes)
{
	return auto_run_routes(
		p, o, x, result, propagate_bounds, kernels, routes, NULL);
}


int simplex_auto_run(const struct simplex_Problem *p,
	const struct lp_simplex_Options *o, double *x,
	struct lp_simplex_Result *result, int propagate_bounds,
	const struct simplex_AutoKernels *kernels)
{
	return simplex_auto_run_routes(
		p, o, x, result, propagate_bounds, kernels, NULL);
}


int simplex_auto_solve_problem(const struct simplex_Problem *p,
	const struct lp_simplex_Options *o, double *x,
	struct lp_simplex_Result *result, int propagate_bounds)
{
	struct simplex_AutoKernels kernels;
	struct simplex_AutoRoutes routes;

	kernels.dual = simplex_dual_solve_problem;
	kernels.hybrid = simplex_pan_bounded_solve;
	routes.solve = simplex_auto_route_solve;
	return simplex_auto_run_routes(
		p, o, x, result, propagate_bounds, &kernels, &routes);
}


int simplex_auto_solve_initialized(const struct simplex_Problem *p,
	const struct lp_simplex_Options *o, double *x,
	struct lp_simplex_Result *result, int propagate_bounds,
	const struct simplex_AutoInitializer *initializer)
{
	struct simplex_AutoKernels kernels;
	struct simplex_AutoRoutes routes;

	kernels.dual = simplex_dual_solve_problem;
	kernels.hybrid = simplex_pan_bounded_solve;
	routes.solve = simplex_auto_route_solve;
	return auto_run_routes(p, o, x, result, propagate_bounds, &kernels,
		&routes, initializer);
}


int simplex_auto_run_initialized(const struct simplex_Problem *p,
	const struct lp_simplex_Options *o, double *x,
	struct lp_simplex_Result *result, int propagate_bounds,
	const struct simplex_AutoKernels *kernels,
	const struct simplex_AutoInitializer *initializer)
{
	return auto_run_routes(
		p, o, x, result, propagate_bounds, kernels, NULL, initializer);
}
