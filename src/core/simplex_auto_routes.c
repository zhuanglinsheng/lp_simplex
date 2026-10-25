/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_auto_routes.c
 * @brief Route capability/resource admission and kernel adaptation. No AUTO
 * state.
 */

#include "simplex_stage.h"
#include "simplex_execution.h"
#include "simplex_auto_routes.h"
#include "simplex_cover.h"
#include "policy/resources.h"
#include "simplex_pan_bounded.h"
#include "simplex_pan.h"
#include "simplex_pan_config.h"
#include "simplex_pan_service.h"
#include "simplex_pan_basis.h"
#include "simplex_presolve.h"
#include "policy/controls.h"
#include "simplex_singleton_dual.h"
#include "simplex_tableau_solver.h"
#include "utils.h"
#include <math.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>

static int auto_original_binding_valid(
	const struct simplex_AutoOriginalInitializer *c)
{
	const struct lp_Model *m = c->model, *b = &c->binding;
	return m && m->m == b->m && m->n == b->n && m->nnz == b->nnz &&
	       m->objective_offset == b->objective_offset &&
	       m->objective_sense == b->objective_sense &&
	       m->objective == b->objective && m->bounds == b->bounds &&
	       m->constraints == b->constraints &&
	       m->coefficients == b->coefficients &&
	       m->column_start == b->column_start &&
	       m->row_index == b->row_index && m->value == b->value &&
	       m->row_start == b->row_start &&
	       m->column_index == b->column_index &&
	       m->row_value == b->row_value;
}

/* Only the read-only shell may differ. Every input storage/revision and the
 * live presolve binding must agree before dimensions or mapping are read. */
static int auto_initializer_binding_valid(
	const struct simplex_AutoOriginalInitializer *c,
	const struct simplex_Problem *input)
{
	struct simplex_ProblemIdentity current;
	if (!c || !input || !c->presolve || !c->controls || !c->original ||
		!auto_original_binding_valid(c) ||
		c->presolve->original != c->model || !c->scope.context ||
		c->scope.model <= 0 ||
		c->scope.context != input->execution.context ||
		c->scope.model == input->execution.model ||
		!simplex_problem_identity_unchanged(
			&c->target, c->presolve->reduced))
		return 0;
	current = simplex_problem_identity(input);
	current.input = c->target.input;
	return simplex_problem_identity_equal(&current, &c->target) &&
	       input->matrix.rows == input->rows &&
	       input->matrix.columns == input->columns &&
	       input->columns <= c->model->n && c->presolve->column_map;
}

int simplex_auto_routes_initialize(
	const struct simplex_AutoInitializer *initializer,
	const struct simplex_Problem *input,
	const struct lp_simplex_Options *options, double *point, int *work)
{
	*work = 0;
	return initializer && initializer->run
		       ? initializer->run(initializer->context, input, options,
				 point, work)
		       : 0;
}

/* NNLS provides a candidate, not an LP terminal classification. Capability
 * admission and model adaptation live here; the exact lease comes from AUTO. */
int simplex_auto_routes_initialize_original(void *context,
	const struct simplex_Problem *input, const struct lp_simplex_Options *o,
	double *point, int *work)
{
	const struct simplex_AutoOriginalInitializer *c = context;
	struct simplex_Problem raw;
	struct simplex_PanStandard a;
	struct simplex_PanBasis basis;
	struct simplex_PanConfig config;
	struct simplex_PanServiceContext services;
	struct simplex_PanServiceOutcome outcome;
	size_t expanded_rows;
	double *z;
	int j, status, valid = 0;
	*work = 0;
	if (!auto_initializer_binding_valid(c, input) || point == c->original) {
		if (input && input->report) {
			input->report->stop = SIMPLEX_STOP_NUMERICAL;
			input->report->failure = SIMPLEX_FAILURE_STAGE_CONTRACT;
		}
		return 0;
	}
	if (input->report)
		input->report->stop = SIMPLEX_STOP_NOT_APPLICABLE;
	for (j = 0; j < input->columns; j++)
		if (c->presolve->column_map[j] < 0 ||
			c->presolve->column_map[j] >= c->model->n) {
			if (input->report) {
				input->report->stop = SIMPLEX_STOP_NUMERICAL;
				input->report->failure =
					SIMPLEX_FAILURE_STAGE_CONTRACT;
			}
			return 0;
		}
	expanded_rows = (size_t)input->rows;
	for (j = 0; j < input->columns; j++)
		if (input->bounds[j].b_type == optm_BOUND_T_BS &&
			input->bounds[j].ub > input->bounds[j].lb)
			expanded_rows++;
	if (expanded_rows <= 2 * (size_t)c->model->m ||
		o->iteration_limit < SIMPLEX_AUTO_PROBE_CAP ||
		!simplex_pan_original_form_smaller(c->model, input))
		return 0;
	if (simplex_problem_from_model(&raw, c->model)) {
		if (input->report)
			input->report->stop = SIMPLEX_STOP_RESOURCE;
		return 0;
	}
	raw.execution = c->scope;
	raw.execution.history = input->execution.history;
	simplex_execution_inherit_parent(&raw.execution, &input->execution);
	raw.report = input->report;
	if (!simplex_auto_route_allowed(&raw, SIMPLEX_AUTO_CLASSIC)) {
		simplex_problem_destroy(&raw);
		return 0;
	}
	if (simplex_pan_standard_create(&a, &raw)) {
		if (input->report)
			input->report->stop = SIMPLEX_STOP_RESOURCE;
		simplex_problem_destroy(&raw);
		return 0;
	}
	simplex_pan_service_context_init(&services, &raw);
	simplex_pan_service_bind(&a, &services);
	config = *c->controls;
	config.model_form = SIMPLEX_PAN_STANDARD;
	config.direction_method = SIMPLEX_PAN_BASIS_DIRECTION;
	config.oracle = SIMPLEX_PAN_NO_ORACLE;
	if (simplex_pan_basis_create_configured(
		    &basis, &a, o->pivot_tolerance, &config)) {
		if (input->report)
			input->report->stop = SIMPLEX_STOP_RESOURCE;
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
		(void)simplex_pan_service_nnls(&basis, z, o->primal_tolerance,
			o->dual_tolerance, o->iteration_limit, work, &status,
			&outcome);
		simplex_pan_standard_recover(&a, z, c->original);
		if (outcome.valid && auto_initializer_binding_valid(c, input) &&
			simplex_problem_primal_infeasibility(
				&raw, c->original) <= o->primal_tolerance) {
			for (j = 0; j < input->columns; j++)
				point[j] =
					c->original[c->presolve->column_map[j]];
			valid = simplex_problem_primal_infeasibility(
					input, point) <= o->primal_tolerance;
		}
		if (input->report)
			input->report->stop =
				valid ? SIMPLEX_STOP_HANDOFF : outcome.stop;
		lp_simplex_free(z);
	} else if (input->report)
		input->report->stop = SIMPLEX_STOP_RESOURCE;
	if (getenv(SIMPLEX_ENV_PROFILE))
		fprintf(stderr,
			"auto: original-feasibility expanded-rows=%lu original-rows=%d work=%d accepted=%d\n",
			(unsigned long)expanded_rows, raw.rows, *work, valid);
	simplex_pan_basis_destroy(&basis);
	simplex_pan_standard_destroy(&a);
	simplex_problem_destroy(&raw);
	return valid;
}

/* A workspace admission limit, not a claim about total RSS or sparse fill.
 * All arithmetic is promoted before multiplication and before integer casts. */
int simplex_auto_route_admission(const struct simplex_Problem *p, int route)
{
	long double rows, cols, cells, nonzeros;
	int i, j, positive, negative;

	if (route < 0 || route >= SIMPLEX_AUTO_ROUTE_COUNT)
		return SIMPLEX_AUTO_MODEL_ROUTE;
	if (!p || p->rows <= 0 || p->columns <= 0)
		return SIMPLEX_AUTO_MODEL_SHAPE;
	rows = p->rows;
	cols = p->columns;
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
		return SIMPLEX_AUTO_MODEL_INDEX;
	if (route == SIMPLEX_AUTO_CLASSIC ||
		route == SIMPLEX_AUTO_STANDARD_IMPLICIT ||
		route == SIMPLEX_AUTO_STANDARD_EXPLICIT) {
		/* Free columns duplicate nonzeros; widths add a row, a coefficient
		 * and a slack. Budget materialized storage before performing any
		 * conversion. Sparse QR fill is separate and not predictable here. */
		nonzeros =
			2.L * p->matrix.nonzeros + 2.L * p->columns + p->rows;
		if (nonzeros >= INT_MAX || cols + rows >= INT_MAX / 2)
			return SIMPLEX_AUTO_MODEL_INDEX;
		cells = 8 * (rows + cols + 1) + 4 * nonzeros;
		if (cells * sizeof(double) > SIMPLEX_AUTO_CONVERSION_BYTES)
			return SIMPLEX_AUTO_MODEL_CONVERSION;
	}
	if (route == SIMPLEX_AUTO_TABLEAU || route == SIMPLEX_AUTO_SINGLETON) {
		/* Original dense model + transformed matrix + worst Phase-I table,
		 * plus conservative vector/index storage. Upper bounds add rows,
		 * free variables add columns; neither may be ignored. */
		cells = (long double)p->rows * p->columns + rows * cols +
			(rows + 1) * (cols + 2 * rows + 2) +
			16 * (rows + cols + 1);
		if (cells * sizeof(double) > SIMPLEX_AUTO_DENSE_BYTES)
			return SIMPLEX_AUTO_MODEL_DENSE;
	}
#ifndef LP_SIMPLEX_HAVE_SPQR
	if (route == SIMPLEX_AUTO_HYBRID || route == SIMPLEX_AUTO_EXPLICIT ||
		route == SIMPLEX_AUTO_CLASSIC ||
		route == SIMPLEX_AUTO_STANDARD_IMPLICIT ||
		route == SIMPLEX_AUTO_STANDARD_EXPLICIT) {
		/* Dense QR stores Q, R, factors and reconstruction scratch. */
		if (8 * (rows + 1) * (rows + 1) * sizeof(double) >
			SIMPLEX_AUTO_DENSE_BYTES)
			return SIMPLEX_AUTO_MODEL_DENSE_QR;
	}
#endif
	if (route != SIMPLEX_AUTO_SINGLETON)
		return SIMPLEX_AUTO_MODEL_READY;
	/* The transformation is admissible only with a paired positive-cost
	 * singleton for each equality and nonnegative, lower-only variables. */
	for (j = 0; j < p->columns; j++)
		if (p->bounds[j].b_type != optm_BOUND_T_LO ||
			p->bounds[j].lb != 0.)
			return SIMPLEX_AUTO_MODEL_SINGLETON_BOUNDS;
	for (i = 0; i < p->rows; i++) {
		if (p->row_type[i] != optm_CONS_T_EQ)
			return SIMPLEX_AUTO_MODEL_SINGLETON_ROWS;
		positive = negative = 0;
		/* CSR already indexes the columns incident to this row. Scanning all
		 * columns for every row made this admission check O(rows * columns)
		 * even when the model has only O(rows + columns) nonzeros. */
		for (j = p->matrix.row_start[i]; j < p->matrix.row_start[i + 1];
			j++) {
			int column = p->matrix.column_index[j];
			int k = p->matrix.column_start[column];

			if (p->matrix.column_start[column + 1] - k != 1 ||
				p->matrix.row_index[k] != i ||
				p->objective[column] <= 0.)
				continue;
			positive += p->matrix.value[k] > 0.;
			negative += p->matrix.value[k] < 0.;
		}
		if (positive != 1 || negative != 1)
			return SIMPLEX_AUTO_MODEL_SINGLETON_PAIRS;
	}
	return SIMPLEX_AUTO_MODEL_READY;
}

int simplex_auto_route_allowed(const struct simplex_Problem *p, int route)
{
	return simplex_auto_route_admission(p, route) ==
	       SIMPLEX_AUTO_MODEL_READY;
}

const char *simplex_auto_route_admission_name(int admission)
{
	static const char *const names[] = {"unchecked", "ready",
		"unknown-route", "model-shape", "index-capacity",
		"standard-form-storage", "dense-model-storage",
		"dense-qr-storage", "singleton-bounds", "singleton-rows",
		"singleton-pairs"};
	return admission >= 0 && (size_t)admission <
					 sizeof(names) / sizeof(names[0])
		       ? names[admission]
		       : "unknown";
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

/* Invocation output is never a continuation of a previous route's result. */
static void auto_route_reset(
	const struct simplex_Problem *p, struct lp_simplex_Result *r)
{
	lp_simplex_memset(r, 0, sizeof(*r));
	r->status = lp_simplex_CondUnsatisfied;
	r->objective = r->primal_infeasibility = r->dual_infeasibility = NAN;
	r->dual_residual_space = lp_simplex_RESIDUAL_UNAVAILABLE;
	if (p && p->report)
		simplex_stage_reset(p->report);
}

int simplex_auto_route_solve(int route, const struct simplex_Problem *p,
	const struct lp_simplex_Options *o, const double *seed, double *x,
	struct lp_simplex_Result *r)
{
	struct lp_simplex_Options stage = *o;
	struct simplex_PanConfig config;
	struct lp_Model *m;
	int state, applicable = 0;

	auto_route_reset(p, r);
	if (!simplex_auto_route_allowed(p, route)) {
		r->status = lp_simplex_CondUnsatisfied;
		r->iterations = 0;
		if (p && p->report)
			p->report->stop = SIMPLEX_STOP_NOT_APPLICABLE;
		return lp_simplex_EXIT_FAILURE;
	}
	if (route == SIMPLEX_AUTO_HYBRID || route == SIMPLEX_AUTO_EXPLICIT ||
		route == SIMPLEX_AUTO_CLASSIC ||
		route == SIMPLEX_AUTO_STANDARD_IMPLICIT ||
		route == SIMPLEX_AUTO_STANDARD_EXPLICIT) {
		int native = route == SIMPLEX_AUTO_HYBRID ||
			     route == SIMPLEX_AUTO_EXPLICIT;
		int oracle =
			route == SIMPLEX_AUTO_CLASSIC
				? SIMPLEX_PAN_NO_ORACLE
				: (route == SIMPLEX_AUTO_HYBRID ||
							  route ==
								  SIMPLEX_AUTO_STANDARD_IMPLICIT
						  ? SIMPLEX_PAN_IMPLICIT_ORACLE
						  : SIMPLEX_PAN_EXPLICIT_ORACLE);

		(void)simplex_pan_config_default(&config,
			native ? SIMPLEX_PAN_NATIVE_BOUNDS
			       : SIMPLEX_PAN_STANDARD,
			route == SIMPLEX_AUTO_CLASSIC
				? SIMPLEX_PAN_BASIS_DIRECTION
				: SIMPLEX_PAN_CONE_DIRECTION,
			oracle);
		stage.algorithm = lp_simplex_ALGORITHM_PAN_BDA;
		if (!native)
			stage.pricing = lp_simplex_PRICING_PAN_NORMALIZED;
		return simplex_pan_solve_configured(
			p, &stage, seed, x, r, &config);
	}
	m = auto_dense_model(p);
	if (!m) {
		r->status = lp_simplex_MemoryAllocError;
		if (p->report)
			p->report->stop = SIMPLEX_STOP_RESOURCE;
		return lp_simplex_EXIT_FAILURE;
	}
	if (route == SIMPLEX_AUTO_SINGLETON) {
		stage.algorithm = lp_simplex_ALGORITHM_DUAL_REVISED;
		state = simplex_singleton_dual_solve_reported(
			m, &stage, x, r, &applicable, p->report, &p->execution);
		if (!applicable && r->status != lp_simplex_MemoryAllocError) {
			r->status = lp_simplex_CondUnsatisfied;
			state = lp_simplex_EXIT_FAILURE;
		}
	} else {
		stage.algorithm = lp_simplex_ALGORITHM_TABLEAU;
		stage.pricing = lp_simplex_PRICING_BLAND;
		state = simplex_tableau_solve_model_reported(m, &stage, x,
			&r->objective, &r->status, &r->iterations, p->report,
			&p->execution);
		/* Tableau owns its Phase-II optimality check but does not export
		 * multipliers. Never invent a zero dual diagnostic. */
		r->dual_infeasibility = NAN;
		r->dual_residual_space = lp_simplex_RESIDUAL_UNAVAILABLE;
	}
	lp_model_free(m);
	return state;
}

/* Production binding resolves BDA choices here; no AUTO-mode flag reaches
 * the configured BDA kernel, its PDHG/oracle services or QR/NNLS layers. */
static int auto_native_bda(const struct simplex_Problem *p,
	const struct lp_simplex_Options *o, const double *seed, double *x,
	struct lp_simplex_Result *r)
{
	return simplex_auto_route_solve(SIMPLEX_AUTO_HYBRID, p, o, seed, x, r);
}

void simplex_auto_routes_defaults(
	struct simplex_AutoKernels *kernels, struct simplex_AutoRoutes *routes)
{
	kernels->dual = simplex_dual_solve_problem;
	kernels->hybrid = auto_native_bda;
	routes->solve = simplex_auto_route_solve;
}

int simplex_auto_routes_available(int route,
	const struct simplex_AutoKernels *kernels,
	const struct simplex_AutoRoutes *routes)
{
	if (route == SIMPLEX_AUTO_COVER)
		return 1;
	if (route == SIMPLEX_AUTO_DUAL_EXPLORE ||
		route == SIMPLEX_AUTO_DUAL_FALLBACK)
		return kernels && kernels->dual;
	if (route == SIMPLEX_AUTO_HYBRID)
		return kernels && kernels->hybrid;
	return route >= 0 && route < SIMPLEX_AUTO_ROUTE_COUNT && routes &&
	       routes->solve;
}

int simplex_auto_routes_has_checkpoint(
	const struct simplex_AutoKernels *kernels)
{
	return kernels && kernels->dual == simplex_dual_solve_problem;
}

int simplex_auto_routes_invoke(int route, const struct simplex_Problem *p,
	const struct lp_simplex_Options *o, const double *seed, double *x,
	struct lp_simplex_Result *r, int propagate_bounds,
	const struct simplex_AutoKernels *kernels,
	const struct simplex_AutoRoutes *routes,
	struct simplex_DualCheckpoint *checkpoint)
{
	if (!simplex_auto_routes_available(route, kernels, routes)) {
		r->status = lp_simplex_CondUnsatisfied;
		r->iterations = 0;
		r->objective = r->primal_infeasibility = r->dual_infeasibility =
			NAN;
		r->dual_residual_space = lp_simplex_RESIDUAL_UNAVAILABLE;
		if (p->report) {
			simplex_stage_reset(p->report);
			p->report->stop = SIMPLEX_STOP_NOT_APPLICABLE;
		}
		return lp_simplex_EXIT_FAILURE;
	}
	if (route == SIMPLEX_AUTO_COVER)
		return simplex_cover_solve(p, o, x, NULL, r);
	if (route == SIMPLEX_AUTO_DUAL_EXPLORE ||
		route == SIMPLEX_AUTO_DUAL_FALLBACK) {
		if (checkpoint && checkpoint->basis.index &&
			simplex_auto_routes_has_checkpoint(kernels))
			return simplex_dual_solve_problem_checkpoint(
				p, o, x, NULL, r, propagate_bounds, checkpoint);
		return kernels->dual(p, o, x, NULL, r, propagate_bounds);
	}
	if (route == SIMPLEX_AUTO_HYBRID)
		return kernels->hybrid(p, o, seed, x, r);
	if (route >= 0 && route < SIMPLEX_AUTO_ROUTE_COUNT && routes &&
		routes->solve)
		return routes->solve(route, p, o, seed, x, r);
	r->status = lp_simplex_CondUnsatisfied;
	r->iterations = 0;
	return lp_simplex_EXIT_FAILURE;
}
