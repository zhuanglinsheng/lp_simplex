/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
/* Bounded dual revised simplex over an immutable CSC constraint matrix. */
#include "simplex_dual.h"
#include "simplex_dual_internal.h"
#include "simplex_dual_pricing.h"
#include "linalg.h"
#include "utils.h"

#include <lp_simplex/status.h>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>


#define DUAL_RELATIVE_PIVOT_TOLERANCE 1e-8


static int dual_is_finite_lower(double value)
{
	return value > __lp_simplex_NINF__;
}


static int dual_is_finite_upper(double value)
{
	return value < __lp_simplex_INF__;
}


static int dual_pan_structure_enabled(const struct simplex_DualState *state)
{
	return state->rows >= 4096 && state->structural_basic < state->rows;
}


/* Highly regular column degrees are a cheap symmetry signal.  On these
 * models (assignment/QAP-like formulations in particular), lexicographic
 * near-tie handling can walk very long degenerate orbits; maximum-pivot
 * ordering is the better default. */
static int dual_matrix_has_regular_columns(
		const struct simplex_CscMatrix *matrix, int *modal_degree)
{
	int column, best = 0;
	int *frequency;
	frequency = (int *)lp_simplex_malloc(
		(size_t)(matrix->rows + 1) * sizeof(int));
	if (frequency == NULL)
		return 0;
	lp_simplex_memset(frequency, 0,
		(size_t)(matrix->rows + 1) * sizeof(int));
	for (column = 0; column < matrix->columns; column++) {
		int degree = matrix->column_start[column + 1] -
			matrix->column_start[column];
		frequency[degree]++;
		if (frequency[degree] > best) {
			best = frequency[degree];
			*modal_degree = degree;
		}
	}
	lp_simplex_free(frequency);
	return matrix->columns > 0 &&
		(long)best * 10L >= (long)matrix->columns * 9L;
}


static void dual_change_status(
		struct simplex_DualState *state, const int variable,
		const unsigned char new_status)
{
	simplex_degeneracy_update_status(&state->degeneracy, variable,
		state->status[variable], new_status);
	state->status[variable] = new_status;
}


static void dual_change_basis(
		struct simplex_DualState *state, const int position,
		const int new_variable)
{
	simplex_degeneracy_update_basis(&state->degeneracy, position,
		state->basis[position], new_variable);
	state->basis[position] = new_variable;
}


static void dual_destroy(struct simplex_DualState *state)
{
	simplex_dual_feasibility_destroy(&state->feasibility);
	simplex_basis_destroy(&state->factor);
	simplex_csc_destroy(&state->matrix);
	lp_simplex_free(state->basis);
	lp_simplex_free(state->position);
	lp_simplex_free(state->status);
	lp_simplex_free(state->lower);
	lp_simplex_free(state->upper);
	lp_simplex_free(state->cost);
	lp_simplex_free(state->value);
	lp_simplex_free(state->basic_value);
	lp_simplex_free(state->basic_lower);
	lp_simplex_free(state->basic_upper);
	lp_simplex_free(state->reduced);
	lp_simplex_free(state->pi);
	lp_simplex_free(state->rho);
	simplex_dual_pricing_destroy(state);
	simplex_sparse_vector_destroy(&state->packed_direction);
	lp_simplex_free(state->alpha);
	lp_simplex_free(state->breakpoint);
	lp_simplex_free(state->direction);
	lp_simplex_free(state->work);
	lp_simplex_free(state->edge_weight);
	simplex_sparse_vector_destroy(&state->flip_rhs);
}


static int dual_allocate(
		struct simplex_DualState *state, const struct simplex_Problem *problem,
		const struct lp_simplex_Options *options)
{
	int variables = problem->columns + problem->rows;
	lp_simplex_memset(state, 0, sizeof(*state));
	state->rows = problem->rows;
	state->structural = problem->columns;
	state->variables = variables;
	state->options = options;
	state->profile_enabled = getenv("LP_SIMPLEX_PROFILE") != NULL;
	state->profile_rejected_relative = 0;
	state->profile_rejected_ftran = 0;
	state->profile_bound_flips = 0;
	state->profile_flip_batches = 0;
	state->pan_enabled = getenv("LP_SIMPLEX_DISABLE_PAN") == NULL;
	state->matrix = problem->matrix;
	state->matrix.owns_storage = 0;
	state->modal_column_degree = 0;
	state->regular_columns = dual_matrix_has_regular_columns(
		&state->matrix, &state->modal_column_degree);
	state->stable_candidate_order = state->rows >= 2048 &&
		!state->regular_columns;
	state->basis = (int *)lp_simplex_malloc((size_t)problem->rows * sizeof(int));
	state->position = (int *)lp_simplex_malloc((size_t)variables * sizeof(int));
	state->status = (unsigned char *)lp_simplex_malloc(
		(size_t)variables * sizeof(unsigned char));
	state->lower = (double *)lp_simplex_malloc((size_t)variables * sizeof(double));
	state->upper = (double *)lp_simplex_malloc((size_t)variables * sizeof(double));
	state->cost = (double *)lp_simplex_malloc((size_t)variables * sizeof(double));
	state->value = (double *)lp_simplex_malloc((size_t)variables * sizeof(double));
	state->basic_value = (double *)lp_simplex_malloc(
		(size_t)problem->rows * sizeof(double));
	state->basic_lower = (double *)lp_simplex_malloc(
		(size_t)problem->rows * sizeof(double));
	state->basic_upper = (double *)lp_simplex_malloc(
		(size_t)problem->rows * sizeof(double));
	state->reduced = (double *)lp_simplex_malloc((size_t)variables * sizeof(double));
	state->pi = (double *)lp_simplex_malloc((size_t)problem->rows * sizeof(double));
	state->rho = (double *)lp_simplex_malloc((size_t)problem->rows * sizeof(double));
	if (simplex_dual_pricing_create(state) == lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	if (simplex_sparse_vector_create(&state->packed_direction, problem->rows) ==
	    lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	if (simplex_sparse_vector_create(&state->flip_rhs, problem->rows) ==
	    lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	state->alpha = (double *)lp_simplex_malloc((size_t)variables * sizeof(double));
	state->breakpoint = (double *)lp_simplex_malloc(
		(size_t)variables * sizeof(double));
	state->direction = (double *)lp_simplex_malloc((size_t)problem->rows * sizeof(double));
	state->work = (double *)lp_simplex_malloc((size_t)problem->rows * sizeof(double));
	state->edge_weight = (double *)lp_simplex_malloc(
		(size_t)problem->rows * sizeof(double));
	if (simplex_dual_feasibility_create(&state->feasibility, problem->rows) ==
	    lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	if (state->basis == NULL || state->position == NULL || state->status == NULL ||
	    state->lower == NULL || state->upper == NULL || state->cost == NULL ||
	    state->value == NULL || state->basic_value == NULL ||
	    state->basic_lower == NULL || state->basic_upper == NULL ||
	    state->reduced == NULL || state->pi == NULL ||
	    state->rho == NULL || state->alpha == NULL || state->breakpoint == NULL ||
	    state->direction == NULL ||
	    state->work == NULL || state->edge_weight == NULL)
		return lp_simplex_EXIT_FAILURE;
	lp_simplex_memset(state->alpha, 0,
		(size_t)variables * sizeof(double));
	if (simplex_basis_create(&state->factor, &state->matrix, problem->columns,
				 state->basis) == lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	state->factor.allow_sparse_eta = !state->regular_columns ||
		state->modal_column_degree == 1;
	simplex_degeneracy_initialize(&state->degeneracy);
	return lp_simplex_EXIT_SUCCESS;
}


static void dual_set_column_bounds(
		struct simplex_DualState *state, const struct simplex_Problem *problem)
{
	int i, j;
	for (j = 0; j < problem->columns; j++) {
		const struct optm_VariableBound *bound = problem->bounds + j;
		state->lower[j] = __lp_simplex_NINF__;
		state->upper[j] = __lp_simplex_INF__;
		if (bound->b_type == optm_BOUND_T_LO || bound->b_type == optm_BOUND_T_BS)
			state->lower[j] = bound->lb;
		if (bound->b_type == optm_BOUND_T_UP || bound->b_type == optm_BOUND_T_BS)
			state->upper[j] = bound->ub;
		state->cost[j] = problem->objective[j];
	}
	for (i = 0; i < problem->rows; i++) {
		int variable = problem->columns + i;
		int type = problem->row_type[i];
		state->lower[variable] = __lp_simplex_NINF__;
		state->upper[variable] = __lp_simplex_INF__;
		if (type == optm_CONS_T_EQ || type == optm_CONS_T_GE)
			state->lower[variable] = problem->rhs[i];
		if (type == optm_CONS_T_EQ || type == optm_CONS_T_LE)
			state->upper[variable] = problem->rhs[i];
		state->cost[variable] = 0.;
	}
}


/* Tighten structural bounds by interval propagation over the original rows. */
static void dual_propagate_bounds(
		struct simplex_DualState *state, const struct simplex_Problem *problem)
{
	int pass, i, k;
	for (pass = 0; pass < 8; pass++) {
		int changed = 0;
		for (i = 0; i < problem->rows; i++) {
			int type = problem->row_type[i];
			double rhs = problem->rhs[i];
			double finite_minimum = 0., finite_maximum = 0.;
			int minimum_infinite = 0, maximum_infinite = 0;
			for (k = state->matrix.row_start[i];
			     k < state->matrix.row_start[i + 1]; k++) {
				int j = state->matrix.column_index[k];
				double a = state->matrix.row_value[k];
				if ((a > 0. && !dual_is_finite_lower(state->lower[j])) ||
				    (a < 0. && !dual_is_finite_upper(state->upper[j])))
					minimum_infinite++;
				else
					finite_minimum += a * (a > 0.
						? state->lower[j] : state->upper[j]);
				if ((a > 0. && !dual_is_finite_upper(state->upper[j])) ||
				    (a < 0. && !dual_is_finite_lower(state->lower[j])))
					maximum_infinite++;
				else
					finite_maximum += a * (a > 0.
						? state->upper[j] : state->lower[j]);
			}
			for (k = state->matrix.row_start[i];
			     k < state->matrix.row_start[i + 1]; k++) {
				int j = state->matrix.column_index[k];
				double a = state->matrix.row_value[k];
				int own_minimum_infinite, own_maximum_infinite;
				double own_minimum = 0., own_maximum = 0.;
				own_minimum_infinite =
					(a > 0. && !dual_is_finite_lower(state->lower[j])) ||
					(a < 0. && !dual_is_finite_upper(state->upper[j]));
				own_maximum_infinite =
					(a > 0. && !dual_is_finite_upper(state->upper[j])) ||
					(a < 0. && !dual_is_finite_lower(state->lower[j]));
				if (!own_minimum_infinite)
					own_minimum = a * (a > 0.
						? state->lower[j] : state->upper[j]);
				if (!own_maximum_infinite)
					own_maximum = a * (a > 0.
						? state->upper[j] : state->lower[j]);
				if (type != optm_CONS_T_GE &&
				    minimum_infinite - own_minimum_infinite == 0) {
					double bound = (rhs -
						(finite_minimum - own_minimum)) / a;
					if (a > 0. && bound < state->upper[j]) {
						state->upper[j] = bound;
						changed = 1;
					} else if (a < 0. && bound > state->lower[j]) {
						state->lower[j] = bound;
						changed = 1;
					}
				}
				if (type != optm_CONS_T_LE &&
				    maximum_infinite - own_maximum_infinite == 0) {
					double bound = (rhs -
						(finite_maximum - own_maximum)) / a;
					if (a > 0. && bound > state->lower[j]) {
						state->lower[j] = bound;
						changed = 1;
					} else if (a < 0. && bound < state->upper[j]) {
						state->upper[j] = bound;
						changed = 1;
					}
				}
			}
		}
		if (!changed)
			break;
	}
}


static unsigned char dual_status_for_cost(
		const struct simplex_DualState *state, const int variable)
{
	if (state->lower[variable] == state->upper[variable])
		return DUAL_STATUS_FIXED;
	if (state->cost[variable] >= 0. &&
	    dual_is_finite_lower(state->lower[variable]))
		return DUAL_STATUS_LOWER;
	if (state->cost[variable] < 0. &&
	    dual_is_finite_upper(state->upper[variable]))
		return DUAL_STATUS_UPPER;
	if (dual_is_finite_lower(state->lower[variable]))
		return DUAL_STATUS_LOWER;
	if (dual_is_finite_upper(state->upper[variable]))
		return DUAL_STATUS_UPPER;
	return DUAL_STATUS_FREE;
}


static void dual_set_nonbasic_value(
		struct simplex_DualState *state, const int variable)
{
	switch (state->status[variable]) {
	case DUAL_STATUS_LOWER:
	case DUAL_STATUS_FIXED:
		state->value[variable] = state->lower[variable];
		break;
	case DUAL_STATUS_UPPER:
		state->value[variable] = state->upper[variable];
		break;
	case DUAL_STATUS_FREE:
		state->value[variable] = 0.;
		break;
	}
}


static void dual_initialize_basis(
		struct simplex_DualState *state, const struct simplex_Problem *problem)
{
	int i, j;
	for (j = 0; j < state->variables; j++)
		state->position[j] = -1;
	for (j = 0; j < problem->columns; j++) {
		state->status[j] = dual_status_for_cost(state, j);
		dual_set_nonbasic_value(state, j);
	}
	for (i = 0; i < problem->rows; i++) {
		int variable = problem->columns + i;
		state->basis[i] = variable;
		state->position[variable] = i;
		state->status[variable] = DUAL_STATUS_BASIC;
		state->value[variable] = 0.;
	}
}


static double dual_column_dot(
		const struct simplex_DualState *state,
		const int variable, const double *vector)
{
	if (variable >= state->structural)
		return -vector[variable - state->structural];
	return simplex_csc_column_dot(&state->matrix, variable, vector);
}


static int dual_compute_reduced_costs(struct simplex_DualState *state)
{
	int i, j;
	for (i = 0; i < state->rows; i++)
		state->pi[i] = state->cost[state->basis[i]];
	if (simplex_basis_btran(&state->factor, state->pi) == lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	for (j = 0; j < state->variables; j++)
		state->reduced[j] = state->cost[j] -
			dual_column_dot(state, j, state->pi);
	for (i = 0; i < state->rows; i++)
		state->reduced[state->basis[i]] = 0.;
	return lp_simplex_EXIT_SUCCESS;
}


static double dual_variable_infeasibility(
		const struct simplex_DualState *state, const int variable)
{
	double reduced = state->reduced[variable];
	switch (state->status[variable]) {
	case DUAL_STATUS_LOWER:
		return reduced < 0. ? -reduced : 0.;
	case DUAL_STATUS_UPPER:
		return reduced > 0. ? reduced : 0.;
	case DUAL_STATUS_FREE:
		return __lp_simplex_ABS__(reduced);
	default:
		return 0.;
	}
}


static int dual_bound_status_for_reduced(
		const struct simplex_DualState *state, const int variable,
		const double reduced, unsigned char *status)
{
	if (state->lower[variable] == state->upper[variable]) {
		*status = DUAL_STATUS_FIXED;
		return 1;
	}
	if (reduced >= -state->options->dual_tolerance &&
	    dual_is_finite_lower(state->lower[variable])) {
		*status = DUAL_STATUS_LOWER;
		return 1;
	}
	if (reduced <= state->options->dual_tolerance &&
	    dual_is_finite_upper(state->upper[variable])) {
		*status = DUAL_STATUS_UPPER;
		return 1;
	}
	if (__lp_simplex_ABS__(reduced) <= state->options->dual_tolerance &&
	    !dual_is_finite_lower(state->lower[variable]) &&
	    !dual_is_finite_upper(state->upper[variable])) {
		*status = DUAL_STATUS_FREE;
		return 1;
	}
	return 0;
}


/* Build a dual-feasible crash basis without changing the immutable CSC model. */
static int dual_crash(struct simplex_DualState *state)
{
	int attempt, limit = 4 * state->variables + state->rows;
	for (attempt = 0; attempt < limit; attempt++) {
		int i, p = -1, q = -1, rejected = 0;
		unsigned char leaving_status = DUAL_STATUS_FREE;
		if (dual_compute_reduced_costs(state) == lp_simplex_EXIT_FAILURE)
			return lp_simplex_EXIT_FAILURE;
		if (attempt == 0 && state->profile_enabled) {
			int infeasible = 0;
			for (i = 0; i < state->variables; i++)
				if (state->position[i] < 0 &&
				    dual_variable_infeasibility(state, i) >
				    state->options->dual_tolerance)
					infeasible++;
			fprintf(stderr,
				"dual profile: crash initial infeasible variables=%d\n",
				infeasible);
		}
		lp_simplex_memset(state->candidate_sign, 0,
			(size_t)state->variables * sizeof(int));
		for (;;) {
			double best_pivot = 0.;
			q = -1;
			for (i = 0; i < state->variables; i++) {
				double violation;
				if (state->position[i] >= 0 || state->candidate_sign[i])
					continue;
				violation = dual_variable_infeasibility(state, i);
				if (violation > state->options->dual_tolerance) {
					q = i;
					break;
				}
			}
			if (q < 0) {
				if (rejected == 0)
					return lp_simplex_EXIT_SUCCESS;
				if (state->profile_enabled)
					fprintf(stderr,
						"dual profile: crash stalled attempt=%d rejected=%d\n",
						attempt, rejected);
				return lp_simplex_EXIT_FAILURE;
			}
			simplex_csc_column_to_dense(&state->matrix, state->structural,
						    q, state->direction);
			if (simplex_basis_ftran(&state->factor, state->direction) ==
			    lp_simplex_EXIT_FAILURE)
				return lp_simplex_EXIT_FAILURE;
			p = -1;
			for (i = 0; i < state->rows; i++) {
				int leaving = state->basis[i];
				double pivot = state->direction[i];
				double leaving_reduced;
				unsigned char proposed;
				if (__lp_simplex_ABS__(pivot) <=
				    state->options->pivot_tolerance)
					continue;
				leaving_reduced = -state->reduced[q] / pivot;
				if (!dual_bound_status_for_reduced(state, leaving,
							   leaving_reduced, &proposed))
					continue;
				if (__lp_simplex_ABS__(pivot) > best_pivot) {
					best_pivot = __lp_simplex_ABS__(pivot);
					p = i;
					leaving_status = proposed;
				}
			}
			if (p >= 0)
				break;
			state->candidate_sign[q] = 1;
			rejected++;
		}
		{
			int leaving = state->basis[p];
			state->structural_basic +=
				(q < state->structural) - (leaving < state->structural);
			state->position[leaving] = -1;
			dual_change_status(state, leaving, leaving_status);
			dual_set_nonbasic_value(state, leaving);
		}
		dual_change_basis(state, p, q);
		state->position[q] = p;
		dual_change_status(state, q, DUAL_STATUS_BASIC);
		{
			int update = simplex_basis_update(&state->factor, p,
						  state->direction);
			if (update == 1) {
				if (simplex_basis_factorize(&state->factor) ==
				    lp_simplex_EXIT_FAILURE)
					return lp_simplex_EXIT_FAILURE;
			} else if (update == lp_simplex_EXIT_FAILURE) {
				return lp_simplex_EXIT_FAILURE;
			}
		}
	}
	if (state->profile_enabled)
		fprintf(stderr, "dual profile: crash iteration limit=%d\n", limit);
	return lp_simplex_EXIT_FAILURE;
}


static int dual_compute_primal_values(struct simplex_DualState *state)
{
	int i, j;
	lp_simplex_memset(state->work, 0, (size_t)state->rows * sizeof(double));
	for (j = 0; j < state->variables; j++) {
		if (state->position[j] >= 0)
			continue;
		dual_set_nonbasic_value(state, j);
		if (state->value[j] != 0.)
			simplex_csc_column_axpy(&state->matrix, state->structural,
						  j, -state->value[j], state->work);
	}
	if (simplex_basis_ftran(&state->factor, state->work) ==
	    lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	for (i = 0; i < state->rows; i++) {
		state->value[state->basis[i]] = state->basic_value[i] = state->work[i];
		state->basic_lower[i] = state->lower[state->basis[i]];
		state->basic_upper[i] = state->upper[state->basis[i]];
	}
	return lp_simplex_EXIT_SUCCESS;
}


/* Controller-level reinversion transaction.  All exact-state consumers are
 * refreshed together so the factor, primal/dual values and hypersparse
 * infeasibility queues can never describe different bases. */
static int dual_reinvert(struct simplex_DualState *state)
{
	if (simplex_basis_factorize(&state->factor) == lp_simplex_EXIT_FAILURE ||
	    dual_compute_primal_values(state) == lp_simplex_EXIT_FAILURE ||
	    dual_compute_reduced_costs(state) == lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	state->reinversions++;
	simplex_dual_feasibility_rebuild(state);
	return lp_simplex_EXIT_SUCCESS;
}


static int dual_update_edge_weights(
		struct simplex_DualState *state, const int p)
{
	int k;
	double pivot = state->direction[p];
	double old_p = state->edge_weight[p];
	if (__lp_simplex_ABS__(pivot) <= state->options->pivot_tolerance)
		return lp_simplex_EXIT_FAILURE;
	for (k = 0; k < state->packed_direction.count; k++) {
		int i = state->packed_direction.index[k];
		double ratio;
		double updated;
		if (i == p)
			continue;
		ratio = state->packed_direction.value[k] / pivot;
		updated = state->edge_weight[i] - 2. * ratio * state->work[i] +
			ratio * ratio * old_p;
		state->edge_weight[i] = __lp_simplex_MAX__(updated, 1e-12);
	}
	state->edge_weight[p] = __lp_simplex_MAX__(old_p / (pivot * pivot), 1e-12);
	return lp_simplex_EXIT_SUCCESS;
}


static int dual_opposite_bound_crossed(
		const struct simplex_DualState *state, const int variable,
		const int sign, const double new_value, double *movement,
		unsigned char *new_status)
{
	if (sign > 0 && dual_is_finite_upper(state->upper[variable]) &&
	    new_value > state->upper[variable] + state->options->primal_tolerance) {
		*movement = state->upper[variable] - state->value[variable];
		*new_status = DUAL_STATUS_UPPER;
		return 1;
	}
	if (sign < 0 && dual_is_finite_lower(state->lower[variable]) &&
	    new_value < state->lower[variable] - state->options->primal_tolerance) {
		*movement = state->lower[variable] - state->value[variable];
		*new_status = DUAL_STATUS_LOWER;
		return 1;
	}
	return 0;
}


static int dual_accumulate_bound_flip(
		struct simplex_DualState *state, const int variable,
		const double movement)
{
	int k;
	if (variable >= state->structural)
		return simplex_sparse_vector_add(&state->flip_rhs,
			variable - state->structural, -movement);
	for (k = state->matrix.column_start[variable];
	     k < state->matrix.column_start[variable + 1]; k++)
		if (simplex_sparse_vector_add(&state->flip_rhs,
				state->matrix.row_index[k],
				movement * state->matrix.value[k]) ==
		    lp_simplex_EXIT_FAILURE)
			return lp_simplex_EXIT_FAILURE;
	return lp_simplex_EXIT_SUCCESS;
}


static int dual_flush_bound_flips(
		struct simplex_DualState *state, const int leaving_position)
{
	int i;
	if (state->flip_rhs.count == 0)
		return lp_simplex_EXIT_SUCCESS;
	lp_simplex_memcpy(state->direction, state->flip_rhs.dense,
		(size_t)state->rows * sizeof(double));
	if (simplex_basis_ftran(&state->factor, state->direction) ==
	    lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	simplex_sparse_vector_pack(
		&state->packed_direction, state->direction, 0.);
	for (i = 0; i < state->packed_direction.count; i++) {
		int row = state->packed_direction.index[i];
		/* The leaving-row value is maintained incrementally from alpha so
		 * subsequent flips can be classified before this batched FTRAN. */
		if (row != leaving_position)
			state->basic_value[row] -=
				state->packed_direction.value[i];
	}
	simplex_dual_feasibility_update_packed(
		state, &state->packed_direction);
	simplex_sparse_vector_clear(&state->flip_rhs);
	state->profile_flip_batches++;
	return lp_simplex_EXIT_SUCCESS;
}


static double dual_max_dual_infeasibility(const struct simplex_DualState *state)
{
	int j;
	double maximum = 0.;
	for (j = 0; j < state->variables; j++) {
		double violation;
		if (state->position[j] >= 0)
			continue;
		violation = dual_variable_infeasibility(state, j);
		if (violation > maximum)
			maximum = violation;
	}
	return maximum;
}


static void dual_update_dual_values(
		struct simplex_DualState *state, const int kappa,
		const double theta, const int leaving, const int entering)
{
	int i;
	double tau = -kappa * theta;
	if (tau != 0.) {
		for (i = 0; i < state->rows; i++)
			state->pi[i] += tau * state->rho[i];
		if (state->packed_alpha_valid) {
			for (i = 0; i < state->packed_alpha.count; i++) {
				int j = state->packed_alpha.index[i];
				/* Reduced costs are maintained only for the active
				 * nonbasic set.  Basic entries must remain exactly zero. */
				if (state->nonbasic_slot[j] >= 0)
					state->reduced[j] -=
						tau * state->packed_alpha.value[i];
			}
		} else {
			lp_simplex_linalg_daxpy(state->variables, -tau,
				state->alpha, 1, state->reduced, 1);
		}
	}
	state->reduced[entering] = 0.;
	state->reduced[leaving] = -tau;
}


static int dual_run(struct simplex_DualState *state, int *status)
{
	int terminal_repairs = 0;
	while (state->iterations < state->options->iteration_limit) {
		int i, p, q, kappa, leaving, sign, restart_iteration = 0;
		int compact_was_active;
		int recomputed_dual = 0;
		int bound_flips = 0;
		int weight_work_ready = 0;
		int rejected_relative = 0;
		double relative_pivot_tolerance = DUAL_RELATIVE_PIVOT_TOLERANCE;
		double target, primal_infeasibility, theta = 0.;
		double delta, new_value, movement = 0.;
		unsigned char flipped_status = DUAL_STATUS_FREE;
		simplex_sparse_vector_clear(&state->flip_rhs);
		p = simplex_dual_feasibility_choose(state, &target, &kappa,
					&primal_infeasibility,
					state->pan_enabled &&
					simplex_degeneracy_should_probe(&state->degeneracy) &&
					dual_pan_structure_enabled(state));
		if (p < 0) {
			if (dual_reinvert(state) == lp_simplex_EXIT_FAILURE) {
				*status = lp_simplex_Singularity;
				return lp_simplex_EXIT_FAILURE;
			}
			p = simplex_dual_feasibility_choose(state, &target, &kappa,
						&primal_infeasibility,
						state->pan_enabled &&
						simplex_degeneracy_should_probe(&state->degeneracy) &&
						dual_pan_structure_enabled(state));
			if (p < 0) {
				double dual_error = dual_max_dual_infeasibility(state);
				if (dual_error <= state->options->dual_tolerance) {
					*status = lp_simplex_Success;
					return lp_simplex_EXIT_SUCCESS;
				}
				if (dual_error <= 10. * state->options->dual_tolerance &&
				    terminal_repairs < 2 &&
				    dual_crash(state) == lp_simplex_EXIT_SUCCESS &&
				    dual_reinvert(state) == lp_simplex_EXIT_SUCCESS) {
					terminal_repairs++;
					for (i = 0; i < state->rows; i++)
						state->edge_weight[i] = 1.;
					simplex_dual_feasibility_rebuild(state);
					lp_simplex_memset(state->alpha, 0,
						(size_t)state->variables * sizeof(double));
					lp_simplex_memset(state->candidate_sign, 0,
						(size_t)state->variables * sizeof(int));
					state->candidate_count = 0;
					simplex_dual_nonbasic_initialize(state);
					state->structural_basic = 0;
					for (i = 0; i < state->rows; i++)
						if (state->basis[i] < state->structural)
							state->structural_basic++;
					continue;
				}
				*status = lp_simplex_PrecisionError;
				return lp_simplex_EXIT_FAILURE;
			}
			continue;
		}
		lp_simplex_memset(state->rho, 0,
			(size_t)state->rows * sizeof(double));
		state->rho[p] = 1.;
		if (simplex_basis_btran(&state->factor, state->rho) ==
		    lp_simplex_EXIT_FAILURE) {
			*status = lp_simplex_Singularity;
			return lp_simplex_EXIT_FAILURE;
		}
		{
			double exact_weight = 0.;
			for (i = 0; i < state->rows; i++)
				exact_weight += state->rho[i] * state->rho[i];
			state->edge_weight[p] =
				__lp_simplex_MAX__(exact_weight, 1e-12);
		}
		if (simplex_dual_prepare_ratio_test(state, kappa) == 0) {
			if (!simplex_dual_nonbasic_is_consistent(state)) {
				if (state->profile_enabled)
					fprintf(stderr,
						"dual profile: repairing nonbasic index set\n");
				simplex_dual_nonbasic_initialize(state);
				if (simplex_dual_prepare_ratio_test(state, kappa) != 0)
					continue;
			}
			if (state->factor.update_count != 0 &&
			    dual_reinvert(state) == lp_simplex_EXIT_SUCCESS)
				continue;
			*status = lp_simplex_Infeasibility;
			return lp_simplex_EXIT_FAILURE;
		}
		q = -1;
		for (;;) {
			double direction_maximum = 0.;
			q = simplex_dual_next_ratio_candidate(state, &theta,
				state->pan_enabled &&
				simplex_degeneracy_should_probe(&state->degeneracy),
				state->basis[p]);
			state->ratio_minimum_valid = 0;
			if (q < 0) {
				if (bound_flips > 0) {
					if (dual_flush_bound_flips(state, p) ==
					    lp_simplex_EXIT_FAILURE) {
						*status = lp_simplex_Singularity;
						return lp_simplex_EXIT_FAILURE;
					}
					restart_iteration = 1;
					break;
				}
				if (rejected_relative > 0 &&
				    relative_pivot_tolerance > 1e-14) {
					relative_pivot_tolerance *= 1e-4;
					rejected_relative = 0;
					simplex_dual_prepare_ratio_test(state, kappa);
					continue;
				}
				if (state->factor.update_count != 0 &&
				    dual_reinvert(state) == lp_simplex_EXIT_SUCCESS) {
					restart_iteration = 1;
					break;
				}
				*status = lp_simplex_Infeasibility;
				return lp_simplex_EXIT_FAILURE;
			}
			/* BFRT classifies a bound flip from the tableau-row coefficient
			 * alpha.  Do this before FTRAN and aggregate all flipped columns
			 * into one RHS; only an actual entering column needs its individual
			 * direction and stability validation. */
			if (!state->regular_columns &&
			    __lp_simplex_ABS__(state->alpha[q]) >
			    state->options->pivot_tolerance) {
				double tentative_delta =
					(state->basic_value[p] - target) / state->alpha[q];
				double tentative_value = state->value[q] + tentative_delta;
				sign = state->candidate_sign[q];
				if (sign * tentative_delta >=
					    -state->options->primal_tolerance &&
				    state->status[q] != DUAL_STATUS_FREE &&
				    dual_opposite_bound_crossed(state, q, sign,
					    tentative_value, &movement, &flipped_status) &&
				    movement != 0.) {
					if (dual_accumulate_bound_flip(state, q, movement) ==
					    lp_simplex_EXIT_FAILURE) {
						*status = lp_simplex_PrecisionError;
						return lp_simplex_EXIT_FAILURE;
					}
					state->basic_value[p] -= state->alpha[q] * movement;
					state->value[q] += movement;
					dual_change_status(state, q, flipped_status);
					state->candidate_sign[q] = 0;
					bound_flips++;
					state->profile_bound_flips++;
					movement = 0.;
					continue;
				}
				movement = 0.;
			}
			if (dual_flush_bound_flips(state, p) ==
			    lp_simplex_EXIT_FAILURE) {
				*status = lp_simplex_Singularity;
				return lp_simplex_EXIT_FAILURE;
			}
			simplex_csc_column_to_dense(&state->matrix, state->structural,
						    q, state->direction);
			if (!weight_work_ready) {
				lp_simplex_memcpy(state->work, state->rho,
					(size_t)state->rows * sizeof(double));
				if (simplex_basis_ftran_pair(&state->factor,
						state->direction, state->work) ==
				    lp_simplex_EXIT_FAILURE) {
					state->profile_rejected_ftran++;
					state->candidate_sign[q] = 0;
					continue;
				}
				weight_work_ready = 1;
			} else if (simplex_basis_ftran(&state->factor, state->direction) ==
				   lp_simplex_EXIT_FAILURE) {
				state->profile_rejected_ftran++;
				state->candidate_sign[q] = 0;
				continue;
			}
			simplex_sparse_vector_pack(
				&state->packed_direction, state->direction, 0.);
			for (i = 0; i < state->packed_direction.count; i++)
				direction_maximum = __lp_simplex_MAX__(direction_maximum,
					__lp_simplex_ABS__(
						state->packed_direction.value[i]));
			if (__lp_simplex_ABS__(state->direction[p]) <=
			    state->options->pivot_tolerance ||
			    __lp_simplex_ABS__(state->direction[p]) <
			    relative_pivot_tolerance * direction_maximum) {
				state->candidate_sign[q] = 0;
				rejected_relative++;
				state->profile_rejected_relative++;
				continue;
			}
			if (__lp_simplex_ABS__(state->direction[p] - state->alpha[q]) >
			    100. * state->options->pivot_tolerance *
			    __lp_simplex_MAX__(1.,
				__lp_simplex_ABS__(state->alpha[q]))) {
				if (state->factor.update_count != 0 &&
				    dual_reinvert(state) == lp_simplex_EXIT_SUCCESS) {
					restart_iteration = 1;
					break;
				}
				*status = lp_simplex_PrecisionError;
				return lp_simplex_EXIT_FAILURE;
			}
			delta = (state->basic_value[p] - target) /
				state->direction[p];
			sign = state->candidate_sign[q];
			new_value = state->value[q] + delta;
			if (sign * delta >= -state->options->primal_tolerance &&
			    (state->status[q] == DUAL_STATUS_FREE ||
			     !dual_opposite_bound_crossed(state, q, sign, new_value,
							 &movement, &flipped_status)))
				break;
			if (state->status[q] == DUAL_STATUS_FREE || movement == 0.) {
				state->candidate_sign[q] = 0;
				continue;
			}
			for (i = 0; i < state->packed_direction.count; i++) {
				int row = state->packed_direction.index[i];
				state->basic_value[row] -=
					state->packed_direction.value[i] * movement;
			}
			simplex_dual_feasibility_update_packed(
				state, &state->packed_direction);
			state->value[q] += movement;
			dual_change_status(state, q, flipped_status);
			state->candidate_sign[q] = 0;
			bound_flips++;
			state->profile_bound_flips++;
			movement = 0.;
		}
		if (restart_iteration)
			continue;
		for (i = 0; i < state->packed_direction.count; i++) {
			int row = state->packed_direction.index[i];
			if (row != p)
				state->basic_value[row] -=
					state->packed_direction.value[i] * delta;
		}
		if (dual_update_edge_weights(state, p) == lp_simplex_EXIT_FAILURE) {
			*status = lp_simplex_PrecisionError;
			return lp_simplex_EXIT_FAILURE;
		}
		{
			leaving = state->basis[p];
			dual_update_dual_values(state, kappa, theta, leaving, q);
			simplex_dual_nonbasic_remove(state, q);
			state->structural_basic +=
				(q < state->structural) - (leaving < state->structural);
			state->value[leaving] = target;
			dual_change_status(state, leaving, kappa > 0
				? DUAL_STATUS_LOWER : DUAL_STATUS_UPPER);
			if (state->lower[leaving] == state->upper[leaving])
				dual_change_status(state, leaving, DUAL_STATUS_FIXED);
			else
				simplex_dual_nonbasic_add(state, leaving);
			state->position[leaving] = -1;
			dual_change_basis(state, p, q);
			state->basic_value[p] = new_value;
			state->basic_lower[p] = state->lower[q];
			state->basic_upper[p] = state->upper[q];
			state->position[q] = p;
			dual_change_status(state, q, DUAL_STATUS_BASIC);
			state->alpha[q] = 0.;
			state->value[q] = new_value;
		}
		simplex_dual_feasibility_update_packed(
			state, &state->packed_direction);
		state->iterations++;
		compact_was_active = state->factor.compact_active;
		simplex_degeneracy_observe(&state->degeneracy,
			theta * primal_infeasibility,
			10. * state->options->dual_tolerance,
			state->factor.compact_ever_active);
		if (state->pan_enabled &&
		    simplex_degeneracy_should_probe(&state->degeneracy))
			state->factor.compact_requested = 1;
		if (simplex_degeneracy_should_probe(&state->degeneracy))
			simplex_degeneracy_record_state(&state->degeneracy,
				state->basis, state->rows,
				state->status, state->variables);
		{
			int update;
			int compact_desired;
			compact_desired = state->factor.compact_requested &&
				state->rows >= 4096 &&
				state->structural_basic * 3 <= state->rows;
			if (!compact_was_active && compact_desired)
				update = 1;
			else
				update = simplex_basis_update(&state->factor, p,
							state->direction);
			if (update == 1) {
				update = dual_reinvert(state);
				if (update == lp_simplex_EXIT_SUCCESS) {
					recomputed_dual = 1;
				}
			}
			if (update == lp_simplex_EXIT_FAILURE) {
				*status = lp_simplex_Singularity;
				return lp_simplex_EXIT_FAILURE;
			}
		}
		if (recomputed_dual && dual_max_dual_infeasibility(state) >
		    10. * state->options->dual_tolerance) {
			if (dual_reinvert(state) == lp_simplex_EXIT_FAILURE) {
				*status = lp_simplex_PrecisionError;
				return lp_simplex_EXIT_FAILURE;
			}
			if (dual_max_dual_infeasibility(state) >
			    10. * state->options->dual_tolerance) {
				if (dual_crash(state) == lp_simplex_EXIT_FAILURE ||
				    dual_reinvert(state) == lp_simplex_EXIT_FAILURE ||
				    dual_max_dual_infeasibility(state) >
				    10. * state->options->dual_tolerance) {
					*status = lp_simplex_PrecisionError;
					return lp_simplex_EXIT_FAILURE;
				}
				for (i = 0; i < state->rows; i++)
					state->edge_weight[i] = 1.;
				simplex_dual_feasibility_rebuild(state);
				lp_simplex_memset(state->alpha, 0,
					(size_t)state->variables * sizeof(double));
				lp_simplex_memset(state->candidate_sign, 0,
					(size_t)state->variables * sizeof(int));
				state->candidate_count = 0;
				simplex_dual_nonbasic_initialize(state);
				state->structural_basic = 0;
				for (i = 0; i < state->rows; i++)
					if (state->basis[i] < state->structural)
						state->structural_basic++;
			}
		}
	}
	*status = lp_simplex_ExceedIterLimit;
	return lp_simplex_EXIT_FAILURE;
}


static double dual_primal_infeasibility(
		const struct simplex_DualState *state)
{
	int i, j;
	double maximum = 0.;
	for (j = 0; j < state->variables; j++) {
		double violation = 0.;
		if (state->value[j] < state->lower[j])
			violation = state->lower[j] - state->value[j];
		else if (state->value[j] > state->upper[j])
			violation = state->value[j] - state->upper[j];
		if (violation > maximum)
			maximum = violation;
	}
	lp_simplex_memset(state->work, 0,
		(size_t)state->rows * sizeof(double));
	for (j = 0; j < state->structural; j++) {
		if (state->value[j] != 0.)
			simplex_csc_column_axpy(&state->matrix, state->structural,
						  j, state->value[j], state->work);
	}
	for (i = 0; i < state->rows; i++) {
		double residual = __lp_simplex_ABS__(
			state->work[i] - state->value[state->structural + i]);
		if (residual > maximum)
			maximum = residual;
	}
	return maximum;
}


int simplex_dual_solve_problem(
		const struct simplex_Problem *problem,
		const struct lp_simplex_Options *options,
		double *x, double *row_dual, struct lp_simplex_Result *result,
		const int propagate_bounds)
{
	struct simplex_DualState state;
	int j, status = lp_simplex_CondUnsatisfied;
	int solve_state;
	double objective = 0.;
	clock_t profile_started = 0;
	if (getenv("LP_SIMPLEX_PROFILE") != NULL)
		profile_started = clock();
	if (dual_allocate(&state, problem, options) == lp_simplex_EXIT_FAILURE) {
		dual_destroy(&state);
		result->status = lp_simplex_MemoryAllocError;
		return lp_simplex_EXIT_FAILURE;
	}
	dual_set_column_bounds(&state, problem);
	if (propagate_bounds)
		dual_propagate_bounds(&state, problem);
	dual_initialize_basis(&state, problem);
	{
		const char *stage = NULL;
		if (simplex_basis_factorize(&state.factor) == lp_simplex_EXIT_FAILURE)
			stage = "initial factorization";
		else if (dual_crash(&state) == lp_simplex_EXIT_FAILURE)
			stage = "dual crash";
		else if (simplex_basis_factorize(&state.factor) ==
			 lp_simplex_EXIT_FAILURE)
			stage = "crash-basis factorization";
		else if (dual_compute_primal_values(&state) == lp_simplex_EXIT_FAILURE)
			stage = "initial primal values";
		else if (dual_compute_reduced_costs(&state) == lp_simplex_EXIT_FAILURE)
			stage = "initial reduced costs";
		else
			for (j = 0; j < state.rows; j++)
				state.edge_weight[j] = 1.;
		if (stage != NULL) {
			if (state.profile_enabled)
				fprintf(stderr, "dual profile: initialization failed at %s\n",
					stage);
			result->status = lp_simplex_PrecisionError;
			dual_destroy(&state);
			return lp_simplex_EXIT_FAILURE;
		}
	}
	simplex_dual_feasibility_rebuild(&state);
	lp_simplex_memset(state.alpha, 0,
		(size_t)state.variables * sizeof(double));
	lp_simplex_memset(state.candidate_sign, 0,
		(size_t)state.variables * sizeof(int));
	state.candidate_count = 0;
	simplex_dual_nonbasic_initialize(&state);
	state.structural_basic = 0;
	for (j = 0; j < state.rows; j++)
		if (state.basis[j] < state.structural)
			state.structural_basic++;
	solve_state = dual_run(&state, &status);
	for (j = 0; j < state.rows; j++)
		state.value[state.basis[j]] = state.basic_value[j];
	for (j = 0; j < problem->columns; j++) {
		x[j] = state.value[j];
		objective += problem->objective[j] * x[j];
	}
	result->status = status;
	result->iterations = state.iterations;
	result->objective = objective;
	result->primal_infeasibility = dual_primal_infeasibility(&state);
	result->dual_infeasibility = dual_max_dual_infeasibility(&state);
	if (row_dual != NULL)
		for (j = 0; j < problem->rows; j++)
			row_dual[j] = state.pi[j];
	if (state.profile_enabled) {
		int active_rank = 0;
		for (j = 0; j < state.rows; j++) {
			int variable = state.basis[j];
			if (state.basic_value[j] > state.lower[variable] +
			    state.options->primal_tolerance &&
			    state.basic_value[j] < state.upper[variable] -
			    state.options->primal_tolerance)
				active_rank++;
		}
		state.profile_total_seconds =
			(double)(clock() - profile_started) / (double)CLOCKS_PER_SEC;
		fprintf(stderr,
			"dual profile: total=%.6f factor=%.6f/%ld ftran=%.6f/%ld "
			"btran=%.6f/%ld ratio=%.6f residual=%.6f "
			"pan_degenerate=%ld pan_activations=%ld pan_probes=%ld "
			"pan_rank=%d/%d factor_core=%d/%d compact=%ld[%d,%d] "
			"eta_nnz=%ld/%ld cert=%ld/%ld refine=%ld/%ld "
			"reject=%ld/%ld flips=%ld/%ld\n",
			state.profile_total_seconds,
			state.factor.profile_factor_seconds,
			state.factor.profile_factor_calls,
			state.factor.profile_ftran_seconds,
			state.factor.profile_ftran_calls,
			state.factor.profile_btran_seconds,
			state.factor.profile_btran_calls,
			state.profile_ratio_seconds,
			state.profile_total_seconds -
			state.factor.profile_factor_seconds -
			state.factor.profile_ftran_seconds -
			state.factor.profile_btran_seconds -
			state.profile_ratio_seconds,
			state.degeneracy.degenerate_pivots,
			state.degeneracy.activations,
			state.degeneracy.probes, active_rank, state.rows,
			state.factor.factor_size, state.rows,
			state.factor.profile_compact_calls,
			state.factor.profile_compact_min,
			state.factor.profile_compact_max,
			state.factor.profile_eta_nonzeros,
			state.factor.profile_eta_slots,
			state.factor.profile_compact_ftran_validations,
			state.factor.profile_compact_btran_validations,
			state.factor.profile_compact_ftran_refinements,
			state.factor.profile_compact_btran_refinements,
			state.profile_rejected_relative,
			state.profile_rejected_ftran,
			state.profile_bound_flips,
			state.profile_flip_batches);
	}
	dual_destroy(&state);
	return solve_state;
}


int simplex_dual_solve(
		const struct lp_Model *model,
		const struct lp_simplex_Options *options,
		double *x, double *row_dual, struct lp_simplex_Result *result,
		const int propagate_bounds)
{
	struct simplex_Problem problem;
	int state;
	if (simplex_problem_from_model(&problem, model) ==
	    lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	state = simplex_dual_solve_problem(&problem, options, x, row_dual,
		result, propagate_bounds);
	simplex_problem_destroy(&problem);
	return state;
}
