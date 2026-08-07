/* Bounded dual revised simplex over an immutable CSC constraint matrix. */
#include "simplex_dual.h"
#include "simplex_basis.h"
#include "simplex_csc.h"
#include "utils.h"
#include <lp_simplex/status.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define DUAL_STATUS_BASIC 0
#define DUAL_STATUS_LOWER 1
#define DUAL_STATUS_UPPER 2
#define DUAL_STATUS_FIXED 3
#define DUAL_STATUS_FREE  4
#define DUAL_RELATIVE_PIVOT_TOLERANCE 1e-8


struct simplex_DualState {
	int rows;
	int structural;
	int variables;
	int iterations;
	int reinversions;
	const struct lp_simplex_Options *options;
	struct simplex_CscMatrix matrix;
	struct simplex_Basis factor;
	int *basis;
	int *position;
	unsigned char *status;
	double *lower;
	double *upper;
	double *cost;
	double *value;
	double *reduced;
	double *pi;
	double *rho;
	double *alpha;
	double *breakpoint;
	double *direction;
	double *work;
	double *edge_weight;
	int *candidate_sign;
	int profile_enabled;
	double profile_ratio_seconds;
	double profile_total_seconds;
};


static int dual_is_finite_lower(double value)
{
	return value > __lp_simplex_NINF__;
}


static int dual_is_finite_upper(double value)
{
	return value < __lp_simplex_INF__;
}


static void dual_destroy(struct simplex_DualState *state)
{
	simplex_basis_destroy(&state->factor);
	simplex_csc_destroy(&state->matrix);
	lp_simplex_free(state->basis);
	lp_simplex_free(state->position);
	lp_simplex_free(state->status);
	lp_simplex_free(state->lower);
	lp_simplex_free(state->upper);
	lp_simplex_free(state->cost);
	lp_simplex_free(state->value);
	lp_simplex_free(state->reduced);
	lp_simplex_free(state->pi);
	lp_simplex_free(state->rho);
	lp_simplex_free(state->alpha);
	lp_simplex_free(state->breakpoint);
	lp_simplex_free(state->direction);
	lp_simplex_free(state->work);
	lp_simplex_free(state->edge_weight);
	lp_simplex_free(state->candidate_sign);
}


static int dual_allocate(
		struct simplex_DualState *state, const struct lp_Model *model,
		const struct lp_simplex_Options *options)
{
	int variables = model->n + model->m;
	lp_simplex_memset(state, 0, sizeof(*state));
	state->rows = model->m;
	state->structural = model->n;
	state->variables = variables;
	state->options = options;
	state->profile_enabled = getenv("LP_SIMPLEX_PROFILE") != NULL;
	if (simplex_csc_from_model(model, &state->matrix) == lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	state->basis = (int *)lp_simplex_malloc((size_t)model->m * sizeof(int));
	state->position = (int *)lp_simplex_malloc((size_t)variables * sizeof(int));
	state->status = (unsigned char *)lp_simplex_malloc(
		(size_t)variables * sizeof(unsigned char));
	state->lower = (double *)lp_simplex_malloc((size_t)variables * sizeof(double));
	state->upper = (double *)lp_simplex_malloc((size_t)variables * sizeof(double));
	state->cost = (double *)lp_simplex_malloc((size_t)variables * sizeof(double));
	state->value = (double *)lp_simplex_malloc((size_t)variables * sizeof(double));
	state->reduced = (double *)lp_simplex_malloc((size_t)variables * sizeof(double));
	state->pi = (double *)lp_simplex_malloc((size_t)model->m * sizeof(double));
	state->rho = (double *)lp_simplex_malloc((size_t)model->m * sizeof(double));
	state->alpha = (double *)lp_simplex_malloc((size_t)variables * sizeof(double));
	state->breakpoint = (double *)lp_simplex_malloc(
		(size_t)variables * sizeof(double));
	state->direction = (double *)lp_simplex_malloc((size_t)model->m * sizeof(double));
	state->work = (double *)lp_simplex_malloc((size_t)model->m * sizeof(double));
	state->edge_weight = (double *)lp_simplex_malloc(
		(size_t)model->m * sizeof(double));
	state->candidate_sign = (int *)lp_simplex_malloc((size_t)variables * sizeof(int));
	if (state->basis == NULL || state->position == NULL || state->status == NULL ||
	    state->lower == NULL || state->upper == NULL || state->cost == NULL ||
	    state->value == NULL || state->reduced == NULL || state->pi == NULL ||
	    state->rho == NULL || state->alpha == NULL || state->breakpoint == NULL ||
	    state->direction == NULL ||
	    state->work == NULL || state->edge_weight == NULL ||
	    state->candidate_sign == NULL)
		return lp_simplex_EXIT_FAILURE;
	if (simplex_basis_create(&state->factor, &state->matrix, model->n,
				 state->basis) == lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	return lp_simplex_EXIT_SUCCESS;
}


static void dual_set_column_bounds(
		struct simplex_DualState *state, const struct lp_Model *model)
{
	int i, j;
	for (j = 0; j < model->n; j++) {
		const struct optm_VariableBound *bound = model->bounds + j;
		state->lower[j] = __lp_simplex_NINF__;
		state->upper[j] = __lp_simplex_INF__;
		if (bound->b_type == optm_BOUND_T_LO || bound->b_type == optm_BOUND_T_BS)
			state->lower[j] = bound->lb;
		if (bound->b_type == optm_BOUND_T_UP || bound->b_type == optm_BOUND_T_BS)
			state->upper[j] = bound->ub;
		state->cost[j] = model->objective[j];
	}
	for (i = 0; i < model->m; i++) {
		int variable = model->n + i;
		const struct optm_LinearConstraint *constraint = model->constraints + i;
		state->lower[variable] = __lp_simplex_NINF__;
		state->upper[variable] = __lp_simplex_INF__;
		if (constraint->type == optm_CONS_T_EQ ||
		    constraint->type == optm_CONS_T_GE)
			state->lower[variable] = constraint->rhs;
		if (constraint->type == optm_CONS_T_EQ ||
		    constraint->type == optm_CONS_T_LE)
			state->upper[variable] = constraint->rhs;
		state->cost[variable] = 0.;
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
		struct simplex_DualState *state, const struct lp_Model *model)
{
	int i, j;
	for (j = 0; j < state->variables; j++)
		state->position[j] = -1;
	for (j = 0; j < model->n; j++) {
		state->status[j] = dual_status_for_cost(state, j);
		dual_set_nonbasic_value(state, j);
	}
	for (i = 0; i < model->m; i++) {
		int variable = model->n + i;
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
		int i, p = -1, q = -1;
		double largest = state->options->dual_tolerance;
		double best_pivot = 0.;
		unsigned char leaving_status = DUAL_STATUS_FREE;
		if (dual_compute_reduced_costs(state) == lp_simplex_EXIT_FAILURE)
			return lp_simplex_EXIT_FAILURE;
		for (i = 0; i < state->variables; i++) {
			double violation;
			if (state->position[i] >= 0)
				continue;
			violation = dual_variable_infeasibility(state, i);
			if (violation > largest) {
				largest = violation;
				q = i;
			}
		}
		if (q < 0)
			return lp_simplex_EXIT_SUCCESS;
		simplex_csc_column_to_dense(&state->matrix, state->structural,
					    q, state->direction);
		if (simplex_basis_ftran(&state->factor, state->direction) ==
		    lp_simplex_EXIT_FAILURE)
			return lp_simplex_EXIT_FAILURE;
		for (i = 0; i < state->rows; i++) {
			int leaving = state->basis[i];
			double pivot = state->direction[i];
			double leaving_reduced;
			unsigned char proposed;
			if (__lp_simplex_ABS__(pivot) <= state->options->pivot_tolerance)
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
		if (p < 0) {
			return lp_simplex_EXIT_FAILURE;
		}
		state->position[state->basis[p]] = -1;
		state->status[state->basis[p]] = leaving_status;
		dual_set_nonbasic_value(state, state->basis[p]);
		state->basis[p] = q;
		state->position[q] = p;
		state->status[q] = DUAL_STATUS_BASIC;
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
	for (i = 0; i < state->rows; i++)
		state->value[state->basis[i]] = state->work[i];
	return lp_simplex_EXIT_SUCCESS;
}


static int dual_choose_leaving(
		struct simplex_DualState *state,
		double *target, int *kappa, double *maximum)
{
	int i, p = -1;
	double largest = state->options->primal_tolerance;
	double best_score = 0.;
	for (i = 0; i < state->rows; i++) {
		int variable = state->basis[i];
		double violation = 0.;
		double score;
		if (state->value[variable] < state->lower[variable] -
		    state->options->primal_tolerance)
			violation = state->lower[variable] - state->value[variable];
		else if (state->value[variable] > state->upper[variable] +
			 state->options->primal_tolerance)
			violation = state->value[variable] - state->upper[variable];
		if (violation <= state->options->primal_tolerance)
			continue;
		score = violation * violation / state->edge_weight[i];
		if (p < 0 || score > best_score) {
			best_score = score;
			largest = violation;
			p = i;
		}
	}
	*maximum = largest;
	if (p < 0)
		return -1;
	if (state->value[state->basis[p]] < state->lower[state->basis[p]]) {
		*target = state->lower[state->basis[p]];
		*kappa = 1;
	} else {
		*target = state->upper[state->basis[p]];
		*kappa = -1;
	}
	return p;
}


static int dual_compute_edge_weights(struct simplex_DualState *state)
{
	int i, k;
	if (state->factor.update_count == 0)
		return simplex_basis_edge_weights(&state->factor, state->edge_weight);
	for (i = 0; i < state->rows; i++) {
		double weight = 0.;
		lp_simplex_memset(state->work, 0,
			(size_t)state->rows * sizeof(double));
		state->work[i] = 1.;
		if (simplex_basis_btran(&state->factor, state->work) ==
		    lp_simplex_EXIT_FAILURE)
			return lp_simplex_EXIT_FAILURE;
		for (k = 0; k < state->rows; k++)
			weight += state->work[k] * state->work[k];
		state->edge_weight[i] = __lp_simplex_MAX__(weight, 1e-12);
	}
	return lp_simplex_EXIT_SUCCESS;
}


static int dual_update_edge_weights(
		struct simplex_DualState *state, const int p)
{
	int i;
	double pivot = state->direction[p];
	double old_p = state->edge_weight[p];
	if (__lp_simplex_ABS__(pivot) <= state->options->pivot_tolerance)
		return lp_simplex_EXIT_FAILURE;
	lp_simplex_memcpy(state->work, state->rho,
		(size_t)state->rows * sizeof(double));
	if (simplex_basis_ftran(&state->factor, state->work) ==
	    lp_simplex_EXIT_FAILURE)
		return lp_simplex_EXIT_FAILURE;
	for (i = 0; i < state->rows; i++) {
		double ratio;
		double updated;
		if (i == p)
			continue;
		ratio = state->direction[i] / pivot;
		updated = state->edge_weight[i] - 2. * ratio * state->work[i] +
			ratio * ratio * old_p;
		state->edge_weight[i] = __lp_simplex_MAX__(updated, 1e-12);
	}
	state->edge_weight[p] = __lp_simplex_MAX__(old_p / (pivot * pivot), 1e-12);
	return lp_simplex_EXIT_SUCCESS;
}


static int dual_sign_for_candidate(
		const struct simplex_DualState *state,
		const int variable, const int kappa, const double alpha)
{
	if (state->status[variable] == DUAL_STATUS_LOWER)
		return 1;
	if (state->status[variable] == DUAL_STATUS_UPPER)
		return -1;
	if (state->status[variable] == DUAL_STATUS_FREE)
		return kappa * alpha > 0. ? -1 : 1;
	return 0;
}


static int dual_prepare_ratio_test(
		struct simplex_DualState *state, const int kappa)
{
	int candidates = 0, j;
	clock_t started = 0;
	if (state->profile_enabled)
		started = clock();
	for (j = 0; j < state->variables; j++) {
		int sign;
		double denominator, signed_reduced, theta;
		state->candidate_sign[j] = 0;
		state->breakpoint[j] = __lp_simplex_INF__;
		if (state->position[j] >= 0 || state->status[j] == DUAL_STATUS_FIXED)
			continue;
		state->alpha[j] = dual_column_dot(state, j, state->rho);
		sign = dual_sign_for_candidate(state, j, kappa, state->alpha[j]);
		state->candidate_sign[j] = sign;
		denominator = -kappa * sign * state->alpha[j];
		if (denominator <= state->options->pivot_tolerance)
			continue;
		signed_reduced = sign * state->reduced[j];
		if (signed_reduced < 0. &&
		    signed_reduced >= -state->options->dual_tolerance)
			signed_reduced = 0.;
		if (signed_reduced < 0.)
			continue;
		theta = signed_reduced / denominator;
		state->breakpoint[j] = theta;
		candidates++;
	}
	if (state->profile_enabled)
		state->profile_ratio_seconds +=
			(double)(clock() - started) / (double)CLOCKS_PER_SEC;
	return candidates;
}


static int dual_next_ratio_candidate(
		struct simplex_DualState *state, double *chosen_theta)
{
	int j, q = -1;
	double minimum = __lp_simplex_INF__;
	double best_pivot = 0.;
	for (j = 0; j < state->variables; j++) {
		if (state->candidate_sign[j] != 0 &&
		    state->breakpoint[j] < minimum)
			minimum = state->breakpoint[j];
	}
	if (minimum == __lp_simplex_INF__)
		return -1;
	for (j = 0; j < state->variables; j++) {
		double theta, relaxed;
		if (state->candidate_sign[j] == 0)
			continue;
		theta = state->breakpoint[j];
		/* Keep the exact minimum ratio as the dual-feasibility boundary.
		 * Pivot magnitude only breaks numerically indistinguishable ties. */
		relaxed = minimum + 1e-12 * (1. + __lp_simplex_ABS__(minimum));
		if (theta <= relaxed &&
		    __lp_simplex_ABS__(state->alpha[j]) > best_pivot) {
			best_pivot = __lp_simplex_ABS__(state->alpha[j]);
			q = j;
			*chosen_theta = theta;
		}
	}
	return q;
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
	int i, j;
	double tau = -kappa * theta;
	for (i = 0; i < state->rows; i++)
		state->pi[i] += tau * state->rho[i];
	for (j = 0; j < state->variables; j++) {
		if (state->position[j] < 0)
			state->reduced[j] -= tau * state->alpha[j];
	}
	state->reduced[entering] = 0.;
	state->reduced[leaving] = -tau;
}


static int dual_run(struct simplex_DualState *state, int *status)
{
	while (state->iterations < state->options->iteration_limit) {
		int i, p, q, kappa, sign, restart_iteration = 0;
		int bound_flips = 0;
		int rejected_relative = 0;
		double relative_pivot_tolerance = DUAL_RELATIVE_PIVOT_TOLERANCE;
		double target, primal_infeasibility, theta = 0.;
		double delta, new_value, movement = 0.;
		unsigned char flipped_status = DUAL_STATUS_FREE;
		p = dual_choose_leaving(state, &target, &kappa,
					&primal_infeasibility);
		if (p < 0) {
			if (simplex_basis_factorize(&state->factor) ==
			    lp_simplex_EXIT_FAILURE ||
			    dual_compute_primal_values(state) == lp_simplex_EXIT_FAILURE ||
			    dual_compute_reduced_costs(state) == lp_simplex_EXIT_FAILURE) {
				*status = lp_simplex_Singularity;
				return lp_simplex_EXIT_FAILURE;
			}
			p = dual_choose_leaving(state, &target, &kappa,
						&primal_infeasibility);
			if (p < 0 && dual_max_dual_infeasibility(state) <=
			    state->options->dual_tolerance) {
				*status = lp_simplex_Success;
				return lp_simplex_EXIT_SUCCESS;
			}
			if (p < 0) {
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
		lp_simplex_memset(state->alpha, 0,
			(size_t)state->variables * sizeof(double));
		lp_simplex_memset(state->candidate_sign, 0,
			(size_t)state->variables * sizeof(int));
		if (dual_prepare_ratio_test(state, kappa) == 0) {
			if (state->factor.update_count != 0 &&
			    simplex_basis_factorize(&state->factor) ==
			    lp_simplex_EXIT_SUCCESS &&
			    dual_compute_primal_values(state) ==
			    lp_simplex_EXIT_SUCCESS &&
			    dual_compute_reduced_costs(state) ==
			    lp_simplex_EXIT_SUCCESS)
				continue;
			*status = lp_simplex_Infeasibility;
			return lp_simplex_EXIT_FAILURE;
		}
		q = -1;
		for (;;) {
			double direction_maximum = 0.;
			q = dual_next_ratio_candidate(state, &theta);
			if (q < 0) {
				if (bound_flips > 0) {
					restart_iteration = 1;
					break;
				}
				if (rejected_relative > 0 &&
				    relative_pivot_tolerance > 1e-14) {
					relative_pivot_tolerance *= 1e-4;
					rejected_relative = 0;
					dual_prepare_ratio_test(state, kappa);
					continue;
				}
				if (state->factor.update_count != 0 &&
				    simplex_basis_factorize(&state->factor) ==
				    lp_simplex_EXIT_SUCCESS &&
				    dual_compute_primal_values(state) ==
				    lp_simplex_EXIT_SUCCESS &&
				    dual_compute_reduced_costs(state) ==
				    lp_simplex_EXIT_SUCCESS) {
					restart_iteration = 1;
					break;
				}
				*status = lp_simplex_Infeasibility;
				return lp_simplex_EXIT_FAILURE;
			}
			simplex_csc_column_to_dense(&state->matrix, state->structural,
						    q, state->direction);
			if (simplex_basis_ftran(&state->factor, state->direction) ==
			    lp_simplex_EXIT_FAILURE) {
				state->candidate_sign[q] = 0;
				continue;
			}
			for (i = 0; i < state->rows; i++)
				direction_maximum = __lp_simplex_MAX__(direction_maximum,
					__lp_simplex_ABS__(state->direction[i]));
			if (__lp_simplex_ABS__(state->direction[p]) <=
			    state->options->pivot_tolerance ||
			    __lp_simplex_ABS__(state->direction[p]) <
			    relative_pivot_tolerance * direction_maximum) {
				state->candidate_sign[q] = 0;
				rejected_relative++;
				continue;
			}
			if (__lp_simplex_ABS__(state->direction[p] - state->alpha[q]) >
			    100. * state->options->pivot_tolerance *
			    __lp_simplex_MAX__(1.,
				__lp_simplex_ABS__(state->alpha[q]))) {
				if (state->factor.update_count != 0 &&
				    simplex_basis_factorize(&state->factor) ==
				    lp_simplex_EXIT_SUCCESS &&
				    dual_compute_primal_values(state) ==
				    lp_simplex_EXIT_SUCCESS &&
				    dual_compute_reduced_costs(state) ==
				    lp_simplex_EXIT_SUCCESS) {
					restart_iteration = 1;
					break;
				}
				*status = lp_simplex_PrecisionError;
				return lp_simplex_EXIT_FAILURE;
			}
			delta = (state->value[state->basis[p]] - target) /
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
			for (i = 0; i < state->rows; i++)
				state->value[state->basis[i]] -=
					state->direction[i] * movement;
			state->value[q] += movement;
			state->status[q] = flipped_status;
			state->candidate_sign[q] = 0;
			bound_flips++;
			movement = 0.;
		}
		if (restart_iteration)
			continue;
		for (i = 0; i < state->rows; i++) {
			if (i != p)
				state->value[state->basis[i]] -=
					state->direction[i] * delta;
		}
		if (dual_update_edge_weights(state, p) == lp_simplex_EXIT_FAILURE) {
			*status = lp_simplex_PrecisionError;
			return lp_simplex_EXIT_FAILURE;
		}
		{
			int leaving = state->basis[p];
			dual_update_dual_values(state, kappa, theta, leaving, q);
			state->value[leaving] = target;
			state->status[leaving] = kappa > 0
				? DUAL_STATUS_LOWER : DUAL_STATUS_UPPER;
			if (state->lower[leaving] == state->upper[leaving])
				state->status[leaving] = DUAL_STATUS_FIXED;
			state->position[leaving] = -1;
			state->basis[p] = q;
			state->position[q] = p;
			state->status[q] = DUAL_STATUS_BASIC;
			state->value[q] = new_value;
		}
		state->iterations++;
		{
			int update = simplex_basis_update(&state->factor, p,
						  state->direction);
			if (update == 1) {
				update = simplex_basis_factorize(&state->factor);
				state->reinversions++;
				if (update == lp_simplex_EXIT_SUCCESS)
					update = dual_compute_primal_values(state);
				if (update == lp_simplex_EXIT_SUCCESS)
					update = dual_compute_reduced_costs(state);
			}
			if (update == lp_simplex_EXIT_FAILURE) {
				*status = lp_simplex_Singularity;
				return lp_simplex_EXIT_FAILURE;
			}
		}
		if (dual_max_dual_infeasibility(state) >
		    10. * state->options->dual_tolerance) {
			if (simplex_basis_factorize(&state->factor) ==
			    lp_simplex_EXIT_FAILURE ||
			    dual_compute_primal_values(state) == lp_simplex_EXIT_FAILURE ||
			    dual_compute_reduced_costs(state) == lp_simplex_EXIT_FAILURE) {
				*status = lp_simplex_PrecisionError;
				return lp_simplex_EXIT_FAILURE;
			}
			if (dual_max_dual_infeasibility(state) >
			    10. * state->options->dual_tolerance) {
				if (dual_crash(state) == lp_simplex_EXIT_FAILURE ||
				    simplex_basis_factorize(&state->factor) ==
				    lp_simplex_EXIT_FAILURE ||
				    dual_compute_primal_values(state) ==
				    lp_simplex_EXIT_FAILURE ||
				    dual_compute_reduced_costs(state) ==
				    lp_simplex_EXIT_FAILURE ||
				    dual_max_dual_infeasibility(state) >
				    10. * state->options->dual_tolerance) {
					*status = lp_simplex_PrecisionError;
					return lp_simplex_EXIT_FAILURE;
				}
				for (i = 0; i < state->rows; i++)
					state->edge_weight[i] = 1.;
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


int simplex_dual_solve(
		const struct lp_Model *model,
		const struct lp_simplex_Options *options,
		double *x, double *row_dual, struct lp_simplex_Result *result)
{
	struct simplex_DualState state;
	int j, status = lp_simplex_CondUnsatisfied;
	int solve_state;
	double objective = 0.;
	clock_t profile_started = 0;
	if (getenv("LP_SIMPLEX_PROFILE") != NULL)
		profile_started = clock();
	if (dual_allocate(&state, model, options) == lp_simplex_EXIT_FAILURE) {
		dual_destroy(&state);
		result->status = lp_simplex_MemoryAllocError;
		return lp_simplex_EXIT_FAILURE;
	}
	dual_set_column_bounds(&state, model);
	dual_initialize_basis(&state, model);
	if (simplex_basis_factorize(&state.factor) == lp_simplex_EXIT_FAILURE ||
	    dual_crash(&state) == lp_simplex_EXIT_FAILURE ||
	    simplex_basis_factorize(&state.factor) == lp_simplex_EXIT_FAILURE ||
	    dual_compute_primal_values(&state) == lp_simplex_EXIT_FAILURE ||
	    dual_compute_reduced_costs(&state) == lp_simplex_EXIT_FAILURE ||
	    dual_compute_edge_weights(&state) == lp_simplex_EXIT_FAILURE) {
		result->status = lp_simplex_PrecisionError;
		dual_destroy(&state);
		return lp_simplex_EXIT_FAILURE;
	}
	solve_state = dual_run(&state, &status);
	for (j = 0; j < model->n; j++) {
		x[j] = state.value[j];
		objective += model->objective[j] * x[j];
	}
	result->status = status;
	result->iterations = state.iterations;
	result->objective = objective;
	result->primal_infeasibility = dual_primal_infeasibility(&state);
	result->dual_infeasibility = dual_max_dual_infeasibility(&state);
	if (row_dual != NULL)
		for (j = 0; j < model->m; j++)
			row_dual[j] = state.pi[j];
	if (state.profile_enabled) {
		state.profile_total_seconds =
			(double)(clock() - profile_started) / (double)CLOCKS_PER_SEC;
		fprintf(stderr,
			"dual profile: total=%.6f factor=%.6f/%ld ftran=%.6f/%ld "
			"btran=%.6f/%ld ratio=%.6f residual=%.6f\n",
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
			state.profile_ratio_seconds);
	}
	dual_destroy(&state);
	return solve_state;
}
