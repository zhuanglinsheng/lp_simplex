/**
 * @file simplex_model_audit.c
 * @brief Original public-model audits, independent of algorithm and
 * coordinates.
 */

#include "simplex_model_audit.h"
#include "utils.h"
#include <math.h>

double simplex_primal_infeasibility(
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

double simplex_original_objective(const struct lp_Model *model, const double *x)
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

int simplex_prefer_point(const struct lp_Model *model,
	const struct lp_simplex_Options *o, const double *candidate,
	const double *saved)
{
	double primal = simplex_primal_infeasibility(model, candidate);
	double cost = model->objective_sense *
		      simplex_original_objective(model, candidate);
	double saved_primal, saved_cost;

	if (!isfinite(primal) || primal > o->primal_tolerance ||
		!isfinite(cost))
		return 0;
	saved_primal = simplex_primal_infeasibility(model, saved);
	saved_cost = model->objective_sense *
		     simplex_original_objective(model, saved);
	return !isfinite(saved_primal) || saved_primal > o->primal_tolerance ||
	       !isfinite(saved_cost) || cost < saved_cost;
}

/**
 * Original-model gate in front of every final public result.
 *
 * A solver may only claim `lp_simplex_Success` when the vector it hands back
 * is finite and feasible for the caller's own model.  Without this check the
 * presolve path can report success while its reduced problem overflowed
 * (objective inf/NaN) and the tableau path reports success while its
 * recovered vector violates the constraints outright.
 */
int simplex_finalize_result(const struct lp_Model *model,
	const struct lp_simplex_Options *options, const double *x,
	struct lp_simplex_Result *result, int state)
{
	result->primal_infeasibility = simplex_primal_infeasibility(model, x);

	result->objective = simplex_original_objective(model, x);
	if (result->dual_residual_space == lp_simplex_RESIDUAL_UNAVAILABLE)
		result->dual_infeasibility = NAN;
	/* A tolerance-feasible counterexample invalidates an infeasibility
	 * claim independently of objective overflow or dual diagnostics. */
	if (result->status == lp_simplex_Infeasibility &&
		isfinite(result->primal_infeasibility) &&
		result->primal_infeasibility <= options->primal_tolerance)
		result->status = lp_simplex_PrecisionError;
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
