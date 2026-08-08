/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
/* Row activity classification and implied-bound propagation. */
#include "simplex_presolve_activity.h"
#include "simplex_presolve_rules.h"
#include "utils.h"

#include <math.h>


static int presolve_apply_implied_lower(
		struct optm_VariableBound *bound, double candidate,
		const double tolerance)
{
	double margin;
	if (!isfinite(candidate))
		return 0;
	margin = tolerance * (1. + __lp_simplex_ABS__(candidate));
	candidate -= margin;
	if (simplex_presolve_has_upper(bound) && candidate > bound->ub + margin)
		return -1;
	if (simplex_presolve_has_upper(bound) && candidate > bound->ub)
		candidate = bound->ub;
	if (simplex_presolve_has_lower(bound) &&
	    candidate <= bound->lb + margin)
		return 0;
	simplex_presolve_set_lower(bound, candidate);
	return 1;
}


static int presolve_apply_implied_upper(
		struct optm_VariableBound *bound, double candidate,
		const double tolerance)
{
	double margin;
	if (!isfinite(candidate))
		return 0;
	margin = tolerance * (1. + __lp_simplex_ABS__(candidate));
	candidate += margin;
	if (simplex_presolve_has_lower(bound) && candidate < bound->lb - margin)
		return -1;
	if (simplex_presolve_has_lower(bound) && candidate < bound->lb)
		candidate = bound->lb;
	if (simplex_presolve_has_upper(bound) &&
	    candidate >= bound->ub - margin)
		return 0;
	simplex_presolve_set_upper(bound, candidate);
	return 1;
}


/* Analyze one sparse row and propagate bounds from its residual activity.
 * Return -1 for infeasible, 1 for redundant, and 0 when the row must stay. */
struct presolve_RowActivity {
	long double minimum;
	long double maximum;
	long double magnitude;
	int minimum_infinite;
	int maximum_infinite;
	int start;
	int end;
	int type;
};


static double presolve_row_entry(
		const struct lp_Model *model, const int row,
		const int position, int *column)
{
	if (model->row_start != NULL) {
		*column = model->column_index[position];
		return model->row_value[position];
	}
	*column = position;
	return model->constraints[row].coef[position];
}


static void presolve_compute_row_activity(
		const struct lp_Model *model, const int row,
		const unsigned char *column_active,
		const struct optm_VariableBound *bounds,
		struct presolve_RowActivity *activity)
{
	int column, position;
	lp_simplex_memset(activity, 0, sizeof(*activity));
	activity->type = model->constraints[row].type;
	activity->start = model->row_start != NULL ? model->row_start[row] : 0;
	activity->end = model->row_start != NULL
		? model->row_start[row + 1] : model->n;
	for (position = activity->start; position < activity->end; position++) {
		double coefficient = presolve_row_entry(
			model, row, position, &column);
		const struct optm_VariableBound *bound;
		double lower, upper;
		if (!column_active[column] || coefficient == 0.)
			continue;
		bound = bounds + column;
		if (coefficient > 0.) {
			if (simplex_presolve_has_lower(bound)) {
				lower = coefficient * bound->lb;
				activity->minimum += lower;
				activity->magnitude += fabsl((long double)lower);
			} else activity->minimum_infinite++;
			if (simplex_presolve_has_upper(bound)) {
				upper = coefficient * bound->ub;
				activity->maximum += upper;
				activity->magnitude += fabsl((long double)upper);
			} else activity->maximum_infinite++;
		} else {
			if (simplex_presolve_has_upper(bound)) {
				lower = coefficient * bound->ub;
				activity->minimum += lower;
				activity->magnitude += fabsl((long double)lower);
			} else activity->minimum_infinite++;
			if (simplex_presolve_has_lower(bound)) {
				upper = coefficient * bound->lb;
				activity->maximum += upper;
				activity->magnitude += fabsl((long double)upper);
			} else activity->maximum_infinite++;
		}
	}
}


static int presolve_classify_row_activity(
		const struct presolve_RowActivity *activity,
		const double rhs, const double tolerance,
		int *force_minimum, int *force_maximum)
{
	long double margin = (long double)tolerance *
		(1. + fabsl((long double)rhs) + activity->magnitude);
	*force_minimum = 0;
	*force_maximum = 0;
	if ((activity->type == optm_CONS_T_GE ||
	     activity->type == optm_CONS_T_EQ) &&
	    activity->maximum_infinite == 0 &&
	    activity->maximum < (long double)rhs - margin)
		return -1;
	if ((activity->type == optm_CONS_T_LE ||
	     activity->type == optm_CONS_T_EQ) &&
	    activity->minimum_infinite == 0 &&
	    activity->minimum > (long double)rhs + margin)
		return -1;
	if (activity->type == optm_CONS_T_GE &&
	    activity->minimum_infinite == 0 &&
	    activity->minimum >= (long double)rhs + margin)
		return 1;
	if (activity->type == optm_CONS_T_LE &&
	    activity->maximum_infinite == 0 &&
	    activity->maximum <= (long double)rhs - margin)
		return 1;
	*force_minimum =
		(activity->type == optm_CONS_T_LE ||
		 activity->type == optm_CONS_T_EQ) &&
		activity->minimum_infinite == 0 &&
		activity->minimum == (long double)rhs;
	*force_maximum =
		(activity->type == optm_CONS_T_GE ||
		 activity->type == optm_CONS_T_EQ) &&
		activity->maximum_infinite == 0 &&
		activity->maximum == (long double)rhs;
	return 0;
}


static int presolve_force_row_bounds(
		const struct lp_Model *model, const int row,
		const unsigned char *column_active,
		struct optm_VariableBound *bounds,
		const struct presolve_RowActivity *activity,
		const int force_minimum)
{
	int column, position, tightened = 0;
	for (position = activity->start; position < activity->end; position++) {
		double coefficient = presolve_row_entry(
			model, row, position, &column);
		struct optm_VariableBound *bound;
		double value;
		if (!column_active[column] || coefficient == 0.)
			continue;
		bound = bounds + column;
		if (bound->b_type == optm_BOUND_T_BS && bound->lb == bound->ub)
			continue;
		if (force_minimum)
			value = coefficient > 0. ? bound->lb : bound->ub;
		else
			value = coefficient > 0. ? bound->ub : bound->lb;
		simplex_presolve_set_lower(bound, value);
		simplex_presolve_set_upper(bound, value);
		tightened++;
	}
	return tightened;
}


static int presolve_propagate_row_bounds(
		const struct lp_Model *model, const int row,
		const unsigned char *column_active,
		struct optm_VariableBound *bounds,
		const struct presolve_RowActivity *activity,
		const double rhs, const double propagation_tolerance,
		int *tightened)
{
	int column, position;
	for (position = activity->start; position < activity->end; position++) {
		double coefficient = presolve_row_entry(
			model, row, position, &column);
		struct optm_VariableBound *bound;
		int own_minimum_infinite, own_maximum_infinite, state;
		long double own_minimum = 0., own_maximum = 0., residual;
		if (!column_active[column] || coefficient == 0.)
			continue;
		bound = bounds + column;
		own_minimum_infinite = coefficient > 0.
			? !simplex_presolve_has_lower(bound)
			: !simplex_presolve_has_upper(bound);
		own_maximum_infinite = coefficient > 0.
			? !simplex_presolve_has_upper(bound)
			: !simplex_presolve_has_lower(bound);
		if (!own_minimum_infinite)
			own_minimum = (long double)coefficient *
				(coefficient > 0. ? bound->lb : bound->ub);
		if (!own_maximum_infinite)
			own_maximum = (long double)coefficient *
				(coefficient > 0. ? bound->ub : bound->lb);
		if ((activity->type == optm_CONS_T_GE ||
		     activity->type == optm_CONS_T_EQ) &&
		    (activity->maximum_infinite == 0 ||
		     (activity->maximum_infinite == 1 && own_maximum_infinite))) {
			residual = activity->maximum - own_maximum;
			if (coefficient > 0.)
				state = presolve_apply_implied_lower(bound,
					(double)(((long double)rhs - residual) / coefficient),
					propagation_tolerance);
			else
				state = presolve_apply_implied_upper(bound,
					(double)(((long double)rhs - residual) / coefficient),
					propagation_tolerance);
			if (state < 0)
				return -1;
			*tightened += state;
		}
		if ((activity->type == optm_CONS_T_LE ||
		     activity->type == optm_CONS_T_EQ) &&
		    (activity->minimum_infinite == 0 ||
		     (activity->minimum_infinite == 1 && own_minimum_infinite))) {
			residual = activity->minimum - own_minimum;
			if (coefficient > 0.)
				state = presolve_apply_implied_upper(bound,
					(double)(((long double)rhs - residual) / coefficient),
					propagation_tolerance);
			else
				state = presolve_apply_implied_lower(bound,
					(double)(((long double)rhs - residual) / coefficient),
					propagation_tolerance);
			if (state < 0)
				return -1;
			*tightened += state;
		}
	}
	return 0;
}


/* Analyze one row and propagate bounds from its residual activity.
 * Return -1 for infeasible, 1 for redundant, and 0 when the row must stay. */
int simplex_presolve_analyze_row(
		const struct lp_Model *model, const int row,
		const unsigned char *column_active,
		struct optm_VariableBound *bounds,
		const double rhs, const double tolerance,
		const double propagation_tolerance, const int propagate,
		int *tightened, int *forcing)
{
	struct presolve_RowActivity activity;
	int classification, force_minimum, force_maximum;
	*tightened = 0;
	*forcing = 0;
	presolve_compute_row_activity(
		model, row, column_active, bounds, &activity);
	classification = presolve_classify_row_activity(
		&activity, rhs, tolerance, &force_minimum, &force_maximum);
	if (classification != 0)
		return classification;
	if (force_minimum || force_maximum) {
		*tightened = presolve_force_row_bounds(model, row, column_active,
			bounds, &activity, force_minimum);
		*forcing = *tightened != 0;
		return 0;
	}
	return propagate ? presolve_propagate_row_bounds(model, row,
		column_active, bounds, &activity, rhs, propagation_tolerance,
		tightened) : 0;
}
