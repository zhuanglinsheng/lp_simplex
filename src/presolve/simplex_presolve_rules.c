/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#include "simplex_presolve_rules.h"
#include "utils.h"

#include <lp_simplex/status.h>


int simplex_presolve_has_lower(const struct optm_VariableBound *bound)
{
	return bound->b_type == optm_BOUND_T_LO ||
		bound->b_type == optm_BOUND_T_BS;
}


int simplex_presolve_has_upper(const struct optm_VariableBound *bound)
{
	return bound->b_type == optm_BOUND_T_UP ||
		bound->b_type == optm_BOUND_T_BS;
}


void simplex_presolve_set_lower(
		struct optm_VariableBound *bound, const double value)
{
	bound->lb = value;
	bound->b_type = simplex_presolve_has_upper(bound)
		? optm_BOUND_T_BS : optm_BOUND_T_LO;
}


void simplex_presolve_set_upper(
		struct optm_VariableBound *bound, const double value)
{
	bound->ub = value;
	bound->b_type = simplex_presolve_has_lower(bound)
		? optm_BOUND_T_BS : optm_BOUND_T_UP;
}


int simplex_presolve_empty_row_feasible(
		const int row_type, const double rhs, const double tolerance)
{
	if (row_type == optm_CONS_T_EQ)
		return __lp_simplex_ABS__(rhs) <= tolerance;
	if (row_type == optm_CONS_T_GE)
		return rhs <= tolerance;
	return rhs >= -tolerance;
}


int simplex_presolve_choose_empty_column(
		const double cost, const struct optm_VariableBound *bound,
		double *value)
{
	if (cost > 0.) {
		if (!simplex_presolve_has_lower(bound))
			return lp_simplex_Unboundedness;
		*value = bound->lb;
	} else if (cost < 0.) {
		if (!simplex_presolve_has_upper(bound))
			return lp_simplex_Unboundedness;
		*value = bound->ub;
	} else if (simplex_presolve_has_lower(bound) && bound->lb > 0.) {
		*value = bound->lb;
	} else if (simplex_presolve_has_upper(bound) && bound->ub < 0.) {
		*value = bound->ub;
	} else {
		*value = 0.;
	}
	return lp_simplex_Success;
}


int simplex_presolve_tighten_singleton(
		struct optm_VariableBound *bound, const int row_type,
		const double coefficient, const double rhs,
		const double tolerance, int *tightened)
{
	double implied = rhs / coefficient;
	double margin = tolerance * (1. + __lp_simplex_ABS__(implied));
	*tightened = 0;
	if (row_type == optm_CONS_T_EQ) {
		if ((simplex_presolve_has_lower(bound) &&
		     implied < bound->lb - margin) ||
		    (simplex_presolve_has_upper(bound) &&
		     implied > bound->ub + margin))
			return lp_simplex_Infeasibility;
		if (!(bound->b_type == optm_BOUND_T_BS &&
		      bound->lb == implied && bound->ub == implied))
			*tightened = 1;
		simplex_presolve_set_lower(bound, implied);
		simplex_presolve_set_upper(bound, implied);
		return lp_simplex_Success;
	}
	if ((row_type == optm_CONS_T_GE && coefficient > 0.) ||
	    (row_type == optm_CONS_T_LE && coefficient < 0.)) {
		if (simplex_presolve_has_upper(bound) &&
		    implied > bound->ub + margin)
			return lp_simplex_Infeasibility;
		if (!simplex_presolve_has_lower(bound) || implied > bound->lb) {
			simplex_presolve_set_lower(bound, implied);
			*tightened = 1;
		}
	} else {
		if (simplex_presolve_has_lower(bound) &&
		    implied < bound->lb - margin)
			return lp_simplex_Infeasibility;
		if (!simplex_presolve_has_upper(bound) || implied < bound->ub) {
			simplex_presolve_set_upper(bound, implied);
			*tightened = 1;
		}
	}
	return lp_simplex_Success;
}
