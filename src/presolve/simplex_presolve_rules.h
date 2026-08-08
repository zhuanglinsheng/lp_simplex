/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#ifndef LP_SIMPLEX_PRESOLVE_RULES_INTERNAL_H
#define LP_SIMPLEX_PRESOLVE_RULES_INTERNAL_H

#include <lp_simplex/model.h>


int simplex_presolve_has_lower(const struct optm_VariableBound *bound);

int simplex_presolve_has_upper(const struct optm_VariableBound *bound);

void simplex_presolve_set_lower(
		struct optm_VariableBound *bound, double value);

void simplex_presolve_set_upper(
		struct optm_VariableBound *bound, double value);

int simplex_presolve_empty_row_feasible(
		int row_type, double rhs, double tolerance);

int simplex_presolve_choose_empty_column(
		double cost, const struct optm_VariableBound *bound, double *value);

int simplex_presolve_tighten_singleton(
		struct optm_VariableBound *bound, int row_type,
		double coefficient, double rhs, double tolerance, int *tightened);

#endif
