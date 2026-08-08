/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#ifndef LP_SIMPLEX_PRESOLVE_ACTIVITY_INTERNAL_H
#define LP_SIMPLEX_PRESOLVE_ACTIVITY_INTERNAL_H

#include <lp_simplex/model.h>


/* Return -1 for infeasible, 1 for redundant, and 0 when the row must stay. */
int simplex_presolve_analyze_row(
		const struct lp_Model *model, int row,
		const unsigned char *column_active,
		struct optm_VariableBound *bounds,
		double rhs, double tolerance,
		double propagation_tolerance, int propagate,
		int *tightened, int *forcing);

#endif
