/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#ifndef LP_SIMPLEX_PRESOLVE_DOMINANCE_H
#define LP_SIMPLEX_PRESOLVE_DOMINANCE_H
#include <lp_simplex/model.h>
#include <stddef.h>

/* Optional exact-comparison reductions. Rows are deactivated immediately;
 * columns are only fixed, then removed/journalled by the ordinary queue.
 * work is shared across all fixed-point passes; exhaustion loses reductions,
 * never correctness. Allocation failure also leaves the ordinary rules usable. */
void simplex_presolve_dominance(const struct lp_Model *model,
	struct optm_VariableBound *bounds, const double *rhs,
	unsigned char *row_active, const unsigned char *column_active,
	const int *row_degree, int *column_degree, size_t *work,
	int *removed_rows, int *fixed_columns);
#endif
