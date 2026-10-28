/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_presolve_dependencies.h
 * @brief Optional exact signed-unit conservation dependencies.
 */
#ifndef LP_SIMPLEX_PRESOLVE_DEPENDENCIES_H
#define LP_SIMPLEX_PRESOLVE_DEPENDENCIES_H
#include <lp_simplex/model.h>
#include <stddef.h>

/* Exact integer witnesses precede modular independence selection. All accepted
 * witnesses refer to the same unchanged input; publish deletions together.
 * Resource or arithmetic refusal loses optional reductions, never classifies
 * the LP. Work counts index visits and modular multiply-adds, not solve work. */
int simplex_presolve_dependencies(const struct lp_Model *model,
	const double *rhs, unsigned char *row_active,
	const unsigned char *column_active, int *column_degree,
	size_t work_limit, size_t byte_limit);
#endif
