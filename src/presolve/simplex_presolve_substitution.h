/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_presolve_substitution.h
 * @brief Sparse substitution and dependent-row presolve reductions.
 */

#ifndef LP_SIMPLEX_PRESOLVE_SUBSTITUTION_INTERNAL_H
#define LP_SIMPLEX_PRESOLVE_SUBSTITUTION_INTERNAL_H

struct simplex_Presolve;

int simplex_presolve_substitute_doubletons(
	struct simplex_Presolve *presolve, double tolerance);

#endif
