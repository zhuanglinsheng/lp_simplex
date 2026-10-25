/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_presolve_record.h
 * @brief Record reversible presolve operations for postsolve.
 */

#ifndef LP_SIMPLEX_PRESOLVE_RECORD_INTERNAL_H
#define LP_SIMPLEX_PRESOLVE_RECORD_INTERNAL_H

#include "simplex_presolve.h"

/* Reserve one expression before changing its model. Failure retains all
 * completed expressions and term counts; physical roots/capacities may grow.
 * A nonempty expression becomes visible to postsolve only after its final
 * add_term. Fill exactly term_count terms before beginning another record.
 * A zero-term expression is complete on return. Adding terms cannot allocate. */
struct simplex_PresolveSubstitution *simplex_presolve_record_begin(
	struct simplex_Presolve *presolve, int column, double constant,
	int term_count);

void simplex_presolve_record_add_term(
	struct simplex_Presolve *presolve, int column, double multiplier);

#endif
