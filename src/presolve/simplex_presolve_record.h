/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#ifndef LP_SIMPLEX_PRESOLVE_RECORD_INTERNAL_H
#define LP_SIMPLEX_PRESOLVE_RECORD_INTERNAL_H

#include "simplex_presolve.h"


struct simplex_PresolveSubstitution *simplex_presolve_record_begin(
		struct simplex_Presolve *presolve, int column,
		double constant, int term_count);

void simplex_presolve_record_add_term(
		struct simplex_Presolve *presolve, int column, double multiplier);

#endif
