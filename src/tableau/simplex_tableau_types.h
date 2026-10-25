/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_tableau_types.h
 * @brief Tableau representation and variable transformation state.
 */

#ifndef LP_SIMPLEX_TABLEAU_TYPES_INTERNAL_H
#define LP_SIMPLEX_TABLEAU_TYPES_INTERNAL_H

/* Tableau-local phase identities, not Dual phases or additional algorithms. */
enum simplex_TableauPhase {
	SIMPLEX_TABLEAU_INPUT = 0,
	SIMPLEX_TABLEAU_PHASE_ONE,
	SIMPLEX_TABLEAU_ARTIFICIAL_CLEANUP,
	SIMPLEX_TABLEAU_PHASE_TWO
};

/* Raw invocation diagnostics. Phase work is disjoint and sums to iterations.
 * Point availability means a fully initialized output, not an input audit.
 * Tableau exports no reusable basis, factor, multiplier or lower bound. */
struct simplex_TableauReport {
	int available, phase, phase_one_work, cleanup_work, phase_two_work;
	int point_available;
	/* Model ids share the caller's live solve context; zero is unknown.
	 * A bounds conversion is a representation, never another algorithm. */
	int input_model, standard_model, standard_rows, standard_columns;
	int converted, recovered_point;
};
#endif
