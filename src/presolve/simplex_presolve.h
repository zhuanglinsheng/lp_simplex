/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_presolve.h
 * @brief Presolve reductions, owned mapping metadata and algebraic point lift.
 */

#ifndef LP_SIMPLEX_PRESOLVE_INTERNAL_H
#define LP_SIMPLEX_PRESOLVE_INTERNAL_H

#include "simplex_problem.h"

#include <lp_simplex/model.h>

struct simplex_PresolveSubstitution {
	int column;
	double constant;
	int term_start;
	int term_count;
};

struct simplex_PresolveJournal {
	struct simplex_PresolveSubstitution *record;
	int count;
	int capacity;
	int *term_column;
	double *term_multiplier;
	int term_count;
	int term_capacity;
};

struct simplex_PresolveStats {
	int removed_columns;
	int removed_rows;
	int fixed_columns;
	int empty_columns;
	int singleton_columns;
	int singleton_inequality_columns;
	int singleton_rows;
	int redundant_rows;
	int duplicate_rows;
	int dependent_rows;
	int dominated_rows;
	int dominated_columns;
	int forcing_rows;
	int forced_columns;
	int doubleton_rows;
	int singleton_column_rows;
	int implied_free_columns;
	int singleton_projection_columns;
	int singleton_column_candidates;
	int free_singleton_columns;
	int equality_singleton_columns;
	int outward_bound_propagation;
	int tightened_bounds;
	int queue_epochs;
};

struct simplex_Presolve {
	const struct lp_Model *original;
	struct simplex_Problem *reduced;
	/* Owned until destroy. column_map maps live reduced columns only while
	 * status is CondUnsatisfied and reduced is present; terminal reductions
	 * release it. row_map is retained for original-model progress bounds. */
	int *column_map;
	int *row_map;
	double *eliminated_value;
	unsigned char *eliminated;
	struct simplex_PresolveJournal journal;
	struct simplex_PresolveStats stats;
	int status;
};

/* Exact dependency reduction remains opt-in pending performance acceptance. */
int simplex_presolve_run_configured(struct simplex_Presolve *presolve,
	const struct lp_Model *model, double tolerance, int dependencies);
int simplex_presolve_run(struct simplex_Presolve *presolve,
	const struct lp_Model *model, double tolerance);

/* Read-only algebraic lift, shared by final output and sampled progress.
 * Replays the complete substitution journal on each call; does not audit,
 * solve, consume iteration budget or destroy metadata. original_x must hold
 * original->n values and must not alias reduced_x. Terminal reductions may
 * pass NULL reduced_x. The owner must outlive every progress bridge. */
void simplex_presolve_postsolve(const struct simplex_Presolve *presolve,
	const double *reduced_x, double *original_x);

void simplex_presolve_destroy(struct simplex_Presolve *presolve);

void simplex_presolve_print_profile(
	const struct simplex_Presolve *presolve, double seconds);

#endif
