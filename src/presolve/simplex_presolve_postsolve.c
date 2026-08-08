/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#include "simplex_presolve.h"
#include "utils.h"


void simplex_presolve_postsolve(
		const struct simplex_Presolve *presolve,
		const double *reduced_x, double *original_x)
{
	int j;
	for (j = 0; j < presolve->original->n; j++)
		original_x[j] = presolve->eliminated[j]
			? presolve->eliminated_value[j] : 0.;
	if (reduced_x != NULL && presolve->reduced != NULL)
		for (j = 0; j < presolve->reduced->columns; j++)
			original_x[presolve->column_map[j]] = reduced_x[j];
	/* Substitutions may depend on columns eliminated later.  Replaying the
	 * journal backwards restores those dependencies before their consumers. */
	for (j = presolve->journal.count - 1; j >= 0; j--) {
		const struct simplex_PresolveSubstitution *record =
			presolve->journal.record + j;
		double value = record->constant;
		int k;
		for (k = 0; k < record->term_count; k++) {
			int position = record->term_start + k;
			value += presolve->journal.term_multiplier[position] *
				original_x[presolve->journal.term_column[position]];
		}
		original_x[record->column] = value;
	}
}


void simplex_presolve_destroy(struct simplex_Presolve *presolve)
{
	if (presolve == NULL)
		return;
	simplex_problem_free(presolve->reduced);
	lp_simplex_free(presolve->column_map);
	lp_simplex_free(presolve->row_map);
	lp_simplex_free(presolve->eliminated_value);
	lp_simplex_free(presolve->eliminated);
	lp_simplex_free(presolve->journal.record);
	lp_simplex_free(presolve->journal.term_column);
	lp_simplex_free(presolve->journal.term_multiplier);
	presolve->reduced = NULL;
	presolve->column_map = NULL;
	presolve->row_map = NULL;
	presolve->eliminated_value = NULL;
	presolve->eliminated = NULL;
	presolve->journal.record = NULL;
	presolve->journal.term_column = NULL;
	presolve->journal.term_multiplier = NULL;
	presolve->journal.count = 0;
	presolve->journal.capacity = 0;
	presolve->journal.term_count = 0;
	presolve->journal.term_capacity = 0;
}
