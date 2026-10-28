/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_presolve_postsolve.c
 * @brief Recover original variables after presolve reductions.
 */

#include "simplex_presolve.h"
#include "utils.h"

#include <math.h>

void simplex_presolve_postsolve(const struct simplex_Presolve *presolve,
	const double *reduced_x, double *original_x)
{
	int j;

	for (j = 0; j < presolve->original->n; j++)
		original_x[j] = presolve->eliminated[j]
					? presolve->eliminated_value[j]
					: 0.;
	if (reduced_x != NULL && presolve->reduced != NULL)
		for (j = 0; j < presolve->reduced->columns; j++)
			original_x[presolve->column_map[j]] = reduced_x[j];
	/* Substitutions may depend on columns eliminated later.  Replaying the
	 * journal backwards restores those dependencies before their consumers. */
	for (j = presolve->journal.count - 1; j >= 0; j--) {
		const struct simplex_PresolveSubstitution *record =
			presolve->journal.record + j;
		double value = record->constant, correction = 0.;
		int k;

		for (k = 0; k < record->term_count; k++) {
			int position = record->term_start + k;

			double multiplier =
				presolve->journal.term_multiplier[position];
			double coordinate = original_x[presolve->journal
					.term_column[position]];
			double product = multiplier * coordinate;
			double next = value + product;

			/* A short aggregate can still cancel a large RHS. Keep both
			 * addition and multiplication errors while replaying the exact
			 * same reverse dependencies; no coordinate is clipped or dropped. */
			correction +=
				(fabs(value) >= fabs(product)
						? (value - next) + product
						: (product - next) + value) +
				fma(multiplier, coordinate, -product);
			value = next;
		}
		original_x[record->column] = value + correction;
	}
}
