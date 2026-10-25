/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_presolve_record.c
 * @brief Record reversible presolve operations for postsolve.
 */

#include "simplex_presolve_record.h"
#include "simplex_presolve_growth.h"
#include "utils.h"

#include <lp_simplex/status.h>

struct simplex_PresolveSubstitution *simplex_presolve_record_begin(
	struct simplex_Presolve *presolve, const int column,
	const double constant, const int term_count)
{
	struct simplex_PresolveJournal *journal = &presolve->journal;
	struct simplex_PresolveSubstitution *record;
	int record_capacity, term_capacity;
	size_t record_bytes, term_bytes, column_bytes;
	size_t width =
		sizeof(int) > sizeof(double) ? sizeof(int) : sizeof(double);

	/* Validate both requests before allocating or publishing either record. */
	if (!simplex_presolve_growth(journal->capacity, journal->count, 1, 32,
		    sizeof(*record), (size_t)-1, &record_capacity,
		    &record_bytes) ||
		!simplex_presolve_growth(journal->term_capacity,
			journal->term_count, term_count, 64, width, (size_t)-1,
			&term_capacity, &term_bytes))
		return NULL;
	if (record_capacity != journal->capacity) {
		record = (struct simplex_PresolveSubstitution *)
			lp_simplex_realloc(journal->record, record_bytes);
		if (!record)
			return NULL;
		journal->record = record;
		journal->capacity = record_capacity;
	}
	if (term_capacity != journal->term_capacity) {
		int *columns;
		double *multipliers;

		/* The wider-element check above proves both array products fit. */
		column_bytes = (size_t)term_capacity * sizeof(int);
		columns = (int *)lp_simplex_realloc(
			journal->term_column, column_bytes);
		if (!columns)
			return NULL;
		journal->term_column = columns;
		multipliers =
			(double *)lp_simplex_realloc(journal->term_multiplier,
				(size_t)term_capacity * sizeof(double));
		if (!multipliers)
			return NULL;
		journal->term_multiplier = multipliers;
		journal->term_capacity = term_capacity;
	}
	record = journal->record + journal->count;
	record->column = column;
	record->constant = constant;
	record->term_start = journal->term_count;
	record->term_count = term_count;
	if (!term_count)
		journal->count++;
	return record;
}

void simplex_presolve_record_add_term(struct simplex_Presolve *presolve,
	const int column, const double multiplier)
{
	struct simplex_PresolveJournal *journal = &presolve->journal;
	const struct simplex_PresolveSubstitution *record =
		journal->record + journal->count;
	int position = journal->term_count++;

	journal->term_column[position] = column;
	journal->term_multiplier[position] = multiplier;
	/* Publish only after every reserved term is present; no allocation here. */
	if (journal->term_count == record->term_start + record->term_count)
		journal->count++;
}
