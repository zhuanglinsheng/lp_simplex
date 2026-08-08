/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#include "simplex_presolve_record.h"
#include "utils.h"

#include <lp_simplex/status.h>

#include <limits.h>


static int presolve_grow_capacity(
		const int current, const int required, const int initial)
{
	int capacity = current == 0 ? initial : current;
	if (required < 0)
		return -1;
	while (capacity < required) {
		if (capacity > INT_MAX / 2)
			return -1;
		capacity *= 2;
	}
	return capacity;
}


static int presolve_record_reserve_terms(
		struct simplex_Presolve *presolve, const int additional)
{
	int required, capacity;
	int *columns;
	double *multipliers;
	if (additional < 0 || presolve->journal.term_count >
	    INT_MAX - additional)
		return lp_simplex_EXIT_FAILURE;
	required = presolve->journal.term_count + additional;
	if (required <= presolve->journal.term_capacity)
		return lp_simplex_EXIT_SUCCESS;
	capacity = presolve_grow_capacity(
		presolve->journal.term_capacity, required, 64);
	if (capacity < 0)
		return lp_simplex_EXIT_FAILURE;
	columns = (int *)lp_simplex_realloc(presolve->journal.term_column,
		(size_t)capacity * sizeof(int));
	if (columns == NULL)
		return lp_simplex_EXIT_FAILURE;
	presolve->journal.term_column = columns;
	multipliers = (double *)lp_simplex_realloc(
		presolve->journal.term_multiplier,
		(size_t)capacity * sizeof(double));
	if (multipliers == NULL)
		return lp_simplex_EXIT_FAILURE;
	presolve->journal.term_multiplier = multipliers;
	presolve->journal.term_capacity = capacity;
	return lp_simplex_EXIT_SUCCESS;
}


struct simplex_PresolveSubstitution *simplex_presolve_record_begin(
		struct simplex_Presolve *presolve, const int column,
		const double constant, const int term_count)
{
	struct simplex_PresolveSubstitution *record;
	int capacity;
	if (presolve->journal.count == presolve->journal.capacity) {
		capacity = presolve_grow_capacity(presolve->journal.capacity,
			presolve->journal.count + 1, 32);
		if (capacity < 0)
			return NULL;
		record = (struct simplex_PresolveSubstitution *)lp_simplex_realloc(
			presolve->journal.record,
			(size_t)capacity * sizeof(*record));
		if (record == NULL)
			return NULL;
		presolve->journal.record = record;
		presolve->journal.capacity = capacity;
	}
	if (presolve_record_reserve_terms(presolve, term_count) ==
	    lp_simplex_EXIT_FAILURE)
		return NULL;
	record = presolve->journal.record + presolve->journal.count++;
	record->column = column;
	record->constant = constant;
	record->term_start = presolve->journal.term_count;
	record->term_count = term_count;
	return record;
}


void simplex_presolve_record_add_term(
		struct simplex_Presolve *presolve, const int column,
		const double multiplier)
{
	int position = presolve->journal.term_count++;
	presolve->journal.term_column[position] = column;
	presolve->journal.term_multiplier[position] = multiplier;
}
