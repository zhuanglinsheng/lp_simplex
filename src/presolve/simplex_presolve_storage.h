/**
 * @file simplex_presolve_storage.h
 * @brief Presolve startup array sizes and bounded optional search work.
 */

#ifndef LP_SIMPLEX_PRESOLVE_STORAGE_INTERNAL_H
#define LP_SIMPLEX_PRESOLVE_STORAGE_INTERNAL_H

#include "simplex_storage.h"
#include "policy/search.h"
#include <lp_simplex/model.h>

struct simplex_PresolveStorage {
	size_t column_int, row_int, column_double, row_double;
	size_t bounds, column_flags, row_flags;
};

/* Separate allocations: each array must be representable. This is not a
 * total memory quota and does not validate matrix indices or coefficients. */
static int simplex_presolve_storage(int rows, int columns, size_t limit,
	struct simplex_PresolveStorage *storage)
{
	struct simplex_PresolveStorage next = {0};

	*storage = next;
	if (rows < 0 || columns < 0 ||
		!simplex_storage_bytes((size_t)columns, sizeof(int), limit,
			&next.column_int) ||
		!simplex_storage_bytes(
			(size_t)rows, sizeof(int), limit, &next.row_int) ||
		!simplex_storage_bytes((size_t)columns, sizeof(double), limit,
			&next.column_double) ||
		!simplex_storage_bytes((size_t)rows, sizeof(double), limit,
			&next.row_double) ||
		!simplex_storage_bytes((size_t)columns,
			sizeof(struct optm_VariableBound), limit,
			&next.bounds) ||
		!simplex_storage_bytes((size_t)columns, sizeof(unsigned char),
			limit, &next.column_flags) ||
		!simplex_storage_bytes((size_t)rows, sizeof(unsigned char),
			limit, &next.row_flags))
		return 0;
	*storage = next;
	return 1;
}

/* An unrepresentable optional allowance disables this search. It cannot
 * wrap into an unrelated budget or establish an LP termination condition. */
static size_t simplex_presolve_search_work(
	int rows, int columns, int nonzeros, size_t limit)
{
	size_t count, work;

	if (rows < 0 || columns < 0 || nonzeros < 0 || (size_t)rows > limit ||
		(size_t)columns > limit - (size_t)rows)
		return 0;
	count = (size_t)rows + (size_t)columns;
	if ((size_t)nonzeros > limit - count)
		return 0;
	count += (size_t)nonzeros;
	if (!simplex_storage_bytes(count,
		    SIMPLEX_PRESOLVE_DOMINANCE_WORK_FACTOR, limit, &work))
		return 0;
	return work;
}

#endif
