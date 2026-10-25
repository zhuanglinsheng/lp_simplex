/**
 * @file simplex_dual_storage.h
 * @brief Dual variable indexing and owned startup workspace sizes.
 */

#ifndef LP_SIMPLEX_DUAL_STORAGE_INTERNAL_H
#define LP_SIMPLEX_DUAL_STORAGE_INTERNAL_H

#include "simplex_storage.h"
#include <limits.h>

struct simplex_DualStorage {
	int variables;
	size_t index, status, column, row, multiplier, frequency;
};

/* The structural/logical variable ids remain int. Flattened workspace
 * offsets use size_t; this check adds neither a quota nor a total-RSS bound. */
static int simplex_dual_storage(int rows, int columns, size_t limit,
	struct simplex_DualStorage *storage)
{
	struct simplex_DualStorage next = {0};
	*storage = next;
	if (rows < 0 || columns < 0 || columns > INT_MAX - rows)
		return 0;
	next.variables = rows + columns;
	if (!simplex_storage_bytes((size_t)rows + next.variables, sizeof(int),
		    limit, &next.index) ||
		!simplex_storage_bytes((size_t)next.variables,
			sizeof(unsigned char), limit, &next.status) ||
		!simplex_storage_bytes((size_t)next.variables,
			7 * sizeof(double), limit, &next.column) ||
		!simplex_storage_bytes(
			(size_t)rows, 9 * sizeof(double), limit, &next.row) ||
		!simplex_storage_bytes((size_t)rows + 1, sizeof(double), limit,
			&next.multiplier) ||
		!simplex_storage_bytes(
			(size_t)rows + 1, sizeof(int), limit, &next.frequency))
		return 0;
	*storage = next;
	return 1;
}

#endif
