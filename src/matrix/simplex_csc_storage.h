/**
 * @file simplex_csc_storage.h
 * @brief Checked owned CSC/CSR array sizes, independent of solver policy.
 */

#ifndef LP_SIMPLEX_CSC_STORAGE_INTERNAL_H
#define LP_SIMPLEX_CSC_STORAGE_INTERNAL_H

#include "simplex_storage.h"

struct simplex_CscStorage {
	size_t columns, rows, indices, values, cursor;
};

/* Zero-dimensional auxiliary representations are allowed. This checks byte
 * representation, not coordinate validity or an algorithm's variable ids. */
static int simplex_csc_storage(int rows, int columns, int nonzeros,
	size_t limit, struct simplex_CscStorage *storage)
{
	struct simplex_CscStorage next = {0};
	*storage = next;
	if (rows < 0 || columns < 0 || nonzeros < 0 ||
		!simplex_storage_bytes((size_t)columns + 1, sizeof(int), limit,
			&next.columns) ||
		!simplex_storage_bytes(
			(size_t)rows + 1, sizeof(int), limit, &next.rows) ||
		!simplex_storage_bytes(
			(size_t)nonzeros, sizeof(int), limit, &next.indices) ||
		!simplex_storage_bytes((size_t)nonzeros, sizeof(double), limit,
			&next.values) ||
		!simplex_storage_bytes(
			(size_t)rows, sizeof(int), limit, &next.cursor))
		return 0;
	*storage = next;
	return 1;
}

#endif
