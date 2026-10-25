/**
 * @file simplex_presolve_growth.h
 * @brief Checked capacity growth for presolve records and mutable supports.
 */

#ifndef LP_SIMPLEX_PRESOLVE_GROWTH_INTERNAL_H
#define LP_SIMPLEX_PRESOLVE_GROWTH_INTERNAL_H

#include "simplex_storage.h"
#include <limits.h>

/* The logical count is independent of physical roots after partial realloc.
 * Neither published count nor capacity is advanced by a refused request. */
static int simplex_presolve_growth(int current, int used, int additional,
	int initial, size_t width, size_t limit, int *capacity, size_t *bytes)
{
	size_t maximum;
	int required, next;

	*capacity = 0;
	*bytes = 0;
	if (!width || current < 0 || used < 0 || used > current ||
		additional < 0 || used > INT_MAX - additional || initial <= 0)
		return 0;
	maximum = limit / width;
	if (maximum > (size_t)INT_MAX)
		maximum = (size_t)INT_MAX;
	required = used + additional;
	if ((size_t)current > maximum || (size_t)required > maximum)
		return 0;
	next = current;
	if (next < required) {
		if (!next)
			next = (size_t)initial > maximum ? required : initial;
		while (next < required) {
			if ((size_t)next > maximum / 2) {
				next = required;
				break;
			}
			next *= 2;
		}
	}
	if (!simplex_storage_bytes((size_t)next, width, limit, bytes))
		return 0;
	*capacity = next;
	return 1;
}

#endif
