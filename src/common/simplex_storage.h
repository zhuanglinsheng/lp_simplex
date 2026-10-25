/**
 * @file simplex_storage.h
 * @brief Checked array byte counts, independent of numerical or memory policy.
 */

#ifndef LP_SIMPLEX_STORAGE_INTERNAL_H
#define LP_SIMPLEX_STORAGE_INTERNAL_H

#include <stddef.h>

/* Use the host size_t maximum in production. The explicit limit permits
 * narrow-address-space arithmetic tests without allocating large buffers. */
static int simplex_storage_bytes(
	size_t count, size_t width, size_t limit, size_t *bytes)
{
	*bytes = 0;
	if (!width || count > limit / width)
		return 0;
	*bytes = count * width;
	return 1;
}

#endif
