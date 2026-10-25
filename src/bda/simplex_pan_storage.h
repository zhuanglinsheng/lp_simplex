/**
 * @file simplex_pan_storage.h
 * @brief Representation storage checks; these are not a total-memory quota.
 */

#ifndef LP_SIMPLEX_PAN_STORAGE_INTERNAL_H
#define LP_SIMPLEX_PAN_STORAGE_INTERNAL_H
#include <limits.h>
#include <stddef.h>

static int pan_storage_add(int *count, int amount)
{
	if (amount < 0 || *count < 0 || *count > INT_MAX - amount)
		return 0;
	*count += amount;
	return 1;
}

/* Native arrays include one spare element. Standard form only pads column
 * starts, and keeps a dummy column vector when all variables are fixed.
 * The explicit limit permits tests of a narrower address space. */
static int pan_storage_fits(int rows, int columns, int originals, int nonzeros,
	int native, size_t limit)
{
	size_t r, c, n, o;
	if (rows < 0 || columns < 0 || originals < 0 || nonzeros < 0)
		return 0;
	r = (size_t)rows + (native != 0);
	c = native ? (size_t)columns + 1 : (size_t)(columns > 0 ? columns : 1);
	n = (size_t)nonzeros + (native != 0);
	o = (size_t)originals;
	return r <= limit / sizeof(double) && c <= limit / sizeof(double) &&
	       c <= limit / sizeof(int) && n <= limit / sizeof(double) &&
	       n <= limit / sizeof(int) && o <= limit / sizeof(double) &&
	       (size_t)columns + 1 <= limit / sizeof(int);
}
#endif
