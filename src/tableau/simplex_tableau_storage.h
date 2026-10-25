/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_tableau_storage.h
 * @brief Checked Tableau layout and allocation sizes.
 */

#ifndef LP_SIMPLEX_TABLEAU_STORAGE_INTERNAL_H
#define LP_SIMPLEX_TABLEAU_STORAGE_INTERNAL_H
#include <limits.h>
#include <stddef.h>

/* Tableau/transform row offsets are int. A size_t allocation alone cannot
 * make a larger flattened index representable. No total-RSS policy here. */
static int simplex_tableau_storage_size(
	int rows, int columns, size_t width, size_t *bytes)
{
	size_t cells;
	*bytes = 0;
	if (rows <= 0 || columns <= 0 || !width || columns > INT_MAX / rows)
		return 0;
	cells = (size_t)rows * (size_t)columns;
	if (cells > (size_t)-1 / width)
		return 0;
	*bytes = cells * width;
	return 1;
}
#endif
