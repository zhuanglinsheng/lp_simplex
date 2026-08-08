/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#include "utils.h"

#include <stdlib.h>
#include <string.h>


int is_in_arri(const int idx, const int *idxset, const int len)
{
	int i;
	for (i = 0; i < len; i++) {
		if (idx == idxset[i])
			return 1;
	}
	return 0;
}

int maxabs_arri(const int *arr, const int len, const int inc)
{
	int i, value, maximum = 0;
	for (i = 0; i < len; i += inc) {
		value = arr[i] < 0 ? -arr[i] : arr[i];
		if (value > maximum)
			maximum = value;
	}
	return maximum;
}

void *lp_simplex_malloc(size_t size)
{
	return malloc(size);
}

void *lp_simplex_realloc(void *ptr, size_t size)
{
	return realloc(ptr, size);
}

void lp_simplex_free(void *ptr)
{
	free(ptr);
}

void *lp_simplex_memset(void *str, int c, size_t n)
{
	return memset(str, c, n);
}

void *lp_simplex_memcpy(void *dest, const void *src, size_t n)
{
	return memcpy(dest, src, n);
}

int lp_simplex_memcmp(const void *left, const void *right, size_t n)
{
	return memcmp(left, right, n);
}

size_t lp_simplex_strcspn(const char *left, const char *right)
{
	return strcspn(left, right);
}

size_t lp_simplex_strlen(const char *text)
{
	return strlen(text);
}

double lp_simplex_atof(const char *text)
{
	return atof(text);
}
