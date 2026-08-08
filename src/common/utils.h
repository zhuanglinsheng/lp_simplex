/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#ifndef LP_SIMPLEX_UTILS_INTERNAL_H
#define LP_SIMPLEX_UTILS_INTERNAL_H

#include <assert.h>
#include <stddef.h>


#define __lp_simplex_ABS__(x) ((x) >= 0 ? (x) : (-(x)))
#define __lp_simplex_MAX__(x, y) ((x) >= (y) ? (x) : (y))
#define __lp_simplex_MIN__(x, y) ((x) <= (y) ? (x) : (y))
#define __lp_simplex_INF__ (1. / 0.)
#define __lp_simplex_NINF__ (-1. / 0.)


int is_in_arri(int idx, const int *idxset, int len);

int maxabs_arri(const int *arr, int len, int inc);

void *lp_simplex_malloc(size_t size);
void *lp_simplex_realloc(void *ptr, size_t size);
void lp_simplex_free(void *ptr);
void *lp_simplex_memset(void *str, int c, size_t n);
void *lp_simplex_memcpy(void *dest, const void *src, size_t n);
int lp_simplex_memcmp(const void *left, const void *right, size_t n);
size_t lp_simplex_strcspn(const char *left, const char *right);
size_t lp_simplex_strlen(const char *text);
double lp_simplex_atof(const char *text);

#endif
