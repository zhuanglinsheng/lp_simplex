/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#include "utils.h"
#include <stdio.h>
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

int argmaxabs_arrd(const double *arr, const int len, const int inc)
{
	int i, result = 0;
	double maximum = 0.;
	for (i = 0; i < len; i += inc) {
		double value = arr[i] < 0. ? -arr[i] : arr[i];
		if (value > maximum) {
			maximum = value;
			result = i;
		}
	}
	return result;
}

double maxabs_arrd(const double *arr, const int len, const int inc)
{
	int index = argmaxabs_arrd(arr, len, inc);
	double value = arr[index];
	return value < 0. ? -value : value;
}

double maxabs_arrd_gap(
		const double *left, const double *right, const int len, const int inc)
{
	int i;
	double maximum = 0.;
	for (i = 0; i < len; i += inc) {
		double value = left[i] - right[i];
		if (value < 0.)
			value = -value;
		if (value > maximum)
			maximum = value;
	}
	return maximum;
}

void lp_simplex_prt_arri(const int *arr, const int len, const int inc)
{
	int i;
	for (i = 0; i < len; i += inc)
		printf(i + inc < len ? "%d, " : "%d\n", arr[i]);
}

void lp_simplex_prt_arrl(const long *arr, const int len, const int inc)
{
	int i;
	for (i = 0; i < len; i += inc)
		printf(i + inc < len ? "%ld, " : "%ld\n", arr[i]);
}

void lp_simplex_prt_arrd(
		const double *arr, const int len, const int inc, const int scientific)
{
	int i;
	const char *middle = scientific ? "%e, " : "%f, ";
	const char *last = scientific ? "%e\n" : "%f\n";
	for (i = 0; i < len; i += inc)
		printf(i + inc < len ? middle : last, arr[i]);
}

void lp_simplex_prt_arrld(
		const long double *arr, const int len, const int inc, const int scientific)
{
	int i;
	const char *middle = scientific ? "%Le, " : "%Lf, ";
	const char *last = scientific ? "%Le\n" : "%Lf\n";
	for (i = 0; i < len; i += inc)
		printf(i + inc < len ? middle : last, arr[i]);
}

void lp_simplex_prt_matd(
		const double *matrix, const int ld, const int rows, const int columns)
{
	int i, j;
	for (i = 0; i < rows; i++) {
		for (j = 0; j < columns; j++)
			printf(j + 1 < columns ? "%e, " : "%e", matrix[j + i * ld]);
		printf("\n");
	}
}

void *lp_simplex_malloc(size_t size) { return malloc(size); }
void lp_simplex_free(void *ptr) { free(ptr); }
void *lp_simplex_memset(void *str, int c, size_t n) { return memset(str, c, n); }
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
size_t lp_simplex_strlen(const char *text) { return strlen(text); }
double lp_simplex_atof(const char *text) { return atof(text); }
