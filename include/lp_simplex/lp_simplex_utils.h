/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#ifndef __lp_simplex_BASIC_H__
#define __lp_simplex_BASIC_H__

#include <assert.h>
#include <stddef.h>

#ifdef __cpluscplus
extern "C" {
#endif /* __cplusplus */

#ifdef USE_BLAS
extern void dscal_(int *n, double *alpha, double *x, int *incx);
extern void daxpy_(int *n, double *alpha, double *x, int *incx, double *y, int *incy);
#endif

#ifdef LAPACK_USE
extern void dlarfg_(int *n, double *alpha, double *x, int *incx, double *tau);
#endif

void lp_simplex_linalg_daxpy(const int n, const double a, const double *x, const int incx, double *y, const int incy);
void lp_simplex_linalg_dscal(const int n, const double x, double *arr, const int inc);
void lp_simplex_linalg_dlarfg(const int n, double *alpha, double *x, const int incx, double *tau);

int is_in_arri(const int idx, const int *idxset, const int len);
int maxabs_arri(const int *arr, const int len, const int inc);

int argmaxabs_arrd(const double *arr, const int len, const int inc);
double maxabs_arrd(const double *arr, const int len, const int inc);
double maxabs_arrd_gap(const double *arr1, const double *arr2, const int len, const int inc);

void lp_simplex_prt_arri(const int *arr, const int len, const int inc);
void lp_simplex_prt_arrl(const long *arr, const int len, const int inc);
void lp_simplex_prt_arrd(const double *arr, const int len, const int inc, const int sci);
void lp_simplex_prt_arrld(const long double *arr, const int len, const int inc, const int sci);
void lp_simplex_prt_matd(const double *mat, const int ld, const int nrow, const int ncol);

void *lp_simplex_malloc(size_t size);
void lp_simplex_free(void *ptr);

void *lp_simplex_memset(void *str, int c, size_t n);
void *lp_simplex_memcpy(void *dest, const void *src, size_t n);
int lp_simplex_memcmp(const void *str1, const void *str2, size_t n);

size_t lp_simplex_strcspn(const char *str1, const char *str2);
size_t lp_simplex_strlen(const char *str);
double lp_simplex_atof(const char* str);

#ifdef __cpluscplus
}
#endif /* __cpluscplus */

#endif
