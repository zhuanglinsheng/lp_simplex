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
int argmaxabs_arrd(const double *arr, int len, int inc);
double maxabs_arrd(const double *arr, int len, int inc);
double maxabs_arrd_gap(const double *arr1, const double *arr2, int len, int inc);

void lp_simplex_prt_arri(const int *arr, int len, int inc);
void lp_simplex_prt_arrl(const long *arr, int len, int inc);
void lp_simplex_prt_arrd(const double *arr, int len, int inc, int sci);
void lp_simplex_prt_arrld(const long double *arr, int len, int inc, int sci);
void lp_simplex_prt_matd(const double *mat, int ld, int nrow, int ncol);

void *lp_simplex_malloc(size_t size);
void lp_simplex_free(void *ptr);
void *lp_simplex_memset(void *str, int c, size_t n);
void *lp_simplex_memcpy(void *dest, const void *src, size_t n);
int lp_simplex_memcmp(const void *left, const void *right, size_t n);
size_t lp_simplex_strcspn(const char *left, const char *right);
size_t lp_simplex_strlen(const char *text);
double lp_simplex_atof(const char *text);

#endif
