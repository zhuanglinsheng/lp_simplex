#ifndef LP_SIMPLEX_LINALG_INTERNAL_H
#define LP_SIMPLEX_LINALG_INTERNAL_H

void lp_simplex_linalg_daxpy(int n, double a, double *x, int incx,
			     double *y, int incy);
void lp_simplex_linalg_dscal(int n, double x, double *arr, int inc);
void lp_simplex_linalg_dlarfg(int n, double *alpha, double *x, int incx,
			      double *tau);

#endif
