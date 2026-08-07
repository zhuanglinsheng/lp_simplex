#ifndef LP_SIMPLEX_PHASE_H
#define LP_SIMPLEX_PHASE_H

#include <lp_simplex/model.h>

int simplex_solve_standard(const double *objective,
			   const struct optm_LinearConstraint *constraints,
			   int m, int n, const char *criteria,
			   int iteration_limit, double *x, double *value,
			   int *status);

#endif
