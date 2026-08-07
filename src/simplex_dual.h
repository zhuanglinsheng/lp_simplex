#ifndef LP_SIMPLEX_DUAL_INTERNAL_H
#define LP_SIMPLEX_DUAL_INTERNAL_H

#include <lp_simplex/model.h>
#include <lp_simplex/solve.h>

int simplex_dual_solve(const struct lp_Model *model,
		const struct lp_simplex_Options *options,
		double *x, double *row_dual, struct lp_simplex_Result *result);

#endif
