/**
 * @file simplex_model_audit.h
 * @brief Original public-model numerical checks shared by assembly and
 * recovery.
 */

#ifndef LP_SIMPLEX_MODEL_AUDIT_INTERNAL_H
#define LP_SIMPLEX_MODEL_AUDIT_INTERNAL_H
#include <lp_simplex/solve.h>

struct lp_Model;
double simplex_primal_infeasibility(const struct lp_Model *, const double *);
double simplex_original_objective(const struct lp_Model *, const double *);
int simplex_prefer_point(const struct lp_Model *,
	const struct lp_simplex_Options *, const double *, const double *);
int simplex_finalize_result(const struct lp_Model *,
	const struct lp_simplex_Options *, const double *,
	struct lp_simplex_Result *, int);
#endif
