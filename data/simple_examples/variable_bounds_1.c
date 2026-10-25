/**
 * @file variable_bounds_1.c
 * @brief variable_bounds.mps: An upper-only X1 and two variables with two-sided
 * bounds.
 *
 * Manual modeling version: construct and solve the LP through the C API.
 * This standalone example needs only the public lp_simplex API.
 */

#include <lp_simplex/lp_simplex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static struct lp_Model *create_model(void)
{
	static const double a[1][3] = {{1., 0., 0.}};
	static const double c[3] = {1., 1., -1.};
	static const double rhs[1] = {-4.};
	static const int type[1] = {optm_CONS_T_GE};
	struct lp_Model *model = lp_model_create(1, 3);
	int i, j;

	if (model == NULL)
		return NULL;
	/* Column j represents the same variable in C and in the MPS file. */
	memcpy(model->bounds[0].name, "X1", sizeof("X1"));
	memcpy(model->bounds[1].name, "X2", sizeof("X2"));
	memcpy(model->bounds[2].name, "X3", sizeof("X3"));
	for (j = 0; j < 3; j++)
		model->objective[j] = c[j];
	for (i = 0; i < 1; i++) {
		model->constraints[i].type = type[i];
		model->constraints[i].rhs = rhs[i];
		for (j = 0; j < 3; j++)
			model->constraints[i].coef[j] = a[i][j];
	}
	/* X1: lower=None, upper=2.0. */
	model->bounds[0].b_type = optm_BOUND_T_UP;
	model->bounds[0].ub = 2.;
	/* X2: lower=3.0, upper=5.0. */
	model->bounds[1].b_type = optm_BOUND_T_BS;
	model->bounds[1].lb = 3.;
	model->bounds[1].ub = 5.;
	/* X3: lower=-2.0, upper=4.0. */
	model->bounds[2].b_type = optm_BOUND_T_BS;
	model->bounds[2].lb = -2.;
	model->bounds[2].ub = 4.;
	if (lp_model_build_sparse(model) != lp_simplex_EXIT_SUCCESS) {
		lp_model_free(model);
		return NULL;
	}
	return model;
}

int main(void)
{
	const double expected = -5.;
	struct lp_Model *model;
	struct lp_simplex_Options options;
	struct lp_simplex_Result result;
	double *x;
	int state, j, passed;

	model = create_model();
	if (model == NULL) {
		fprintf(stderr, "Cannot create or read model\n");
		return 1;
	}
	x = (double *)calloc((size_t)model->n, sizeof(*x));
	if (x == NULL) {
		lp_model_free(model);
		return 1;
	}
	lp_simplex_default_options(&options, lp_simplex_ALGORITHM_AUTO);
	state = lp_simplex_solve(model, &options, x, &result);
	passed = state == lp_simplex_EXIT_SUCCESS &&
		 result.status == lp_simplex_Success;
	if (passed) {
		printf("objective=%.12g primal=%.6g iterations=%d\n",
			result.objective, result.primal_infeasibility,
			result.iterations);
		for (j = 0; j < model->n; j++)
			printf("%s=%.12g\n", model->bounds[j].name, x[j]);
		passed = result.primal_infeasibility <=
				 options.primal_tolerance &&
			 fabs(result.objective - expected) <=
				 1e-7 * (1. + fabs(expected));
		if (!passed)
			fprintf(stderr,
				"Example validation failed; expected objective %.12g\n",
				expected);
	} else {
		fprintf(stderr,
			"Solve failed: state=%d status=%d iterations=%d\n",
			state, result.status, result.iterations);
	}
	free(x);
	lp_model_free(model);
	return passed ? 0 : 1;
}
