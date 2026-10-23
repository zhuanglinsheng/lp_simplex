/*
 * production_planning.mps: Production planning: minimize negative profit, min -3 X1 - 2 X2.
 * Manual modeling version: construct and solve the LP through the C API.
 * This standalone example needs only the public lp_simplex API.
 */
#include <lp_simplex/lp_simplex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static struct lp_Model *create_model(void)
{
	static const double a[4][2] = {{1., 1.}, {3., 1.}, {1., 0.}, {0., 1.}};
	static const double c[2] = {-3., -2.};
	static const double rhs[4] = {9., 18., 7., 6.};
	static const int type[4] = {
		optm_CONS_T_LE, optm_CONS_T_LE, optm_CONS_T_LE, optm_CONS_T_LE};
	struct lp_Model *model = lp_model_create(4, 2);
	int i, j;

	if (model == NULL)
		return NULL;
	/* Column j represents the same variable in C and in the MPS file. */
	sprintf(model->bounds[0].name, "%s", "X1");
	sprintf(model->bounds[1].name, "%s", "X2");
	for (j = 0; j < 2; j++)
		model->objective[j] = c[j];
	for (i = 0; i < 4; i++) {
		model->constraints[i].type = type[i];
		model->constraints[i].rhs = rhs[i];
		for (j = 0; j < 2; j++)
			model->constraints[i].coef[j] = a[i][j];
	}

	if (lp_model_build_sparse(model) != lp_simplex_EXIT_SUCCESS) {
		lp_model_free(model);
		return NULL;
	}
	return model;
}

int main(void)
{
	const double expected = -22.5;
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
