/*
 * MPS input version of variable_bounds.mps.
 * Read the paired model file and solve it through the public C API.
 */
#include <lp_simplex/lp_simplex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv)
{
	const double expected = -5.;
	struct lp_Model *model;
	struct lp_simplex_Options options;
	struct lp_simplex_Result result;
	double *x;
	int state, j, passed;

	if (argc != 2) {
		fprintf(stderr, "usage: %s MODEL.mps\n", argv[0]);
		return 2;
	}
	model = lp_read_mps(argv[1]);
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
