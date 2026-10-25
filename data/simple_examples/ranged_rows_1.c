/**
 * @file ranged_rows_1.c
 * @brief ranged_rows.mps: Four ranged rows expanded into eight ordinary
 * inequalities.
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
	static const double a[8][1] = {
		{1.}, {1.}, {1.}, {1.}, {1.}, {1.}, {1.}, {1.}};
	static const double c[1] = {-1.};
	static const double rhs[8] = {5., 10., 2., 6., 4., 6., 6., 8.};
	static const int type[8] = {optm_CONS_T_GE, optm_CONS_T_LE,
		optm_CONS_T_GE, optm_CONS_T_LE, optm_CONS_T_GE, optm_CONS_T_LE,
		optm_CONS_T_GE, optm_CONS_T_LE};
	struct lp_Model *model = lp_model_create(8, 1);
	int i, j;

	if (model == NULL)
		return NULL;
	/* RANGES: RL in [5,10], RG in [2,6], EP in [4,6], EN in [6,8].
	 * The API represents each interval with one >= and one <= row. */
	/* Column j represents the same variable in C and in the MPS file. */
	memcpy(model->bounds[0].name, "X", sizeof("X"));
	for (j = 0; j < 1; j++)
		model->objective[j] = c[j];
	for (i = 0; i < 8; i++) {
		model->constraints[i].type = type[i];
		model->constraints[i].rhs = rhs[i];
		for (j = 0; j < 1; j++)
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
	const double expected = -6.;
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
