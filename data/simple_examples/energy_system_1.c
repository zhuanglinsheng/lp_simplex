/**
 * @file energy_system_1.c
 * @brief energy_system.mps: A 27-row, 32-variable energy-system model, stored
 * as nonzero entries.
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
	static const struct {
		int row, column;
		double value;
	} entries[] = {
		{0, 0, -1.}, /* C1, X1 */
		{0, 1, 1.}, /* C1, X2 */
		{0, 2, 1.}, /* C1, X3 */
		{1, 0, -1.06}, /* C2, X1 */
		{1, 3, 1.}, /* C2, X4 */
		{2, 0, 1.}, /* C3, X1 */
		{3, 1, -1.}, /* C4, X2 */
		{3, 12, 1.4}, /* C4, X13 */
		{4, 4, -1.}, /* C5, X5 */
		{4, 5, -1.}, /* C5, X6 */
		{4, 6, -1.}, /* C5, X7 */
		{4, 7, -1.}, /* C5, X8 */
		{4, 12, 1.}, /* C5, X13 */
		{4, 13, 1.}, /* C5, X14 */
		{5, 4, -1.06}, /* C6, X5 */
		{5, 5, -1.06}, /* C6, X6 */
		{5, 6, -0.96}, /* C6, X7 */
		{5, 7, -0.86}, /* C6, X8 */
		{5, 14, 1.}, /* C6, X15 */
		{6, 4, 1.}, /* C7, X5 */
		{6, 8, -1.}, /* C7, X9 */
		{7, 5, 1.}, /* C8, X6 */
		{7, 9, -1.}, /* C8, X10 */
		{8, 6, 1.}, /* C9, X7 */
		{8, 10, -1.}, /* C9, X11 */
		{9, 7, 1.}, /* C10, X8 */
		{9, 11, -1.}, /* C10, X12 */
		{10, 15, -1.}, /* C11, X16 */
		{10, 16, 1.}, /* C11, X17 */
		{10, 17, 1.}, /* C11, X18 */
		{10, 18, 1.}, /* C11, X19 */
		{11, 15, -0.43}, /* C12, X16 */
		{11, 19, 1.}, /* C12, X20 */
		{12, 15, 1.}, /* C13, X16 */
		{13, 16, -1.}, /* C14, X17 */
		{13, 28, 1.4}, /* C14, X29 */
		{14, 20, -0.43}, /* C15, X21 */
		{14, 21, -0.43}, /* C15, X22 */
		{14, 22, -0.39}, /* C15, X23 */
		{14, 23, -0.37}, /* C15, X24 */
		{14, 30, 1.}, /* C15, X31 */
		{15, 20, 1.}, /* C16, X21 */
		{15, 21, 1.}, /* C16, X22 */
		{15, 22, 1.}, /* C16, X23 */
		{15, 23, 1.}, /* C16, X24 */
		{15, 28, -1.}, /* C16, X29 */
		{15, 29, 1.}, /* C16, X30 */
		{15, 31, 1.}, /* C16, X32 */
		{16, 20, 1.}, /* C17, X21 */
		{16, 24, -1.}, /* C17, X25 */
		{17, 21, 1.}, /* C18, X22 */
		{17, 25, -1.}, /* C18, X26 */
		{18, 22, 1.}, /* C19, X23 */
		{18, 26, -1.}, /* C19, X27 */
		{19, 23, 1.}, /* C20, X24 */
		{19, 27, -1.}, /* C20, X28 */
		{20, 8, 2.364}, /* C21, X9 */
		{20, 9, 2.386}, /* C21, X10 */
		{20, 10, 2.408}, /* C21, X11 */
		{20, 11, 2.429}, /* C21, X12 */
		{20, 18, -1.}, /* C21, X19 */
		{20, 24, 2.191}, /* C21, X25 */
		{20, 25, 2.219}, /* C21, X26 */
		{20, 26, 2.249}, /* C21, X27 */
		{20, 27, 2.279}, /* C21, X28 */
		{21, 2, -1.}, /* C22, X3 */
		{21, 15, 0.109}, /* C22, X16 */
		{22, 13, -1.}, /* C23, X14 */
		{22, 20, 0.109}, /* C23, X21 */
		{22, 21, 0.108}, /* C23, X22 */
		{22, 22, 0.108}, /* C23, X23 */
		{22, 23, 0.107}, /* C23, X24 */
		{23, 0, 0.301}, /* C24, X1 */
		{23, 17, -1.}, /* C24, X18 */
		{24, 4, 0.301}, /* C25, X5 */
		{24, 5, 0.313}, /* C25, X6 */
		{24, 6, 0.313}, /* C25, X7 */
		{24, 7, 0.326}, /* C25, X8 */
		{24, 29, -1.}, /* C25, X30 */
		{25, 3, 1.}, /* C26, X4 */
		{25, 19, 1.}, /* C26, X20 */
		{26, 14, 1.}, /* C27, X15 */
		{26, 30, 1.}, /* C27, X31 */
	};

	static const double c[32] = {0., -0.4, 0., 0., 0., 0., 0., 0., 0., 0.,
		0., 0., -0.32, 0., 0., 0., -0.6, 0., 0., 0., 0., 0., 0., 0., 0.,
		0., 0., 0., -0.48, 0., 0., 10.};
	static const double rhs[27] = {0., 0., 80., 0., 0., 0., 80., 0., 0., 0.,
		0., 0., 500., 0., 0., 0., 44., 500., 0., 0., 0., 0., 0., 0., 0.,
		310., 300.};
	static const int type[27] = {optm_CONS_T_EQ, optm_CONS_T_EQ,
		optm_CONS_T_LE, optm_CONS_T_LE, optm_CONS_T_EQ, optm_CONS_T_EQ,
		optm_CONS_T_LE, optm_CONS_T_LE, optm_CONS_T_LE, optm_CONS_T_LE,
		optm_CONS_T_EQ, optm_CONS_T_EQ, optm_CONS_T_LE, optm_CONS_T_LE,
		optm_CONS_T_EQ, optm_CONS_T_EQ, optm_CONS_T_LE, optm_CONS_T_LE,
		optm_CONS_T_LE, optm_CONS_T_LE, optm_CONS_T_LE, optm_CONS_T_LE,
		optm_CONS_T_LE, optm_CONS_T_LE, optm_CONS_T_LE, optm_CONS_T_LE,
		optm_CONS_T_LE};
	struct lp_Model *model = lp_model_create(27, 32);
	int i, j, k;

	if (model == NULL)
		return NULL;
	/* Column j represents the same variable in C and in the MPS file. */
	memcpy(model->bounds[0].name, "X1", sizeof("X1"));
	memcpy(model->bounds[1].name, "X2", sizeof("X2"));
	memcpy(model->bounds[2].name, "X3", sizeof("X3"));
	memcpy(model->bounds[3].name, "X4", sizeof("X4"));
	memcpy(model->bounds[4].name, "X5", sizeof("X5"));
	memcpy(model->bounds[5].name, "X6", sizeof("X6"));
	memcpy(model->bounds[6].name, "X7", sizeof("X7"));
	memcpy(model->bounds[7].name, "X8", sizeof("X8"));
	memcpy(model->bounds[8].name, "X9", sizeof("X9"));
	memcpy(model->bounds[9].name, "X10", sizeof("X10"));
	memcpy(model->bounds[10].name, "X11", sizeof("X11"));
	memcpy(model->bounds[11].name, "X12", sizeof("X12"));
	memcpy(model->bounds[12].name, "X13", sizeof("X13"));
	memcpy(model->bounds[13].name, "X14", sizeof("X14"));
	memcpy(model->bounds[14].name, "X15", sizeof("X15"));
	memcpy(model->bounds[15].name, "X16", sizeof("X16"));
	memcpy(model->bounds[16].name, "X17", sizeof("X17"));
	memcpy(model->bounds[17].name, "X18", sizeof("X18"));
	memcpy(model->bounds[18].name, "X19", sizeof("X19"));
	memcpy(model->bounds[19].name, "X20", sizeof("X20"));
	memcpy(model->bounds[20].name, "X21", sizeof("X21"));
	memcpy(model->bounds[21].name, "X22", sizeof("X22"));
	memcpy(model->bounds[22].name, "X23", sizeof("X23"));
	memcpy(model->bounds[23].name, "X24", sizeof("X24"));
	memcpy(model->bounds[24].name, "X25", sizeof("X25"));
	memcpy(model->bounds[25].name, "X26", sizeof("X26"));
	memcpy(model->bounds[26].name, "X27", sizeof("X27"));
	memcpy(model->bounds[27].name, "X28", sizeof("X28"));
	memcpy(model->bounds[28].name, "X29", sizeof("X29"));
	memcpy(model->bounds[29].name, "X30", sizeof("X30"));
	memcpy(model->bounds[30].name, "X31", sizeof("X31"));
	memcpy(model->bounds[31].name, "X32", sizeof("X32"));
	for (j = 0; j < 32; j++)
		model->objective[j] = c[j];
	for (i = 0; i < 27; i++) {
		model->constraints[i].type = type[i];
		model->constraints[i].rhs = rhs[i];
	}
	for (k = 0; k < (int)(sizeof(entries) / sizeof(entries[0])); k++)
		model->constraints[entries[k].row].coef[entries[k].column] =
			entries[k].value;

	if (lp_model_build_sparse(model) != lp_simplex_EXIT_SUCCESS) {
		lp_model_free(model);
		return NULL;
	}
	return model;
}

int main(void)
{
	const double expected = -464.7531429;
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
