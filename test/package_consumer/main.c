/**
 * @file main.c
 * @brief Public headers and installed target only; no source/build-tree
 * includes.
 */

#include <lp_simplex/lp_simplex.h>
#include <stdio.h>
#include <string.h>

static void progress(void *context, const struct lp_simplex_Progress *event)
{
	int *calls = (int *)context;
	if (event->stage)
		(*calls)++;
}

int main(void)
{
	const enum lp_simplex_Algorithm algorithms[] = {
		lp_simplex_ALGORITHM_AUTO, lp_simplex_ALGORITHM_DUAL_REVISED,
		lp_simplex_ALGORITHM_TABLEAU, lp_simplex_ALGORITHM_PAN_BDA};
	struct lp_Model *m = lp_model_create(1, 2);
	struct lp_simplex_Options o;
	struct lp_simplex_Result r;
	double x[2];
	int i, calls;
	if (strcmp(lp_simplex_version(), LP_SIMPLEX_VERSION) != 0) {
		lp_model_free(m);
		return 3;
	}
	if (!m)
		return 1;
	m->constraints[0].type = optm_CONS_T_GE;
	m->constraints[0].rhs = 1.;
	m->constraints[0].coef[0] = m->constraints[0].coef[1] = 1.;
	m->objective[0] = 1.;
	m->objective[1] = 2.;
	for (i = 0; i < 4; i++) {
		lp_simplex_default_options(&o, algorithms[i]);
		o.iteration_limit = 100;
		o.progress_callback = progress;
		o.progress_context = &calls;
		calls = 0;
		if (lp_simplex_solve(m, &o, x, &r) ||
			r.status != lp_simplex_Success || !calls ||
			!(r.objective > .999999 && r.objective < 1.000001) ||
			!(x[0] > .999999 && x[0] < 1.000001) ||
			!(x[1] > -.000001 && x[1] < .000001)) {
			fprintf(stderr,
				"installed algorithm %d failed: status=%d calls=%d\n",
				i, r.status, calls);
			lp_model_free(m);
			return 2;
		}
	}
	lp_model_free(m);
	puts("installed public C API: AUTO/Dual/Tableau/BDA and progress passed");
	return 0;
}
