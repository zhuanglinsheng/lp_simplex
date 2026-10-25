/**
 * @file test_dual_recovery_restart.c
 * @brief Numerical restart restores the original objective and respects budgets.
 */

#include "simplex_basis.h"
#include "simplex_dual_kernel.h"
#include "simplex_problem.h"
#include "simplex_stage.h"
#include <lp_simplex/lp_simplex.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

static int failures_left, factor_calls;

static int failing_factorize(struct simplex_Basis *basis)
{
	factor_calls++;
	if (failures_left > 0) {
		failures_left--;
		return lp_simplex_EXIT_FAILURE;
	}
	return simplex_basis_factorize(basis);
}

/* Inject only a numerical factorization failure. The actual restart, cost
 * restoration, pivot accounting and terminal checks run without replacements. */
#define simplex_basis_factorize failing_factorize
#define simplex_dual_kernel_solve recovery_kernel_solve
#define simplex_dual_continue recovery_continue
#define simplex_dual_primal_cleanup recovery_primal_cleanup
#include "../src/dual/simplex_dual_kernel.c"
#undef simplex_basis_factorize

static int run_case(int failures, int budget)
{
	struct lp_Model *model = lp_model_create(1, 2);
	struct simplex_Problem problem;
	struct simplex_StageReport report;
	struct lp_simplex_Options options;
	struct lp_simplex_Result result;
	double x[2], dual[1];
	int code, valid;

	if (!model)
		return 0;
	model->constraints[0].type = optm_CONS_T_GE;
	model->constraints[0].rhs = 1.;
	model->constraints[0].coef[0] = model->constraints[0].coef[1] = 1.;
	model->objective[0] = 1.;
	model->objective[1] = 2.;
	if (simplex_problem_from_model(&problem, model)) {
		lp_model_free(model);
		return 0;
	}
	simplex_stage_reset(&report);
	problem.report = &report;
	lp_simplex_default_options(&options, lp_simplex_ALGORITHM_DUAL_REVISED);
	options.iteration_limit = budget;
	memset(&result, 0, sizeof(result));
	failures_left = failures;
	factor_calls = 0;
	code = recovery_kernel_solve(&problem, &options, x, dual, &result, 0,
		SIMPLEX_BASIS_BACKEND_INTERNAL, 0, NULL);
	valid = !failures_left && factor_calls >= 2 &&
		result.iterations <= budget;
	if (failures == 3)
		valid = valid && code == lp_simplex_EXIT_FAILURE &&
			result.status == lp_simplex_PrecisionError &&
			factor_calls == 3 && result.iterations == 0;
	else if (!budget)
		valid = valid && code == lp_simplex_EXIT_FAILURE &&
			result.status == lp_simplex_ExceedIterLimit &&
			result.iterations == 0;
	else
		valid = valid && code == lp_simplex_EXIT_SUCCESS &&
			result.status == lp_simplex_Success &&
			!report.perturbed_costs && report.cost_epoch >= 3 &&
			fabs(x[0] - 1.) < options.primal_tolerance &&
			fabs(x[1]) < options.primal_tolerance &&
			fabs(result.objective - 1.) < options.dual_tolerance &&
			simplex_problem_dual_infeasibility(&problem, x, dual,
				options.primal_tolerance,
				-1) <= options.dual_tolerance;
	if (!valid)
		fprintf(stderr,
			"restart failures=%d budget=%d code=%d status=%d work=%d factors=%d\n",
			failures, budget, code, result.status,
			result.iterations, factor_calls);
	simplex_problem_destroy(&problem);
	lp_model_free(model);
	return valid;
}

static int row_allowance_reset(void)
{
	struct lp_Model *model = lp_model_create(1, 1);
	struct simplex_Problem problem;
	struct simplex_DualState state;
	struct lp_simplex_Options options;
	int valid;

	if (!model)
		return 0;
	model->objective[0] = 1.;
	model->constraints[0].type = optm_CONS_T_GE;
	model->constraints[0].rhs = model->constraints[0].coef[0] = 1.;
	if (simplex_problem_from_model(&problem, model)) {
		lp_model_free(model);
		return 0;
	}
	lp_simplex_default_options(&options, lp_simplex_ALGORITHM_DUAL_REVISED);
	if (simplex_dual_state_create(&state, &problem, &options,
		    SIMPLEX_BASIS_BACKEND_INTERNAL)) {
		simplex_problem_destroy(&problem);
		lp_model_free(model);
		return 0;
	}
	simplex_dual_initialize_bounds(&state, &problem, 0);
	memcpy(state.original_cost, state.cost,
		(size_t)state.variables * sizeof(double));
	state.feasibility_tolerance[0] = 1.;
	state.iterations = 7;
	valid = dual_initialize_stable_phase_one(&state, &problem, 1) == 0 &&
		state.feasibility_tolerance[0] == 0. && state.iterations == 7;
	simplex_dual_state_destroy(&state);
	simplex_problem_destroy(&problem);
	lp_model_free(model);
	return valid;
}

int main(void)
{
	return run_case(1, 20) && run_case(1, 0) && run_case(3, 20) &&
			       row_allowance_reset()
		       ? 0
		       : 1;
}
