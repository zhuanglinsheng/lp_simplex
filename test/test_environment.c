/**
 * @file test_environment.c
 * @brief Regression tests for environment.
 */

#define _POSIX_C_SOURCE 200809L
#include "simplex_solve_context.h"
#include "simplex_execution.h"
#include "simplex_dual_internal.h"
#include "simplex_dual.h"
#include "simplex_dual_state.h"
#include "policy/controls.h"
#include <lp_simplex/lp_simplex.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

#ifdef _WIN32
#define setenv(name, value, replace) _putenv_s((name), (value))
#define unsetenv(name) _putenv_s((name), "")
#endif
#define CHECK(c) \
	do { \
		if (!(c)) { \
			fprintf(stderr, "environment line %d: %s\n", __LINE__, \
				#c); \
			return 1; \
		} \
	} while (0)
static const char *const names[] = {SIMPLEX_ENV_DISABLE_PRESOLVE,
	SIMPLEX_ENV_DISABLE_SCALING, SIMPLEX_ENV_STABLE_KLU,
	SIMPLEX_ENV_CLEANUP_DEVEX, SIMPLEX_ENV_DISABLE_BOUND_PERTURB,
	SIMPLEX_ENV_DISABLE_PAN, SIMPLEX_ENV_FORCE_BOUND_PERTURB,
	SIMPLEX_ENV_INITIAL_COST_PERTURB, SIMPLEX_ENV_PRIMAL_BOUND_PERTURB,
	SIMPLEX_ENV_SINGLE_PERTURB_EPOCH, SIMPLEX_ENV_PAN_AFFINE_MULTIPLIER,
	SIMPLEX_ENV_PAN_BLOCK_CRASH, SIMPLEX_ENV_PAN_CENTER_EPOCH,
	SIMPLEX_ENV_PAN_COMPATIBLE_PRICING, SIMPLEX_ENV_PAN_COMPLETE_START,
	SIMPLEX_ENV_PAN_CONE, SIMPLEX_ENV_PAN_CONE_BATCH,
	SIMPLEX_ENV_PAN_DISABLE_SINGLETON_CRASH, SIMPLEX_ENV_PAN_DUAL_CENTER,
	SIMPLEX_ENV_PAN_DUAL_POLISH, SIMPLEX_ENV_PAN_DUAL_START,
	SIMPLEX_ENV_PAN_FACE_GUARD, SIMPLEX_ENV_PAN_FIRST_ORDER,
	SIMPLEX_ENV_PAN_FULL_REFINEMENT, SIMPLEX_ENV_PAN_IMPLICIT_BULK,
	SIMPLEX_ENV_PAN_IMPLICIT_UPDATES, SIMPLEX_ENV_PAN_NATIVE_ADAPTIVE,
	SIMPLEX_ENV_PAN_NATIVE_BOUNDS, SIMPLEX_ENV_PAN_NATIVE_PDHG,
	SIMPLEX_ENV_PAN_NATIVE_PREDICT, SIMPLEX_ENV_PAN_NEWTON_QR,
	SIMPLEX_ENV_PAN_NNLS_BLOCK, SIMPLEX_ENV_PAN_PDHG_RESTORE,
	SIMPLEX_ENV_PAN_PDHG_WEIGHT, SIMPLEX_ENV_PAN_STABILIZE};

static void clear_controls(void)
{
	size_t i;
	for (i = 0; i < sizeof(names) / sizeof(names[0]); i++)
		unsetenv(names[i]);
}

static void poison(void)
{
	size_t i;
	for (i = 0; i < sizeof(names) / sizeof(names[0]); i++)
		setenv(names[i], "0", 1);
	setenv(SIMPLEX_ENV_PAN_CONE, "implicit", 1);
	setenv(SIMPLEX_ENV_PAN_NATIVE_PDHG, "7", 1);
	setenv(SIMPLEX_ENV_PAN_DUAL_POLISH, "11", 1);
	setenv(SIMPLEX_ENV_PAN_FIRST_ORDER, "9", 1);
	setenv(SIMPLEX_ENV_PAN_NNLS_BLOCK, "3", 1);
	setenv(SIMPLEX_ENV_PAN_PDHG_WEIGHT, "2", 1);
}

static struct lp_Model *fixture(void)
{
	struct lp_Model *m = lp_model_create(2, 3);
	if (!m)
		return NULL;
	m->constraints[0].type = optm_CONS_T_EQ;
	m->constraints[0].rhs = 1.;
	m->constraints[0].coef[0] = m->constraints[0].coef[1] = 1.;
	m->constraints[1].type = optm_CONS_T_GE;
	m->constraints[1].rhs = 1.;
	m->constraints[1].coef[1] = m->constraints[1].coef[2] = 1.;
	m->bounds[0].b_type = optm_BOUND_T_FR;
	m->objective[0] = 1.;
	m->objective[1] = 2.;
	m->objective[2] = .5;
	return m;
}

static int same(double a, double b)
{
	return (isnan(a) && isnan(b)) || a == b;
}

static int production(void)
{
	struct lp_Model *m = fixture();
	struct lp_simplex_Options o;
	struct lp_simplex_Result base, r;
	struct simplex_Execution before, after;
	double x[3], expected[3];
	int a, p, j, kind, problem, code;
	int base_count[SIMPLEX_EXECUTION_PAN_FACE_OPPORTUNITY + 1],
		count[SIMPLEX_EXECUTION_PAN_FACE_OPPORTUNITY + 1];
	CHECK(m);
	for (problem = 0; problem < 3; problem++) {
		m->constraints[0].coef[0] = 1.;
		m->constraints[1].coef[0] = 0.;
		m->constraints[1].coef[2] = 1.;
		m->constraints[1].rhs = 1.;
		m->constraints[1].type = optm_CONS_T_GE;
		if (problem == 1) {
			m->constraints[1].coef[0] = 1.;
			m->constraints[1].coef[2] = 0.;
			m->constraints[1].rhs = 2.;
			m->constraints[1].type = optm_CONS_T_EQ;
		}
		if (problem == 2)
			m->constraints[0].coef[0] = 0.;
		for (a = 0; a <= lp_simplex_ALGORITHM_AUTO; a++)
			for (p = 0; p < 2; p++) {
				clear_controls();
				lp_simplex_default_options(&o, a);
				o.presolve = p;
				code = simplex_solve_reported(
					m, &o, expected, &base, &before);
				CHECK(problem ? code != 0 : code == 0);
				/* This boundary test also covers honest nonterminal numerical exits.
				 * It must not turn a known contradictory LP into a success/other proof. */
				CHECK(base.status ==
						(problem == 0 ? lp_simplex_Success
							: problem == 1
								? lp_simplex_Infeasibility
								: lp_simplex_Unboundedness) ||
					(problem &&
						base.status ==
							lp_simplex_PrecisionError));
				if (!problem)
					CHECK(base.objective == 1.5);
				memset(base_count, 0, sizeof(base_count));
				memset(count, 0, sizeof(count));
				for (j = 0; j < before.count; j++) {
					kind = before.events[j].kind;
					CHECK(kind > 0 &&
						kind < (int)(sizeof(base_count) /
							       sizeof(base_count
									       [0])));
					base_count[kind]++;
				}
				poison();
				CHECK(simplex_solve_reported(
					      m, &o, x, &r, &after) == code);
				CHECK(r.status == base.status &&
					r.iterations == base.iterations);
				CHECK(same(r.objective, base.objective) &&
					same(r.primal_infeasibility,
						base.primal_infeasibility));
				CHECK(same(r.dual_infeasibility,
					      base.dual_infeasibility) &&
					r.dual_residual_space ==
						base.dual_residual_space);
				for (j = 0; j < 3; j++)
					CHECK(same(x[j], expected[j]));
				CHECK(!after.events[0]
						.experimental_environment);
				for (j = 0; j < after.count; j++) {
					kind = after.events[j].kind;
					CHECK(kind > 0 &&
						kind < (int)(sizeof(count) /
							       sizeof(count[0])));
					count[kind]++;
				}
				CHECK(!memcmp(
					base_count, count, sizeof(count)));
			}
	}
	clear_controls();
	lp_model_free(m);
	return 0;
}

static void change_environment(
	void *context, const struct lp_simplex_Progress *progress)
{
	int *called = (int *)context;
	(void)progress;
	if (!*called) {
		clear_controls();
		*called = 1;
	}
}

static int experimental_snapshot(void)
{
	struct lp_Model *m = fixture();
	struct lp_simplex_Options o;
	struct lp_simplex_Result first, r;
	struct simplex_Execution h;
	double x[3];
	int algorithm, i, called;
	CHECK(m);
	for (algorithm = lp_simplex_ALGORITHM_TABLEAU;
		algorithm <= lp_simplex_ALGORITHM_AUTO; algorithm++) {
		poison();
		lp_simplex_default_options(&o, algorithm);
		o.iteration_limit = 2000;
		CHECK(!simplex_solve_reported_experimental(
			m, &o, x, &first, &h));
		CHECK(first.objective == 1.5 &&
			h.events[0].experimental_environment);
		for (i = 0; i < h.count; i++)
			CHECK(h.events[i].kind != SIMPLEX_EXECUTION_PRESOLVE &&
				h.events[i].kind != SIMPLEX_EXECUTION_SCALING);
		if (algorithm == lp_simplex_ALGORITHM_PAN_BDA)
			CHECK(first.dual_residual_space ==
				lp_simplex_RESIDUAL_ORIGINAL);
		poison();
		called = 0;
		o.progress_callback = change_environment;
		o.progress_context = &called;
		CHECK(!simplex_solve_reported_experimental(m, &o, x, &r, &h));
		CHECK(called && r.objective == first.objective &&
			r.iterations == first.iterations);
		for (i = 0; i < h.count; i++)
			CHECK(h.events[i].kind != SIMPLEX_EXECUTION_PRESOLVE &&
				h.events[i].kind != SIMPLEX_EXECUTION_SCALING);
	}
	clear_controls();
	lp_model_free(m);
	return 0;
}

static int kernel_values(void)
{
	struct lp_Model *m = fixture();
	struct simplex_Problem p;
	struct simplex_SolveContext context;
	struct simplex_DualState state;
	struct simplex_RunConfig clean, legacy;
	struct lp_simplex_Options o;
	struct simplex_PanConfig pan;
	unsigned long signature;
	CHECK(m && !simplex_problem_from_model(&p, m));
	poison();
	simplex_run_config_init(&clean, 0);
	simplex_run_config_init(&legacy, 1);
	CHECK(!clean.experimental_environment && !clean.disable_presolve &&
		!clean.disable_scaling);
	CHECK(!clean.dual.force_bounds && clean.dual.disable_bounds &&
		clean.dual.pan_enabled && clean.dual.multi_epoch);
	CHECK(legacy.experimental_environment && legacy.disable_presolve &&
		legacy.disable_scaling);
	CHECK(legacy.dual.force_bounds && legacy.dual.disable_bounds &&
		!legacy.dual.pan_enabled && !legacy.dual.multi_epoch);
	CHECK(legacy.dual.initial_cost_perturb && legacy.dual.cleanup_devex &&
		legacy.dual.stable_klu);
	CHECK(legacy.pan.model_form == SIMPLEX_PAN_NATIVE_BOUNDS &&
		legacy.pan.pdhg_iterations == 7 &&
		legacy.pan.nnls_block_limit == 3);
	simplex_solve_context_init(&context, NULL, &p.execution);
	simplex_solve_context_configure(&context, 1);
	signature = simplex_run_dual_signature(&context.controls.dual);
	clear_controls();
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_DUAL_REVISED);
	CHECK(!simplex_dual_state_create(
		&state, &p, &o, SIMPLEX_BASIS_BACKEND_DEFAULT));
	CHECK(simplex_run_dual_signature(&state.controls) == signature &&
		state.recovery.force_bounds && !state.pan_enabled);
	simplex_dual_state_destroy(&state);
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_PAN_BDA);
	simplex_pan_config_resolve(&pan, &p, &o);
	CHECK(pan.model_form == SIMPLEX_PAN_NATIVE_BOUNDS &&
		pan.pdhg_iterations == 7 && pan.nnls_block_limit == 3);
	simplex_problem_destroy(&p);
	lp_model_free(m);
	simplex_solve_context_destroy(&context);
	return 0;
}

static int scaling_isolation(void)
{
	struct lp_Model *m = lp_model_create(1, 1);
	struct lp_simplex_Options o;
	struct lp_simplex_Result r;
	double x[1];
	CHECK(m);
	clear_controls();
	m->constraints[0].type = optm_CONS_T_EQ;
	m->constraints[0].rhs = 1.;
	m->constraints[0].coef[0] = 1e-12;
	m->bounds[0].b_type = optm_BOUND_T_FR;
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_DUAL_REVISED);
	o.presolve = 0;
	o.iteration_limit = 100;
	setenv(SIMPLEX_ENV_DISABLE_SCALING, "0", 1);
	CHECK(!lp_simplex_solve(m, &o, x, &r));
	CHECK(r.status == lp_simplex_Success &&
		fabs(x[0] / 1e12 - 1.) < 1e-12 &&
		r.primal_infeasibility <= o.primal_tolerance);
	CHECK(lp_simplex_solve_experimental(m, &o, x, &r));
	CHECK(r.status == lp_simplex_PrecisionError ||
		r.status == lp_simplex_Singularity);
	clear_controls();
	lp_model_free(m);
	return 0;
}

static void toggle_scaling(void *context, const struct lp_simplex_Progress *p)
{
	int *called = (int *)context;
	(void)p;
	if (!*called) {
		if (getenv(SIMPLEX_ENV_DISABLE_SCALING))
			unsetenv(SIMPLEX_ENV_DISABLE_SCALING);
		else
			setenv(SIMPLEX_ENV_DISABLE_SCALING, "0", 1);
		*called = 1;
	}
}

/* The dense Dual specialization must retain the public root's configuration
 * and model namespace instead of entering the legacy private environment API. */
static int singleton_root(void)
{
	struct lp_Model *m = lp_model_create(1, 3);
	struct lp_simplex_Options o;
	struct lp_simplex_Result r, base;
	struct simplex_Execution h;
	double x[3], expected[3];
	int experiment, disabled, mutate, called, state, i, outer, model,
		scaling;
	CHECK(m);
	m->constraints[0].type = optm_CONS_T_EQ;
	m->constraints[0].rhs = 2.;
	m->constraints[0].coef[0] = 1.;
	m->constraints[0].coef[1] = -1.;
	m->constraints[0].coef[2] = -1e-12;
	m->objective[0] = m->objective[1] = 1.;
	m->objective[2] = -5e-13;
	for (experiment = 0; experiment < 2; experiment++)
		for (disabled = 0; disabled < 2; disabled++)
			for (mutate = 0; mutate < 2; mutate++) {
				clear_controls();
				if (disabled)
					setenv(SIMPLEX_ENV_DISABLE_SCALING, "0",
						1);
				lp_simplex_default_options(
					&o, lp_simplex_ALGORITHM_DUAL_REVISED);
				o.presolve = 0;
				o.iteration_limit = 100;
				called = 0;
				o.progress_callback =
					mutate ? toggle_scaling : NULL;
				o.progress_context = &called;
				state = experiment
						? simplex_solve_reported_experimental(
							  m, &o, x, &r, &h)
						: simplex_solve_reported(
							  m, &o, x, &r, &h);
				CHECK(!state &&
					r.status == lp_simplex_Success &&
					r.objective == 2.);
				CHECK(called == mutate &&
					h.events[0].experimental_environment ==
						experiment);
				outer = model = scaling = 0;
				for (i = 0; i < h.count; i++) {
					const struct simplex_ExecutionEvent *e =
						h.events + i;
					CHECK(e->complete &&
						e->context ==
							h.events[0].context);
					if (e->kind ==
						SIMPLEX_EXECUTION_SINGLETON) {
						outer = i + 1;
						CHECK(e->model == 1 &&
							e->parent == 1 &&
							e->grant == 100 &&
							e->singleton
								.input_point_checked &&
							e->work ==
								r.iterations);
					}
					if (e->kind ==
						SIMPLEX_EXECUTION_SINGLETON_MODEL) {
						model = i + 1;
						CHECK(outer && e->model > 1 &&
							e->source_model == 1 &&
							e->parent == outer);
					}
					if (e->kind ==
						SIMPLEX_EXECUTION_SCALING) {
						scaling++;
						CHECK(model &&
							e->source_model ==
								h.events[model -
									 1]
									.model &&
							e->parent == model);
					}
				}
				if (!outer || !model)
					fprintf(stderr,
						"public Singleton lost root: outer=%d model=%d scaling=%d\n",
						outer, model, scaling);
				CHECK(outer && model &&
					scaling == !(experiment && disabled));
				if (!mutate) {
					base = r;
					memcpy(expected, x, sizeof(x));
				} else
					CHECK(r.status == base.status &&
						r.iterations ==
							base.iterations &&
						r.objective == base.objective &&
						!memcmp(expected, x,
							sizeof(x)));
			}
	clear_controls();
	lp_model_free(m);
	return 0;
}

int main(void)
{
	if (singleton_root() || production() || experimental_snapshot() ||
		kernel_values() || scaling_isolation())
		return 1;
	clear_controls();
	puts("production ignores numerical environment; explicit experiments use a solve-local snapshot");
	return 0;
}
