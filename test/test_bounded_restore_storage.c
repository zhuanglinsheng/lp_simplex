/**
 * @file test_bounded_restore_storage.c
 * @brief Fault-inject restore workspaces, including the borrowed face hint.
 */

#include "simplex_pan_bounded.h"
#include "simplex_pan_config.h"
#include <lp_simplex/lp_simplex.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

static void *live[32];
static int requests, fail_at, active, invalid_free;

void restore_test_destroy(struct simplex_PanBounded *b);
int restore_test_repair_dual(const struct simplex_PanBounded *b,
	const double *x, double *y, double tolerance, int limit,
	int *iterations, double *violation);

void *restore_test_malloc(size_t bytes)
{
	void *p;
	int i;
	requests++;
	if (requests == fail_at)
		return NULL;
	p = malloc(bytes ? bytes : 1);
	if (!p)
		return NULL;
	for (i = 0; i < 32; i++)
		if (!live[i]) {
			live[i] = p;
			active++;
			return p;
		}
	free(p);
	return NULL;
}

void restore_test_free(void *p)
{
	int i;
	if (!p)
		return;
	for (i = 0; i < 32; i++)
		if (live[i] == p) {
			live[i] = NULL;
			active--;
			free(p);
			return;
		}
	invalid_free++;
}

#define lp_simplex_malloc restore_test_malloc
#define lp_simplex_free restore_test_free
#define simplex_pan_bounded_create restore_test_create
#define simplex_pan_bounded_destroy restore_test_destroy
#define simplex_pan_bounded_primal restore_test_primal
#define simplex_pan_bounded_encode restore_test_encode
#define simplex_pan_bounded_check_dual restore_test_check_dual
#define simplex_pan_bounded_repair_dual restore_test_repair_dual
#define simplex_pan_bounded_repair_dual_sparse restore_test_repair_dual_sparse
#define simplex_pan_bounded_restore_configured_reported restore_test_run
#define simplex_pan_bounded_solve_certified_configured restore_test_solve
#include "../src/bda/simplex_pan_bounded.c"
#undef lp_simplex_malloc
#undef lp_simplex_free
#undef simplex_pan_bounded_create
#undef simplex_pan_bounded_destroy
#undef simplex_pan_bounded_primal

static int restore_case(int scenario)
{
	struct lp_Model *m = lp_model_create(2, 3);
	struct simplex_Problem p;
	struct simplex_PanBounded b;
	struct simplex_PanConfig controls;
	struct lp_simplex_Options options;
	double saved[3], x[3], residual[2], work[2];
	int j, fault, iterations, result, good = 1;
	if (!m)
		return 0;
	m->constraints[0].type = m->constraints[1].type = optm_CONS_T_EQ;
	m->constraints[0].rhs = .5;
	m->constraints[1].rhs = 1.;
	m->constraints[0].coef[0] = 1.;
	m->constraints[1].coef[1] = m->constraints[1].coef[2] = 1.;
	for (j = 0; j < 3; j++) {
		m->bounds[j].b_type = optm_BOUND_T_BS;
		m->bounds[j].ub = scenario == 2 && j ? (j == 1 ? .4 : .7) : 1.;
		m->objective[j] = j ? 1. : 0.;
		saved[j] = j == 0 && scenario != 0 ? .5 : 0.;
	}
	if (simplex_problem_from_model(&p, m) ||
		simplex_pan_bounded_create(&b, &p))
		return 0;
	lp_simplex_default_options(&options, lp_simplex_ALGORITHM_PAN_BDA);
	options.iteration_limit = 128;
	simplex_pan_config_resolve(&controls, b.problem, &options);
	/* Includes known-empty and nonempty free faces, successful cone repair,
	 * failed one-round rollback, and every root allocation in both methods.
	 * Refusing an affine workspace may legitimately recover through the cone. */
	for (fault = 0; fault < 12; fault++) {
		memcpy(x, saved, sizeof(x));
		requests = active = invalid_free = 0;
		fail_at = fault;
		iterations = 0;
		result = restore_test_run(
			&b, x, &options, &iterations, 2., &controls, 1, NULL);
		good = good && !active && !invalid_free && iterations >= 0 &&
		       iterations <= options.iteration_limit;
		if (result == 1)
			good = good && simplex_pan_bounded_primal(
					       &b, x, residual, work) <= 1e-7;
		else
			good = good && memcmp(x, saved, sizeof(x)) == 0;
		if (!good) {
			fprintf(stderr,
				"restore storage scenario=%d fault=%d result=%d live=%d invalid-free=%d\n",
				scenario, fault, result, active, invalid_free);
			break;
		}
	}
	simplex_pan_bounded_destroy(&b);
	simplex_problem_destroy(&p);
	lp_model_free(m);
	return good;
}

static int recovery_reservation(void)
{
	const int remaining[8] = {0, 1, 2, 127, 2048, 20480, 300000, INT_MAX};
	int i, j;
	for (i = 0; i < 8; i++)
		for (j = 0; j < 8; j++) {
			int grant =
				native_pdhg_budget(remaining[i], remaining[j]);
			int reserve = SIMPLEX_PAN_DUAL_POLISH_WORK_CAP +
				      SIMPLEX_PAN_RESTORE_WORK_CAP;
			if (reserve > remaining[i] / 2)
				reserve = remaining[i] / 2;
			if (grant < 0 || grant > remaining[i] ||
				grant > remaining[j] ||
				remaining[i] - grant < reserve)
				return 0;
		}
	return native_pdhg_budget(300000, 500000) > 150000 &&
	       native_pdhg_budget(-1, INT_MAX) == 0 &&
	       native_pdhg_budget(INT_MAX, -1) == 0;
}

static int iterative_multiplier(void)
{
	const double scales[3] = {1e-9, 1., 1e9};
	int scale, fault, good = 1;
	for (scale = 0; scale < 3; scale++) {
		struct lp_Model *m = lp_model_create(3, 2);
		struct simplex_Problem p;
		struct simplex_PanBounded b;
		double x[2] = {1., 1.}, y[3], saved[3], violation;
		int i, j, used, result;
		if (!m)
			return 0;
		for (i = 0; i < 3; i++) {
			m->constraints[i].type = optm_CONS_T_EQ;
			m->constraints[i].rhs = i < 2 ? 3. * scales[scale] : 0.;
			for (j = 0; j < 2; j++)
				m->constraints[i].coef[j] =
					i < 2 ? (j + 1.) * scales[scale] : 0.;
		}
		for (j = 0; j < 2; j++) {
			m->bounds[j].b_type = optm_BOUND_T_FR;
			m->objective[j] = j + 1.;
		}
		if (simplex_problem_from_model(&p, m) ||
			simplex_pan_bounded_create(&b, &p))
			return 0;
		saved[0] = .125 / scales[scale];
		saved[1] = .25 / scales[scale];
		saved[2] = .75;
		memcpy(y, saved, sizeof(y));
		requests = active = invalid_free = fail_at = used = 0;
		result = native_dual_iterative_repair(
			&b, x, y, 1e-7, 4, &used, &violation);
		/* Dependent equations and an unconstrained multiplier: minimum-norm
		 * correction preserves the null-space component, not just the gap. */
		good = good && result && used > 0 && used <= 4 &&
		       fabs(y[0] * scales[scale] + .5625) < 1e-12 &&
		       fabs(y[1] * scales[scale] + .4375) < 1e-12 &&
		       y[2] == saved[2] && !active && !invalid_free;
		/* The public-to-BDA transaction must retain exact rollback under
		 * refusal, allocation failure and inconsistent stationarity. */
		for (fault = 0; fault < 4; fault++) {
			memcpy(y, saved, sizeof(y));
			requests = active = invalid_free = used = 0;
			fail_at = fault;
			b.a.objective[1] = 3.;
			result = restore_test_repair_dual_sparse(
				&b, x, y, 1e-7, 4, &used, &violation);
			good = good && !result && used <= 4 &&
			       memcmp(y, saved, sizeof(y)) == 0 && !active &&
			       !invalid_free;
		}
		fail_at = 0;
		simplex_pan_bounded_destroy(&b);
		simplex_problem_destroy(&p);
		lp_model_free(m);
	}
	return good;
}

static int iterative_upper_bound(void)
{
	struct lp_Model *m = lp_model_create(1, 3);
	struct simplex_Problem p;
	struct simplex_PanBounded b;
	double x[3] = {0., 1., 5.}, y[1] = {2.}, violation;
	int j, used = INT_MAX - 1, result, good;
	if (!m)
		return 0;
	m->constraints[0].type = optm_CONS_T_EQ;
	m->constraints[0].rhs = 1.;
	m->constraints[0].coef[0] = m->constraints[0].coef[1] = 1.;
	m->objective[1] = -1.;
	for (j = 0; j < 3; j++) {
		m->bounds[j].b_type = optm_BOUND_T_BS;
		m->bounds[j].ub = j == 2 ? 5. : 1.;
	}
	m->bounds[2].lb = 5.;
	m->objective[2] = 100000.;
	if (simplex_problem_from_model(&p, m) ||
		simplex_pan_bounded_create(&b, &p))
		return 0;
	requests = active = invalid_free = fail_at = 0;
	result = restore_test_repair_dual_sparse(
		&b, x, y, 1e-7, INT_MAX, &used, &violation);
	/* Upper-bound directions reverse sign; fixed columns must not enter
	 * stationarity. One remaining work unit is sufficient and cannot wrap. */
	good = result && used == INT_MAX && y[0] == 1. && !active &&
	       !invalid_free &&
	       restore_test_check_dual(&b, x, y, 1e-7, &violation);
	y[0] = 2.;
	result = restore_test_repair_dual_sparse(
		&b, x, y, 1e-7, INT_MAX, &used, &violation);
	good = good && !result && used == INT_MAX && y[0] == 2. && !active &&
	       !invalid_free;
	simplex_pan_bounded_destroy(&b);
	simplex_problem_destroy(&p);
	lp_model_free(m);
	return good;
}

int main(void)
{
	int scenario;
	if (!recovery_reservation() || !iterative_multiplier() ||
		!iterative_upper_bound())
		return 1;
	for (scenario = 0; scenario < 3; scenario++)
		if (!restore_case(scenario))
			return 1;
	puts("Restore workspaces preserve ownership, rollback and shared work");
	return 0;
}
