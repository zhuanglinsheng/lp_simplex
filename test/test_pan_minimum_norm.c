/**
 * @file test_pan_minimum_norm.c
 * @brief Numerical backend equations and null-space orthogonality, not LP
 * heuristics.
 */

#include "simplex_pan_basis.h"
#include "simplex_pan_config.h"
#include "simplex_pan.h"
#include "simplex_stage.h"
#include <lp_simplex/lp_simplex.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#define CHECK(c) \
	do { \
		if (!(c)) { \
			fprintf(stderr, "minimum norm line %d: %s\n", \
				__LINE__, #c); \
			return 1; \
		} \
	} while (0)

static int rectangular(void)
{
	struct simplex_PanConfig basis_controls;
	struct simplex_PanStandard a;
	struct simplex_PanBasis b;
	int start[3] = {0, 1, 4}, row[4] = {0, 0, 1, 2};
	double value[4], normalized[4], norms[2], right[2], saved[2], y[3],
		error;
	double epsilon = ldexp(1., -18), first, second;
	int exponent, order, i, selected[2];

	simplex_pan_config_init(&basis_controls, NULL);
	memset(&a, 0, sizeof(a));
	a.rows = 3;
	a.columns = 2;
	a.nonzeros = 4;
	a.column_start = start;
	a.row_index = row;
	a.value = value;
	a.normalized_value = normalized;
	a.column_norm = norms;
	for (exponent = -40; exponent <= 40; exponent += 10)
		for (order = 0; order < 2; order++) {
			first = ldexp(1., exponent);
			second = ldexp(order ? -1. : 1., -exponent);
			value[0] = first;
			value[1] = second;
			value[2] = value[3] = epsilon * second;
			norms[0] = fabs(first);
			norms[1] = hypot(second, hypot(value[2], value[3]));
			for (i = 0; i < 4; i++)
				normalized[i] = value[i] / norms[i ? 1 : 0];
			CHECK(!simplex_pan_basis_create_configured(
				&b, &a, 1e-12, &basis_controls));
			CHECK(!simplex_pan_basis_add(&b, order));
			CHECK(!simplex_pan_basis_add(&b, 1 - order));
			CHECK(!simplex_pan_basis_rebuild(&b));
			CHECK(b.count == 2);
			for (i = 0; i < 2; i++) {
				selected[i] = b.column[i];
				right[i] = b.column[i] == 0
						   ? first
						   : second * (1. + epsilon);
			}
			memcpy(saved, right, sizeof(right));
			CHECK(!simplex_pan_basis_minimum_norm_qr(
				&b, right, y, &error));
			CHECK(isfinite(error) && b.minimum_norm_qr_calls == 1);
			CHECK(!memcmp(saved, right, sizeof(right)) &&
				b.count == 2 && b.column[0] == selected[0] &&
				b.column[1] == selected[1]);
			/* Exact constraints imply y0=1 and y1+y2=1. Orthogonality to
			 * null(B^T)=span(0,1,-1) proves the unique minimum norm. */
			CHECK(fabs(y[0] - 1.) < 1e-8 &&
				fabs(y[1] + y[2] - 1.) < 1e-7);
			CHECK(fabs(y[1] - y[2]) < 1e-8);
			CHECK(fabs(y[1] - .5) < 1e-7 && fabs(y[2] - .5) < 1e-7);
			for (i = 0; i < 2; i++) {
				int c = b.column[i], k;
				long double dot = 0.;
				for (k = start[c]; k < start[c + 1]; k++)
					dot += (long double)value[k] *
					       y[row[k]];
				CHECK(fabsl(dot - right[i]) <
					1e-12L * (1. + fabs(right[i])));
			}
			/* Nonfinite original equations cannot be certified by a cached
			 * factor or by a norm routine that silently turns NaN into zero. */
			value[0] = NAN;
			y[0] = y[1] = y[2] = 123.;
			CHECK(simplex_pan_basis_minimum_norm_qr(
				      &b, right, y, &error) != 0);
			CHECK(isinf(error) && y[0] == 123. && y[1] == 123. &&
				y[2] == 123.);
			value[0] = first;
			right[0] = NAN;
			y[0] = y[1] = y[2] = 123.;
			CHECK(simplex_pan_basis_minimum_norm_qr(
				      &b, right, y, &error) != 0);
			CHECK(isinf(error) && y[0] == 123. && y[1] == 123. &&
				y[2] == 123.);
			simplex_pan_basis_destroy(&b);
		}
	return 0;
}

static int rejected_and_empty(void)
{
	struct simplex_PanConfig basis_controls;
	struct simplex_PanStandard a;
	struct simplex_PanBasis b;
	int start[3] = {0, 1, 2}, row[2] = {0, 0};
	double value[2] = {1., 1.}, norms[2] = {1., 1.};
	double right[2] = {1., 2.}, y[3] = {123., 123., 123.}, error;
	double invalid[3] = {NAN, 0., 0.};

	simplex_pan_config_init(&basis_controls, NULL);
	memset(&a, 0, sizeof(a));
	a.rows = 3;
	a.columns = 2;
	a.nonzeros = 2;
	a.column_start = start;
	a.row_index = row;
	a.value = a.normalized_value = value;
	a.column_norm = norms;
	CHECK(!simplex_pan_basis_create_configured(
		&b, &a, 1e-12, &basis_controls));
	CHECK(!simplex_pan_basis_least_squares(&b, invalid, y, &error));
	CHECK(isinf(error));
	CHECK(!simplex_pan_basis_minimum_norm_qr(&b, NULL, y, &error));
	CHECK(error == 0. && y[0] == 0. && y[1] == 0. && y[2] == 0.);
	CHECK(!simplex_pan_basis_add(&b, 0) && !simplex_pan_basis_add(&b, 1));
	y[0] = y[1] = y[2] = 123.;
	/* A singular/unfactored working set must not publish a candidate. */
	CHECK(simplex_pan_basis_minimum_norm_qr(&b, right, y, &error) != 0);
	CHECK(isinf(error) && y[0] == 123. && y[1] == 123. && y[2] == 123.);
	simplex_pan_basis_destroy(&b);
	/* Metadata-only overlarge request: no matrix storage may be accessed. */
	memset(&b, 0, sizeof(b));
	a.rows = INT_MAX;
	b.standard = &a;
	b.count = 1;
	CHECK(simplex_pan_basis_minimum_norm_qr(&b, right, y, &error) != 0);
	CHECK(isinf(error) && y[0] == 123.);
	return 0;
}

static int multiplier_recovery(void)
{
	struct lp_Model *m = lp_model_create(3, 3);
	struct simplex_Problem p;
	struct simplex_PanConfig config;
	struct simplex_StageReport report;
	struct lp_simplex_Options o;
	struct lp_simplex_Result r;
	double x[3], scales[3] = {1048576., .0009765625, 32.};
	int i;

	CHECK(m);
	for (i = 0; i < 3; i++) {
		m->constraints[i].type = optm_CONS_T_GE;
		m->constraints[i].rhs = scales[i];
		m->constraints[i].coef[i] = scales[i];
		m->constraints[i].coef[(i + 1) % 3] = scales[i];
		m->bounds[i].b_type = optm_BOUND_T_BS;
		m->bounds[i].ub = 1.;
		m->objective[i] = 4096.;
	}
	CHECK(!simplex_problem_from_model(&p, m));
	simplex_stage_reset(&report);
	p.report = &report;
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_PAN_BDA);
	CHECK(simplex_pan_config_default(&config, SIMPLEX_PAN_STANDARD,
		SIMPLEX_PAN_BASIS_DIRECTION, SIMPLEX_PAN_NO_ORACLE));
	CHECK(!simplex_pan_solve_configured(&p, &o, NULL, x, &r, &config));
	CHECK(r.status == lp_simplex_Success &&
		fabs(r.objective - 6144.) < 1e-6);
	CHECK(simplex_problem_primal_infeasibility(&p, x) <=
		o.primal_tolerance);
	CHECK(report.work.minimum_norm_qr_available);
#ifdef TEST_DENSE_BACKEND
	CHECK(report.work.minimum_norm_qr_calls > 0);
#endif
	printf("scaled multiplier recovery work=%d qr=%ld objective=%.12g\n",
		r.iterations, report.work.minimum_norm_qr_calls, r.objective);
	simplex_problem_destroy(&p);
	lp_model_free(m);
	return 0;
}

int main(void)
{
	return rectangular() || rejected_and_empty() || multiplier_recovery();
}
