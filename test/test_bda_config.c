/**
 * @file test_bda_config.c
 * @brief BDA model/implementation controls are independent of algorithm
 * identity.
 */

#include "simplex_auto.h"
#include "simplex_pan.h"
#include "simplex_pan_basis.h"
#include "simplex_pan_bounded.h"
#include "simplex_pan_config.h"
#include <lp_simplex/lp_simplex.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#define CHECK(c) \
	do { \
		if (!(c)) { \
			fprintf(stderr, "BDA config line %d: %s\n", __LINE__, \
				#c); \
			return 1; \
		} \
	} while (0)

static int configurations(void)
{
	struct simplex_PanConfig c;
	int model, method, oracle, valid;
	for (model = -1; model <= 2; model++)
		for (method = -1; method <= 2; method++)
			for (oracle = -1; oracle <= 3; oracle++) {
				valid = (model == SIMPLEX_PAN_STANDARD ||
						model ==
							SIMPLEX_PAN_NATIVE_BOUNDS) &&
					((model == SIMPLEX_PAN_STANDARD &&
						 method ==
							 SIMPLEX_PAN_BASIS_DIRECTION &&
						 oracle ==
							 SIMPLEX_PAN_NO_ORACLE) ||
						(method == SIMPLEX_PAN_CONE_DIRECTION &&
							(oracle == SIMPLEX_PAN_EXPLICIT_ORACLE ||
								oracle ==
									SIMPLEX_PAN_IMPLICIT_ORACLE)));
				CHECK(simplex_pan_config_default(&c, model,
					      method, oracle) == valid);
				CHECK(simplex_pan_config_valid(&c) == valid);
				CHECK(!c.allow_oracle_recovery);
				CHECK(c.certificate_dual_repair);
			}
	CHECK(simplex_pan_config_default(&c, SIMPLEX_PAN_NATIVE_BOUNDS,
		SIMPLEX_PAN_CONE_DIRECTION, SIMPLEX_PAN_IMPLICIT_ORACLE));
	c.pdhg_weight_set = 1;
	c.pdhg_weight = NAN;
	CHECK(!simplex_pan_config_valid(&c));
	c.pdhg_weight_set = 0;
	c.first_order_iterations = -1;
	CHECK(!simplex_pan_config_valid(&c));
	return 0;
}

static int solve_plans(int mixed)
{
	struct lp_Model *m = lp_model_create(1, 2);
	struct simplex_Problem p;
	struct simplex_PanStandard a;
	struct simplex_PanBasis basis;
	struct simplex_PanConfig c, saved;
	struct lp_simplex_Options o, options_saved;
	struct lp_simplex_Result r;
	double seed[2], x[2];
	int route, native, oracle;

	CHECK(m);
	m->constraints[0].type = optm_CONS_T_EQ;
	m->constraints[0].rhs = 1.;
	m->constraints[0].coef[0] = m->constraints[0].coef[1] = 1.;
	m->objective[1] = 1.;
	seed[0] = 0.;
	seed[1] = 1.;
	if (mixed) {
		m->bounds[0].b_type = optm_BOUND_T_FR;
		m->bounds[1].b_type = optm_BOUND_T_BS;
		m->bounds[1].lb = 0.;
		m->bounds[1].ub = 1.;
		m->objective[0] = 1.;
		m->objective[1] = 0.;
		seed[0] = 1.;
		seed[1] = 0.;
	}
	CHECK(!simplex_problem_from_model(&p, m));
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_PAN_BDA);
	o.iteration_limit = 4096;
	options_saved = o;
	for (route = SIMPLEX_AUTO_HYBRID;
		route <= SIMPLEX_AUTO_STANDARD_EXPLICIT; route++) {
		native = route == SIMPLEX_AUTO_HYBRID ||
			 route == SIMPLEX_AUTO_EXPLICIT;
		oracle =
			route == SIMPLEX_AUTO_CLASSIC
				? SIMPLEX_PAN_NO_ORACLE
				: (route == SIMPLEX_AUTO_HYBRID ||
							  route ==
								  SIMPLEX_AUTO_STANDARD_IMPLICIT
						  ? SIMPLEX_PAN_IMPLICIT_ORACLE
						  : SIMPLEX_PAN_EXPLICIT_ORACLE);
		CHECK(simplex_pan_config_default(&c,
			native ? SIMPLEX_PAN_NATIVE_BOUNDS
			       : SIMPLEX_PAN_STANDARD,
			route == SIMPLEX_AUTO_CLASSIC
				? SIMPLEX_PAN_BASIS_DIRECTION
				: SIMPLEX_PAN_CONE_DIRECTION,
			oracle));
		saved = c;
		memset(&r, 0, sizeof(r));
		CHECK(!simplex_pan_solve_configured(&p, &o, seed, x, &r, &c));
		CHECK(r.status == lp_simplex_Success &&
			fabs(r.objective) <= 1e-7);
		CHECK(simplex_problem_primal_infeasibility(&p, x) <=
			o.primal_tolerance);
		CHECK(r.iterations >= 0 && r.iterations <= o.iteration_limit);
		CHECK(!memcmp(&c, &saved, sizeof(c)) &&
			!memcmp(&o, &options_saved, sizeof(o)));
		printf("configured mixed=%d plan=%d work=%d objective=%.12g\n",
			mixed, route, r.iterations, r.objective);
	}
	/* QR/NNLS construction consumes the supplied controls, even when the
	 * process has conflicting legacy experimental settings. */
	CHECK(!simplex_pan_standard_create(&a, &p));
	c.implicit_updates = 0;
	c.nnls_block_limit = 2;
	c.full_refinement = 1;
	CHECK(!simplex_pan_basis_create_configured(
		&basis, &a, o.pivot_tolerance, &c));
	CHECK(!basis.prefer_updates && basis.nnls_block_limit == 2 &&
		basis.full_refinement);
	simplex_pan_basis_destroy(&basis);
	simplex_pan_standard_destroy(&a);
	/* Reject inconsistent identity/configuration before touching the point. */
	x[0] = x[1] = 123.;
	c.model_form = -1;
	CHECK(simplex_pan_solve_configured(&p, &o, seed, x, &r, &c) != 0);
	CHECK(r.status == lp_simplex_CondUnsatisfied && !r.iterations &&
		x[0] == 123.);
	CHECK(simplex_pan_config_default(&c, SIMPLEX_PAN_STANDARD,
		SIMPLEX_PAN_BASIS_DIRECTION, SIMPLEX_PAN_NO_ORACLE));
	o.algorithm = lp_simplex_ALGORITHM_AUTO;
	CHECK(simplex_pan_solve_configured(&p, &o, seed, x, &r, &c) != 0);
	CHECK(r.status == lp_simplex_CondUnsatisfied && x[0] == 123.);
	/* Removed private route aliases are also invalid public algorithms. */
	for (route = -102; route <= -100; route++) {
		o.algorithm = route;
		CHECK(lp_simplex_solve(m, &o, x, &r) != 0);
		CHECK(r.status == lp_simplex_CondUnsatisfied && !r.iterations &&
			x[0] == 123.);
	}
	simplex_problem_destroy(&p);
	lp_model_free(m);
	return 0;
}

int main(void)
{
	return configurations() || solve_plans(0) || solve_plans(1);
}
