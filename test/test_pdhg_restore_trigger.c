/**
 * @file test_pdhg_restore_trigger.c
 * @brief Exercise the PDHG repair producer with a short test-only stall window.
 */

#include "simplex_pan_service.h"
#include "simplex_pan_bounded.h"
#include "simplex_pan_config.h"
#include "simplex_pan_pdhg.h"
#include "simplex_stage.h"
#include "utils.h"
#include "policy/pan.h"
#include <lp_simplex/lp_simplex.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

static int calls, invalid;
static double last_primal;

static int producer_restore(const struct simplex_PanBounded *b, double *point,
	const struct lp_simplex_Options *o, int *iterations,
	long double ceiling, const struct simplex_PanConfig *cfg,
	const struct simplex_PanServiceTrigger *trigger,
	struct simplex_PanServiceOutcome *out);

int producer_pdhg_configured(const struct simplex_PanBounded *b, double *point,
	double *multiplier, const struct lp_simplex_Options *o, int *iterations,
	int budget, const struct simplex_PanConfig *cfg);

/* Only observation timing changes. The real eligibility predicate, PDHG
 * iterations, trigger producer and service wrapper are compiled/executed. */
#undef SIMPLEX_PDHG_AUDIT_PERIOD
#define SIMPLEX_PDHG_AUDIT_PERIOD 1
#undef SIMPLEX_PDHG_OBSERVE_WORK
#define SIMPLEX_PDHG_OBSERVE_WORK 0
#undef SIMPLEX_PDHG_STAGNANT_WORK
#define SIMPLEX_PDHG_STAGNANT_WORK 0
#define simplex_pan_pdhg_progress_reset producer_progress_reset
#define simplex_pan_pdhg_repair_progress producer_repair_progress
#define simplex_pan_pdhg_progress_stalled producer_progress_stalled
#define simplex_pan_pdhg_progress_stalled_pair producer_progress_stalled_pair
#define simplex_pan_pdhg_handoff_gap producer_handoff_gap
#define simplex_pan_pdhg_columns_destroy producer_columns_destroy
#define simplex_pan_pdhg_columns_create producer_columns_create
#define simplex_pan_pdhg_residual producer_residual
#define simplex_pan_bounded_pdhg_configured_reported producer_pdhg_reported
#define simplex_pan_bounded_pdhg producer_pdhg
#define simplex_pan_bounded_pdhg_configured producer_pdhg_configured
#define simplex_pan_service_restore producer_restore
#include "../src/bda/simplex_pan_pdhg.c"
#undef simplex_pan_service_restore

static int producer_restore(const struct simplex_PanBounded *b, double *point,
	const struct lp_simplex_Options *o, int *iterations,
	long double ceiling, const struct simplex_PanConfig *cfg,
	const struct simplex_PanServiceTrigger *trigger,
	struct simplex_PanServiceOutcome *out)
{
	struct lp_simplex_Options exhausted = *o;
	struct simplex_PanServiceOutcome actual;
	double residual[1], work[1];
	double primal = simplex_pan_bounded_primal(b, point, residual, work);
	int code;

	if (!trigger || trigger->kind != SIMPLEX_PAN_TRIGGER_PDHG_STAGNATION ||
		trigger->work != *iterations || trigger->primal != primal ||
		(calls && trigger->previous_primal != last_primal) ||
		(!calls && !isinf(trigger->previous_primal)) ||
		!producer_repair_progress(
			primal, trigger->previous_primal, o->primal_tolerance))
		invalid = 1;
	last_primal = primal;
	calls++;
	/* Force a failed zero-work repair through the real service boundary.
	 * Subsequent attempts must still improve the producer's residual baseline. */
	exhausted.iteration_limit = *iterations;
	code = simplex_pan_service_restore(b, point, &exhausted, iterations,
		ceiling, cfg, trigger, &actual);
	if (code == 1 || !actual.valid || actual.output)
		invalid = 1;
	if (out)
		*out = actual;
	return code;
}

/* A tiny coordinate is removable only when its aggregate row activity is
 * harmless. The exact same coordinate with a large coefficient is essential.
 * Neither a face adjustment nor a feasible point alone grants KKT status. */
static int face_certificate(void)
{
	int essential;
	for (essential = 0; essential < 2; essential++) {
		struct lp_Model *m = lp_model_create(1, 2);
		struct simplex_Problem p;
		struct simplex_PanBounded b;
		struct lp_simplex_Options o;
		double point[2] = {1e-10, 1. - 1e-10}, saved[2];
		double multiplier[1] = {0.}, physical[1], r[1], work[1];
		int improvements = 0, iterations = 17, changed, certified;
		if (!m)
			return 0;
		m->constraints[0].type = optm_CONS_T_EQ;
		m->constraints[0].rhs = essential ? 2. : 1.;
		m->constraints[0].coef[0] = essential ? 1e10 : 1.;
		m->constraints[0].coef[1] = 1.;
		m->objective[0] = 1.;
		if (essential)
			point[1] = 1.;
		if (simplex_problem_from_model(&p, m) ||
			simplex_pan_bounded_create(&b, &p))
			return 0;
		lp_simplex_default_options(&o, lp_simplex_ALGORITHM_PAN_BDA);
		o.iteration_limit = 128;
		certified = pdhg_audit(&b, point, multiplier, 1., point,
			multiplier, &o, physical, r, work, &improvements, NULL);
		if (certified)
			return 0;
		changed = simplex_pan_service_face(NULL, &b, point, multiplier,
			&o, &iterations, saved, r, work,
			SIMPLEX_PAN_FACE_PDHG_IMPROVEMENT, 1, 0, NULL, NULL);
		certified = pdhg_audit(&b, point, multiplier, 1., point,
			multiplier, &o, physical, r, work, &improvements, NULL);
		if (iterations != 17 ||
			(essential ? (changed != 0 || certified ||
					     point[0] != 1e-10)
				   : (changed <= 0 || !certified ||
					     point[0] != 0.)))
			return 0;
		simplex_pan_bounded_destroy(&b);
		simplex_problem_destroy(&p);
		lp_model_free(m);
	}
	return 1;
}

static int cancelling_face_certificate(void)
{
	struct lp_Model *m = lp_model_create(1, 3);
	struct simplex_Problem p;
	struct simplex_PanBounded b;
	struct lp_simplex_Options o;
	double point[3] = {1e-6, 1e-6, 100.}, saved[3];
	double multiplier[1] = {0.}, physical[1], r[1], work[1];
	int improvements = 0, iterations = 17, changed, certified;
	if (!m)
		return 0;
	m->constraints[0].type = optm_CONS_T_EQ;
	m->constraints[0].rhs = 100.;
	m->constraints[0].coef[0] = m->constraints[0].coef[2] = 1.;
	m->constraints[0].coef[1] = -1.;
	m->objective[0] = m->objective[1] = 1.;
	m->objective[2] = -10000.;
	m->bounds[2].b_type = optm_BOUND_T_BS;
	m->bounds[2].ub = 100.;
	if (simplex_problem_from_model(&p, m) ||
		simplex_pan_bounded_create(&b, &p))
		return 0;
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_PAN_BDA);
	o.iteration_limit = 128;
	if (pdhg_audit(&b, point, multiplier, 1., point, multiplier, &o,
		    physical, r, work, &improvements, NULL))
		return 0;
	changed = simplex_pan_service_face(NULL, &b, point, multiplier, &o,
		&iterations, saved, r, work, SIMPLEX_PAN_FACE_PDHG_IMPROVEMENT,
		1, 0, NULL, NULL);
	certified = pdhg_audit(&b, point, multiplier, 1., point, multiplier, &o,
		physical, r, work, &improvements, NULL);
	simplex_pan_bounded_destroy(&b);
	simplex_problem_destroy(&p);
	lp_model_free(m);
	return changed == 2 && certified && iterations == 17 &&
	       point[0] == 0. && point[1] == 0. && point[2] == 100.;
}

static int fresh_feasible_certificate(void)
{
	struct lp_Model *m = lp_model_create(1, 2);
	struct simplex_Problem p;
	struct simplex_PanBounded b;
	struct lp_simplex_Options o;
	double candidate[2] = {1., 0.}, old[2] = {1. - 3e-8, 1e-8};
	double point[2], multiplier[1] = {0.}, seed[1] = {0.};
	double physical[1], r[1], work[1];
	int improved = 0, used = 0, certified;
	if (!m)
		return 0;
	m->constraints[0].type = optm_CONS_T_EQ;
	m->constraints[0].rhs = 1.;
	m->constraints[0].coef[0] = m->constraints[0].coef[1] = 1.;
	m->objective[0] = 1.;
	m->objective[1] = 2.;
	m->bounds[0].b_type = optm_BOUND_T_BS;
	m->bounds[0].ub = 1.;
	if (simplex_problem_from_model(&p, m) ||
		simplex_pan_bounded_create(&b, &p))
		return 0;
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_PAN_BDA);
	memcpy(point, old, sizeof(point));
	/* The incumbent is tolerance-feasible and cheaper, but its two interior
	 * variables impose incompatible stationarity equations. Zero-work refusal
	 * must retain it. A fresh exact vertex can replace it only with full KKT. */
	certified = pdhg_refine_candidate(&b, candidate, seed, 1., point,
		multiplier, &o, 0, &used, physical, r, work, &improved);
	if (certified || memcmp(point, old, sizeof(point)) || used)
		return 0;
	certified = pdhg_refine_candidate(&b, candidate, seed, 1., point,
		multiplier, &o, 16, &used, physical, r, work, &improved);
	simplex_pan_bounded_destroy(&b);
	simplex_problem_destroy(&p);
	lp_model_free(m);
	return certified && improved == 1 && used > 0 && used <= 16 &&
	       point[0] == 1. && point[1] == 0. && candidate[0] == 1. &&
	       candidate[1] == 0.;
}

int main(void)
{
	struct lp_Model *m = lp_model_create(1, 2);
	struct simplex_Problem p;
	struct simplex_PanBounded b;
	struct simplex_PanConfig cfg;
	struct simplex_PanServiceContext c;
	struct simplex_PanServiceFrame outer;
	struct simplex_SolveContext context;
	struct simplex_Execution history;
	struct simplex_StageReport report;
	struct lp_simplex_Options o;
	double point[3], multiplier[1], saved[3];
	int mode, iterations, saved_work = 0, saved_calls = 0, stop, output,
			      code, i, triggers;

	if (!face_certificate() || !cancelling_face_certificate() ||
		!fresh_feasible_certificate())
		return 1;
	if (!m)
		return 1;
	m->constraints[0].type = optm_CONS_T_GE;
	m->constraints[0].rhs = 1.;
	m->constraints[0].coef[0] = m->constraints[0].coef[1] = 1.;
	m->objective[1] = 1.;
	if (simplex_problem_from_model(&p, m) ||
		simplex_pan_bounded_create(&b, &p))
		return 2;
	p.report = &report;
	if (!simplex_pan_config_default(&cfg, SIMPLEX_PAN_NATIVE_BOUNDS,
		    SIMPLEX_PAN_CONE_DIRECTION, SIMPLEX_PAN_IMPLICIT_ORACLE))
		return 3;
	cfg.restore = 1;
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_PAN_BDA);
	o.iteration_limit = 128;
	for (mode = 0; mode < 3; mode++) {
		simplex_stage_reset(&report);
		simplex_solve_context_init(&context, &history, &p.execution);
		if (mode == 1)
			p.execution.history = NULL;
		simplex_pan_service_context_init(&c, &p);
		simplex_pan_service_bind(&b.a, &c);
		if (mode == 2)
			history.count = SIMPLEX_EXECUTION_CAPACITY;
		iterations = calls = invalid = 0;
		last_primal = HUGE_VAL;
		point[0] = point[1] = .5;
		point[2] = multiplier[0] = 0.;
		if (!simplex_pan_service_begin(&outer, &b.a,
			    SIMPLEX_PAN_SERVICE_PDHG, 0, o.iteration_limit,
			    o.iteration_limit))
			return 4;
		code = producer_pdhg_reported(&b, point, multiplier, &o,
			&iterations, o.iteration_limit, &cfg, &stop, &output);
		if (!calls || invalid || c.failed ||
			!simplex_pan_service_finish(&outer, &iterations, code,
				code == 1 ? lp_simplex_Success
					  : lp_simplex_CondUnsatisfied,
				stop, output, NULL))
			return 5;
		if (!mode) {
			saved_work = iterations;
			saved_calls = calls;
			memcpy(saved, point, sizeof(saved));
			triggers = 0;
			for (i = 0; i < history.count; i++)
				triggers += history.events[i].kind ==
					    SIMPLEX_EXECUTION_PAN_TRIGGER;
			if (triggers != calls)
				return 6;
		} else if (iterations != saved_work || calls != saved_calls ||
			   memcmp(saved, point, sizeof(saved)))
			return 7;
		simplex_solve_context_destroy(&context);
	}
	printf("PDHG producer: %d failed repairs, %d charged work; journal-independent\n",
		saved_calls, saved_work);
	simplex_pan_bounded_destroy(&b);
	simplex_problem_destroy(&p);
	lp_model_free(m);
	return 0;
}
