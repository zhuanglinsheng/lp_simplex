/**
 * @file test_pan_services.c
 * @brief Phase leases, identity, typed handoffs and causal nesting within BDA.
 */

#include "simplex_solve_context.h"
#include "simplex_execution.h"
#include "simplex_pan_service.h"
#include "simplex_pan_config.h"
#include "simplex_pan_basis.h"
#include "simplex_pan_bounded.h"
#include "simplex_stage.h"
#include <lp_simplex/lp_simplex.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#define CHECK(c) \
	do { \
		if (!(c)) { \
			fprintf(stderr, "PAN service line %d: %s\n", __LINE__, \
				#c); \
			return 1; \
		} \
	} while (0)

static int leases(void)
{
	struct lp_Model *m = lp_model_create(1, 2);
	struct simplex_Problem p;
	struct simplex_PanStandard a;
	struct simplex_PanServiceContext c;
	struct simplex_PanServiceFrame f, child;
	struct simplex_PanServiceOutcome out;
	struct simplex_StageReport report;
	struct simplex_SolveContext context;
	struct simplex_Execution execution;
	int diagnostics, work, i, outer, raw;
	CHECK(m);
	m->constraints[0].type = optm_CONS_T_EQ;
	m->constraints[0].rhs = 1.;
	m->constraints[0].coef[0] = m->constraints[0].coef[1] = 1.;
	CHECK(!simplex_problem_from_model(&p, m));
	p.report = &report;
	CHECK(!simplex_pan_standard_create(&a, &p));
	for (diagnostics = 0; diagnostics < 3; diagnostics++) {
		simplex_stage_reset(&report);
		simplex_solve_context_init(&context, &execution, &p.execution);
		if (diagnostics == 0)
			p.execution.history = NULL;
		if (diagnostics == 2)
			for (i = 0; i < SIMPLEX_EXECUTION_CAPACITY; i++)
				(void)simplex_execution_begin(&p.execution,
					SIMPLEX_EXECUTION_SOLVE, 0, 0, 0);
		simplex_pan_service_context_init(&c, &p);
		simplex_pan_service_bind(&a, &c);
		work = 5;
		CHECK(simplex_pan_service_begin(
			&f, &a, SIMPLEX_PAN_SERVICE_NNLS, work, 12, 3));
		CHECK(f.outcome.grant == 3);
		work = 7;
		CHECK(simplex_pan_service_finish(&f, &work, 0,
			lp_simplex_CondUnsatisfied, SIMPLEX_STOP_UNKNOWN,
			SIMPLEX_PAN_OUTPUT_NONE, &out));
		CHECK(work == 7 &&
			out.stop ==
				SIMPLEX_STOP_UNKNOWN); /* No inferred yield. */
		CHECK(report.services.calls[SIMPLEX_PAN_SERVICE_NNLS] == 1);
		/* A diagnostic summary is deliberately poisoned; it is never runtime state. */
		report.services.output[SIMPLEX_PAN_SERVICE_NNLS] =
			SIMPLEX_PAN_OUTPUT_NATIVE_KKT;
		work = 10;
		CHECK(simplex_pan_service_begin(
			&f, &a, SIMPLEX_PAN_SERVICE_PDHG, work, 30, 7));
		outer = f.outcome.event;
		CHECK(simplex_pan_service_begin(&child, &a,
			SIMPLEX_PAN_SERVICE_RESTORE, work, 17, 100));
		CHECK(child.outcome.grant == 7);
		if (diagnostics == 1)
			CHECK(execution.events[child.outcome.event - 1]
					.parent == outer);
		work += 2;
		CHECK(simplex_pan_service_finish(&child, &work, 1,
			lp_simplex_CondUnsatisfied, SIMPLEX_STOP_HANDOFF,
			SIMPLEX_PAN_OUTPUT_NATIVE_FEASIBLE_POINT, &out));
		CHECK(c.parent == outer);
		work++;
		CHECK(simplex_pan_service_finish(&f, &work, 0,
			lp_simplex_CondUnsatisfied, SIMPLEX_STOP_HANDOFF,
			SIMPLEX_PAN_OUTPUT_NATIVE_FEASIBLE_POINT, &out));
		CHECK(work == 13 && c.parent == p.execution.parent);
		CHECK(report.services.work[SIMPLEX_PAN_SERVICE_PDHG] == 3 &&
			report.services.work[SIMPLEX_PAN_SERVICE_RESTORE] ==
				2); /* Inclusive. */
		CHECK(!c.failed && out.valid);
		/* Corrupt cumulative work is retained raw but charged to the lease end. */
		for (i = 0; i < 3; i++) {
			simplex_pan_service_context_init(&c, &p);
			simplex_pan_service_bind(&a, &c);
			work = 5;
			CHECK(simplex_pan_service_begin(
				&f, &a, SIMPLEX_PAN_SERVICE_NNLS, work, 12, 3));
			raw = i == 0 ? 4 : i == 1 ? 9 : INT_MIN;
			work = raw;
			CHECK(!simplex_pan_service_finish(&f, &work, 0,
				lp_simplex_CondUnsatisfied,
				SIMPLEX_STOP_HANDOFF,
				SIMPLEX_PAN_OUTPUT_NONNEGATIVE_EQUALITY_POINT,
				&out));
			CHECK(work == 8 && out.reported_end == raw &&
				!out.valid && c.failed);
			CHECK(out.raw_output ==
					SIMPLEX_PAN_OUTPUT_NONNEGATIVE_EQUALITY_POINT &&
				!out.output);
			CHECK(!simplex_pan_service_begin(
				&f, &a, SIMPLEX_PAN_SERVICE_PDHG, work, 12, 4));
			CHECK(!simplex_pan_service_finish(&f, &work, 0,
				lp_simplex_CondUnsatisfied,
				SIMPLEX_STOP_UNKNOWN, SIMPLEX_PAN_OUTPUT_NONE,
				&out));
			CHECK(work == 8 && out.grant == 0);
		}
		/* A role cannot be promoted to KKT, even with well-formed work. */
		simplex_pan_service_context_init(&c, &p);
		simplex_pan_service_bind(&a, &c);
		work = 0;
		CHECK(simplex_pan_service_begin(
			&f, &a, SIMPLEX_PAN_SERVICE_NNLS, work, 10, 10));
		CHECK(!simplex_pan_service_finish(&f, &work, 0,
			lp_simplex_Success, SIMPLEX_STOP_OPTIMAL,
			SIMPLEX_PAN_OUTPUT_NATIVE_KKT, &out));
		CHECK(!out.output && out.raw_status == lp_simplex_Success);
		/* Both immutable source identity and represented-input generation matter. */
		for (i = 0; i < 2; i++) {
			simplex_pan_service_context_init(&c, &p);
			simplex_pan_service_bind(&a, &c);
			work = 1;
			CHECK(simplex_pan_service_begin(
				&f, &a, SIMPLEX_PAN_SERVICE_NNLS, work, 10, 2));
			if (i == 0)
				simplex_problem_mark_changed(&p);
			else
				simplex_pan_service_derive(&a, &a);
			work++;
			CHECK(!simplex_pan_service_finish(&f, &work, 0,
				lp_simplex_CondUnsatisfied,
				SIMPLEX_STOP_HANDOFF,
				SIMPLEX_PAN_OUTPUT_NONNEGATIVE_EQUALITY_POINT,
				&out));
			CHECK(work == 2 && !out.output && c.failed);
		}
		simplex_solve_context_destroy(&context);
	}
	/* Local-only contexts must still reject a sticky contract failure. */
	simplex_pan_service_context_init(&c, NULL);
	CHECK(simplex_pan_service_context_valid(&c));
	c.failed = 1;
	CHECK(!simplex_pan_service_context_valid(&c));
	simplex_pan_service_context_init(&c, &p);
	c.representations = ULONG_MAX;
	simplex_pan_service_bind(&a, &c);
	CHECK(!simplex_pan_service_context_valid(&c));
	simplex_pan_standard_destroy(&a);
	simplex_problem_destroy(&p);
	lp_model_free(m);
	return 0;
}

static int roles(void)
{
	struct lp_Model *m = lp_model_create(1, 2);
	struct simplex_Problem p;
	struct simplex_PanBounded b;
	struct simplex_PanBasis basis;
	struct simplex_PanStandard aux;
	struct simplex_PanServiceContext c;
	struct simplex_PanServiceOutcome out;
	struct simplex_SolveContext context;
	struct simplex_Execution execution;
	struct simplex_StageReport report;
	struct simplex_PanConfig cfg;
	struct lp_simplex_Options o;
	struct lp_simplex_Result result;
	double point[3], y[1] = {0.}, root[2], violation = HUGE_VAL;
	int work, status, code, before, i, services = 0;
	int starts[2] = {0, 1}, rows[1] = {0};
	unsigned char free_mask[1] = {0};
	double entry[1] = {1.}, norm[1] = {1.}, rhs[2] = {0., -1.},
	       cost[1] = {0.}, d[1], residual[2];
	CHECK(m);
	m->constraints[0].type = optm_CONS_T_GE;
	m->constraints[0].rhs = 1.;
	m->constraints[0].coef[0] = m->constraints[0].coef[1] = 1.;
	m->objective[1] = 1.;
	CHECK(!simplex_problem_from_model(&p, m));
	simplex_stage_reset(&report);
	p.report = &report;
	simplex_solve_context_init(&context, &execution, &p.execution);
	simplex_pan_service_context_init(&c, &p);
	CHECK(!simplex_pan_bounded_create(&b, &p));
	simplex_pan_service_bind(&b.a, &c);
	CHECK(simplex_pan_config_default(&cfg, SIMPLEX_PAN_NATIVE_BOUNDS,
		SIMPLEX_PAN_CONE_DIRECTION, SIMPLEX_PAN_IMPLICIT_ORACLE));
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_PAN_BDA);
	o.iteration_limit = 128;
	CHECK(!simplex_pan_basis_create_configured(
		&basis, &b.a, o.pivot_tolerance, &cfg));
	work = 3;
	status = lp_simplex_CondUnsatisfied;
	CHECK(!simplex_pan_service_nnls(&basis, point, o.primal_tolerance,
		o.dual_tolerance, 32, &work, &status, &out));
	CHECK(out.valid &&
		out.output == SIMPLEX_PAN_OUTPUT_NONNEGATIVE_EQUALITY_POINT);
	CHECK(out.stop == SIMPLEX_STOP_HANDOFF &&
		out.raw_status == lp_simplex_CondUnsatisfied);
	simplex_pan_basis_destroy(&basis);
	before = work;
	CHECK(!simplex_pan_service_dual_feasible(
		&c, &p, &o, root, &work, &result, &out));
	CHECK(out.output == SIMPLEX_PAN_OUTPUT_INPUT_FEASIBLE_POINT &&
		!out.representation);
	CHECK(out.columns == p.columns && out.relative_work &&
		result.iterations == work - before);
	CHECK(simplex_problem_primal_infeasibility(&p, root) <=
		o.primal_tolerance);
	before = work;
	CHECK(!simplex_pan_service_pdhg(
		&b, point, y, &o, &work, 0, &cfg, &out));
	CHECK(out.stop == SIMPLEX_STOP_BUDGET && !out.output && work == before);
	code = simplex_pan_service_pdhg(
		&b, point, y, &o, &work, 16, &cfg, &out);
	CHECK(out.valid && work <= before + 16);
	CHECK(code != 1 || out.output == SIMPLEX_PAN_OUTPUT_NATIVE_KKT);
	point[0] = .6;
	point[1] = .3;
	point[2] = 0.;
	CHECK(simplex_pan_service_restore(
		      &b, point, &o, &work, 2., &cfg, NULL, &out) == 1);
	CHECK(out.valid &&
		out.output == SIMPLEX_PAN_OUTPUT_NATIVE_FEASIBLE_POINT);
	cfg.dual_polish_iterations = 0;
	before = work;
	CHECK(!simplex_pan_service_first_order(&b.a, point, y, &o, &work,
		&violation, NULL, NULL, 0, &cfg, &out));
	CHECK(out.stop == SIMPLEX_STOP_BUDGET &&
		out.service == SIMPLEX_PAN_SERVICE_DUAL_REFINE &&
		work == before);
	/* An auxiliary NNLS stationary residual does not classify the root LP. */
	memset(&aux, 0, sizeof(aux));
	aux.rows = 2;
	aux.columns = 1;
	aux.nonzeros = 1;
	aux.column_start = starts;
	aux.row_index = rows;
	aux.value = aux.normalized_value = entry;
	aux.column_norm = norm;
	aux.rhs = rhs;
	aux.objective = cost;
	aux.report = &report;
	simplex_pan_service_bind(&aux, &c);
	CHECK(!simplex_pan_basis_create_configured(
		&basis, &aux, o.pivot_tolerance, &cfg));
	status = lp_simplex_CondUnsatisfied;
	CHECK(simplex_pan_service_cone(&basis, 0, d, free_mask, residual, 1e-9,
		      o.iteration_limit, &work, &status, NULL, NULL,
		      &out) == 1);
	CHECK(out.valid &&
		out.output == SIMPLEX_PAN_OUTPUT_STATIONARY_RESIDUAL &&
		out.stop == SIMPLEX_STOP_HANDOFF);
	CHECK(out.implementation == SIMPLEX_PAN_EXPLICIT_ORACLE &&
		status != lp_simplex_Infeasibility);
	CHECK(simplex_problem_primal_infeasibility(&p, root) <=
			o.primal_tolerance &&
		!c.failed);
	for (i = 0; i < execution.count; i++) {
		struct simplex_ExecutionEvent *e = execution.events + i;
		CHECK(e->complete && e->model > 0 && e->parent <= i);
		if (e->kind == SIMPLEX_EXECUTION_PAN_SERVICE) {
			services++;
			CHECK(e->service_valid &&
				e->source_model == p.execution.model);
			CHECK(e->work >= 0 && e->work <= e->grant);
		}
	}
	CHECK(services >= 6);
	simplex_pan_basis_destroy(&basis);
	simplex_pan_bounded_destroy(&b);
	simplex_problem_destroy(&p);
	lp_model_free(m);
	simplex_solve_context_destroy(&context);
	return 0;
}

/* Identity replacement is checked without dereferencing substituted arrays. */
static int identities(void)
{
	struct lp_Model *m = lp_model_create(1, 1);
	struct simplex_Problem p, saved, copy;
	struct simplex_PanStandard a, original;
	struct simplex_PanBounded b, saved_b;
	struct simplex_PanServiceContext c;
	struct simplex_PanServiceFrame f;
	struct simplex_PanServiceOutcome out;
	struct simplex_SolveContext context;
	struct simplex_Execution history;
	struct simplex_StageReport other_report;
	struct simplex_PanConfig cfg;
	struct lp_simplex_Options o;
	double replacement = 2., point[2] = {1., 0.}, multiplier[1] = {-1.},
	       violation;
	int diagnostics, phase, field, work, i;
	CHECK(m);
	m->constraints[0].type = optm_CONS_T_GE;
	m->constraints[0].rhs = m->constraints[0].coef[0] = m->objective[0] =
		1.;
	CHECK(!simplex_problem_from_model(&p, m));
	CHECK(!simplex_pan_standard_create(&a, &p));
	CHECK(!simplex_pan_bounded_create(&b, &p));
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_PAN_BDA);
	o.iteration_limit = 16;
	CHECK(simplex_pan_config_default(&cfg, SIMPLEX_PAN_NATIVE_BOUNDS,
		SIMPLEX_PAN_CONE_DIRECTION, SIMPLEX_PAN_IMPLICIT_ORACLE));
	cfg.dual_polish_iterations = 8;
	for (diagnostics = 0; diagnostics < 3; diagnostics++) {
		simplex_solve_context_init(
			&context, diagnostics ? &history : NULL, &p.execution);
		if (diagnostics == 2)
			for (i = 0; i < SIMPLEX_EXECUTION_CAPACITY; i++)
				(void)simplex_execution_begin(&p.execution,
					SIMPLEX_EXECUTION_SOLVE, 0, 0, 0);
		saved = p;
		for (phase = 0; phase < 2; phase++)
			for (field = 0; field < 19; field++) {
				p = saved;
				copy = p;
				simplex_pan_service_context_init(&c, &p);
				simplex_pan_service_bind(&a, &c);
				work = 2;
				if (phase)
					CHECK(simplex_pan_service_begin(&f, &a,
						SIMPLEX_PAN_SERVICE_DUAL_REFINE,
						work, 10, 3));
				switch (field) {
				case 0:
					p.rhs++;
					break;
				case 1:
					p.objective++;
					break;
				case 2:
					p.bounds++;
					break;
				case 3:
					p.row_type++;
					break;
				case 4:
					p.matrix.value++;
					break;
				case 5:
					p.matrix.column_start++;
					break;
				case 6:
					p.matrix.row_index++;
					break;
				case 7:
					p.matrix.row_start++;
					break;
				case 8:
					p.matrix.column_index++;
					break;
				case 9:
					p.matrix.row_value++;
					break;
				case 10:
					p.matrix.nonzeros++;
					break;
				case 11:
					p.matrix.rows++;
					break;
				case 12:
					p.matrix.columns++;
					break;
				case 13:
					p.rows++;
					break;
				case 14:
					p.columns++;
					break;
				case 15:
					simplex_problem_mark_changed(&p);
					break;
				case 16:
					p.execution.context = NULL;
					break;
				case 17:
					p.execution.model++;
					break;
				case 18:
					c.input = &copy;
					break;
				}
				CHECK(!simplex_pan_service_context_valid(&c));
				if (!phase)
					CHECK(!simplex_pan_service_begin(&f, &a,
						SIMPLEX_PAN_SERVICE_DUAL_REFINE,
						work, 10, 3));
				else
					work++;
				CHECK(!simplex_pan_service_finish(&f, &work, 1,
					lp_simplex_CondUnsatisfied,
					SIMPLEX_STOP_HANDOFF,
					SIMPLEX_PAN_OUTPUT_CHECKED_DUAL, &out));
				CHECK(c.failed && !out.valid && !out.output &&
					out.raw_output ==
						SIMPLEX_PAN_OUTPUT_CHECKED_DUAL);
				CHECK(work == (phase ? 3 : 2) &&
					out.grant == (phase ? 3 : 0));
			}
		p = saved;
		for (field = 0; field < 20; field++) {
			simplex_pan_service_context_init(&c, &p);
			simplex_pan_service_bind(&a, &c);
			original = a;
			work = 2;
			CHECK(simplex_pan_service_begin(&f, &a,
				SIMPLEX_PAN_SERVICE_DUAL_REFINE, work, 10, 3));
			switch (field) {
			case 0:
				a.value++;
				break;
			case 1:
				a.normalized_value++;
				break;
			case 2:
				a.column_norm++;
				break;
			case 3:
				a.rhs++;
				break;
			case 4:
				a.objective++;
				break;
			case 5:
				a.column_start++;
				break;
			case 6:
				a.row_index++;
				break;
			case 7:
				a.original_column++;
				break;
			case 8:
				a.original_scale++;
				break;
			case 9:
				a.original_shift++;
				break;
			case 10:
				a.original_columns++;
				break;
			case 11:
				a.objective_offset++;
				break;
			case 12:
				a.audit_objective = a.objective;
				break;
			case 13:
				a.audit_objective_offset++;
				break;
			case 14:
				a.nonzeros++;
				break;
			case 15:
				a.rows++;
				break;
			case 16:
				a.columns++;
				break;
			case 17:
				a.service_model++;
				break;
			case 18:
				a.service_revision++;
				break;
			case 19:
				a.services = NULL;
				break;
			}
			work++;
			CHECK(!simplex_pan_service_finish(&f, &work, 1,
				lp_simplex_CondUnsatisfied,
				SIMPLEX_STOP_HANDOFF,
				SIMPLEX_PAN_OUTPUT_CHECKED_DUAL, &out));
			CHECK(c.failed && work == 3 && !out.output);
			a = original;
		}
		for (field = 0; field < 3; field++) {
			simplex_pan_service_context_init(&c, &p);
			simplex_pan_service_bind(&b.a, &c);
			saved_b = b;
			copy = p;
			work = 2;
			CHECK(simplex_pan_service_begin(&f, &b.a,
				SIMPLEX_PAN_SERVICE_PDHG, work, 10, 3));
			/* The wrapper captures this same borrowed native envelope. */
			f.bounded = &b;
			f.bounded_input = b.problem;
			f.lower = b.lower;
			f.upper = b.upper;
			if (field == 0)
				b.lower++;
			else if (field == 1)
				b.upper++;
			else
				b.problem = &copy;
			work++;
			CHECK(!simplex_pan_service_finish(&f, &work, 1,
				lp_simplex_CondUnsatisfied,
				SIMPLEX_STOP_OPTIMAL,
				SIMPLEX_PAN_OUTPUT_NATIVE_KKT, &out));
			CHECK(work == 3 && !out.output && c.failed);
			b = saved_b;
		}
		/* Reinitializing a context cannot replace an in-flight input identity. */
		simplex_pan_service_context_init(&c, &p);
		simplex_pan_service_bind(&a, &c);
		copy = p;
		work = 2;
		CHECK(simplex_pan_service_begin(
			&f, &a, SIMPLEX_PAN_SERVICE_DUAL_REFINE, work, 10, 3));
		simplex_pan_service_context_init(&c, &copy);
		work++;
		CHECK(!simplex_pan_service_finish(&f, &work, 1,
			lp_simplex_CondUnsatisfied, SIMPLEX_STOP_HANDOFF,
			SIMPLEX_PAN_OUTPUT_CHECKED_DUAL, &out));
		/* Diagnostic redirection is not a numerical model replacement. */
		simplex_pan_service_context_init(&c, &p);
		simplex_pan_service_bind(&a, &c);
		work = 2;
		original = a;
		CHECK(simplex_pan_service_begin(
			&f, &a, SIMPLEX_PAN_SERVICE_DUAL_REFINE, work, 10, 3));
		p.execution.history = NULL;
		p.report = &other_report;
		a.report = &other_report;
		work++;
		CHECK(simplex_pan_service_finish(&f, &work, 1,
			lp_simplex_CondUnsatisfied, SIMPLEX_STOP_HANDOFF,
			SIMPLEX_PAN_OUTPUT_CHECKED_DUAL, &out));
		a = original;
		p = saved;
		/* Real dual refinement must not certify the old RHS after rebinding. */
		simplex_pan_service_context_init(&c, &p);
		simplex_pan_service_bind(&b.a, &c);
		p.rhs = &replacement;
		work = 0;
		multiplier[0] = -1.;
		CHECK(!simplex_pan_service_first_order(&b.a, point, multiplier,
			&o, &work, &violation, NULL, NULL, 1, &cfg, &out));
		CHECK(!out.valid && !out.output && !work &&
			multiplier[0] == -1. &&
			simplex_problem_primal_infeasibility(&p, point) == 1.);
		p = saved;
		copy = p;
		/* Native bounds belong to the same source object as the phase context. */
		simplex_pan_service_context_init(&c, &copy);
		simplex_pan_service_bind(&b.a, &c);
		work = 0;
		CHECK(!simplex_pan_service_pdhg(
			&b, point, multiplier, &o, &work, 8, &cfg, &out));
		CHECK(!out.valid && !out.output && !work && point[0] == 1. &&
			multiplier[0] == -1.);
		simplex_solve_context_destroy(&context);
	}
	p = saved;
	simplex_pan_bounded_destroy(&b);
	simplex_pan_standard_destroy(&a);
	simplex_problem_destroy(&p);
	lp_model_free(m);
	return 0;
}

static int event_lifetime(void)
{
	struct lp_Model *model = lp_model_create(1, 2);
	struct simplex_Problem p;
	struct simplex_PanStandard standard;
	struct simplex_PanServiceContext services;
	struct simplex_PanServiceFrame frame;
	struct simplex_PanServiceOutcome outcome;
	struct simplex_SolveContext context;
	struct simplex_Execution history, replacement;
	int role, variant, work, id, count, j;
	CHECK(model);
	model->constraints[0].type = optm_CONS_T_EQ;
	model->constraints[0].rhs = 1.;
	model->constraints[0].coef[0] = model->constraints[0].coef[1] = 1.;
	CHECK(!simplex_problem_from_model(&p, model));
	CHECK(!simplex_pan_standard_create(&standard, &p));
	for (role = 0; role < SIMPLEX_PAN_SERVICE_COUNT; role++)
		for (variant = 0; variant < 2; variant++) {
			simplex_solve_context_init(
				&context, &history, &p.execution);
			simplex_pan_service_context_init(&services, &p);
			simplex_pan_service_bind(&standard, &services);
			CHECK(simplex_pan_service_begin(
				&frame, &standard, role, 0, 3, 3));
			id = frame.outcome.event;
			CHECK(id > 0);
			if (variant == 0) {
				count = history.count;
				memset(&history, 0, sizeof(history));
				for (j = 0; j < count; j++)
					(void)simplex_execution_begin(
						&p.execution,
						SIMPLEX_EXECUTION_PAN_SERVICE,
						role, 0, 0);
			} else {
				replacement = history;
				services.scope.history = &replacement;
			}
			work = 1;
			CHECK(simplex_pan_service_finish(&frame, &work, 0,
				lp_simplex_CondUnsatisfied,
				SIMPLEX_STOP_HANDOFF, SIMPLEX_PAN_OUTPUT_NONE,
				&outcome));
			CHECK(work == 1 && outcome.valid && !outcome.event &&
				!services.failed);
			CHECK(!services.scope.history->events[id - 1].complete);
			simplex_solve_context_destroy(&context);
		}
	simplex_pan_standard_destroy(&standard);
	simplex_problem_destroy(&p);
	lp_model_free(model);
	return 0;
}

static int repair_trigger(void)
{
	struct lp_Model *m = lp_model_create(1, 2);
	struct simplex_Problem p;
	struct simplex_PanBounded b;
	struct simplex_PanServiceContext c;
	struct simplex_PanServiceFrame outer;
	struct simplex_PanServiceOutcome out;
	struct simplex_PanServiceTrigger trigger;
	struct simplex_SolveContext context;
	struct simplex_Execution h;
	struct simplex_StageReport report;
	struct simplex_PanConfig cfg;
	struct lp_simplex_Options o;
	double point[3], saved[3], residual[1], scratch[1];
	int mode, work, saved_work = 0, id, cause;

	CHECK(m);
	m->constraints[0].type = optm_CONS_T_GE;
	m->constraints[0].rhs = 1.;
	m->constraints[0].coef[0] = m->constraints[0].coef[1] = 1.;
	m->objective[1] = 1.;
	CHECK(!simplex_problem_from_model(&p, m));
	p.report = &report;
	CHECK(!simplex_pan_bounded_create(&b, &p));
	CHECK(simplex_pan_config_default(&cfg, SIMPLEX_PAN_NATIVE_BOUNDS,
		SIMPLEX_PAN_CONE_DIRECTION, SIMPLEX_PAN_IMPLICIT_ORACLE));
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_PAN_BDA);
	o.iteration_limit = 128;
	for (mode = 0; mode < 6; mode++) {
		simplex_stage_reset(&report);
		simplex_solve_context_init(&context, &h, &p.execution);
		if (mode == 1)
			p.execution.history = NULL;
		simplex_pan_service_context_init(&c, &p);
		simplex_pan_service_bind(&b.a, &c);
		work = 17;
		CHECK(simplex_pan_service_begin(&outer, &b.a,
			SIMPLEX_PAN_SERVICE_PDHG, work, o.iteration_limit,
			100));
		if (mode == 2)
			h.count = SIMPLEX_EXECUTION_CAPACITY;
		if (mode == 3)
			h.count = SIMPLEX_EXECUTION_CAPACITY - 1;
		if (mode == 4)
			h.count = SIMPLEX_EXECUTION_CAPACITY - 2;
		if (mode == 5)
			context.execution_serial = ULONG_MAX;
		point[0] = .6;
		point[1] = .3;
		point[2] = 0.;
		trigger.kind = SIMPLEX_PAN_TRIGGER_PDHG_STAGNATION;
		trigger.work = work;
		trigger.epoch = 7;
		trigger.primal = simplex_pan_bounded_primal(
			&b, point, residual, scratch);
		trigger.previous_primal = .4;
		CHECK(simplex_pan_service_restore(&b, point, &o, &work, 2.,
			      &cfg, &trigger, &out) == 1);
		CHECK(out.valid && !c.failed &&
			out.trigger.kind == trigger.kind);
		CHECK(out.trigger.primal == trigger.primal &&
			out.trigger.previous_primal == .4 &&
			out.trigger.epoch == 7);
		CHECK(report.services.calls[SIMPLEX_PAN_SERVICE_RESTORE] == 1);
		if (mode == 0) {
			saved_work = work;
			memcpy(saved, point, sizeof(saved));
		} else
			CHECK(work == saved_work &&
				!memcmp(point, saved, sizeof(saved)));
		id = out.event;
		if (mode == 0 || mode == 4) {
			CHECK(id > 0 && h.events[id - 1].reason ==
						SIMPLEX_REASON_PDHG_STAGNATION);
			cause = h.events[id - 1].cause;
			CHECK(cause > 0 && cause != id);
			CHECK(h.events[cause - 1].kind ==
				SIMPLEX_EXECUTION_PAN_TRIGGER);
			CHECK(h.events[cause - 1].complete &&
				!h.events[cause - 1].result_available);
			CHECK(!h.events[cause - 1].work &&
				h.events[cause - 1].observation_work == 17);
			CHECK(h.events[cause - 1].primal.value ==
				trigger.primal);
			CHECK(h.events[cause - 1].primal.epoch == 7);
			CHECK(h.events[cause - 1].primal.space ==
				SIMPLEX_SPACE_WORKING);
			CHECK(h.events[id - 1].parent == outer.outcome.event &&
				h.events[cause - 1].parent ==
					outer.outcome.event);
		} else
			CHECK(!id);
		CHECK(simplex_pan_service_finish(&outer, &work, 0,
			lp_simplex_CondUnsatisfied, SIMPLEX_STOP_HANDOFF,
			SIMPLEX_PAN_OUTPUT_NONE, NULL));
		simplex_solve_context_destroy(&context);
	}
	simplex_pan_bounded_destroy(&b);
	simplex_problem_destroy(&p);
	lp_model_free(m);
	return 0;
}

static int stale_trigger(void)
{
	struct lp_Model *m = lp_model_create(1, 1);
	struct simplex_Problem p;
	struct simplex_PanStandard a;
	struct simplex_PanServiceContext c;
	struct simplex_PanServiceFrame frame;
	struct simplex_PanServiceOutcome out;
	struct simplex_SolveContext context;
	struct simplex_Execution h;
	struct simplex_ExecutionScope scope;
	int mutation, cause, work;

	CHECK(m);
	CHECK(!simplex_problem_from_model(&p, m));
	memset(&a, 0, sizeof(a));
	for (mutation = 0; mutation < 2; mutation++) {
		simplex_solve_context_init(&context, &h, &p.execution);
		simplex_pan_service_context_init(&c, &p);
		simplex_pan_service_bind(&a, &c);
		scope = c.scope;
		scope.model = a.service_model;
		cause = simplex_execution_begin(&scope,
			SIMPLEX_EXECUTION_PAN_TRIGGER,
			SIMPLEX_PAN_TRIGGER_PDHG_STAGNATION,
			SIMPLEX_REASON_PDHG_STAGNATION, 0);
		simplex_execution_end(&h, cause, NULL, NULL);
		CHECK(simplex_pan_service_begin(
			&frame, &a, SIMPLEX_PAN_SERVICE_RESTORE, 0, 1, 1));
		frame.trigger_actor = simplex_execution_ref(&scope, cause);
		frame.outcome.trigger.kind =
			SIMPLEX_PAN_TRIGGER_PDHG_STAGNATION;
		if (mutation)
			h.events[cause - 1].serial++;
		work = 0;
		CHECK(simplex_pan_service_finish(&frame, &work, 0,
			lp_simplex_CondUnsatisfied, SIMPLEX_STOP_HANDOFF,
			SIMPLEX_PAN_OUTPUT_NONE, &out));
		CHECK(!work && !c.failed && out.valid);
		CHECK(h.events[out.event - 1].cause == (mutation ? 0 : cause));
		CHECK(h.events[out.event - 1].reason ==
			SIMPLEX_REASON_PDHG_STAGNATION);
		simplex_solve_context_destroy(&context);
	}
	simplex_problem_destroy(&p);
	lp_model_free(m);
	return 0;
}

int main(void)
{
	if (leases() || roles() || repair_trigger() || stale_trigger() ||
		identities() || event_lifetime())
		return 1;
	puts("BDA phase leases, input identities, roles, nesting and journal independence passed");
	return 0;
}
