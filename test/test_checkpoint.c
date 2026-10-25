/**
 * @file test_checkpoint.c
 * @brief Checkpoint ownership and causal basis recovery, without instance
 * names.
 */

#include "simplex_stage.h"
#include "simplex_solve_context.h"
#include "simplex_execution.h"
#include "simplex_dual.h"
#include "simplex_scaling.h"
#include <lp_simplex/lp_simplex.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#define CHECK(c) \
	do { \
		if (!(c)) { \
			fprintf(stderr, "checkpoint line %d: %s\n", __LINE__, \
				#c); \
			return 1; \
		} \
	} while (0)

static int retry_contract(void)
{
	struct simplex_DualBasisSeed seed;
	memset(&seed, 0, sizeof(seed));
	CHECK(!simplex_dual_seed_retry_allowed(NULL));
	CHECK(!simplex_dual_seed_retry_allowed(&seed));
	seed.import_attempted = 1;
	CHECK(!simplex_dual_seed_retry_allowed(
		&seed)); /* Unknown is not permission. */
	seed.rejection = SIMPLEX_BASIS_LAYOUT_REJECTED;
	CHECK(!simplex_dual_seed_retry_allowed(&seed));
	seed.rejection = SIMPLEX_BASIS_FEASIBILITY_REJECTED;
	CHECK(!simplex_dual_seed_retry_allowed(&seed));
	seed.rejection = SIMPLEX_BASIS_MODEL_REJECTED;
	CHECK(!simplex_dual_seed_retry_allowed(&seed));
	seed.rejection = SIMPLEX_BASIS_FACTORIZATION_REJECTED;
	CHECK(simplex_dual_seed_retry_allowed(&seed));
	seed.rejection = SIMPLEX_BASIS_RESIDUAL_REJECTED;
	CHECK(simplex_dual_seed_retry_allowed(&seed));
	seed.rejection = SIMPLEX_BASIS_NOT_REJECTED;
	seed.accepted = 1;
	CHECK(simplex_dual_seed_retry_allowed(&seed));
	return 0;
}

static int binding(int scaled)
{
	struct lp_Model *m = lp_model_create(1, 2);
	struct simplex_Problem p, other;
	struct simplex_DualCheckpoint checkpoint;
	struct simplex_StageReport report;
	struct simplex_SolveContext context;
	struct simplex_Execution history;
	struct simplex_Scaling scaling;
	struct simplex_Problem *transformed;
	struct lp_simplex_Options o;
	struct lp_simplex_Result r;
	int index[1], saved_index[1];
	unsigned char status[3], saved_status[3];
	double x[2], scale = scaled ? 1048576. : 1.;
	CHECK(m);
	m->constraints[0].type = optm_CONS_T_EQ;
	m->constraints[0].rhs = scale;
	m->constraints[0].coef[0] = m->constraints[0].coef[1] = scale;
	m->objective[1] = 1.;
	CHECK(!simplex_problem_from_model(&p, m));
	CHECK(!simplex_problem_from_model(&other, m));
	simplex_solve_context_init(&context, &history, &p.execution);
	p.report = &report;
	other.report = &report;
	simplex_dual_checkpoint_init(&checkpoint, &p, index, status);
	CHECK(simplex_dual_checkpoint_matches(&checkpoint, &p));
	CHECK(!simplex_dual_checkpoint_matches(&checkpoint, &other));
	transformed = simplex_scaling_create_problem(&p, &scaling);
	CHECK(transformed);
	CHECK(!simplex_dual_checkpoint_matches(&checkpoint, transformed));
	simplex_problem_free(transformed);
	simplex_scaling_destroy(&scaling);
	lp_simplex_default_options(&o, lp_simplex_ALGORITHM_AUTO);
	o.iteration_limit = 0;
	CHECK(simplex_dual_solve_problem_checkpoint(
		&p, &o, x, NULL, &r, 0, &checkpoint));
	CHECK(checkpoint.basis.valid && r.status == lp_simplex_ExceedIterLimit);
	/* A different route parent does not change the immutable model. Internal
	 * scaling also gets a new identity on each invocation, outside this match. */
	simplex_execution_set_parent(&p.execution, history.count);
	CHECK(p.execution.parent > 0 &&
		simplex_dual_checkpoint_matches(&checkpoint, &p));
	o.iteration_limit = 100;
	CHECK(!simplex_dual_solve_problem_checkpoint(
		&p, &o, x, NULL, &r, 0, &checkpoint));
	CHECK(checkpoint.basis.accepted && report.basis_import_attempted);
	CHECK(!report.basis_rejection && !report.factor_available);
	/* Export provenance does not depend on optional stage reports. */
	p.report = NULL;
	CHECK(!simplex_dual_solve_problem_checkpoint(
		&p, &o, x, NULL, &r, 0, &checkpoint));
	CHECK(checkpoint.basis.exported &&
		simplex_dual_seed_source_event(
			&checkpoint.basis, &p.execution) > 0);
	p.report = &report;
	/* Diagnostic sinks are not immutable model identity. */
	p.execution.history = NULL;
	CHECK(simplex_dual_checkpoint_matches(&checkpoint, &p));
	CHECK(!simplex_dual_solve_problem_checkpoint(
		&p, &o, x, NULL, &r, 0, &checkpoint));
	CHECK(checkpoint.basis.accepted && report.basis_import_attempted);
	p.execution.history = &history;
	p.execution.context = NULL;
	CHECK(!simplex_dual_checkpoint_matches(&checkpoint, &p));
	p.execution.context = &context;
	checkpoint.basis.columns--;
	CHECK(simplex_dual_checkpoint_rejection(&checkpoint, &p) ==
		SIMPLEX_BASIS_LAYOUT_REJECTED);
	checkpoint.basis.columns++;
	p.execution.model = 0;
	CHECK(simplex_dual_checkpoint_rejection(&checkpoint, &p) ==
		SIMPLEX_BASIS_MODEL_REJECTED);
	p.execution.model = checkpoint.basis.identity.model;
	memcpy(saved_index, index, sizeof(index));
	memcpy(saved_status, status, sizeof(status));
	/* Matching dimensions do not permit another model to overwrite storage. */
	CHECK(!simplex_dual_solve_problem_checkpoint(
		&other, &o, x, NULL, &r, 0, &checkpoint));
	CHECK(r.status == lp_simplex_Success && !checkpoint.basis.valid);
	CHECK(!checkpoint.basis.import_attempted &&
		checkpoint.basis.rejection == SIMPLEX_BASIS_MODEL_REJECTED);
	CHECK(report.basis_rejection == SIMPLEX_BASIS_MODEL_REJECTED);
	CHECK(!memcmp(saved_index, index, sizeof(index)));
	CHECK(!memcmp(saved_status, status, sizeof(status)));
	/* Rebinding is explicit and discards old feasibility claims. */
	simplex_dual_checkpoint_init(&checkpoint, &other, index, status);
	CHECK(!checkpoint.basis.valid &&
		simplex_dual_checkpoint_matches(&checkpoint, &other));
	CHECK(!simplex_dual_solve_problem_checkpoint(
		&other, &o, x, NULL, &r, 0, &checkpoint));
	simplex_problem_mark_changed(&other);
	other.rhs[0] = 2. * scale;
	CHECK(!simplex_dual_checkpoint_matches(&checkpoint, &other));
	CHECK(!simplex_dual_solve_problem_checkpoint(
		&other, &o, x, NULL, &r, 0, &checkpoint));
	CHECK(!checkpoint.basis.accepted && !checkpoint.basis.valid);
	CHECK(fabs(x[0] + x[1] - 2.) < 1e-7);
	other.revision = ULONG_MAX - 1;
	simplex_problem_mark_changed(&other);
	simplex_dual_checkpoint_init(&checkpoint, &other, index, status);
	CHECK(!simplex_dual_checkpoint_matches(&checkpoint, &other));
	simplex_problem_mark_changed(&other);
	CHECK(other.revision == ULONG_MAX);
	simplex_problem_destroy(&p);
	simplex_problem_destroy(&other);
	lp_model_free(m);
	simplex_solve_context_destroy(&context);
	return 0;
}

int main(void)
{
	return retry_contract() || binding(0) || binding(1);
}
