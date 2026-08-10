/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#ifndef LP_SIMPLEX_DUAL_FEASIBILITY_INTERNAL_H
#define LP_SIMPLEX_DUAL_FEASIBILITY_INTERNAL_H


struct simplex_DualState;
struct simplex_SparseVector;

/* Indexed primal-infeasibility queues.  The second queue lets Pan's
 * structure-preserving mode select a structural basic row without scanning
 * the basis. */
struct simplex_DualFeasibility {
	int rows;
	int *heap;
	int *slot;
	int count;
	int *structural_heap;
	int *structural_slot;
	int structural_count;
	double *score;
	double *merit;
	double total_merit;
	long rebuilds;
	long incremental_updates;
};


int simplex_dual_feasibility_create(
		struct simplex_DualFeasibility *feasibility, int rows);

void simplex_dual_feasibility_destroy(
		struct simplex_DualFeasibility *feasibility);

void simplex_dual_feasibility_rebuild(struct simplex_DualState *state);

void simplex_dual_feasibility_update(struct simplex_DualState *state, int row);

void simplex_dual_feasibility_update_packed(
		struct simplex_DualState *state,
		const struct simplex_SparseVector *changed_rows);

int simplex_dual_feasibility_choose(
		struct simplex_DualState *state, double *target, int *kappa,
		double *maximum, int prefer_structural, int lexicographic,
		int deferred_row);

double simplex_dual_feasibility_merit(
		const struct simplex_DualState *state);

#endif
