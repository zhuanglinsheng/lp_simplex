/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_pan_pdhg.h
 * @brief Internal short-column plan and compensated PDHG residual kernel.
 */

#ifndef LP_SIMPLEX_PAN_PDHG_INTERNAL_H
#define LP_SIMPLEX_PAN_PDHG_INTERNAL_H

#include "simplex_pan_standard.h"

/* Per-call progress history. A metric change must reset the history rather
 * than compare residuals measured with different PDHG weights. */
struct simplex_PdhgProgress {
	int initialized, start, progress_at;
	double best;
};

void simplex_pan_pdhg_progress_reset(struct simplex_PdhgProgress *history);
/* Necessary residual gate for a new repair attempt within one PDHG call.
 * Stagnation, objective ceiling and remaining work are checked by the caller. */
int simplex_pan_pdhg_repair_progress(
	double primal, double previous_primal, double primal_tolerance);
/* Observe the better current/average merit. Stagnation permits intervention,
 * but never establishes termination, feasibility or optimality. */
int simplex_pan_pdhg_progress_stalled(
	struct simplex_PdhgProgress *history, int work, double merit);
struct simplex_PanBounded;
double simplex_pan_pdhg_handoff_gap(const struct simplex_PanBounded *b,
	const double *x, const double *y, double scale);

enum simplex_PanPdhgColumnKind {
	PDHG_COLUMN_SINGLE = 1,
	PDHG_COLUMN_PAIR = 2,
	PDHG_COLUMN_NEGATIVE_FIRST = 4,
	PDHG_COLUMN_NEGATIVE_SECOND = 8
};

/* Owned metadata for immutable signed-unit columns. A zero kind selects
 * the generic kernel; allocation failure leaves the entire plan empty. */
struct simplex_PanPdhgColumns {
	int *indices, *first, *second;
	unsigned char *kind;
	int compiled;
};

void simplex_pan_pdhg_columns_create(const struct simplex_PanStandard *a,
	struct simplex_PanPdhgColumns *plan);

void simplex_pan_pdhg_columns_destroy(struct simplex_PanPdhgColumns *plan);

/* Form b - A*(2*next - previous), preserving the generic kernel's entry
 * order and compensated summation. Nonfinite rows are left for the caller
 * to reject; this iteration kernel does not certify a solution. */
void simplex_pan_pdhg_residual(const struct simplex_PanStandard *a,
	const struct simplex_PanPdhgColumns *plan, const double *next,
	const double *previous, double *residual, double *correction);

#endif
