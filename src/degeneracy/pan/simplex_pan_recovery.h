/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_pan_recovery.h
 * @brief Pan support recovery and bounded original-model corrections.
 */

#ifndef LP_SIMPLEX_PAN_RECOVERY_INTERNAL_H
#define LP_SIMPLEX_PAN_RECOVERY_INTERNAL_H
#include "simplex_pan_standard.h"

#include <lp_simplex/solve.h>


/* One instance per Phase II, shared by both backends and all center epochs. */
struct simplex_PanRecovery {
	int low_progress;
	int have_compressed;
	unsigned long compressed_point[2];
	long directions, accepted, rejected, compressions, compression_skips;
	long restorations, restored;
	double direction_seconds, recovery_seconds;
};

/* Suppress another support recovery at the identical physical point. */
int simplex_pan_recovery_compress(
	struct simplex_PanRecovery *state, const double *point, int columns);

/* Bounded correction on independent positive columns of the original A.
 * Returns 1 on audited improvement, 0 when no correction is needed, -1 when
 * unsuccessful. Failure restores the complete input. Uses shared budget.
 * ceiling is the ORIGINAL cost without its constant offset. */
int simplex_pan_restore(const struct simplex_PanStandard *a, double *point,
	const struct lp_simplex_Options *options, long double ceiling,
	int *iterations);
#endif
