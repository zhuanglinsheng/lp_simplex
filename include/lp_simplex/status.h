/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file status.h
 * @brief Detailed solver termination statuses and API return codes.
 */

#ifndef LP_SIMPLEX_STATUS_H
#define LP_SIMPLEX_STATUS_H

/**
 * Detailed termination reason, reported in lp_simplex_Result.status.
 * These values are distinct from the function return codes below. Only
 * lp_simplex_Success represents an optimal solution; budget or numerical
 * failures do not establish infeasibility or unboundedness.
 */
enum lp_Status {
	/** Optimal solve accepted by the final original-model feasibility check. */
	lp_simplex_Success = 0,

	/** Required allocation failed. */
	lp_simplex_MemoryAllocError = 1,

	/** Invalid input/options or unsupported solver preconditions. */
	lp_simplex_CondUnsatisfied = 2,

	/** Cumulative work budget exhausted; a candidate point may be available. */
	lp_simplex_ExceedIterLimit = 3,

	/** Singular basis or numerical factorisation failure. */
	lp_simplex_Singularity = 4,

	/** Backend could not handle the constraint system; not an infeasibility
	 * proof. */
	lp_simplex_OverDetermination = 5,

	/** Solver determined that the objective is unbounded. */
	lp_simplex_Unboundedness = 6,

	/** Solver determined that no feasible point exists. */
	lp_simplex_Infeasibility = 7,

	/** Degeneracy-related termination without an accepted optimum. */
	lp_simplex_Degeneracy = 8,

	/** Numerical precision or final result validation failed. */
	lp_simplex_PrecisionError = 9
};

/** API failure return; consult result.status for the detailed reason. */
#define lp_simplex_EXIT_FAILURE (-1)

/** API success return; result.status is lp_simplex_Success. */
#define lp_simplex_EXIT_SUCCESS 0

#endif
