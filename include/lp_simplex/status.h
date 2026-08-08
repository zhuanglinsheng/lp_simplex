/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#ifndef LP_SIMPLEX_STATUS_H
#define LP_SIMPLEX_STATUS_H


enum lp_Status {
	lp_simplex_Success           = 0,
	lp_simplex_MemoryAllocError  = 1,
	lp_simplex_CondUnsatisfied   = 2,
	lp_simplex_ExceedIterLimit   = 3,
	lp_simplex_Singularity       = 4,
	lp_simplex_OverDetermination = 5,
	lp_simplex_Unboundedness     = 6,
	lp_simplex_Infeasibility     = 7,
	lp_simplex_Degeneracy        = 8,
	lp_simplex_PrecisionError    = 9
};

#define lp_simplex_EXIT_FAILURE (-1)
#define lp_simplex_EXIT_SUCCESS 0

#endif
