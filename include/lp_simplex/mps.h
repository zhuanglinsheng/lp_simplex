/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file mps.h
 * @brief Fixed-column MPS parsing and model construction.
 */

#ifndef LP_SIMPLEX_MPS_H
#define LP_SIMPLEX_MPS_H

#include "model.h"


#ifdef __cplusplus
extern "C" {
#endif

/**
 * Read a fixed-column MPS file.
 *
 * Supported sections are `ROWS`, `COLUMNS`, `RHS`, `RANGES`, `BOUNDS`,
 * `OBJSENSE` and `ENDATA`.  `OBJSENSE MAXIMIZE` is honoured: the stored
 * objective is negated and `objective_sense` is set to -1, so the reported
 * objective is the maximum.  The sense may appear on the OBJSENSE line or
 * the next non-comment, non-blank line; invalid or missing senses are rejected.
 * Comment lines (`*` or `$` in the first
 * non-blank position) and blank lines are ignored anywhere in the file.
 *
 * Integer declarations (`INTORG`/`INTEND`, `BV`, `LI`, `UI`) are read and
 * recorded in `bounds[].v_type`, but this is an LP solver: lp_simplex_solve()
 * rejects such a model instead of solving its relaxation.  Free-format MPS,
 * `SOS`, `QUADOBJ` and multiple RHS/RANGES/BOUNDS vectors are not supported.
 *
 * The returned model is owned by the caller and must be released with
 * lp_model_free().  NULL indicates an I/O error, a malformed file (a message
 * is printed to stderr), or unsupported input.
 */
struct lp_Model *lp_read_mps(const char *path);

#ifdef __cplusplus
}
#endif

#endif
