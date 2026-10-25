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
 * @brief Read an MPS file into a caller-owned sparse model ready for solving.
 *
 * The file must describe at least one constraint and one variable. Names use
 * the fixed-column, eight-character MPS layout; numeric tokens retain long
 * decimal values rather than being truncated to twelve characters. Returned
 * names may retain fixed-column padding. Default objective sense is minimise.
 * Required numeric fields must contain complete finite tokens. Overflow,
 * underflow to zero, missing numbers and nonfinite ranged-row endpoints are
 * rejected; representable subnormal numbers and long decimals are retained.
 *
 * Supported sections are `ROWS`, `COLUMNS`, `RHS`, `RANGES`, `BOUNDS`,
 * `OBJSENSE` and `ENDATA`.  `OBJSENSE MAXIMIZE` is honoured: the stored
 * objective is negated and `objective_sense` is set to -1, so the reported
 * objective is the maximum.  The sense may appear on the OBJSENSE line or
 * the next non-comment, non-blank line; invalid or missing senses are rejected.
 * Comment lines (`*` or `$` in the first non-blank position) and blank lines
 * are ignored anywhere in the file. Ranged rows are expanded into lower/upper
 * inequality pairs, so model->m can exceed the number of input constraint rows.
 *
 * Integer declarations (`INTORG`/`INTEND`, `BV`, `LI`, `UI`) are read and
 * recorded in `bounds[].v_type`, but this is an LP solver: lp_simplex_solve()
 * rejects such a model instead of solving its relaxation.  Free-format MPS,
 * `SOS`, `QUADOBJ` and multiple RHS/RANGES/BOUNDS vectors are not supported.
 *
 * Both CSC and CSR are populated; coefficients and constraints[i].coef are
 * NULL. Do not call lp_model_build_sparse() on this sparse-only model. To edit
 * its matrix, maintain both sparse representations consistently. Sparse input
 * does not imply that every algorithm avoids dense allocations (see TABLEAU).
 *
 * @param path Non-NULL, NUL-terminated filename of a readable MPS text file.
 * @return Caller-owned model, released with lp_model_free(), or NULL on I/O,
 * allocation, malformed-input or unsupported-input failure. Parsing errors
 * may print diagnostics to stderr; no separate parser error object is returned.
 * @see lp_simplex_solve(), lp_model_free()
 */
struct lp_Model *lp_read_mps(const char *path);

#ifdef __cplusplus
}
#endif

#endif
