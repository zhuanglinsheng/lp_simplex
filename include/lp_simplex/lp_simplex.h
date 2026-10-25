/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file lp_simplex.h
 * @brief Public entry point for the LP model, MPS reader and simplex solver
 *        APIs.
 *
 * Include <lp_simplex/lp_simplex.h> and link against lp_simplex. Declarations
 * support C and C++ callers. Only continuous linear programs are supported.
 *
 * Typical lifecycle:
 * 1. Create a model with lp_model_create() / lp_model_create_sparse(), or read
 *    one with lp_read_mps(). Check for NULL before accessing it.
 * 2. Fill the model. For dense manual input, call lp_model_build_sparse()
 *    after editing matrix entries and check for a zero return value. MPS
 *    input already has both sparse views populated.
 * 3. Call lp_simplex_default_options(), then override selected options.
 * 4. Allocate model->n doubles and call lp_simplex_solve(). Check both its
 *    return value and result.status before using the point as an optimum.
 * 5. Release the solution buffer using its matching allocator, and release
 *    the model with lp_model_free(). Options/results need no library cleanup.
 *
 * See model.h for storage and objective-sense conventions, solve.h for option
 * defaults and diagnostic semantics, and status.h for termination reasons.
 */

#ifndef LP_SIMPLEX_H
#define LP_SIMPLEX_H

#include "model.h"
#include "mps.h"
#include "solve.h"
#include "status.h"
#include <lp_simplex/version.h>

#endif
