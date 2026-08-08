/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
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
 * The returned model is owned by the caller and must be released with
 * lp_model_free().  NULL indicates an I/O error or unsupported input.
 */
struct lp_Model *lp_read_mps(const char *path);

#ifdef __cplusplus
}
#endif

#endif
