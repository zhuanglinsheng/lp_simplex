/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file version.c
 * @brief Report the loaded library version.
 */

#include <lp_simplex/version.h>

const char *lp_simplex_version(void)
{
	return LP_SIMPLEX_VERSION;
}
