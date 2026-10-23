/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_pan_config.c
 * @brief Resolve Pan numerical controls at the solver entry boundary.
 */

#include "simplex_pan_config.h"
#include "simplex_pan_basis.h"
#include "utils.h"

#include <stdlib.h>
#include <string.h>


static int setting_integer(const char *name)
{
	const char *value = getenv(name);

	return value ? atoi(value) : 0;
}


void simplex_pan_config_init(
	struct simplex_PanConfig *c, const struct lp_simplex_Options *o)
{
	const char *value;

	lp_simplex_memset(c, 0, sizeof(*c));
	c->nnls_block_limit = 1;
	c->singleton_crash = 1;
	if (o && (o->algorithm == lp_simplex_ALGORITHM_AUTO ||
			 o->algorithm == SIMPLEX_PAN_AUTO_EXPLICIT ||
			 o->algorithm == SIMPLEX_PAN_AUTO_STANDARD_IMPLICIT ||
			 o->algorithm == SIMPLEX_PAN_AUTO_STANDARD_EXPLICIT)) {
		c->implicit_cone = o->algorithm == lp_simplex_ALGORITHM_AUTO;
		c->implicit_updates = 1;
		c->pdhg_iterations =
			o->algorithm == lp_simplex_ALGORITHM_AUTO ? 500000 : 0;
		c->dual_polish_iterations = 8192;
		c->restore = 1;
		c->adaptive = 1;
		return;
	}
	value = getenv("LP_SIMPLEX_PAN_CONE");
	c->implicit_cone = value && strcmp(value, "implicit") == 0;
	c->implicit_updates = getenv("LP_SIMPLEX_PAN_IMPLICIT_UPDATES") != NULL;
	c->implicit_bulk = getenv("LP_SIMPLEX_PAN_IMPLICIT_BULK") != NULL;
	c->pdhg_iterations = setting_integer("LP_SIMPLEX_PAN_NATIVE_PDHG");
	c->restore = getenv("LP_SIMPLEX_PAN_PDHG_RESTORE") != NULL;
	c->adaptive =
		c->restore || getenv("LP_SIMPLEX_PAN_NATIVE_ADAPTIVE") != NULL;
	c->predict = getenv("LP_SIMPLEX_PAN_NATIVE_PREDICT") != NULL;
	c->dual_start = getenv("LP_SIMPLEX_PAN_DUAL_START") != NULL;
	c->dual_polish_iterations =
		setting_integer("LP_SIMPLEX_PAN_DUAL_POLISH");
	c->first_order_iterations =
		setting_integer("LP_SIMPLEX_PAN_FIRST_ORDER");
	c->newton_qr = getenv("LP_SIMPLEX_PAN_NEWTON_QR") != NULL;
	c->dual_center = getenv("LP_SIMPLEX_PAN_DUAL_CENTER") != NULL;
	value = getenv("LP_SIMPLEX_PAN_NNLS_BLOCK");
	if (value) {
		c->nnls_block_limit = atoi(value);
		if (c->nnls_block_limit < 1)
			c->nnls_block_limit = 1;
		if (c->nnls_block_limit > 32)
			c->nnls_block_limit = 32;
	}
	c->block_crash = getenv("LP_SIMPLEX_PAN_BLOCK_CRASH") != NULL;
	c->singleton_crash =
		getenv("LP_SIMPLEX_PAN_DISABLE_SINGLETON_CRASH") == NULL;
	c->full_refinement = getenv("LP_SIMPLEX_PAN_FULL_REFINEMENT") != NULL;
	value = getenv("LP_SIMPLEX_PAN_PDHG_WEIGHT");
	c->pdhg_weight_set = value != NULL;
	c->pdhg_weight = value ? atof(value) : 0.;
}


void simplex_pan_basis_configure(
	struct simplex_PanBasis *b, const struct simplex_PanConfig *c)
{
	b->prefer_updates = c->implicit_updates;

	b->implicit_bulk = c->implicit_bulk;
	b->nnls_block_limit = c->nnls_block_limit;
	b->block_crash = c->block_crash;
	b->singleton_crash = c->singleton_crash;
	b->full_refinement = c->full_refinement;
}


/* Resolve legacy standard-form dispatch at its entry, without changing the
 * process environment. AUTO routes ignore unrelated manual experiments. */
const char *simplex_pan_setting(
	const struct lp_simplex_Options *o, const char *name)
{
	if (o->algorithm == SIMPLEX_PAN_AUTO_STANDARD_IMPLICIT ||
		o->algorithm == SIMPLEX_PAN_AUTO_STANDARD_EXPLICIT) {
		if (strcmp(name, "LP_SIMPLEX_PAN_CONE") == 0)
			return o->algorithm == SIMPLEX_PAN_AUTO_STANDARD_IMPLICIT
				       ? "implicit"
				       : "1";
		return NULL;
	}
	return o->algorithm == lp_simplex_ALGORITHM_AUTO ? NULL : getenv(name);
}
