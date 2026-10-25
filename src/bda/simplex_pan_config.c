/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_pan_config.c
 * @brief Resolve Pan numerical controls at the solver entry boundary.
 */

#include "simplex_solve_context.h"
#include "policy/pan.h"
#include "policy/controls.h"
#include "simplex_pan_config.h"
#include "simplex_pan_basis.h"
#include "utils.h"

#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <math.h>
#include "simplex_problem.h"

void simplex_pan_config_clean(struct simplex_PanConfig *c)
{
	lp_simplex_memset(c, 0, sizeof(*c));
	c->nnls_block_limit = c->singleton_crash = 1;
	c->allow_oracle_recovery = 1;
	c->model_form = SIMPLEX_PAN_STANDARD;
	c->direction_method = SIMPLEX_PAN_BASIS_DIRECTION;
	c->oracle = SIMPLEX_PAN_NO_ORACLE;
}

void simplex_pan_config_resolve(struct simplex_PanConfig *c,
	const struct simplex_Problem *p, const struct lp_simplex_Options *o)
{
	if (p && p->execution.context && p->execution.context->configured &&
		o && o->algorithm != lp_simplex_ALGORITHM_AUTO)
		*c = p->execution.context->controls.pan;
	else
		simplex_pan_config_init(c, o);
}

int simplex_pan_integer_setting_value(const char *value)
{
	char *end;
	long parsed;

	if (!value)
		return 0;
	errno = 0;
	parsed = strtol(value, &end, 10);
	if (end == value)
		return 0;
	while (isspace((unsigned char)*end))
		end++;
	if (*end || parsed < 0)
		return 0;
	return errno == ERANGE || parsed > INT_MAX ? INT_MAX : (int)parsed;
}

static int setting_integer(const char *name)
{
	return simplex_pan_integer_setting_value(getenv(name));
}

int simplex_pan_config_valid(const struct simplex_PanConfig *c)
{
	if (!c || (c->model_form != SIMPLEX_PAN_STANDARD &&
			  c->model_form != SIMPLEX_PAN_NATIVE_BOUNDS))
		return 0;
	if (c->direction_method == SIMPLEX_PAN_BASIS_DIRECTION) {
		if (c->model_form != SIMPLEX_PAN_STANDARD ||
			c->oracle != SIMPLEX_PAN_NO_ORACLE)
			return 0;
	} else if (c->direction_method != SIMPLEX_PAN_CONE_DIRECTION ||
		   (c->oracle != SIMPLEX_PAN_EXPLICIT_ORACLE &&
			   c->oracle != SIMPLEX_PAN_IMPLICIT_ORACLE))
		return 0;
	if (c->pdhg_weight_set &&
		(!isfinite(c->pdhg_weight) ||
			c->pdhg_weight < SIMPLEX_PDHG_MANUAL_WEIGHT_MIN ||
			c->pdhg_weight > SIMPLEX_PDHG_MANUAL_WEIGHT_MAX))
		return 0;
	return c->pdhg_iterations >= 0 && c->dual_polish_iterations >= 0 &&
	       c->first_order_iterations >= 0 && c->center_epoch_limit >= 0 &&
	       c->nnls_block_limit >= 1 &&
	       c->nnls_block_limit <= SIMPLEX_PAN_NNLS_BLOCK_CAP;
}

int simplex_pan_config_default(struct simplex_PanConfig *c, int model_form,
	int direction_method, int oracle)
{
	lp_simplex_memset(c, 0, sizeof(*c));
	c->model_form = model_form;
	c->direction_method = direction_method;
	c->oracle = oracle;
	c->nnls_block_limit = c->singleton_crash = c->implicit_updates = 1;
	c->pdhg_iterations =
		model_form == SIMPLEX_PAN_NATIVE_BOUNDS &&
				oracle == SIMPLEX_PAN_IMPLICIT_ORACLE
			? SIMPLEX_PAN_PDHG_WORK_CAP
			: 0;
	c->dual_polish_iterations = SIMPLEX_PAN_DUAL_POLISH_WORK_CAP;
	c->restore = c->adaptive = 1;
	c->certificate_dual_repair = 1;
	return simplex_pan_config_valid(c);
}

void simplex_pan_config_init(
	struct simplex_PanConfig *c, const struct lp_simplex_Options *o)
{
	const char *value;

	lp_simplex_memset(c, 0, sizeof(*c));
	c->nnls_block_limit = 1;
	c->singleton_crash = 1;
	if (o && o->algorithm == lp_simplex_ALGORITHM_AUTO) {
		(void)simplex_pan_config_default(c, SIMPLEX_PAN_NATIVE_BOUNDS,
			SIMPLEX_PAN_CONE_DIRECTION,
			SIMPLEX_PAN_IMPLICIT_ORACLE);
		return;
	}
	c->model_form = getenv(SIMPLEX_ENV_PAN_NATIVE_BOUNDS)
				? SIMPLEX_PAN_NATIVE_BOUNDS
				: SIMPLEX_PAN_STANDARD;
	value = getenv(SIMPLEX_ENV_PAN_CONE);
	c->direction_method =
		value || c->model_form == SIMPLEX_PAN_NATIVE_BOUNDS
			? SIMPLEX_PAN_CONE_DIRECTION
			: SIMPLEX_PAN_BASIS_DIRECTION;
	c->oracle = c->direction_method == SIMPLEX_PAN_BASIS_DIRECTION
			    ? SIMPLEX_PAN_NO_ORACLE
			    : (value && strcmp(value, "implicit") == 0
					      ? SIMPLEX_PAN_IMPLICIT_ORACLE
					      : SIMPLEX_PAN_EXPLICIT_ORACLE);
	c->allow_oracle_recovery = 1;
	c->compatible_pricing =
		getenv(SIMPLEX_ENV_PAN_COMPATIBLE_PRICING) != NULL;
	c->stabilize = getenv(SIMPLEX_ENV_PAN_STABILIZE) != NULL;
	c->affine_multiplier =
		getenv(SIMPLEX_ENV_PAN_AFFINE_MULTIPLIER) != NULL;
	c->complete_start = getenv(SIMPLEX_ENV_PAN_COMPLETE_START) != NULL;
	c->face_guard = getenv(SIMPLEX_ENV_PAN_FACE_GUARD) != NULL;
	c->cone_batch = getenv(SIMPLEX_ENV_PAN_CONE_BATCH) != NULL;
	c->first_order_enabled = getenv(SIMPLEX_ENV_PAN_FIRST_ORDER) != NULL;
	c->dual_polish_enabled = getenv(SIMPLEX_ENV_PAN_DUAL_POLISH) != NULL;
	c->center_epoch_enabled = getenv(SIMPLEX_ENV_PAN_CENTER_EPOCH) != NULL;
	c->center_epoch_limit = setting_integer(SIMPLEX_ENV_PAN_CENTER_EPOCH);
	c->implicit_updates = getenv(SIMPLEX_ENV_PAN_IMPLICIT_UPDATES) != NULL;
	c->implicit_bulk = getenv(SIMPLEX_ENV_PAN_IMPLICIT_BULK) != NULL;
	c->pdhg_iterations = setting_integer(SIMPLEX_ENV_PAN_NATIVE_PDHG);
	c->restore = getenv(SIMPLEX_ENV_PAN_PDHG_RESTORE) != NULL;
	c->adaptive =
		c->restore || getenv(SIMPLEX_ENV_PAN_NATIVE_ADAPTIVE) != NULL;
	c->predict = getenv(SIMPLEX_ENV_PAN_NATIVE_PREDICT) != NULL;
	c->dual_start = getenv(SIMPLEX_ENV_PAN_DUAL_START) != NULL;
	c->dual_polish_iterations =
		setting_integer(SIMPLEX_ENV_PAN_DUAL_POLISH);
	c->first_order_iterations =
		setting_integer(SIMPLEX_ENV_PAN_FIRST_ORDER);
	c->newton_qr = getenv(SIMPLEX_ENV_PAN_NEWTON_QR) != NULL;
	c->dual_center = getenv(SIMPLEX_ENV_PAN_DUAL_CENTER) != NULL;
	value = getenv(SIMPLEX_ENV_PAN_NNLS_BLOCK);
	if (value) {
		c->nnls_block_limit = simplex_pan_integer_setting_value(value);
		if (c->nnls_block_limit < 1)
			c->nnls_block_limit = 1;
		if (c->nnls_block_limit > SIMPLEX_PAN_NNLS_BLOCK_CAP)
			c->nnls_block_limit = SIMPLEX_PAN_NNLS_BLOCK_CAP;
	}
	c->block_crash = getenv(SIMPLEX_ENV_PAN_BLOCK_CRASH) != NULL;
	c->singleton_crash =
		getenv(SIMPLEX_ENV_PAN_DISABLE_SINGLETON_CRASH) == NULL;
	c->full_refinement = getenv(SIMPLEX_ENV_PAN_FULL_REFINEMENT) != NULL;
	value = getenv(SIMPLEX_ENV_PAN_PDHG_WEIGHT);
	if (value) {
		char *end;
		double parsed;
		errno = 0;
		parsed = strtod(value, &end);
		if (end == value)
			return;
		while (isspace((unsigned char)*end))
			end++;
		if (!*end && errno != ERANGE && isfinite(parsed) &&
			parsed >= SIMPLEX_PDHG_MANUAL_WEIGHT_MIN &&
			parsed <= SIMPLEX_PDHG_MANUAL_WEIGHT_MAX) {
			c->pdhg_weight_set = 1;
			c->pdhg_weight = parsed;
		}
	}
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
