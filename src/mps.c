/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#include <lp_simplex/mps.h>
#include <lp_simplex/status.h>
#include "utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>



static void file_readline(FILE *f, char *line, const int n)
{
	lp_simplex_memset(line, '\0', n);
	fgets(line, n, f);
	line[lp_simplex_strcspn(line, "\r\n")] = '\0';
}


static int change_sect_code(const char *line, int *sect_code)
{
	int old_code = *sect_code;

	if (lp_simplex_memcmp(line, "ROWS", 4) == 0)
		*sect_code = 1;
	if (lp_simplex_memcmp(line, "COLUMNS", 7) == 0)
		*sect_code = 2;
	if (lp_simplex_memcmp(line, "RHS", 3) == 0)
		*sect_code = 3;
	if (lp_simplex_memcmp(line, "RANGES", 6) == 0)
		*sect_code = 4;
	if (lp_simplex_memcmp(line, "BOUNDS", 6) == 0)
		*sect_code = 5;
	if (lp_simplex_memcmp(line, "ENDATA", 6) == 0)
		*sect_code = 9;
	if (old_code == *sect_code)
		return 0;
	else
		return 1;
}


static int get_mps_info(const char *file, int *nrow, int *nvar)
{
	char line[128];
	char last_var[8];
	int sect_code = 0;

	FILE *f = fopen(file, "r");

	if (f == NULL) {
		printf("Cannot open file: \"%s\"\n", file);
		return lp_simplex_EXIT_FAILURE;
	}
	*nrow = 0;
	*nvar = 0;
	lp_simplex_memset(last_var, '\0', 8);
LOOP:
	file_readline(f, line, 128);

	if (change_sect_code(line, &sect_code))
		goto LOOP;
	switch (sect_code) {
	case 1:
		(*nrow)++;
		break;
	case 2:
		if (lp_simplex_memcmp(last_var, line + 4, 8) != 0) {
			(*nvar)++;
			lp_simplex_memcpy(last_var, line + 4, 8);
		}
		break;
	default:
		break;
	}
	if (feof(f))
		goto END;
	goto LOOP;
END:
	fclose(f);
	return lp_simplex_EXIT_SUCCESS;
}


static double get_filed_1_value(const char *line)
{
	char value_str[13];
	double value;

	lp_simplex_memset(value_str, '\0', 13);
	lp_simplex_memcpy(value_str, line + 24, 12);
	value = lp_simplex_atof(value_str);
	return value;
}


static double get_field_2_value(const char *line)
{
	char value_str[13];
	double value;

	lp_simplex_memset(value_str, '\0', 13);
	lp_simplex_memcpy(value_str, line + 49, 12);
	value = lp_simplex_atof(value_str);
	return value;
}


/* MPS names occupy eight columns and may be padded with spaces or NULs. */
static int mps_name_equal(const char *left, const char *right)
{
	int i;

	for (i = 0; i < 8; i++) {
		char a = left[i] == '\0' ? ' ' : left[i];
		char b = right[i] == '\0' ? ' ' : right[i];
		if (a != b)
			return 0;
	}
	return 1;
}


static void fill_model_coef(
		struct lp_Model *model, const double value,
		const char *field_name, const int nvars)
{
	int i, m = model->m;

	for (i = 0; i < m; i++) {
		char *tmp = model->constraints[i].name;

		if (mps_name_equal(field_name, tmp)) {
			model->constraints[i].coef[nvars - 1] = value;
			break;
		}
	}
}


static void fill_columns_to_model(
		struct lp_Model *model,
		const char *obj_name, const char *field_name,
		const double value, const int nvars)
{
	if (mps_name_equal(field_name, obj_name))
		model->objective[nvars - 1] = value;
	else
		fill_model_coef(model, value, field_name, nvars);
}


static void fill_model_rhs(
		struct lp_Model *model,
		const char *field_name, const double value)
{
	int i, m = model->m;

	for (i = 0; i < m; i++) {
		char *tmp = model->constraints[i].name;

		if (mps_name_equal(field_name, tmp)) {
			model->constraints[i].rhs = value;
			break;
		}
	}
}


static struct optm_VariableBound *find_model_bound(
		struct lp_Model *model, const char *name)
{
	int i;

	for (i = 0; i < model->n; i++) {
		if (mps_name_equal(name, model->bounds[i].name))
			return model->bounds + i;
	}
	return NULL;
}


static int fill_model_bound(struct lp_Model *model, const char *line)
{
	struct optm_VariableBound *bound = find_model_bound(model, line + 14);
	double value = get_filed_1_value(line);

	if (bound == NULL)
		return lp_simplex_EXIT_FAILURE;
	if (lp_simplex_memcmp(line + 1, "FR", 2) == 0) {
		bound->lb = __lp_simplex_NINF__;
		bound->ub = __lp_simplex_INF__;
		bound->b_type = optm_BOUND_T_FR;
	} else if (lp_simplex_memcmp(line + 1, "MI", 2) == 0) {
		bound->lb = __lp_simplex_NINF__;
		bound->b_type = bound->ub < __lp_simplex_INF__ ? optm_BOUND_T_UP : optm_BOUND_T_FR;
	} else if (lp_simplex_memcmp(line + 1, "PL", 2) == 0) {
		bound->ub = __lp_simplex_INF__;
		bound->b_type = bound->lb > __lp_simplex_NINF__ ? optm_BOUND_T_LO : optm_BOUND_T_FR;
	} else if (lp_simplex_memcmp(line + 1, "FX", 2) == 0) {
		bound->lb = value;
		bound->ub = value;
		bound->b_type = optm_BOUND_T_BS;
	} else if (lp_simplex_memcmp(line + 1, "BV", 2) == 0) {
		bound->lb = 0.;
		bound->ub = 1.;
		bound->b_type = optm_BOUND_T_BS;
		bound->v_type = optm_VAR_T_BIN;
	} else if (lp_simplex_memcmp(line + 1, "LO", 2) == 0 ||
		   lp_simplex_memcmp(line + 1, "LI", 2) == 0) {
		bound->lb = value;
		bound->b_type = bound->ub < __lp_simplex_INF__ ? optm_BOUND_T_BS : optm_BOUND_T_LO;
		if (line[1] == 'L' && line[2] == 'I')
			bound->v_type = optm_VAR_T_INT;
	} else if (lp_simplex_memcmp(line + 1, "UP", 2) == 0 ||
		   lp_simplex_memcmp(line + 1, "UI", 2) == 0) {
		bound->ub = value;
		bound->b_type = bound->lb > __lp_simplex_NINF__ ? optm_BOUND_T_BS : optm_BOUND_T_UP;
		if (line[1] == 'U' && line[2] == 'I')
			bound->v_type = optm_VAR_T_INT;
	} else {
		return lp_simplex_EXIT_FAILURE;
	}
	return lp_simplex_EXIT_SUCCESS;
}


static int fill_model(const char *file, struct lp_Model *model)
{
	char line[128];
	char obj_name[9];
	char *last_name = model->bounds->name;
	int sect_code = 0;
	int ncons = 0, nvars = 0;
	double value;

	FILE *f = fopen(file, "r");

	if (f == NULL) {
		printf("Cannot open file: \"%s\"\n", file);
		return lp_simplex_EXIT_FAILURE;
	}
	lp_simplex_memset(line, '\0', 128);
	lp_simplex_memset(obj_name, '\0', 9);
	lp_simplex_memset(last_name, '\0', 9);
LOOP:
	file_readline(f, line, 128);

	if (change_sect_code(line, &sect_code))
		goto LOOP;
	switch (sect_code) {
	case 1:  /* ROWS */
		switch (line[1]) {
		case 'N':
			lp_simplex_memset(obj_name, '\0', 8);
			lp_simplex_memcpy(obj_name, line + 4, 8);
			break;
		case 'L':
			lp_simplex_memcpy(model->constraints[ncons].name, line + 4, 8);
			model->constraints[ncons].type = optm_CONS_T_LE;
			ncons++;
			break;
		case 'G':
			lp_simplex_memcpy(model->constraints[ncons].name, line + 4, 8);
			model->constraints[ncons].type = optm_CONS_T_GE;
			ncons++;
			break;
		case 'E':
			lp_simplex_memcpy(model->constraints[ncons].name, line + 4, 8);
			model->constraints[ncons].type = optm_CONS_T_EQ;
			ncons++;
			break;
		default:
			break;
		}
		break;
	case 2:  /* COLUMNS */
		if (lp_simplex_memcmp(last_name, line + 4, 8) != 0) {
			lp_simplex_memset(model->bounds[nvars].name, '\0',
					  sizeof(model->bounds[nvars].name));
			lp_simplex_memcpy(model->bounds[nvars].name, line + 4, 8);
			last_name = model->bounds[nvars].name;
			nvars++;
		}
		value = get_filed_1_value(line);
		fill_columns_to_model(model, obj_name, line + 14, value, nvars);
		if (lp_simplex_strlen(line) < 40)
			goto LOOP;
		value = get_field_2_value(line);
		fill_columns_to_model(model, obj_name, line + 39, value, nvars);
		break;
	case 3:  /* RHS */
		value = get_filed_1_value(line);
		fill_model_rhs(model, line + 14, value);
		if (lp_simplex_strlen(line) < 40)
			goto LOOP;
		value = get_field_2_value(line);
		fill_model_rhs(model, line + 39, value);
		break;
	case 4:  /* RANGES are not representable by the public model structure. */
		fclose(f);
		return lp_simplex_EXIT_FAILURE;
	case 5:  /* BOUNDS */
		if (fill_model_bound(model, line) == lp_simplex_EXIT_FAILURE) {
			fclose(f);
			return lp_simplex_EXIT_FAILURE;
		}
		break;
	default:
		break;
	}
	if (feof(f))
		goto END;
	goto LOOP;
END:
	fclose(f);
	return lp_simplex_EXIT_SUCCESS;
}


struct lp_Model *lp_read_mps(const char *file)
{
	struct lp_Model *model;
	int m, n;  /* number of constraints and variables */
	int n_sect_row = 0, n_sect_columns = 0;

	if (get_mps_info(file, &n_sect_row, &n_sect_columns) == lp_simplex_EXIT_FAILURE)
		return NULL;
	m = n_sect_row - 1;  /* the objective is also counted */
	n = n_sect_columns;
	model = lp_model_create(m, n);

	if (model == NULL)
		return NULL;
	if (fill_model(file, model) == lp_simplex_EXIT_FAILURE) {
		lp_model_free(model);
		return NULL;
	}
	return model;
}
