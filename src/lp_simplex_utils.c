/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#include <lp_simplex/lp_simplex.h>
#include <lp_simplex/lp_simplex_utils.h>
#include <assert.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void lp_simplex_linalg_daxpy(const int n, const double a, const double *x, const int incx, double *y, const int incy)
{
#if USE_BLAS
	daxpy_((int *)&n, (double *)&a, (double *)x, (int *)&incx, y, (int *)&incy);
#else
	int i = 0, j = 0;

	assert(x != NULL);
	assert(y != NULL);

	while (i < n && j < n) {
		y[i] += a * x[j];
		i += incy;
		j += incx;
	}
#endif
}

void lp_simplex_linalg_dscal(const int n, const double x, double *arr, const int inc)
{
#if USE_BLAS
	dscal_(&n, &x, arr, &inc);
#else
	int i;

	assert(arr != NULL);

	for (i = 0; i < n; i += inc)
		arr[i] *= x;
#endif
}

int is_in_arri(const int idx, const int *idxset, const int len)
{
	int i;

	for (i = 0; i < len; i++) {
		if (idx == idxset[i])
			return 1;
	}
	return 0;
}

int maxabs_arri(const int *arr, const int len, const int inc)
{
	int i, ele, maxv = __lp_simplex_NINF__;

	assert(inc >= 1);

	for (i = 0; i < len; i++) {
		ele = __lp_simplex_ABS__(arr[i]);
		if (ele > maxv)
			maxv = ele;
	}
	return maxv;
}

int argmaxabs_arrd(const double *arr, const int len, const int inc)
{
	int i, idx = 0;
	double ele, maxv = __lp_simplex_NINF__;

	assert(inc >= 1);

	for (i = 0; i < len; i += inc) {
		ele = __lp_simplex_ABS__(arr[i]);
		if (ele > maxv) {
			maxv = ele;
			idx = i;
		}
	}
	return idx;
}

double maxabs_arrd(const double *arr, const int len, const int inc)
{
	int i;
	double ele, maxv = __lp_simplex_NINF__;

	assert(inc >= 1);

	for (i = 0; i < len; i += inc) {
		ele = __lp_simplex_ABS__(arr[i]);
		if (ele > maxv)
			maxv = ele;
	}
	return maxv;
}

double maxabs_arrd_gap(const double *arr1, const double *arr2, const int len, const int inc)
{
	int i;
	double ele, maxv = __lp_simplex_NINF__;

	assert(inc >= 1);

	for (i = 0; i < len; i += inc) {
		ele = __lp_simplex_ABS__(arr1[i] - arr2[i]);
		if (ele > maxv)
			maxv = ele;
	}
	return maxv;
}

void lp_simplex_prt_arri(const int *arr, const int len, const int inc)
{
	int j;

	printf("[");
	for (j = 0; j < len; j += inc) {
		printf("%i", arr[j]);
		if (j + inc < len)
			printf(", ");
	}
	printf("]");
}

void lp_simplex_prt_arrl(const long *arr, const int len, const int inc)
{
	int j;

	printf("[");
	for (j = 0; j < len; j += inc) {
		printf("%li", arr[j]);
		if (j + inc < len)
			printf(", ");
	}
	printf("]");
}

void lp_simplex_prt_arrd(const double *arr, const int len, const int inc, const int sci)
{
	int j;
	double ele;

	printf("[");
	for (j = 0; j < len; j += inc) {
		ele = arr[j];
		if (-0. == ele)
			ele = 0.;
		if (sci)
			printf("%e", ele);
		else
			printf("%f", ele);
		if (j + inc < len)
			printf(", ");
	}
	printf("]");
}

void lp_simplex_prt_arrld(const long double *arr, const int len, const int inc, const int sci)
{
	long double ele;
	int j;

	printf("[");
	for (j = 0; j < len; j += inc) {
		ele = arr[j];
		if (-0. == ele)
			ele = 0.;
		if (sci)
			printf("%Le", ele);
		else
			printf("%Lf", ele);
		if (j + inc < len)
			printf(", ");
	}
	printf("]");
}

void lp_simplex_prt_matd(const double *mat, const int ld, const int nrow, const int ncol)
{
	int i, j;

	for (i = 0; i < nrow; i++) {
		for (j = 0; j < ncol; j++) {
			printf("%e", mat[j + i * ld]);
			if (j < ncol - 1)
				printf(", ");
		}
		if(i < nrow - 1)
			printf("\n");
	}
}

void *lp_simplex_malloc(size_t size)
{
	return malloc(size);
}

void lp_simplex_free(void *ptr)
{
	free(ptr);
}

void *lp_simplex_memset(void *str, int c, size_t n)
{
	return memset(str, c, n);
}

void *lp_simplex_memcpy(void *dest, const void *src, size_t n)
{
	return memcpy(dest, src, n);
}

int lp_simplex_memcmp(const void *str1, const void *str2, size_t n)
{
	return memcmp(str1, str2, n);
}

size_t lp_simplex_strcspn(const char *str1, const char *str2)
{
	return strcspn(str1, str2);
}

size_t lp_simplex_strlen(const char *str)
{
	return strlen(str);
}

double lp_simplex_atof(const char* str)
{
	return atof(str);
}


static struct lp_simplex_Model_LP *create_model(const int m, const int n)
{
	double *obj = NULL;
	double *coefficients = NULL;
	struct lp_simplex_LinearConstraint *constraints = NULL;
	struct lp_simplex_VariableBound *bounds = NULL;
	struct lp_simplex_Model_LP *model = NULL;
	int i;

	model = lp_simplex_malloc(sizeof(struct lp_simplex_Model_LP));
	if (model == NULL)
		return NULL;
	obj = lp_simplex_malloc(n * sizeof(double));
	if (obj == NULL) {
		lp_simplex_free(model);
		return NULL;
	}
	coefficients = lp_simplex_malloc(m * n * sizeof(double));
	if (coefficients == NULL) {
		lp_simplex_free(model);
		lp_simplex_free(obj);
		return NULL;
	}
	constraints = lp_simplex_malloc(m * sizeof(struct lp_simplex_LinearConstraint));
	if (constraints == NULL) {
		lp_simplex_free(model);
		lp_simplex_free(obj);
		lp_simplex_free(coefficients);
		return NULL;
	}
	bounds = lp_simplex_malloc(n * sizeof(struct lp_simplex_VariableBound));
	if (bounds == NULL) {
		lp_simplex_free(model);
		lp_simplex_free(obj);
		lp_simplex_free(coefficients);
		lp_simplex_free(constraints);
		return NULL;
	}
	model->m = m;
	model->n = n;
	model->objective = obj;
	model->coefficients = coefficients;
	model->constraints = constraints;
	model->bounds = bounds;
	lp_simplex_memset(coefficients, 0., m * n * sizeof(double));

	for (i = 0; i < m; i++) {
		constraints[i].coef = coefficients + i * n;
		lp_simplex_memset(constraints[i].name, '\0', 16);
	}
	for (i = 0; i < n; i++) {
		lp_simplex_memset(bounds[i].name, '\0', 16);
		bounds[i].lb = 0;
		bounds[i].ub = __lp_simplex_INF__;
		bounds[i].b_type = lp_simplex_BOUND_T_LO;
		bounds[i].v_type = lp_simplex_VAR_T_REAL;
	}
	return model;
}

void lp_simplex_lp_free(struct lp_simplex_Model_LP *model)
{
	if (model == NULL)
		return;
	if (model->coefficients)
		lp_simplex_free(model->coefficients);
	if (model->constraints)
		lp_simplex_free(model->constraints);
	if (model->bounds)
		lp_simplex_free(model->bounds);
	if (model->objective)
		lp_simplex_free(model->objective);
	lp_simplex_free(model);
}

void file_readline(FILE *f, char *line, const int n)
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

	lp_simplex_memset(value_str, '\0', 12);
	lp_simplex_memcpy(value_str, line + 24, 12);
	value = lp_simplex_atof(value_str);
	return value;
}

static double get_field_2_value(const char *line)
{
	char value_str[13];
	double value;

	lp_simplex_memset(value_str, '\0', 12);
	lp_simplex_memcpy(value_str, line + 49, 12);
	value = lp_simplex_atof(value_str);
	return value;
}

static void fill_model_coef(struct lp_simplex_Model_LP *model, const double value, const char *field_name, const int nvars)
{
	int i, m = model->m;

	for (i = 0; i < m; i++) {
		char *tmp = model->constraints[i].name;

		if (lp_simplex_memcmp(field_name, tmp, lp_simplex_strlen(tmp)) == 0) {
			model->constraints[i].coef[nvars - 1] = value;
			break;
		}
	}
}

static void fill_columns_to_model(struct lp_simplex_Model_LP *model, const char *obj_name, const char *field_name,
				  const double value, const int nvars)
{
	if (lp_simplex_memcmp(field_name, obj_name, lp_simplex_strlen(obj_name)) == 0)
		model->objective[nvars - 1] = value;
	else
		fill_model_coef(model, value, field_name, nvars);
}

static void fill_model_rhs(struct lp_simplex_Model_LP *model, const char *field_name, const double value)
{
	int i, m = model->m;

	for (i = 0; i < m; i++) {
		char *tmp = model->constraints[i].name;

		if (lp_simplex_memcmp(field_name, tmp, lp_simplex_strlen(tmp)) == 0) {
			model->constraints[i].rhs = value;
			break;
		}
	}
}

static int fill_model(const char *file, struct lp_simplex_Model_LP *model)
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
	lp_simplex_memset(last_name, '\0', 16);
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
			model->constraints[ncons].type = lp_simplex_CONS_T_LE;
			ncons++;
			break;
		case 'G':
			lp_simplex_memcpy(model->constraints[ncons].name, line + 4, 8);
			model->constraints[ncons].type = lp_simplex_CONS_T_GE;
			ncons++;
			break;
		case 'E':
			lp_simplex_memcpy(model->constraints[ncons].name, line + 4, 8);
			model->constraints[ncons].type = lp_simplex_CONS_T_EQ;
			ncons++;
			break;
		default:
			break;
		}
		break;
	case 2:  /* COLUMNS */
		if (lp_simplex_memcmp(last_name, line + 4, 8) != 0) {
			lp_simplex_memset(model->bounds[nvars].name, '\0', 16);
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

struct lp_simplex_Model_LP *lp_simplex_lp_readmps(const char *file)
{
	struct lp_simplex_Model_LP *model;
	int m, n;  /* number of constraints and variables */
	int n_sect_row = 0, n_sect_columns = 0;

	if (get_mps_info(file, &n_sect_row, &n_sect_columns) == lp_simplex_EXIT_FAILURE)
		return NULL;
	m = n_sect_row - 1;  /* the objective is also counted */
	n = n_sect_columns;
	model = create_model(m, n);

	if (model == NULL)
		return NULL;
	fill_model(file, model);
	return model;
}
