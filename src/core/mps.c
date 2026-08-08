/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#include "utils.h"

#include <lp_simplex/mps.h>
#include <lp_simplex/status.h>

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


static int mps_has_second_field(const char *line);
static unsigned long mps_name_hash(const char *name);


static int get_mps_info(
		const char *file, int *nrow, int *nvar, int *entry_capacity,
		int *range_fields)
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
	*entry_capacity = 0;
	*range_fields = 0;
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
		(*entry_capacity)++;
		if (mps_has_second_field(line))
			(*entry_capacity)++;
		if (lp_simplex_memcmp(last_var, line + 4, 8) != 0) {
			(*nvar)++;
			lp_simplex_memcpy(last_var, line + 4, 8);
		}
		break;
	case 4:
		(*range_fields)++;
		if (mps_has_second_field(line))
			(*range_fields)++;
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


static int mps_has_second_field(const char *line)
{
	int i, fields = 0, inside = 0;
	for (i = 0; line[i] != '\0'; i++) {
		int whitespace = line[i] == ' ' || line[i] == '\t';
		if (!whitespace && !inside)
			fields++;
		inside = !whitespace;
	}
	return fields >= 5;
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


struct mps_RangeScan {
	int original_rows;
	int ranged_rows;
	int duplicate_capacity;
	unsigned char *ranged;
	unsigned char *row_type;
};


/* Identify ranged rows before allocating the final sparse model.  A ranged
 * MPS row is represented internally by two ordinary one-sided constraints,
 * so its matrix coefficients have to be allocated twice. */
static int mps_scan_ranges(
		const char *file, const int original_rows,
		struct mps_RangeScan *scan)
{
	char line[128];
	char (*name)[9] = NULL;
	int *degree = NULL;
	int *slot = NULL;
	FILE *stream = NULL;
	int section = 0, rows = 0, i, table_capacity = 16;
	int state = lp_simplex_EXIT_FAILURE;

	lp_simplex_memset(scan, 0, sizeof(*scan));
	scan->original_rows = original_rows;
	name = (char (*)[9])lp_simplex_malloc(
		(size_t)original_rows * sizeof(*name));
	degree = (int *)lp_simplex_malloc((size_t)original_rows * sizeof(int));
	scan->ranged = (unsigned char *)lp_simplex_malloc(
		(size_t)original_rows * sizeof(unsigned char));
	scan->row_type = (unsigned char *)lp_simplex_malloc(
		(size_t)original_rows * sizeof(unsigned char));
	while (table_capacity < original_rows * 2)
		table_capacity *= 2;
	slot = (int *)lp_simplex_malloc((size_t)table_capacity * sizeof(int));
	if (name == NULL || degree == NULL || scan->ranged == NULL ||
	    scan->row_type == NULL || slot == NULL)
		goto finish;
	lp_simplex_memset(name, 0, (size_t)original_rows * sizeof(*name));
	lp_simplex_memset(degree, 0, (size_t)original_rows * sizeof(int));
	lp_simplex_memset(scan->ranged, 0,
		(size_t)original_rows * sizeof(unsigned char));
	for (i = 0; i < table_capacity; i++)
		slot[i] = -1;
	stream = fopen(file, "r");
	if (stream == NULL)
		goto finish;
	while (!feof(stream)) {
		file_readline(stream, line, (int)sizeof(line));
		if (change_sect_code(line, &section))
			continue;
		if (section == 1) {
			char kind = line[1] == ' ' ? line[2] : line[1];
			unsigned long position;
			if (kind == 'L' || kind == 'G' || kind == 'E') {
				if (rows >= original_rows)
					goto finish;
				lp_simplex_memcpy(name[rows], line + 4, 8);
				scan->row_type[rows] = (unsigned char)(kind == 'L'
					? optm_CONS_T_LE : kind == 'G'
					? optm_CONS_T_GE : optm_CONS_T_EQ);
				position = mps_name_hash(name[rows]) &
					(unsigned long)(table_capacity - 1);
				while (slot[position] >= 0)
					position = (position + 1) &
						(unsigned long)(table_capacity - 1);
				slot[position] = rows++;
			}
		} else if (section == 2 || section == 4) {
			int fields = mps_has_second_field(line) ? 2 : 1;
			int field;
			for (field = 0; field < fields; field++) {
				const char *row_name = line + (field == 0 ? 14 : 39);
				unsigned long position = mps_name_hash(row_name) &
					(unsigned long)(table_capacity - 1);
				i = -1;
				while (slot[position] >= 0) {
					int candidate = slot[position];
					if (mps_name_equal(row_name, name[candidate])) {
						i = candidate;
						break;
					}
					position = (position + 1) &
						(unsigned long)(table_capacity - 1);
				}
				if (i < 0) {
					/* Objective cards in COLUMNS are not row entries. */
					if (section == 4)
						goto finish;
					continue;
				}
				if (section == 2)
					degree[i]++;
				else if (!scan->ranged[i]) {
					scan->ranged[i] = 1;
					scan->ranged_rows++;
				}
			}
		}
	}
	if (rows != original_rows)
		goto finish;
	for (i = 0; i < original_rows; i++)
		if (scan->ranged[i])
			scan->duplicate_capacity += degree[i];
	state = lp_simplex_EXIT_SUCCESS;
finish:
	if (stream != NULL)
		fclose(stream);
	lp_simplex_free(name);
	lp_simplex_free(degree);
	lp_simplex_free(slot);
	if (state == lp_simplex_EXIT_FAILURE) {
		lp_simplex_free(scan->ranged);
		lp_simplex_free(scan->row_type);
		scan->ranged = NULL;
		scan->row_type = NULL;
	}
	return state;
}


static int apply_model_bound(
		struct optm_VariableBound *bound, const char *line)
{
	double value = get_filed_1_value(line);
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






struct mps_NameTable {
	int *slot;
	int capacity;
};


struct mps_SparseEntry {
	int row;
	double value;
};


static int mps_compare_sparse_entry(const void *left, const void *right)
{
	const struct mps_SparseEntry *a =
		(const struct mps_SparseEntry *)left;
	const struct mps_SparseEntry *b =
		(const struct mps_SparseEntry *)right;
	return a->row < b->row ? -1 : a->row != b->row;
}


static int mps_sort_csc_columns(struct lp_Model *model)
{
	struct mps_SparseEntry *buffer;
	int j, k, maximum = 0;
	for (j = 0; j < model->n; j++)
		maximum = __lp_simplex_MAX__(maximum,
			model->column_start[j + 1] - model->column_start[j]);
	if (maximum <= 1)
		return lp_simplex_EXIT_SUCCESS;
	buffer = (struct mps_SparseEntry *)lp_simplex_malloc(
		(size_t)maximum * sizeof(*buffer));
	if (buffer == NULL)
		return lp_simplex_EXIT_FAILURE;
	for (j = 0; j < model->n; j++) {
		int count = model->column_start[j + 1] - model->column_start[j];
		for (k = 0; k < count; k++) {
			int position = model->column_start[j] + k;
			buffer[k].row = model->row_index[position];
			buffer[k].value = model->value[position];
		}
		qsort(buffer, (size_t)count, sizeof(*buffer),
			mps_compare_sparse_entry);
		for (k = 0; k < count; k++) {
			int position = model->column_start[j] + k;
			model->row_index[position] = buffer[k].row;
			model->value[position] = buffer[k].value;
		}
	}
	lp_simplex_free(buffer);
	return lp_simplex_EXIT_SUCCESS;
}


static unsigned long mps_name_hash(const char *name)
{
	unsigned long hash = 2166136261UL;
	int i;
	for (i = 0; i < 8; i++) {
		unsigned char value = (unsigned char)(name[i] == '\0' ? ' ' : name[i]);
		hash = (hash ^ value) * 16777619UL;
	}
	return hash;
}


static int mps_name_table_create(
		struct mps_NameTable *table, const int entries)
{
	int i;
	table->capacity = 16;
	while (table->capacity < entries * 2)
		table->capacity *= 2;
	table->slot = (int *)lp_simplex_malloc(
		(size_t)table->capacity * sizeof(int));
	if (table->slot == NULL)
		return lp_simplex_EXIT_FAILURE;
	for (i = 0; i < table->capacity; i++)
		table->slot[i] = -1;
	return lp_simplex_EXIT_SUCCESS;
}


static void mps_name_table_insert(
		struct mps_NameTable *table, const char *name, const int index)
{
	unsigned long position = mps_name_hash(name) &
		(unsigned long)(table->capacity - 1);
	while (table->slot[position] >= 0)
		position = (position + 1) & (unsigned long)(table->capacity - 1);
	table->slot[position] = index;
}


static int mps_name_table_find(
		const struct mps_NameTable *table, const char *name,
		const struct lp_Model *model, const int rows)
{
	unsigned long position = mps_name_hash(name) &
		(unsigned long)(table->capacity - 1);
	while (table->slot[position] >= 0) {
		int index = table->slot[position];
		const char *stored = rows ? model->constraints[index].name
			: model->bounds[index].name;
		if (mps_name_equal(name, stored))
			return index;
		position = (position + 1) &
			(unsigned long)(table->capacity - 1);
	}
	return -1;
}


static int mps_sparse_add_field(
		struct lp_Model *model, const struct mps_NameTable *rows,
		const char *objective_name, const char *field_name,
		const double value, const int column, const int *range_pair,
		int *nonzeros,
		const int capacity)
{
	int row, targets[2], count, target_index;
	if (mps_name_equal(field_name, objective_name)) {
		model->objective[column] = value;
		return lp_simplex_EXIT_SUCCESS;
	}
	row = mps_name_table_find(rows, field_name, model, 1);
	if (row < 0)
		return lp_simplex_EXIT_FAILURE;
	targets[0] = row;
	count = 1;
	if (range_pair != NULL && range_pair[row] >= 0)
		targets[count++] = range_pair[row];
	for (target_index = 0; target_index < count; target_index++) {
		int target = targets[target_index];
		int k;
		int found = 0;
		/* Repeated (row,column) cards retain last-card-wins semantics. */
		for (k = model->column_start[column]; k < *nonzeros; k++)
			if (model->row_index[k] == target) {
				model->value[k] = value;
				found = 1;
				break;
			}
		if (found)
			continue;
		if (*nonzeros >= capacity)
			return lp_simplex_EXIT_FAILURE;
		model->row_index[*nonzeros] = target;
		model->value[(*nonzeros)++] = value;
	}
	return lp_simplex_EXIT_SUCCESS;
}


static void mps_apply_range(
		struct lp_Model *model, const int row, const int pair,
		const int original_type, const double range)
{
	double rhs = model->constraints[row].rhs;
	double width = __lp_simplex_ABS__(range);
	int type = original_type;

	if (type == optm_CONS_T_LE) {
		model->constraints[row].type = optm_CONS_T_LE;
		model->constraints[row].rhs = rhs;
		model->constraints[pair].type = optm_CONS_T_GE;
		model->constraints[pair].rhs = rhs - width;
	} else if (type == optm_CONS_T_GE) {
		model->constraints[row].type = optm_CONS_T_GE;
		model->constraints[row].rhs = rhs;
		model->constraints[pair].type = optm_CONS_T_LE;
		model->constraints[pair].rhs = rhs + width;
	} else if (range >= 0.) {
		model->constraints[row].type = optm_CONS_T_GE;
		model->constraints[row].rhs = rhs;
		model->constraints[pair].type = optm_CONS_T_LE;
		model->constraints[pair].rhs = rhs + width;
	} else {
		model->constraints[row].type = optm_CONS_T_LE;
		model->constraints[row].rhs = rhs;
		model->constraints[pair].type = optm_CONS_T_GE;
		model->constraints[pair].rhs = rhs - width;
	}
}


static int mps_build_csr(
		struct lp_Model *model, const int nonzeros)
{
	int i, j, k;
	int *cursor = (int *)lp_simplex_malloc(
		(size_t)model->m * sizeof(int));
	if (cursor == NULL)
		return lp_simplex_EXIT_FAILURE;
	lp_simplex_memset(model->row_start, 0,
		(size_t)(model->m + 1) * sizeof(int));
	for (k = 0; k < nonzeros; k++)
		model->row_start[model->row_index[k] + 1]++;
	for (i = 0; i < model->m; i++) {
		model->row_start[i + 1] += model->row_start[i];
		cursor[i] = model->row_start[i];
	}
	for (j = 0; j < model->n; j++)
		for (k = model->column_start[j]; k < model->column_start[j + 1]; k++) {
			i = model->row_index[k];
			model->column_index[cursor[i]] = j;
			model->row_value[cursor[i]++] = model->value[k];
		}
	lp_simplex_free(cursor);
	model->nnz = nonzeros;
	return lp_simplex_EXIT_SUCCESS;
}


static int fill_sparse_model(
		const char *file, struct lp_Model *model, const int entry_capacity,
		const struct mps_RangeScan *range_scan)
{
	char line[128], objective_name[9], last_column[9];
	struct mps_NameTable rows, columns;
	FILE *stream;
	int section = 0, row_count = 0, column = -1, nonzeros = 0;
	int *range_pair = NULL;
	int state = lp_simplex_EXIT_FAILURE;
	rows.slot = columns.slot = NULL;
	range_pair = (int *)lp_simplex_malloc(
		(size_t)range_scan->original_rows * sizeof(int));
	if (range_pair == NULL)
		goto finish;
	{
		int original, duplicate = range_scan->original_rows;
		for (original = 0; original < range_scan->original_rows; original++)
			if (range_scan->ranged[original])
				range_pair[original] = duplicate++;
			else
				range_pair[original] = -1;
	}
	if (mps_name_table_create(&rows, model->m) == lp_simplex_EXIT_FAILURE ||
	    mps_name_table_create(&columns, model->n) == lp_simplex_EXIT_FAILURE)
		goto finish;
	stream = fopen(file, "r");
	if (stream == NULL)
		goto finish;
	lp_simplex_memset(objective_name, 0, sizeof(objective_name));
	lp_simplex_memset(last_column, 0, sizeof(last_column));
	while (!feof(stream)) {
		double value;
		file_readline(stream, line, (int)sizeof(line));
		if (change_sect_code(line, &section))
			continue;
		if (section == 1) {
			char kind = line[1] == ' ' ? line[2] : line[1];
			if (kind == 'N') {
				lp_simplex_memset(objective_name, 0,
					sizeof(objective_name));
				lp_simplex_memcpy(objective_name, line + 4, 8);
			} else if (kind == 'L' || kind == 'G' || kind == 'E') {
				if (row_count >= range_scan->original_rows)
					goto close;
				lp_simplex_memcpy(model->constraints[row_count].name,
					line + 4, 8);
				model->constraints[row_count].type = kind == 'L'
					? optm_CONS_T_LE : kind == 'G'
					? optm_CONS_T_GE : optm_CONS_T_EQ;
				mps_name_table_insert(&rows,
					model->constraints[row_count].name, row_count);
				if (range_pair[row_count] >= 0) {
					int pair = range_pair[row_count];
					lp_simplex_memcpy(model->constraints[pair].name,
						line + 4, 8);
					model->constraints[pair].type =
						model->constraints[row_count].type;
				}
				row_count++;
			}
		} else if (section == 2) {
			if (strstr(line, "'MARKER'") != NULL)
				continue;
			if (!mps_name_equal(last_column, line + 4)) {
				if (column >= 0)
					model->column_start[column + 1] = nonzeros;
				column++;
				if (column >= model->n)
					goto close;
				lp_simplex_memset(model->bounds[column].name, 0,
					sizeof(model->bounds[column].name));
				lp_simplex_memcpy(model->bounds[column].name, line + 4, 8);
				lp_simplex_memcpy(last_column, line + 4, 8);
				mps_name_table_insert(&columns,
					model->bounds[column].name, column);
				model->column_start[column] = nonzeros;
			}
			value = get_filed_1_value(line);
			if (mps_sparse_add_field(model, &rows, objective_name,
					line + 14, value, column, range_pair, &nonzeros,
					entry_capacity) == lp_simplex_EXIT_FAILURE)
				goto close;
			if (mps_has_second_field(line)) {
				value = get_field_2_value(line);
				if (mps_sparse_add_field(model, &rows, objective_name,
						line + 39, value, column, range_pair, &nonzeros,
						entry_capacity) == lp_simplex_EXIT_FAILURE)
					goto close;
			}
		} else if (section == 3) {
			int row;
			if (mps_name_equal(line + 14, objective_name)) {
				model->objective_offset = -get_filed_1_value(line);
			} else {
				row = mps_name_table_find(&rows, line + 14, model, 1);
				if (row < 0)
					goto close;
				model->constraints[row].rhs = get_filed_1_value(line);
				if (range_pair[row] >= 0)
					model->constraints[range_pair[row]].rhs =
						model->constraints[row].rhs;
			}
			if (mps_has_second_field(line)) {
				if (mps_name_equal(line + 39, objective_name)) {
					model->objective_offset = -get_field_2_value(line);
				} else {
					row = mps_name_table_find(&rows, line + 39, model, 1);
					if (row < 0)
						goto close;
					model->constraints[row].rhs = get_field_2_value(line);
					if (range_pair[row] >= 0)
						model->constraints[range_pair[row]].rhs =
							model->constraints[row].rhs;
				}
			}
		} else if (section == 4) {
			int fields = mps_has_second_field(line) ? 2 : 1;
			int field;
			for (field = 0; field < fields; field++) {
				const char *row_name = line + (field == 0 ? 14 : 39);
				int row = mps_name_table_find(&rows, row_name, model, 1);
				double range = field == 0 ? get_filed_1_value(line)
					: get_field_2_value(line);
				if (row < 0 || range_pair[row] < 0)
					goto close;
				mps_apply_range(model, row, range_pair[row],
					range_scan->row_type[row], range);
			}
		} else if (section == 5) {
			int bound_column = mps_name_table_find(
				&columns, line + 14, model, 0);
			if (bound_column < 0 || apply_model_bound(
					model->bounds + bound_column, line) ==
				    lp_simplex_EXIT_FAILURE)
				goto close;
		}
	}
	if (column + 1 != model->n ||
	    row_count != range_scan->original_rows)
		goto close;
	model->column_start[model->n] = nonzeros;
	if (mps_sort_csc_columns(model) == lp_simplex_EXIT_FAILURE)
		goto close;
	state = mps_build_csr(model, nonzeros);
close:
	if (state == lp_simplex_EXIT_FAILURE &&
	    getenv("LP_SIMPLEX_PROFILE") != NULL)
		fprintf(stderr, "mps sparse parse failure: section=%d rows=%d/%d "
			"columns=%d/%d nonzeros=%d/%d line=%.80s\n",
			section, row_count, model->m, column + 1, model->n,
			nonzeros, entry_capacity, line);
	fclose(stream);
finish:
	lp_simplex_free(range_pair);
	lp_simplex_free(rows.slot);
	lp_simplex_free(columns.slot);
	return state;
}


struct lp_Model *lp_read_mps(const char *file)
{
	struct lp_Model *model;
	struct mps_RangeScan range_scan;
	int m, n;  /* number of constraints and variables */
	int n_sect_row = 0, n_sect_columns = 0;
	int entry_capacity = 0, range_fields = 0;

	if (get_mps_info(file, &n_sect_row, &n_sect_columns,
			&entry_capacity, &range_fields) == lp_simplex_EXIT_FAILURE)
		return NULL;
	m = n_sect_row - 1;  /* the objective is also counted */
	n = n_sect_columns;
	if (range_fields > 0) {
		if (mps_scan_ranges(file, m, &range_scan) ==
		    lp_simplex_EXIT_FAILURE)
			return NULL;
	} else {
		lp_simplex_memset(&range_scan, 0, sizeof(range_scan));
		range_scan.original_rows = m;
		range_scan.ranged = (unsigned char *)lp_simplex_malloc(
			(size_t)m * sizeof(unsigned char));
		if (range_scan.ranged == NULL)
			return NULL;
		lp_simplex_memset(range_scan.ranged, 0,
			(size_t)m * sizeof(unsigned char));
	}
	m += range_scan.ranged_rows;
	entry_capacity += range_scan.duplicate_capacity;
	model = lp_model_create_sparse(m, n, entry_capacity);

	if (model == NULL) {
		lp_simplex_free(range_scan.ranged);
		lp_simplex_free(range_scan.row_type);
		return NULL;
	}
	if (fill_sparse_model(file, model, entry_capacity, &range_scan) ==
	    lp_simplex_EXIT_FAILURE) {
		lp_model_free(model);
		lp_simplex_free(range_scan.ranged);
		lp_simplex_free(range_scan.row_type);
		return NULL;
	}
	lp_simplex_free(range_scan.ranged);
	lp_simplex_free(range_scan.row_type);
	return model;
}
