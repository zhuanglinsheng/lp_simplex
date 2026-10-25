/**
 * @file test_sparse_lu_storage.c
 * @brief Sparse LU ownership, failed builds and atomic FT allocation checks.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>

static void *live[4096];
static int active, bad_free, requests, fail_at;

void *lu_storage_malloc(size_t bytes)
{
	void *p;
	int i;
	requests++;
	if (requests == fail_at || bytes > 1024 * 1024)
		return NULL;
	p = malloc(bytes ? bytes : 1);
	if (!p)
		return NULL;
	memset(p, 0xa5, bytes);
	for (i = 0; i < 4096; i++)
		if (!live[i]) {
			live[i] = p;
			active++;
			return p;
		}
	free(p);
	return NULL;
}

void lu_storage_free(void *p)
{
	int i;
	if (!p)
		return;
	for (i = 0; i < 4096; i++)
		if (live[i] == p) {
			live[i] = NULL;
			active--;
			free(p);
			return;
		}
	bad_free++;
}

#define lp_simplex_malloc lu_storage_malloc
#define lp_simplex_free lu_storage_free
#define simplex_sparse_lu_create storage_lu_create
#define simplex_sparse_lu_destroy storage_lu_destroy
#define simplex_sparse_lu_factorize storage_lu_factorize
#define simplex_sparse_lu_factorize_submatrix storage_lu_factorize_submatrix
#define simplex_sparse_lu_solve storage_lu_solve
#define simplex_sparse_lu_solve_pair storage_lu_solve_pair
#define simplex_sparse_lu_ft_update storage_lu_ft_update
#define simplex_sparse_lu_ft_active storage_lu_ft_active
#ifndef LU_STORAGE_SOURCE
#define LU_STORAGE_SOURCE "../src/basis/simplex_sparse_lu.c"
#endif
#include LU_STORAGE_SOURCE
#undef lp_simplex_malloc
#undef lp_simplex_free

#define CHECK(c) \
	do { \
		if (!(c)) { \
			fprintf(stderr, "LU storage line %d: %s\n", __LINE__, \
				#c); \
			return 1; \
		} \
	} while (0)
#ifdef LU_STORAGE_OLD
#define ROW_RESERVE(f, r, n) sparse_row_reserve(r, n)
#define COLUMN_APPEND(f, c, r) sparse_column_rows_append(c, r)
#else
#define ROW_RESERVE(f, r, n) sparse_row_reserve(f, r, n)
#define COLUMN_APPEND(f, c, r) sparse_column_rows_append(f, c, r)
#endif

static int constructors(void)
{
	struct simplex_SparseLu factor;
	int failure, allocations;
	requests = fail_at = 0;
	CHECK(!storage_lu_create(&factor, 8));
	allocations = requests;
	CHECK(allocations == 8);
	storage_lu_destroy(&factor);
	for (failure = 1; failure <= allocations; failure++) {
		requests = 0;
		fail_at = failure;
		CHECK(storage_lu_create(&factor, 8));
#ifndef LU_STORAGE_OLD
		CHECK(factor.storage_failed && !factor.factor_valid);
#endif
		storage_lu_destroy(&factor);
		CHECK(!active && !bad_free);
	}
	fail_at = 0;
	return 0;
}

static int reserves(void)
{
	struct simplex_SparseLu factor;
	struct simplex_SparseRow row;
	struct simplex_SparseColumnRows column;
	int failure, i, *saved_column;
	double *saved_value;
	for (failure = 1; failure <= 2; failure++) {
		memset(&factor, 0, sizeof(factor));
		memset(&row, 0, sizeof(row));
		fail_at = 0;
		CHECK(!ROW_RESERVE(&factor, &row, 4));
		row.count = 4;
		for (i = 0; i < 4; i++) {
			row.column[i] = i;
			row.value[i] = i + .5;
		}
		saved_column = row.column;
		saved_value = row.value;
		requests = 0;
		fail_at = failure;
		CHECK(ROW_RESERVE(&factor, &row, 5));
#ifndef LU_STORAGE_OLD
		CHECK(factor.storage_failed);
#endif
		CHECK(row.column == saved_column && row.value == saved_value &&
			row.capacity == 4 && row.count == 4);
		for (i = 0; i < 4; i++)
			CHECK(row.column[i] == i && row.value[i] == i + .5);
		fail_at = 0;
		CHECK(!ROW_RESERVE(&factor, &row, 5));
		sparse_row_destroy(&row);
		CHECK(!active && !bad_free);
	}
	memset(&factor, 0, sizeof(factor));
	memset(&column, 0, sizeof(column));
	for (i = 0; i < 8; i++)
		CHECK(!COLUMN_APPEND(&factor, &column, i));
	requests = 0;
	fail_at = 1;
	saved_column = column.row;
	CHECK(COLUMN_APPEND(&factor, &column, 8));
	CHECK(column.row == saved_column && column.count == 8 &&
		column.capacity == 8);
	for (i = 0; i < 8; i++)
		CHECK(column.row[i] == i);
	fail_at = 0;
	CHECK(!COLUMN_APPEND(&factor, &column, 8));
	sparse_column_rows_destroy(&column);
	CHECK(!active && !bad_free);
	return 0;
}

struct fixture {
	struct simplex_CscMatrix matrix;
	int starts[9], rows[8], basis[8];
	double values[8], dense[64];
};

static void fixture_init(struct fixture *f, int structural)
{
	int i;
	memset(f, 0, sizeof(*f));
	f->matrix.rows = 8;
	f->matrix.columns = structural ? 8 : 0;
	f->matrix.nonzeros = structural ? 8 : 0;
	f->matrix.column_start = f->starts;
	f->matrix.row_index = f->rows;
	f->matrix.value = f->values;
	for (i = 0; i < 8; i++) {
		f->starts[i] = i;
		f->rows[i] = f->basis[i] = i;
		f->values[i] = 1.;
		f->dense[8 * i + i] = structural ? 1. : -1.;
	}
	f->starts[8] = 8;
	if (!structural)
		f->starts[0] = 0;
}

static int residual(struct simplex_SparseLu *factor, const double *matrix)
{
	double x[8], y[8];
	int i, j;
	for (i = 0; i < 8; i++)
		x[i] = y[i] = i + 1.;
	CHECK(!storage_lu_solve(factor, x, 0));
	CHECK(!storage_lu_solve(factor, y, 1));
	for (i = 0; i < 8; i++) {
		double a = -(i + 1.), b = a;
		for (j = 0; j < 8; j++) {
			a += matrix[i * 8 + j] * x[j];
			b += matrix[j * 8 + i] * y[j];
		}
		CHECK(fabs(a) < 1e-10 && fabs(b) < 1e-10);
	}
	return 0;
}

static int build(
	struct simplex_SparseLu *factor, struct fixture *f, int submatrix)
{
	return submatrix ? storage_lu_factorize_submatrix(
				   factor, &f->matrix, f->basis, f->rows)
			 : storage_lu_factorize(factor, &f->matrix,
				   f->matrix.columns, f->basis);
}

static int builds(int submatrix)
{
	struct fixture f;
	struct simplex_SparseLu factor;
	int allocations, failure, i;
	double first[8], second[8], saved[8];
	fixture_init(&f, submatrix);
	requests = fail_at = 0;
	CHECK(!storage_lu_create(&factor, 8));
	requests = 0;
	CHECK(!build(&factor, &f, submatrix));
	allocations = requests;
	CHECK(allocations > 0 && !residual(&factor, f.dense));
#ifndef LU_STORAGE_OLD
	/* A failed rebuild must revoke an already usable packed factor, too. */
	requests = 0;
	fail_at = 1;
	CHECK(build(&factor, &f, submatrix));
	CHECK(factor.storage_failed && !factor.factor_valid);
	for (i = 0; i < 8; i++)
		first[i] = second[i] = saved[i] = i + 1.;
	CHECK(storage_lu_solve(&factor, first, 0));
	CHECK(storage_lu_solve_pair(&factor, first, second));
	CHECK(!memcmp(first, saved, sizeof(saved)) &&
		!memcmp(second, saved, sizeof(saved)));
	fail_at = 0;
	CHECK(!build(&factor, &f, submatrix));
	CHECK(!residual(&factor, f.dense));
#endif
	storage_lu_destroy(&factor);
	for (failure = 1; failure <= allocations; failure++) {
		fail_at = 0;
		CHECK(!storage_lu_create(&factor, 8));
		requests = 0;
		fail_at = failure;
		CHECK(build(&factor, &f, submatrix));
#ifndef LU_STORAGE_OLD
		CHECK(factor.storage_failed && !factor.factor_valid);
		for (i = 0; i < 8; i++)
			first[i] = second[i] = saved[i] = i + 1.;
		CHECK(storage_lu_solve(&factor, first, 0));
		CHECK(storage_lu_solve_pair(&factor, first, second));
		CHECK(!memcmp(first, saved, sizeof(saved)) &&
			!memcmp(second, saved, sizeof(saved)));
#else
		(void)i;
		(void)first;
		(void)second;
		(void)saved;
#endif
		fail_at = 0;
		CHECK(!build(&factor, &f, submatrix));
		CHECK(!residual(&factor, f.dense));
		storage_lu_destroy(&factor);
		CHECK(!active && !bad_free);
	}
	printf("build submatrix=%d allocations=%d\n", submatrix, allocations);
	return 0;
}

static void proposal(int stage, double *direction)
{
	int i;
	for (i = 0; i < 8; i++)
		direction[i] = stage ? .25 : .5;
	direction[stage] = 1.;
}

static void commit_matrix(double *matrix, int column, const double *direction)
{
	double entering[8];
	int i, j;
	for (i = 0; i < 8; i++) {
		entering[i] = 0.;
		for (j = 0; j < 8; j++)
			entering[i] += matrix[i * 8 + j] * direction[j];
	}
	for (i = 0; i < 8; i++)
		matrix[i * 8 + column] = entering[i];
}

static int ft_case(int stage, int failure, int *allocations)
{
	struct fixture f;
	struct simplex_SparseLu factor;
	double direction[8];
	int count, order[8];
	int sizes[8], columns[64], i;
	double entries[64];
	fixture_init(&f, 0);
	fail_at = 0;
	CHECK(!storage_lu_create(&factor, 8));
	CHECK(!build(&factor, &f, 0));
	if (stage) {
		proposal(0, direction);
		CHECK(!storage_lu_ft_update(&factor, 0, direction));
		commit_matrix(f.dense, 0, direction);
	}
	CHECK(!residual(&factor, f.dense));
	proposal(stage, direction);
	count = factor.ft_update_count;
	memcpy(order, factor.ft_order, sizeof(order));
	for (i = 0; i < 8; i++) {
		sizes[i] = factor.ft_column[i].count;
		memcpy(columns + 8 * i, factor.ft_column[i].column,
			(size_t)sizes[i] * sizeof(int));
		memcpy(entries + 8 * i, factor.ft_column[i].value,
			(size_t)sizes[i] * sizeof(double));
	}
	requests = 0;
	fail_at = failure;
	if (failure) {
		CHECK(storage_lu_ft_update(&factor, stage, direction) ==
			lp_simplex_EXIT_FAILURE);
#ifndef LU_STORAGE_OLD
		CHECK(factor.storage_failed && factor.factor_valid);
#endif
		CHECK(factor.ft_update_count == count &&
			!memcmp(order, factor.ft_order, sizeof(order)));
		for (i = 0; i < 8; i++)
			CHECK(factor.ft_column[i].count == sizes[i] &&
				!memcmp(columns + 8 * i,
					factor.ft_column[i].column,
					(size_t)sizes[i] * sizeof(int)) &&
				!memcmp(entries + 8 * i,
					factor.ft_column[i].value,
					(size_t)sizes[i] * sizeof(double)));
		CHECK(!residual(&factor, f.dense));
		fail_at = 0;
		CHECK(!storage_lu_ft_update(&factor, stage, direction));
	} else {
		CHECK(!storage_lu_ft_update(&factor, stage, direction));
		*allocations = requests;
		CHECK(*allocations > 0);
	}
	commit_matrix(f.dense, stage, direction);
	CHECK(!residual(&factor, f.dense));
	storage_lu_destroy(&factor);
	CHECK(!active && !bad_free);
	return 0;
}

static int ft_failures(void)
{
	int stage, failure, allocations = 0;
	for (stage = 0; stage < 2; stage++) {
		CHECK(!ft_case(stage, 0, &allocations));
		printf("FT stage=%d allocations=%d\n", stage, allocations);
		for (failure = 1; failure <= allocations; failure++)
			CHECK(!ft_case(stage, failure, &allocations));
	}
	return 0;
}

static int extremes(const char *mode)
{
	struct simplex_SparseLu factor;
	struct simplex_SparseRow row, rows[2];
	struct simplex_SparseColumnRows column;
	memset(&factor, 0, sizeof(factor));
	memset(&row, 0, sizeof(row));
	memset(rows, 0, sizeof(rows));
	memset(&column, 0, sizeof(column));
	requests = fail_at = 0;
	if (!strcmp(mode, "row")) {
		row.capacity = INT_MAX - 1;
		CHECK(ROW_RESERVE(&factor, &row, INT_MAX));
	} else if (!strcmp(mode, "column")) {
		column.count = column.capacity = INT_MAX / 2 + 1;
		CHECK(COLUMN_APPEND(&factor, &column, 0));
	} else if (!strcmp(mode, "packed")) {
		factor.dimension = 2;
		factor.row = rows;
		rows[0].count = INT_MAX;
		rows[1].count = 1;
		CHECK(sparse_lu_pack_rows(&factor));
	} else if (!strcmp(mode, "degree"))
		sparse_lu_sift_degree_heap(NULL, INT_MAX, INT_MAX - 1, NULL);
	else
		return 1;
	CHECK(!active && !bad_free);
	return 0;
}

int main(int argc, char **argv)
{
	if (argc > 1) {
		if (!strcmp(argv[1], "constructors"))
			return constructors();
		if (!strcmp(argv[1], "ft"))
			return ft_failures();
		return extremes(argv[1]);
	}
#ifndef LU_STORAGE_OLD
	{
		struct sparse_LuStorage layout;
		CHECK(!sparse_lu_storage(-1, (size_t)-1, &layout));
		CHECK(!sparse_lu_storage(
			INT_MAX, (size_t)4294967295UL, &layout));
		CHECK(sparse_lu_storage(8, (size_t)-1, &layout));
		CHECK(layout.integer == 49 * sizeof(int) &&
			layout.numeric == 40 * sizeof(double));
	}
#endif
	return constructors() || reserves() || builds(0) || builds(1) ||
	       ft_failures() || extremes("row") || extremes("column") ||
	       extremes("packed") || extremes("degree");
}
