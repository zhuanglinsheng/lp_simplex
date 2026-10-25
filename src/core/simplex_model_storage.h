/**
 * @file simplex_model_storage.h
 * @brief Checked sizes for model storage; no algorithm or total-RSS policy.
 */

#ifndef LP_SIMPLEX_MODEL_STORAGE_INTERNAL_H
#define LP_SIMPLEX_MODEL_STORAGE_INTERNAL_H
#include <stddef.h>
#include <lp_simplex/model.h>

struct simplex_ModelStorage {
	size_t objective, coefficients, constraints, bounds;
	size_t column_start, row_start, indices, values, row_cursor;
};

static int simplex_model_array_bytes(
	size_t count, size_t width, size_t limit, size_t *bytes)
{
	*bytes = 0;
	if (!width || count > limit / width)
		return 0;
	*bytes = count * width;
	return 1;
}

/* The explicit limit lets tests check a narrower address space without
 * allocating giant arrays. Production uses the host size_t limit. */
static int simplex_model_storage(int rows, int columns, int nonzeros, int dense,
	size_t limit, struct simplex_ModelStorage *storage)
{
	size_t cells = 0;
	if (rows <= 0 || columns <= 0 || nonzeros < 0)
		return 0;
	if (dense) {
		if ((size_t)rows > limit / (size_t)columns)
			return 0;
		cells = (size_t)rows * (size_t)columns;
	}
	return simplex_model_array_bytes((size_t)columns, sizeof(double), limit,
		       &storage->objective) &&
	       simplex_model_array_bytes(
		       cells, sizeof(double), limit, &storage->coefficients) &&
	       simplex_model_array_bytes((size_t)rows,
		       sizeof(struct optm_LinearConstraint), limit,
		       &storage->constraints) &&
	       simplex_model_array_bytes((size_t)columns,
		       sizeof(struct optm_VariableBound), limit,
		       &storage->bounds) &&
	       simplex_model_array_bytes((size_t)columns + 1, sizeof(int),
		       limit, &storage->column_start) &&
	       simplex_model_array_bytes((size_t)rows + 1, sizeof(int), limit,
		       &storage->row_start) &&
	       simplex_model_array_bytes((size_t)nonzeros, sizeof(int), limit,
		       &storage->indices) &&
	       simplex_model_array_bytes((size_t)nonzeros, sizeof(double),
		       limit, &storage->values) &&
	       simplex_model_array_bytes(
		       (size_t)rows, sizeof(int), limit, &storage->row_cursor);
}
#endif
