#include "simplex_presolve.h"
#include "utils.h"
#include <lp_simplex/status.h>


static int presolve_finite_lower(const struct optm_VariableBound *bound)
{
	return bound->b_type == optm_BOUND_T_LO ||
		bound->b_type == optm_BOUND_T_BS;
}


static int presolve_finite_upper(const struct optm_VariableBound *bound)
{
	return bound->b_type == optm_BOUND_T_UP ||
		bound->b_type == optm_BOUND_T_BS;
}


static int presolve_empty_row_feasible(
		const struct optm_LinearConstraint *constraint,
		const double rhs, const double tolerance)
{
	if (constraint->type == optm_CONS_T_EQ)
		return __lp_simplex_ABS__(rhs) <= tolerance;
	if (constraint->type == optm_CONS_T_GE)
		return rhs <= tolerance;
	return rhs >= -tolerance;
}


static int presolve_choose_empty_column(
		const struct lp_Model *model, const int column, double *value)
{
	const struct optm_VariableBound *bound = model->bounds + column;
	double cost = model->objective[column];
	if (cost > 0.) {
		if (!presolve_finite_lower(bound))
			return lp_simplex_Unboundedness;
		*value = bound->lb;
	} else if (cost < 0.) {
		if (!presolve_finite_upper(bound))
			return lp_simplex_Unboundedness;
		*value = bound->ub;
	} else if (presolve_finite_lower(bound) && bound->lb > 0.) {
		*value = bound->lb;
	} else if (presolve_finite_upper(bound) && bound->ub < 0.) {
		*value = bound->ub;
	} else {
		*value = 0.;
	}
	return lp_simplex_Success;
}


static int presolve_column_empty(
		const struct lp_Model *model, const int column)
{
	int i;
	if (model->column_start != NULL)
		return model->column_start[column] == model->column_start[column + 1];
	for (i = 0; i < model->m; i++)
		if (model->constraints[i].coef[column] != 0.)
			return 0;
	return 1;
}


static void presolve_remove_column_from_rhs(
		const struct lp_Model *model, const int column,
		const double value, double *rhs)
{
	int i, k;
	if (model->column_start != NULL) {
		for (k = model->column_start[column];
		     k < model->column_start[column + 1]; k++)
			rhs[model->row_index[k]] -= model->value[k] * value;
		return;
	}
	for (i = 0; i < model->m; i++)
		rhs[i] -= model->constraints[i].coef[column] * value;
}


static struct lp_Model *presolve_build_sparse_model(
		const struct lp_Model *model, const int kept_rows,
		const int kept_columns, const int *column_map, const int *column_kept,
		const int *row_map, const int *row_kept, const double *adjusted_rhs)
{
	struct lp_Model *reduced;
	int i, j, k, nonzeros = 0, next = 0;
	for (i = 0; i < kept_rows; i++) {
		int original_row = row_map[i];
		if (model->row_start != NULL) {
			for (k = model->row_start[original_row];
			     k < model->row_start[original_row + 1]; k++)
				if (column_kept[model->column_index[k]] >= 0)
					nonzeros++;
		} else {
			for (j = 0; j < model->n; j++)
				if (column_kept[j] >= 0 &&
				    model->constraints[original_row].coef[j] != 0.)
					nonzeros++;
		}
	}
	reduced = (struct lp_Model *)lp_simplex_malloc(sizeof(*reduced));
	if (reduced == NULL)
		return NULL;
	lp_simplex_memset(reduced, 0, sizeof(*reduced));
	reduced->m = kept_rows;
	reduced->n = kept_columns;
	reduced->nnz = nonzeros;
	reduced->objective = (double *)lp_simplex_malloc(
		(size_t)kept_columns * sizeof(double));
	reduced->bounds = (struct optm_VariableBound *)lp_simplex_malloc(
		(size_t)kept_columns * sizeof(*reduced->bounds));
	reduced->constraints = (struct optm_LinearConstraint *)lp_simplex_malloc(
		(size_t)kept_rows * sizeof(*reduced->constraints));
	reduced->column_start = (int *)lp_simplex_malloc(
		(size_t)(kept_columns + 1) * sizeof(int));
	reduced->row_start = (int *)lp_simplex_malloc(
		(size_t)(kept_rows + 1) * sizeof(int));
	reduced->row_index = nonzeros > 0 ? (int *)lp_simplex_malloc(
		(size_t)nonzeros * sizeof(int)) : NULL;
	reduced->value = nonzeros > 0 ? (double *)lp_simplex_malloc(
		(size_t)nonzeros * sizeof(double)) : NULL;
	reduced->column_index = nonzeros > 0 ? (int *)lp_simplex_malloc(
		(size_t)nonzeros * sizeof(int)) : NULL;
	reduced->row_value = nonzeros > 0 ? (double *)lp_simplex_malloc(
		(size_t)nonzeros * sizeof(double)) : NULL;
	if (reduced->objective == NULL || reduced->bounds == NULL ||
	    reduced->constraints == NULL || reduced->column_start == NULL ||
	    reduced->row_start == NULL ||
	    (nonzeros > 0 && (reduced->row_index == NULL ||
	     reduced->value == NULL || reduced->column_index == NULL ||
	     reduced->row_value == NULL))) {
		lp_model_free(reduced);
		return NULL;
	}
	for (j = 0; j < kept_columns; j++) {
		int original_column = column_map[j];
		reduced->objective[j] = model->objective[original_column];
		reduced->bounds[j] = model->bounds[original_column];
		reduced->column_start[j] = next;
		if (model->column_start != NULL) {
			for (k = model->column_start[original_column];
			     k < model->column_start[original_column + 1]; k++) {
				int row = row_kept[model->row_index[k]];
				if (row >= 0) {
					reduced->row_index[next] = row;
					reduced->value[next++] = model->value[k];
				}
			}
		} else {
			for (i = 0; i < kept_rows; i++) {
				double value = model->constraints[row_map[i]].coef[
					original_column];
				if (value != 0.) {
					reduced->row_index[next] = i;
					reduced->value[next++] = value;
				}
			}
		}
	}
	reduced->column_start[kept_columns] = next;
	next = 0;
	for (i = 0; i < kept_rows; i++) {
		int original_row = row_map[i];
		reduced->constraints[i] = model->constraints[original_row];
		reduced->constraints[i].coef = NULL;
		reduced->constraints[i].rhs = adjusted_rhs[original_row];
		reduced->row_start[i] = next;
		if (model->row_start != NULL) {
			for (k = model->row_start[original_row];
			     k < model->row_start[original_row + 1]; k++) {
				int column = column_kept[model->column_index[k]];
				if (column >= 0) {
					reduced->column_index[next] = column;
					reduced->row_value[next++] = model->row_value[k];
				}
			}
		} else {
			for (j = 0; j < kept_columns; j++) {
				double value = model->constraints[original_row].coef[
					column_map[j]];
				if (value != 0.) {
					reduced->column_index[next] = j;
					reduced->row_value[next++] = value;
				}
			}
		}
	}
	reduced->row_start[kept_rows] = next;
	return reduced;
}


int simplex_presolve_run(
		struct simplex_Presolve *presolve, const struct lp_Model *model,
		const double tolerance)
{
	int i, j, k, kept_rows = 0, kept_columns = 0;
	int changed;
	int *column_kept = NULL;
	int *row_kept = NULL;
	double *adjusted_rhs = NULL;
	if (presolve == NULL || model == NULL || model->bounds == NULL)
		return lp_simplex_EXIT_FAILURE;
	lp_simplex_memset(presolve, 0, sizeof(*presolve));
	presolve->original = model;
	presolve->terminal_status = lp_simplex_CondUnsatisfied;
	presolve->column_map = (int *)lp_simplex_malloc(
		(size_t)model->n * sizeof(int));
	presolve->row_map = (int *)lp_simplex_malloc(
		(size_t)model->m * sizeof(int));
	presolve->eliminated_value = (double *)lp_simplex_malloc(
		(size_t)model->n * sizeof(double));
	presolve->eliminated = (unsigned char *)lp_simplex_malloc(
		(size_t)model->n * sizeof(unsigned char));
	column_kept = (int *)lp_simplex_malloc((size_t)model->n * sizeof(int));
	row_kept = (int *)lp_simplex_malloc((size_t)model->m * sizeof(int));
	adjusted_rhs = (double *)lp_simplex_malloc(
		(size_t)model->m * sizeof(double));
	if (presolve->column_map == NULL || presolve->row_map == NULL ||
	    presolve->eliminated_value == NULL || presolve->eliminated == NULL ||
	    column_kept == NULL || row_kept == NULL || adjusted_rhs == NULL)
		goto failure;
	lp_simplex_memset(presolve->eliminated, 0,
		(size_t)model->n * sizeof(unsigned char));
	for (j = 0; j < model->n; j++)
		column_kept[j] = -1;
	for (i = 0; i < model->m; i++)
		row_kept[i] = -1;
	for (i = 0; i < model->m; i++)
		adjusted_rhs[i] = model->constraints[i].rhs;
	for (j = 0; j < model->n; j++) {
		const struct optm_VariableBound *bound = model->bounds + j;
		int fixed = bound->b_type == optm_BOUND_T_BS &&
			bound->lb == bound->ub;
		double value = fixed ? bound->lb : 0.;
		if (!fixed && presolve_column_empty(model, j)) {
			int status = presolve_choose_empty_column(model, j, &value);
			if (status != lp_simplex_Success) {
				presolve->terminal = 1;
				presolve->terminal_status = status;
				goto finish;
			}
			fixed = 1;
			presolve->empty_columns++;
		}
		if (!fixed) {
			column_kept[j] = 1;
			continue;
		}
		column_kept[j] = -1;
		presolve->eliminated[j] = 1;
		presolve->eliminated_value[j] = value;
		presolve->removed_columns++;
		if (!presolve_column_empty(model, j))
			presolve->fixed_columns++;
		presolve_remove_column_from_rhs(model, j, value, adjusted_rhs);
	}
	/* Equality singleton elimination is iterated because fixing one column can
	 * expose another singleton.  The value record is sufficient for exact
	 * primal postsolve; the emptied equality row is removed below. */
	do {
		changed = 0;
		for (i = 0; i < model->m; i++) {
			int entries = 0, singleton = -1;
			double coefficient = 0., value;
			const struct optm_VariableBound *bound;
			if (model->constraints[i].type != optm_CONS_T_EQ)
				continue;
			if (model->row_start != NULL) {
				for (k = model->row_start[i];
				     k < model->row_start[i + 1]; k++) {
					j = model->column_index[k];
					if (column_kept[j] < 0)
						continue;
					entries++;
					singleton = j;
					coefficient = model->row_value[k];
					if (entries > 1)
						break;
				}
			} else {
				for (j = 0; j < model->n; j++) {
					if (column_kept[j] < 0 ||
					    model->constraints[i].coef[j] == 0.)
						continue;
					entries++;
					singleton = j;
					coefficient = model->constraints[i].coef[j];
					if (entries > 1)
						break;
				}
			}
			if (entries != 1 || coefficient == 0.)
				continue;
			value = adjusted_rhs[i] / coefficient;
			bound = model->bounds + singleton;
			if ((presolve_finite_lower(bound) &&
			     value < bound->lb - tolerance) ||
			    (presolve_finite_upper(bound) &&
			     value > bound->ub + tolerance)) {
				presolve->terminal = 1;
				presolve->terminal_status = lp_simplex_Infeasibility;
				goto finish;
			}
			if (presolve_finite_lower(bound) && value < bound->lb)
				value = bound->lb;
			if (presolve_finite_upper(bound) && value > bound->ub)
				value = bound->ub;
			column_kept[singleton] = -1;
			presolve->eliminated[singleton] = 1;
			presolve->eliminated_value[singleton] = value;
			presolve->removed_columns++;
			presolve->singleton_columns++;
			presolve_remove_column_from_rhs(
				model, singleton, value, adjusted_rhs);
			changed = 1;
		}
	} while (changed);
	for (j = 0; j < model->n; j++)
		if (column_kept[j] >= 0) {
			column_kept[j] = kept_columns;
			presolve->column_map[kept_columns++] = j;
		}
	for (i = 0; i < model->m; i++) {
		int entries = 0;
		if (model->row_start != NULL) {
			for (k = model->row_start[i]; k < model->row_start[i + 1]; k++)
				if (column_kept[model->column_index[k]] >= 0)
					entries++;
		} else {
			for (j = 0; j < model->n; j++)
				if (column_kept[j] >= 0 &&
				    model->constraints[i].coef[j] != 0.)
					entries++;
		}
		if (entries == 0) {
			if (!presolve_empty_row_feasible(model->constraints + i,
					adjusted_rhs[i], tolerance)) {
				presolve->terminal = 1;
				presolve->terminal_status = lp_simplex_Infeasibility;
				goto finish;
			}
			row_kept[i] = -1;
			presolve->removed_rows++;
		} else {
			row_kept[i] = kept_rows;
			presolve->row_map[kept_rows++] = i;
		}
	}
	if (kept_columns == 0) {
		presolve->terminal = 1;
		presolve->terminal_status = lp_simplex_Success;
		goto finish;
	}
	if (presolve->removed_columns == 0 && presolve->removed_rows == 0)
		goto finish;
	if (kept_rows == 0) {
		presolve->terminal = 1;
		presolve->terminal_status = lp_simplex_Success;
		goto finish;
	}
	presolve->reduced = presolve_build_sparse_model(model, kept_rows,
		kept_columns, presolve->column_map, column_kept,
		presolve->row_map, row_kept, adjusted_rhs);
	if (presolve->reduced == NULL)
		goto failure;
finish:
	lp_simplex_free(column_kept);
	lp_simplex_free(row_kept);
	lp_simplex_free(adjusted_rhs);
	return lp_simplex_EXIT_SUCCESS;
failure:
	lp_simplex_free(column_kept);
	lp_simplex_free(row_kept);
	lp_simplex_free(adjusted_rhs);
	simplex_presolve_destroy(presolve);
	return lp_simplex_EXIT_FAILURE;
}


void simplex_presolve_postsolve(
		const struct simplex_Presolve *presolve,
		const double *reduced_x, double *original_x)
{
	int j;
	for (j = 0; j < presolve->original->n; j++)
		original_x[j] = presolve->eliminated[j]
			? presolve->eliminated_value[j] : 0.;
	if (reduced_x != NULL && presolve->reduced != NULL)
		for (j = 0; j < presolve->reduced->n; j++)
			original_x[presolve->column_map[j]] = reduced_x[j];
}


void simplex_presolve_destroy(struct simplex_Presolve *presolve)
{
	if (presolve == NULL)
		return;
	lp_model_free(presolve->reduced);
	lp_simplex_free(presolve->column_map);
	lp_simplex_free(presolve->row_map);
	lp_simplex_free(presolve->eliminated_value);
	lp_simplex_free(presolve->eliminated);
	presolve->reduced = NULL;
	presolve->column_map = NULL;
	presolve->row_map = NULL;
	presolve->eliminated_value = NULL;
	presolve->eliminated = NULL;
}
