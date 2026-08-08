/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#ifndef LP_SIMPLEX_BASIS_INTERNAL_H
#define LP_SIMPLEX_BASIS_INTERNAL_H

#include "simplex_csc.h"
#include "simplex_sparse_lu.h"

#ifdef LP_SIMPLEX_HAVE_KLU
#include <klu.h>
#endif


struct simplex_Basis {
	int rows;
	int structural_columns;
	const struct simplex_CscMatrix *matrix;
	int *index;
	int core_size;
	int factor_size;
	int compact_active;
	int compact_ever_active;
	int compact_requested;
	int allow_sparse_eta;
	int *core_basis;
	int *core_position;
	int *core_row;
	int *row_to_core;
	int *logical_position;
	int *base_index;
	double *base_work;
	double *core_work;
	double *refine_work;
	double *correction_work;
#ifdef LP_SIMPLEX_HAVE_KLU
	int *base_column_start;
	int *base_row_index;
	double *base_value;
	klu_symbolic *symbolic;
	klu_numeric *numeric;
	klu_common common;
#else
	struct simplex_SparseLu sparse;
#endif
	double *eta_value;
	double *eta;
	double *eta_pivot_value;
	int *eta_start;
	int *eta_index;
	int *eta_pivot;
	unsigned char *eta_sparse;
	int eta_capacity;
	int update_count;
	int update_limit;
	long update_work;
	double minimum_relative_pivot;
	double last_factor_seconds;
	double eta_apply_seconds;
	int compact_ftran_validation_countdown;
	int compact_btran_validation_countdown;
	int profile_enabled;
	double profile_factor_seconds;
	double profile_ftran_seconds;
	double profile_btran_seconds;
	long profile_factor_calls;
	long profile_ftran_calls;
	long profile_btran_calls;
	long profile_compact_calls;
	int profile_compact_min;
	int profile_compact_max;
	long profile_eta_nonzeros;
	long profile_eta_slots;
	long profile_fill_reinversions;
	long profile_stability_reinversions;
	long profile_compact_ftran_validations;
	long profile_compact_btran_validations;
	long profile_compact_ftran_refinements;
	long profile_compact_btran_refinements;
};


int simplex_basis_create(
		struct simplex_Basis *basis,
		const struct simplex_CscMatrix *matrix,
		int structural_columns, int *index);

void simplex_basis_destroy(struct simplex_Basis *basis);

int simplex_basis_factorize(struct simplex_Basis *basis);

int simplex_basis_ftran(const struct simplex_Basis *basis, double *vector);

int simplex_basis_ftran_pair(
		const struct simplex_Basis *basis,
		double *first, double *second);

int simplex_basis_btran(const struct simplex_Basis *basis, double *vector);

int simplex_basis_update(
		struct simplex_Basis *basis,
		int leaving_position, const double *direction);

#endif
