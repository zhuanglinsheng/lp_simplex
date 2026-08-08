/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#ifndef LP_SIMPLEX_BASIS_INTERNAL_H
#define LP_SIMPLEX_BASIS_INTERNAL_H

#include "simplex_csc.h"


struct simplex_BasisImpl;

struct simplex_Basis {
	struct simplex_BasisImpl *impl;
};


struct simplex_BasisProfile {
	int factor_size;
	double factor_seconds;
	double ftran_seconds;
	double btran_seconds;
	long factor_calls;
	long ftran_calls;
	long btran_calls;
	long compact_calls;
	int compact_min;
	int compact_max;
	long eta_nonzeros;
	long eta_slots;
	long compact_ftran_validations;
	long compact_btran_validations;
	long compact_ftran_refinements;
	long compact_btran_refinements;
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

int simplex_basis_update_packed(
		struct simplex_Basis *basis,
		int leaving_position, const double *direction,
		const int *index, int nonzeros);

void simplex_basis_set_sparse_eta(
		struct simplex_Basis *basis, int allowed);

int simplex_basis_update_count(const struct simplex_Basis *basis);

int simplex_basis_compact_active(const struct simplex_Basis *basis);

int simplex_basis_compact_ever_active(const struct simplex_Basis *basis);

int simplex_basis_compact_requested(const struct simplex_Basis *basis);

void simplex_basis_request_compact(struct simplex_Basis *basis);

void simplex_basis_get_profile(
		const struct simplex_Basis *basis,
		struct simplex_BasisProfile *profile);

#endif
