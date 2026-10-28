/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_klu_sparse.h
 * @brief Optional support-driven solves of immutable KLU factors.
 */

#ifndef LP_SIMPLEX_KLU_SPARSE_H
#define LP_SIMPLEX_KLU_SPARSE_H

struct simplex_KluSparseImpl;
struct simplex_SparseVector;

struct simplex_KluSparse {
	struct simplex_KluSparseImpl *impl;
};

#ifdef LP_SIMPLEX_HAVE_KLU
#include <klu.h>

/* Borrows permutations/scaling from these factors until destroy. Refused
 * optional storage leaves an empty cache; the original KLU factor stays live. */
int simplex_klu_sparse_build(struct simplex_KluSparse *cache,
	klu_symbolic *symbolic, klu_numeric *numeric, klu_common *common);
void simplex_klu_sparse_destroy(struct simplex_KluSparse *cache);

/* 0 success, 1 declined (input unchanged), -1 nonfinite arithmetic (input
 * unchanged). Dense input/support falls back to the original KLU solve. */
int simplex_klu_sparse_solve(
	struct simplex_KluSparse *cache, double *vector, int transpose);
int simplex_klu_sparse_solve_packed(struct simplex_KluSparse *cache,
	struct simplex_SparseVector *vector, int transpose);
#endif

#endif
