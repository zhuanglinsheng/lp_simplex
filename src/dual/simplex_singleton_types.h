/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**
 * @file simplex_singleton_types.h
 * @brief Singleton transformation maps and derived-model state.
 */

#ifndef LP_SIMPLEX_SINGLETON_TYPES_INTERNAL_H
#define LP_SIMPLEX_SINGLETON_TYPES_INTERNAL_H

/* Classification belongs to the input, not automatically to its dual LP. */
enum simplex_SingletonMapping {
	SIMPLEX_SINGLETON_NO_MAPPING = 0,
	SIMPLEX_SINGLETON_OPTIMUM,
	SIMPLEX_SINGLETON_UNBOUNDED,
	SIMPLEX_SINGLETON_REJECTED
};

struct simplex_SingletonReport {
	int available, input_model, auxiliary_model;
	int auxiliary_available, auxiliary_status, mapping;
	int structure_verified, input_point_checked;
};
#endif
