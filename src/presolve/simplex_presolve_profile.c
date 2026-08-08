/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#include "simplex_presolve.h"

#include <lp_simplex/status.h>

#include <stdio.h>


void simplex_presolve_print_profile(
		const struct simplex_Presolve *presolve, const double seconds)
{
	const struct simplex_PresolveStats *stats = &presolve->stats;
	int rows = presolve->reduced != NULL ? presolve->reduced->rows : 0;
	int columns = presolve->reduced != NULL ? presolve->reduced->columns : 0;
	fprintf(stderr, "presolve: time=%.6f seconds\n", seconds);
	if (stats->removed_rows == 0 && stats->removed_columns == 0 &&
	    stats->tightened_bounds == 0)
		return;
	fprintf(stderr, "presolve: removed rows=%d columns=%d "
		"[fixed=%d empty=%d singleton-col=%d singleton-ineq-col=%d "
		"singleton-row=%d redundant-row=%d duplicate-row=%d "
		"forcing-row=%d forced-column=%d doubleton=%d "
		"singleton-column=%d implied-free-column=%d "
		"singleton-projection=%d tightened=%d epochs=%d "
		"singleton-column-candidates=%d free-singleton=%d "
		"equality-singleton=%d exact-propagation=%d] remaining=%d/%d\n",
		stats->removed_rows, stats->removed_columns,
		stats->fixed_columns, stats->empty_columns,
		stats->singleton_columns, stats->singleton_inequality_columns,
		stats->singleton_rows, stats->redundant_rows,
		stats->duplicate_rows, stats->forcing_rows,
		stats->forced_columns, stats->doubleton_rows,
		stats->singleton_column_rows, stats->implied_free_columns,
		stats->singleton_projection_columns, stats->tightened_bounds,
		stats->queue_epochs, stats->singleton_column_candidates,
		stats->free_singleton_columns, stats->equality_singleton_columns,
		stats->exact_bound_propagation, rows, columns);
}
