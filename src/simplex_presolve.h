#ifndef LP_SIMPLEX_PRESOLVE_INTERNAL_H
#define LP_SIMPLEX_PRESOLVE_INTERNAL_H

#include <lp_simplex/model.h>

struct simplex_Presolve {
	const struct lp_Model *original;
	struct lp_Model *reduced;
	int *column_map;
	int *row_map;
	double *eliminated_value;
	unsigned char *eliminated;
	int removed_columns;
	int removed_rows;
	int fixed_columns;
	int empty_columns;
	int singleton_columns;
	int terminal;
	int terminal_status;
};

int simplex_presolve_run(struct simplex_Presolve *presolve,
		const struct lp_Model *model, double tolerance);
void simplex_presolve_postsolve(const struct simplex_Presolve *presolve,
		const double *reduced_x, double *original_x);
void simplex_presolve_destroy(struct simplex_Presolve *presolve);

#endif
