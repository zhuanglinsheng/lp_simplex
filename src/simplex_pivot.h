#ifndef LP_SIMPLEX_PIVOT_H
#define LP_SIMPLEX_PIVOT_H

void simplex_apply_pivot(double *table, int ld, int m, int n, int leaving_row,
			 int entering_column, int normalize, int eliminate_rows,
			 int eliminate_objective);

int simplex_run_pivots(int *iteration, double *table, int ld, int *basis,
		       int m, int n, int real_columns, const char *criteria,
		       int iteration_limit);

#endif
