/*
 * Dynamic sparse QR R-factor for Pan's deficient active basis.
 *
 * The numerical factor is stored as sparse, sorted rows.  Symbolic fill is
 * discovered incrementally when Givens rotations merge two row patterns.
 */
#ifndef LP_SIMPLEX_PAN_DYNAMIC_SPARSE_QR_H
#define LP_SIMPLEX_PAN_DYNAMIC_SPARSE_QR_H


struct pan_DynamicSparseRow {
	int head;
	int tail;
	int count;
};


struct pan_DynamicSparseNode {
	double value;
	union {
		int row;
		int free_next;
	} owner;
	int column;
	int row_previous;
	int row_next;
	int column_previous;
	int column_next;
};


struct pan_DynamicSparseQR {
	struct pan_DynamicSparseRow *row;
	struct pan_DynamicSparseNode *node;
	int *column_head;
	int *column_tail;
	int *column_count;
	int *tag;
	int *union_column;
	int *scratch_stamp;
	double *rotation_first;
	double *rotation_second;
	int *event_row;
	double *event_cosine;
	double *event_sine;
	/* The linked node pool is mutation-oriented.  Repeated triangular solves
	 * use this generation-cached packed row view for locality. */
	int *packed_start;
	int *packed_column;
	double *packed_value;
	double *packed_diagonal;
	int event_count;
	int count;
	int capacity;
	int node_count;
	int node_capacity;
	int free_node;
	int packed_capacity;
	int packed_valid;
	int scratch_epoch;
	long extensions;
	long downdates;
	long rotations;
	long symbolic_updates;
};


int pan_dynamic_sparse_qr_initialize(struct pan_DynamicSparseQR *factor);

void pan_dynamic_sparse_qr_clear(struct pan_DynamicSparseQR *factor);

void pan_dynamic_sparse_qr_destroy(struct pan_DynamicSparseQR *factor);

/* Append an already orthogonalized QR border.  border has factor->count
 * entries and diagonal is the norm of the orthogonal residual. */
int pan_dynamic_sparse_qr_append(
		struct pan_DynamicSparseQR *factor, int tag,
		const double *border, double diagonal, double rank_tolerance);

/* Replace the numerical factor from a column-compressed upper trapezoid.
 * Tags remain in the caller's fixed column order. */
int pan_dynamic_sparse_qr_load_csc(
		struct pan_DynamicSparseQR *factor, const int *tag, int count,
		const int *column_start, const int *row_index,
		const double *value, double rank_tolerance);

/* Delete a principal column from R by removing the column and chasing the
 * resulting bulge with orthogonal Givens rotations. */
int pan_dynamic_sparse_qr_delete(
		struct pan_DynamicSparseQR *factor, int position,
		double rank_tolerance);

/* Move the last appended column left to position using adjacent column swaps
 * and Givens retriangularization. */
int pan_dynamic_sparse_qr_move_last(
		struct pan_DynamicSparseQR *factor, int position,
		double rank_tolerance);

int pan_dynamic_sparse_qr_solve(
		struct pan_DynamicSparseQR *factor, double *right);

int pan_dynamic_sparse_qr_solve_upper(
		struct pan_DynamicSparseQR *factor, double *right);

int pan_dynamic_sparse_qr_solve_upper_transpose(
		struct pan_DynamicSparseQR *factor, double *right);

/* Expensive structural audit used by tests and diagnostics, never by the
 * numerical hot path. */
int pan_dynamic_sparse_qr_validate(
		const struct pan_DynamicSparseQR *factor);

#endif
