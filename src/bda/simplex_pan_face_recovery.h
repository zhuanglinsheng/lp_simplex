/**
 * @file simplex_pan_face_recovery.h
 * @brief Scoped native-face recovery and exact no-change request memory.
 */

#ifndef LP_SIMPLEX_PAN_FACE_RECOVERY_H
#define LP_SIMPLEX_PAN_FACE_RECOVERY_H

#include "simplex_pan_bounded.h"
#include "simplex_problem_identity.h"
#include "simplex_pan_face_types.h"
#include <stddef.h>

/* Fresh zero initialization; destroy before releasing/rebuilding the bounded
 * view. Its matrix, objective, norms, RHS and bounds are immutable while this
 * owner lives. In-place edits require clearing the owner before reuse.
 * No factors or point permissions are retained; counters are output only. */
struct simplex_PanFaceRecovery {
	const struct simplex_PanBounded *model;
	struct simplex_PanStandard representation;
	struct simplex_ProblemIdentity input;
	const double *lower, *upper;
	double *storage;
	size_t columns_bytes, rows_bytes;
	double primal_tolerance, dual_tolerance;
	int known;
	long attempts, skipped, unremembered;
};

/* >0: audited coordinates changed; 0: unchanged, including a skipped exact
 * request; -1: optional sorting storage refused. Failure/skip preserves x/y;
 * saved/residual/work are scratch, never outputs. Iterations are not charged.
 * A NULL owner performs the audited action without remembering a request. */
int simplex_pan_face_recover(struct simplex_PanFaceRecovery *owner,
	const struct simplex_PanBounded *model, double *point,
	const double *multiplier, const struct lp_simplex_Options *options,
	double *saved, double *residual, double *work);
void simplex_pan_face_recovery_destroy(struct simplex_PanFaceRecovery *owner);
int simplex_pan_face_recover_reported(struct simplex_PanFaceRecovery *owner,
	const struct simplex_PanBounded *b, double *point,
	const double *multiplier, const struct lp_simplex_Options *options,
	double *saved, double *residual, double *work,
	struct simplex_PanFaceResult *result);

#endif
