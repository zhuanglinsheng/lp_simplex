/**
 * @file simplex_pan_face_types.h
 * @brief Native-face opportunities and invocation-local diagnostic results.
 */

#ifndef LP_SIMPLEX_PAN_FACE_TYPES_H
#define LP_SIMPLEX_PAN_FACE_TYPES_H

enum simplex_PanFaceOpportunity {
	SIMPLEX_PAN_FACE_UNSPECIFIED = 0,
	SIMPLEX_PAN_FACE_INITIAL,
	SIMPLEX_PAN_FACE_LOW_PROGRESS,
	SIMPLEX_PAN_FACE_AFTER_MULTIPLIER_SEARCH,
	SIMPLEX_PAN_FACE_PDHG_IMPROVEMENT
};

enum simplex_PanFaceAction {
	SIMPLEX_PAN_FACE_NOT_RUN = 0,
	SIMPLEX_PAN_FACE_NO_CHANGE,
	SIMPLEX_PAN_FACE_ADJUSTED,
	SIMPLEX_PAN_FACE_EXACT_SKIP,
	SIMPLEX_PAN_FACE_RESOURCE
};

/* Diagnostics only. An adjustment still needs the native loop's audits;
 * neither this result nor a journal event grants point/KKT permissions. */
struct simplex_PanFaceResult {
	int action, changed, remembered;
};

struct simplex_PanFaceReport {
	int opportunity, steps, low_progress;
	struct simplex_PanFaceResult result;
};

#endif
