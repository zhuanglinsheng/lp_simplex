/**
 * @file resources.h
 * @brief Policy defaults; units and conservative failure contracts:
 *
 * docs/numerical-policy.md.
 */

#ifndef LP_SIMPLEX_POLICY_RESOURCES_H
#define LP_SIMPLEX_POLICY_RESOURCES_H

/* Optional dense multiplier recovery, including LAPACK scratch. */
#define SIMPLEX_PAN_MINIMUM_NORM_QR_BYTES (32UL * 1024UL * 1024UL)
/* Exact failed basis-start snapshots, shared by one solve. */
#define SIMPLEX_RETRY_BASIS_BYTES (4UL * 1024UL * 1024UL)
/* One exact optional no-change face request, local to a native BDA call. */
#define SIMPLEX_PAN_FACE_RETRY_BYTES (4UL * 1024UL * 1024UL)

#define SIMPLEX_AUTO_DENSE_BYTES (32UL * 1024UL * 1024UL)
#define SIMPLEX_AUTO_CONVERSION_BYTES (64UL * 1024UL * 1024UL)
#define COVER_WORKSPACE_BYTES (256.L * 1024.L * 1024.L)
/* Peak owned optional dominance arrays, including temporary index builders. */
#define DOMINANCE_BYTES (32UL * 1024UL * 1024UL)
#define PRESOLVE_PROPAGATION_LIMIT 1e20
#define CERTIFICATE_MAX_RANK 128
#define CERTIFICATE_MAX_BITS 8192
#define SIMPLEX_SCALING_PASSES 4
#define SIMPLEX_SCALING_MIN_EXPONENT (-48)
#define SIMPLEX_SCALING_MAX_EXPONENT 48

#endif
