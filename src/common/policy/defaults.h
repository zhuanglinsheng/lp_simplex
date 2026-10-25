/**
 * @file defaults.h
 * @brief Policy defaults; units and conservative failure contracts:
 *
 * docs/numerical-policy.md.
 */

#ifndef LP_SIMPLEX_POLICY_DEFAULTS_H
#define LP_SIMPLEX_POLICY_DEFAULTS_H

#define SIMPLEX_DEFAULT_ITERATION_LIMIT 100000
#define SIMPLEX_DEFAULT_PRIMAL_TOLERANCE 1e-7
#define SIMPLEX_DEFAULT_DUAL_TOLERANCE 1e-7
#define SIMPLEX_DEFAULT_PIVOT_TOLERANCE 1e-9

/* Deterministic observation polls, independent of callback refresh time. */
#define SIMPLEX_PROGRESS_OBSERVATION_PERIOD 128
/* CPU seconds between progress callback reports; not a solve timeout. */
#define SIMPLEX_DEFAULT_PROGRESS_INTERVAL 5.

#endif
