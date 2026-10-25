# BDA phase contracts

[简体中文](./pan-services.md) | English | [Documentation](./README_en.md) | [AUTO](./auto-policy_en.md)

BDA owns its invocation state. NNLS, PDHG, cones, restoration and dual refinement are internal roles. Standard/native bounds are representations; explicit/implicit are oracle implementations; QR is a numerical backend.

The [service adapter](../src/bda/simplex_pan_service.c) centralizes calls and [types](../src/bda/simplex_pan_service_types.h) define roles/outputs. BDA loops decide when to call; subprocesses own local iteration and numerical checks. AUTO owns cross-algorithm order.

## Services and output permissions

| Service | Input and applicability | Permitted output | Grant and exit |
| --- | --- | --- | --- |
| `nnls-initialize` | Equality representation, live basis and tolerances; no initial point | Nonnegative equality-feasible point in this representation | Parent remainder; handoff, budget, numerical/resource failure or verified current-model infeasibility |
| `dual-feasibility` | Stable input LP and resolved controls; needs an initial point | Input-space feasible point | Remaining work; convert Dual local work to the parent cumulative count |
| `pdhg-advance` | Native point/multipliers, sparse model with nonzero rows | Audited feasible improvement or native KKT | Minimum of request/remainder; budget, stagnation, handoff, certification, numerical/resource exit |
| `cone-direction` | Derived cone, active basis and free mask | Direction proposal, stationary residual or accepted direction | Remaining work; successful roles hand off, preserving actual oracle/failure reason |
| `feasibility-restore` | Standard/native point and objective ceiling | Feasible point in its representation, no multipliers/basis/LP classification | Remaining work; handoff, not-applicable, budget, numerical/resource exit |
| `dual-refine` | Representation, point and multipliers | Checked multipliers | Minimum of `dual_polish_iterations` and remaining work |
| `direction-predict` | Representation, point/multipliers, free mask and optional acceptance function | Checked direction or multipliers | Explicit request or `first_order_iterations`, bounded by remaining work |
| `native-face-recovery` | Native point/multipliers and local request owner; startup, low progress or after multiplier search | Internal face adjustment, no reusable point/KKT authority | Existing zero-iteration charge; adjustment, no change, exact skip or sorting resource refusal |

Approximate points, NNLS feasibility, auxiliary stationary residuals and checked directions do not independently certify input optimality. Ordinary NNLS infeasibility requires an independently verified separator for `Az=b, z>=0`. No entering column, overflow, uncertain signs or certificate allocation failure is no proof. Mixed NNLS stationarity belongs only to the direction role.

Zero grants permit no extra iteration. Already unnecessary standard restoration may return zero-work not-applicable without output; zero-row PDHG is inapplicable. A zero-iteration face check still costs CPU.

## Identity and grants

Each BDA invocation creates a local `PanServiceContext` fixing complete input identity and root scope. Standard/native representations, cones and centered auxiliary models bind or explicitly derive their own identities. Equal dimensions do not imply the same model.

`begin` snapshots matrix, norms, maps, scaling/shifts and objective binding, granting `min(requested, limit-before)`. Native adapters also bind the input object and bound arrays. `finish` checks source/representation identity, basis ownership, cumulative work and output role. Borrowed sources/maps stay alive and immutable; in-place edits must mark revision, and metadata checks do not scan coefficient contents.

Malformed work charges the grant end while preserving the raw report. Identity/permission violations revoke output and set a sticky invocation failure checked again before final commit. Ordinary budget/resource/numerical failure is not automatically a contract violation. Optional snapshot refusal and required workspace refusal differ; see [resource policy](./numerical-policy_en.md).

`Outcome` describes one call, owns no point/multiplier/basis arrays and is no cross-model checkpoint. Dual feasibility uses stable input identity, not the address of a temporary read-only shell, for failed-request comparison.

## Runtime state and history

Contexts/frames enforce protection, Outcome describes output, and StageReport/events are diagnostic. Disabled/full/reset history does not change mathematical eligibility. References spanning callbacks bind context and actor token; stale provenance becomes zero.

Restoration inside PDHG is parented by that PDHG invocation. Both share one cumulative count; restoration is included in the PDHG summary and cannot be summed again. Dual raw local work and adapted cumulative work remain separate.

Triggers cite actual observations: standard trial residual uses `primal-margin`; PDHG stagnation uses `pdhg-stagnation`; face opportunities distinguish startup before the first cone step, low progress and the actual multiplier-search producer. Trigger, work position, metric epoch, point source and numerical output do not come from the latest event. Observations charge no extra work and carry no mathematical result.

## Restoration and repeat admission

| Mechanism | Consumed state and limit |
| --- | --- |
| NNLS | Residual, live basis and blocked set advance together. Factor rejection blocks the column; support changes may reopen it. Successful/block insertions consume work; exhausted correlations are checked under the model's role. |
| Explicit/implicit cone | Support removal, blocked columns, implicit repeated-pattern checks and charged direction steps advance state. Reopening frozen faces, switching first-order/NNLS or tightening tolerance changes the request; inner retries are finite. |
| First-order/CG/Newton | Updates, monotone restarts and CG substeps share the end position. Finite line search is numerical protection, not AUTO reentry. Directions/multipliers remain audited. |
| Standard restoration | Uses a new trial and saved-point cost ceiling. Failure rolls back and switches method, tightens tolerance or exits. Acceptance rechecks residual, cost and nonzero point change. |
| Native compression/certificate repair | Failed compression rolls back, rebuilds a seed and proceeds to the ordinary cone; failed multiplier projection rolls back. No reusable checkpoints are exported. |
| PDHG restoration | Requires a full stagnant window, remaining infeasibility, objective below ceiling and remaining work. Each attempt requires physical residual strictly below half the previous attempt's; success/failure both update the baseline, which ordinary/weight restarts retain. |
| Native face recovery | Only actual no-change requests are remembered. Compare model/representation, point, multipliers and tolerances exactly. Changed requests execute again; missing memory remains unknown. |

Only the latest face request is retained, with a 4 MiB snapshot cap. Sorting refusal is a resource exit; snapshot refusal may preserve an already completed no-change check without remembering it. Destroy the owner before its box representation and clear it before in-place edits. See [numerical policy](./numerical-policy_en.md) for constants/error scales.

## Numerical entries and tests

Full solving uses `simplex_pan_solve_configured`; the native kernel uses `simplex_pan_bounded_solve_certified_configured`. Services call the sole `*_configured_reported` entries for PDHG, restoration and first-order routines. Controls come from the root snapshot or AUTO adapter; private entries have no ABI stability promise.

`test_pan_services` checks roles, nested grants, malformed work, identity and history lifetime. `test_pan_recovery`, `test_pdhg_restore_trigger`, `test_pan_face_recovery` and `test_pan_face_service` exercise actual restoration and failed-state consumption. Mathematical derivations are in the [algorithm collection](./algorithm/README_en.md).
