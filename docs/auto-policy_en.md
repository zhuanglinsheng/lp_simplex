# AUTO scheduling architecture

[简体中文](./auto-policy.md) | English | [Documentation](./README_en.md)

AUTO combines Dual, BDA and Tableau within one cumulative work budget, retains audited candidates and schedules recovery. It is a scheduler. Model representations, direction implementations and LU/QR backends are separate layers.

## Source ownership and dependency direction

| Owner | Responsibility |
| --- | --- |
| [simplex_auto.c](../src/core/simplex_auto.c) | Sole owner of AUTO state advancement, calls, charging, audits, candidate commits and final assembly |
| [simplex_auto_policy.c](../src/core/simplex_auto_policy.c) | Read-only proposals for triggers, routes, grants and priority cursor |
| [simplex_auto_routes.c](../src/core/simplex_auto_routes.c) | Capability/model/resource admission and solver/initializer adaptation |
| `src/runtime` | Canonical models, identity, configuration snapshots, progress bridges, execution history and failed-request ledger |
| `src/dual`, `src/bda`, `src/tableau` | Independent entries, per-call state, iteration, numerical protection and certification |
| `src/presolve` | Reductions, recovery journal and read-only algebraic mapping |
| `src/basis`, `src/matrix`, `src/common` | Factorizations, sparse storage, shared arithmetic and internal policies |

Algorithms and runtime do not include private core headers; only the controller includes `simplex_auto_state.h`. Runtime may use algorithm-owned diagnostic data, configuration resolution and diagnostic names, plus read-only presolve mapping, without entering complete solving. The Dual degeneracy monitor and the BDA solver are distinct mechanisms.

`python3 tools/check_architecture.py` checks these include boundaries, not runtime identity or mathematical validity.

## Entry and call flow

The private entry is `simplex_auto_solve`. Bindings describe kernels, routes and cover capability; the initializer and output history are separate arguments. A NULL binding object selects production defaults. Missing callbacks in an explicit object remain absent. Missing required Dual causes capability refusal before candidate allocation or numerical calls; optional missing capabilities receive no reserved work.

1. Check callable capabilities, model conditions and workspaces; optionally try cover.
2. Grant the initial Dual invocation, allowing live continuation or handoff based on progress.
3. If Dual yields without a feasible point, optionally invoke the original-model NNLS initializer.
4. Advance through the recovery order below, trying each configuration at most once per AUTO call.
5. Give remaining work to final Dual, then assemble the exit and returned point.

| Recovery order | Configuration and purpose | Additional prerequisite |
| --- | --- | --- |
| Singleton | Dualization of a specialized singleton structure | Matching bounds, row types and positive/negative column pairs; an accepted initializer point may supersede it |
| Native implicit BDA | Advance or optimize a feasible point in native bounds | Representation/backend resources available |
| Native explicit BDA | Replace the native implicit oracle | Corresponding implicit invocation explicitly reports oracle failure |
| Standard basis-direction BDA | Use standard form and basis directions | Conversion and QR resources available |
| Standard implicit BDA | Replace direction method | Standard form applicable; may inherit an audited point |
| Standard explicit BDA | Replace the standard implicit oracle | Corresponding implicit invocation explicitly reports oracle failure |
| Tableau | Switch to two-phase simplex | Dense conversion resources available |

These are composite execution configurations, not seven public algorithms. Generic precision errors, exhausted work and factorization failure alone cannot authorize explicit-oracle replacement. See the [capability inventory](./features_en.md).

## Work grants and attempt completion

Work counts algorithmic iterations rather than CPU seconds. For remaining work R, k reserved opportunities and observation quantum Q=1024, the recovery reserve is `min(R/2, k*Q)`; the implementation checks ranges before multiplication.

- Cover receives R/2 only when width and work conditions permit it.
- Initial Dual reserves observation opportunities for available recovery routes and final Dual; without recovery routes it receives all remaining work.
- Initialization receives at most 8Q after reservation.
- Native hybrid reserves later recovery opportunities; other recoveries divide available opportunities, and final Dual receives the remainder.

[Policy source](../src/core/simplex_auto_policy.c) and [numerical policy](./numerical-policy_en.md) define exact formulas. Kernels report local work and the controller charges once. Nested service summaries include child work and must not be added to their parent again.

Every attempt follows one completion protocol: call → charge → independent audit → terminal check → record → commit candidate/evidence. Raw output, charged work and audit verdict remain separate. Negative or excessive work charges the full grant and revokes output permission. Changed model identity rejects handoff before coordinate reads.

## Candidates, checkpoints and termination

| Output | Admission | Permission |
| --- | --- | --- |
| Feasible point | Matching input identity, finite coordinates, independently recomputed feasibility/objective | Retain, compare and seed later stages; not automatically optimal |
| NNLS initializer point | Both original/reduced audits pass, with valid work | Point only; no basis, factor, dual or lower-bound export |
| Dual basis hint | Matching input/context/model/revision/storage and valid indices/statuses | Re-factor, recompute and classify on import; no live factor preservation |
| Row multiplier or lower bound | Matching model, sign, coordinates and valid audit | Only its declared authority; auxiliary diagnostics are not input certificates |

Identity includes the input object, live context, model id, revision, dimensions and mathematical storage pointers. Active models/maps are immutable. In-place edits must mark revision; pointer checks do not hash all coefficients. Restricted-column mapping has its own [handoff contract](./column-handoff_en.md).

Success requires kernel certification and input acceptance. Kernels own infeasibility/unboundedness evidence; a tolerance-feasible candidate defeats a conflicting infeasibility claim. Exhausted routes and exhausted work return precision failure and iteration limit respectively, possibly retaining the best feasible point. Final reporting describes the actual returned point and clears dual diagnostics from replaced coordinates.

## State and causality

`AutoRuntime` owns cumulative work, cursor, candidates, checkpoint and decision evidence. `AutoAttempt` holds one raw output and audit. `AutoHistory` and execution events are diagnostic outputs. Disabled, full or reset journals do not determine numerical admission.

Each switch has separate prerequisite, purpose, cause, point source and basis source. A Dual handoff may motivate hybrid while NNLS supplies its point. Explicit-oracle replacement must cite the corresponding implicit failure rather than the last unrelated action. Admission arrays retain the latest observation, not a complete transition log.

References spanning callbacks bind the live context, actor token, kind and model, and are rechecked before use. Reset/reused slots lose provenance rather than creating numerical authority. Call parent, trigger and numerical producer are distinct relationships.

## Cross-layer recovery chain

| Layer | Trigger and responsibility | Work/result rule |
| --- | --- | --- |
| Dual backend | Default backend stops early numerically; try internal LU once | Same model and remaining grant; mathematical classifications/resource failure do not trigger replacement |
| AUTO | Current route is nonterminal and later capabilities exist | Forward cursor and cumulative grant; independent point/basis/candidate admission |
| Original recovery | Reduced/recovered output fails original acceptance | At most one cold solve with caller algorithm/root configuration; no recursive postsolve |
| Progress archive | A nonoptimal exit has a better original-model feasible point | May replace coordinates; retain charged work and clear incompatible dual diagnostics |

Algebraic recovery only replays the presolve journal; it neither solves, audits nor consumes iteration work. Scaling bridges restore original units; irreversible conversion falls back to the exact borrowed input. Checkpoints and progress bridges bind the actual input, not merely dimensions.

Backend events report local work; `POINT_SELECTION` reports the whole invocation's cumulative work. Selection and attempts share their caller as owner, while retry cause and returned-point source remain separate. Never add the cumulative summary to the attempt sum.

## Comparing failed actions

AUTO's cursor bounds route attempts. The shared ledger separately compares known failed numerical requests exactly: model, backend, controls, tolerances, mode and seeded indices/statuses. A shorter budget prefix may legitimately yield. Resource failures, budget exits and unknown records are not known numerical failures. Full ledgers/refused snapshots lose memory only.

New cover masters must add columns or strictly tighten dual tolerance; a failed master exits immediately. BDA NNLS blocked sets, live support changes, finite inner retries, PDHG residual baselines and exact face requests each consume their own state; see [BDA contracts](./pan-services_en.md).

## Verification entry points

Tests are registered in [test/CMakeLists.txt](../test/CMakeLists.txt). Main entry points include `test_auto_contract`, `test_auto_state`, `test_auto_capabilities`, `test_checkpoint`, `test_retry` and `test_execution`. See [development history](./development/README.md) for recorded evidence and [benchmark instructions](../data/README_en.md) for measurement methods.
