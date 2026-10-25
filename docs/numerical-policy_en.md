# Numerical and resource policies

[简体中文](./numerical-policy.md) | English | [Parameters/results](./parameters_en.md) | [Documentation](./README_en.md)

Public tolerances express user requirements; internal policies govern search, stability and resources. They are not interchangeable. Constants live in [src/common/policy](../src/common/policy); source is authoritative for complete defaults. A policy change must state units, scope, conservative failure behavior and validation rationale.

## Policy categories

| Category | Definition | Meaning |
| --- | --- | --- |
| Production defaults | [defaults.h](../src/common/policy/defaults.h) | Tolerances, work, progress observation and notification |
| LU/QR and refinement | [factorization.h](../src/common/policy/factorization.h) | Pivot admission, update lengths, residuals and refinement limits |
| Pricing | [factorization.h](../src/common/policy/factorization.h) | Ties/audits, protected weights and Devex resets |
| Search and observations | [search.h](../src/common/policy/search.h) | AUTO grants, Dual observation, finite restarts and presolve search |
| BDA/PDHG | [pan.h](../src/common/policy/pan.h) | Steps, stagnation, restoration, directions and local proposals |
| Workspace and certificates | [resources.h](../src/common/policy/resources.h) | Optional storage, scaling and evidence limits |
| Experiment names | [controls.h](../src/common/policy/controls.h) | Names parsed by the explicit experiment entry, not production options |

## Key defaults and units

| Item | Default | Interpretation |
| --- | --- | --- |
| Primal / dual / pivot tolerance | `1e-7 / 1e-7 / 1e-9` | Different units for feasibility, search and pivot/rank admission |
| Work budget | Manual `100000`, AUTO `300000` | Iteration work, not time; public default options are authoritative |
| Progress observation / notification | `128` work / `5` CPU seconds | Deterministic sampling versus UI throttling; no solve timeout |
| Relative stable pivot / absolute floor | `1e-7 / 1e-14` | Weak updates request refactorization, not singularity/infeasibility |
| Internal LU zero-pivot threshold | `1e-13` | Elimination/update protection, distinct from the user pivot tolerance |
| Markowitz / KLU pivot preference | `.1 / .1` | Independent backend choices, not final accuracy |
| Basis / FT update length | `512 / 512` | Finite-chain protection, not an error guarantee |
| Compact residual / refinement / audit period | `1e-11 / 2 / 32` | Internal relative residual and bounded correction; terminal refinement permits `3` steps |
| Initial Dual observation | `8 * 128 = 1024` work | Opportunity to observe continuation/handoff, not a fixed solve deadline |
| Dual crash attempts | `min(INT_MAX, 4*variables+rows)` | Check before arithmetic; initialization does not expand the main grant |
| Dual nonpivot restarts | Saturated `4*rows+64` | Prevent zero-iteration restart loops; separate from pivot work |
| PDHG observation / stagnant window | `128 / 8192` steps | More than `1%` cumulative gain over the historical best in one metric |
| PDHG step contraction | `.99` | Safety margin for the diagonally preconditioned norm bound |
| First-order backtracking / trials / CG cap | `.5 / 20 / 256` | Finite local proposal cost |

These are empirical or numerical policies. Local `16/32/64 * DBL_EPSILON` margins must be read with their absolute sums, dimensions and error scale, not merged into one global epsilon. Absolute LU floors, size gates and some gain floors remain heuristic; replacing them needs conditioning, error or performance evidence.

## Certification and scale

### Tableau

Artificial-basis cleanup checks a nonzero pivot against row scale and pivot tolerance. A redundant row is retained only when all real coefficients and RHS are exactly zero. Unreliable elimination or nonfinite operations return precision failure without erasing an equation or inferring infeasibility. Nonfinite Phase I objective is no infeasibility proof. Bland uses requested tolerance and the smallest eligible index, not decimal threshold schedules.

### BDA and PDHG

Low-progress gain uses `max(64*DBL_EPSILON, units*dual_tolerance) * max(|previous cost without constant|, sum|c_j*x_j|)`, with units `100/10/1000` for ordinary/face-expansion/face-contraction progress. Objective constants do not participate; arbitrary variable transformations need not preserve floating-point trajectories.

PDHG point acceptance checks absolute primal error and row scale: `|b_i-A_i*x| <= primal_tolerance*(|b_i|+sum_j |A_ij|*max(1,|x_j|))`. Changing the weighted merit resets its baseline; ordinary restart retains the historical best. Restoration requires stagnation, remaining grant, objective ceiling and physical residual below half the previous attempt's. See [phase contracts](./pan-services_en.md).

Handoff gap uses the actual bound-constrained dual infimum: a missing bound needed by a nonzero reduced cost makes the infimum infinite; zero reduced cost contributes zero. Relative gap divides by `max(|primal|,|dual|)`, handling equality separately. It is a handoff diagnostic, not a replacement for native KKT.

### Direct QR multiplier recovery

If standard-basis reconstruction still fails its multiplier audit, the [independent QR module](../src/bda/simplex_pan_basis_qr.c) uses DGELSY on original basis equations without normal equations. `RCOND=DBL_EPSILON` resolves arithmetic rank, rather than reapplying user pivot tolerance as a global condition-number cutoff. Full row rank, finite output and original-basis residual are required; reduced costs, feasibility and gap checks remain unchanged. Only success writes multipliers; failure preserves output and does not change the basis or LP grant.

### Input, scaling and presolve

Public input rejects NaN/infinity in active numerical fields and invalid sparse layouts. Unused bound endpoints are ignored. MPS numbers must parse completely; truncation, trailing junk, nonfinite values and range-expansion overflow are rejected rather than replaced with zero.

Parallel-row reduction confirms exact proportionality of raw coefficients; approximate equality cannot authorize contradictory classification. Scaling must be finite and reversible, otherwise original coordinates are retained. Recovery journals prepare before publication and replay dependencies backwards; final output still receives an original-model audit.

## Storage and failure boundaries

| Owner / optional workspace | Limit or commit rule |
| --- | --- |
| AUTO dense model / conversion | `32 / 64 MiB`; failed admission skips the route |
| Cover fixed owned workspace | `256 MiB`, including point scratch/bitmap, not every delegated factorization |
| BDA direct QR | `32 MiB`, including new matrix, indices and LAPACK scratch |
| Exact failed-basis snapshot / native face request | `4 MiB` each; refusal only loses optional memory |
| Dominance indexes | `32 MiB` peak and `64*(nnz+m+n)` search work; refusal disables optional search |
| Infeasibility evidence | rank `128`, bits `8192`; unconfirmed evidence grants no classification |
| Scaling | `4` passes, exponents `[-48,48]`; irreversible conversion falls back |

Check element counts, padding, byte products/sums and backend integer ranges before allocation. Required model/startup refusal returns `MemoryAllocError`; optional refusal revokes that capability rather than inventing numerical or mathematical failure.

- **Model/representation:** sparse caches, substitution records and compaction prepare storage before publishing layouts; failure preserves committed input and cleans scratch.
- **LU/basis:** failed complete factorization revokes numerical validity; resources and singularity are distinct. FT/eta prepare growth before modifying committed factors, preserving prior solvability on refusal. Basis import re-factors and classifies.
- **QR:** dynamic append prepares capacity/temporary rows before damaging old columns, preserving their mathematical usability on append refusal. This does not extend transaction guarantees to every deletion, swap or reload.
- **Presolve:** startup/reduced-output refusal cleans its owners; journal dependencies restore backwards. After dominance fixes a column, the next queue epoch eliminates and records it.

Local limits are not total RSS; work is not CPU time. Dimension/narrow-address probes are not giant LP runs or native 32-bit validation. See [delivery scope](../README_en.md#delivery-and-compatibility-boundary) for platform and deployment boundaries.

## Configuration and reproduction

Production fixes defaults. Only the explicit experiment entry captures numerical environment settings and then uses a root snapshot. Presence flags still treat string `0` as present; integral work settings reject negative/malformed input and saturate positive overflow; weights must parse fully, be finite and in range. See [configuration boundaries](./runtime-configuration_en.md).

Record source, build/backends, public options, input, threads, runs and solution audits when comparing. [Data instructions](../data/README_en.md) own benchmark methodology; [development history](./development/README.md) owns recorded results and fixes.
