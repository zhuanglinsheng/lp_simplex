# Changelog

## 1.1.0 — 2026-10-28

Improve sparse simplex execution, phase-aware AUTO recovery and reproducible
benchmark reporting while preserving the public C/C++ source interfaces.

### Solver improvements

- Retain productive Dual searches using phase-relevant input audits and an
  objective progress signal guarded against cancellation and working-cost drift.
  Numerical failures, cycle protection and cumulative work limits keep priority.
- Carry sparse supports through compact bases, KLU solves and eta updates, with
  existing residual audits and dense fallbacks. Restore cost-based KLU
  refactorization on large bases.
- Audit optimistic steepest-edge weights and correctly track dense pricing
  writes in shared sparse storage.
- Separate presolve bound-rounding allowances from feasibility tolerances,
  preserve small coordinates during postsolve cancellation, and extend
  conservative no-fill signed-unit aggregation.
- Use bounded initial cost separation and limit repeated working-objective
  perturbations, restoring and auditing the original objective before success.
- Distinguish feasibility-only PDHG stagnation from joint stagnation, retain
  bounded recovery and add sparse multiplier repair with rollback.
- Add bounded, witnessed equality-dependency elimination as an explicit
  experiment. It remains disabled by default because reduced dimensions did
  not consistently improve complete solves.

### Benchmark tooling and validation

- Separate report rendering from solver execution. Persist full-precision
  measurements and fingerprints in bilingual Markdown reports; reuse unchanged
  solver measurements independently on the same machine.
- Add explicit `--rerun`, `--refresh-reports` and compatible legacy-report import,
  preserving CLI timing, timeout cleanup and independent original-model audits.
- The recorded macOS arm64 full-backend Release candidate passed all 182 CTest
  cases. Default AUTO passed Netlib 114/114, infeasible Netlib 29/29, PDS 9/9 and
  Railway 7/7 under the existing limits. These measurements describe their
  recorded numerical build and host, not a uniform speedup across cases or
  platforms. Earlier reports are retained as dated Markdown snapshots.

### Compatibility

Public headers are unchanged from 1.0.0. Continue to use matching headers and
libraries and rebuild clients when upgrading; cross-version ABI compatibility
is not promised. Versioned CMake consumers should request
`find_package(lp_simplex 1.1.0 EXACT CONFIG REQUIRED)`.
The published 1.0.0 source and assets remain unchanged.

## 1.0.0 — 2026-10-26

First stable release of lp_simplex, a continuous linear programming solver with
a C90 API, C++ integration and MPS input.

- AUTO orchestration, revised dual simplex, tableau simplex and basis-deficiency
  allowing (BDA/Pan) optimization, with internal PDHG and feasibility recovery.
- Presolve/postsolve, original-model solution auditing, progress callbacks and
  explicit termination/diagnostic results.
- Optional SuiteSparse KLU/SPQR numerical backends and GMP exact terminal
  infeasibility verification; BLAS and LAPACK are required.
- Relocatable CMake installation, public version metadata and runtime version
  query, independent C/C++ consumer verification and platform release archives.
- Native Windows x64 MSVC Release/Debug validation and ZIP packages with bundled
  numerical DLLs; the Windows baseline requires only BLAS/LAPACK.
- Bilingual usage, algorithm, parameter and benchmark documentation.

Native Windows Release/Debug and relocated ZIP consumers have passed validation.
Cross-platform release acceptance is tracked in the
[release rehearsal](https://github.com/zhuanglinsheng/lp_simplex/actions/runs/40226709492).
The frozen source is branch `1.0.0` at annotated tag `v1.0.0`.
This branch is retained unchanged; development continues on `main`, and release
fixes use a new version branch and tag.

### Compatibility and limits

Public interfaces are under `include/lp_simplex/`; `src/` is private.
Use matching headers and rebuild clients when upgrading. This release does not
promise cross-version ABI compatibility. Versioned CMake discovery requires an
exact version. Future source-breaking changes must be documented explicitly.

Only continuous LPs are supported. Callers must inspect the solve return code
and result status before consuming a solution. Work budgets are not CPU timeouts
or total RSS limits. Numerical backend availability affects capabilities.

Linux/macOS archives require compatible host numerical shared libraries; Windows
archives bundle numerical DLLs and require the MSVC runtime. See
[INSTALL.md](INSTALL.md). No new full-suite performance measurement is claimed
for 1.0.0; published benchmark tables retain their historical measurement scope.
