# Changelog

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
- Bilingual usage, algorithm, parameter and benchmark documentation.

### Compatibility and limits

Public interfaces are under `include/lp_simplex/`; `src/` is private.
Use matching headers and rebuild clients when upgrading. This release does not
promise cross-version ABI compatibility. Versioned CMake discovery requires an
exact version. Future source-breaking changes must be documented explicitly.

Only continuous LPs are supported. Callers must inspect the solve return code
and result status before consuming a solution. Work budgets are not CPU timeouts
or total RSS limits. Numerical backend availability affects capabilities.

Precompiled archives require compatible host numerical shared libraries; see
[INSTALL.md](INSTALL.md). No new full-suite performance measurement is claimed
for 1.0.0; published benchmark tables retain their historical measurement scope.
