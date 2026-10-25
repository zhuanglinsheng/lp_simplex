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
