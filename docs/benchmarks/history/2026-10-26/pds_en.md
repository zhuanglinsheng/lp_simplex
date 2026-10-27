> Historical results: 2026-10-26 (Asia/Shanghai), before the AUTO review. These runs include failed or timed-out cases and do not validate the current optimization. Original tables and UTC timestamps are retained. Only Markdown reports are archived; raw run directories mentioned below have been removed. [Current reports](../../README_en.md) · [History](../README_en.md)

# PDS Benchmark

[简体中文](./pds.md) | English | [Project home](../../../../README_en.md) | [Benchmark summary](./README_en.md)

Updated by `data/benchmark.py` after each instance: generate an MPS, run the configured solvers, and clean up the MPS.

## Test and statistical criteria

- Timing is wall-clock time for model loading, solving, and solution export. It excludes MPS generation and deletion, format compatibility processing, shared adapter initialization, and independent audits. CLI timing includes process startup; HiGHS uses a Python binding, so overhead differs on small instances.
- PASS requires an optimal solver status and successful audits of original-model rows, bounds, and objective values. Reference objectives are also checked when available, within tolerance times `1 + abs(reference)`. Missing references are shown as `—`; their absence does not establish independent optimality certification.
- An instance is solved only when every repetition is PASS within the time limit. Take each instance's median time, then an ordinary geometric mean, with no shift or failure penalty. Empty sets are shown as `—`.
- Each solver is compared separately with lp_simplex on their pairwise intersection of solved instances. Geometric mean times are listed as lp_simplex / solver; sets may differ between rows. Relative time is lp_simplex time divided by the solver time; above 1 means this solver is faster. The lp_simplex row omits self-comparisons.
- Successes count instances; timeouts and other unsuccessful results count individual runs. Statistics apply only to the specified run and case selection. In-progress results do not establish full-suite conclusions.
- NOT_RUN / PENDING / INCOMPLETE mean never tested / not yet tested in this run / repetitions incomplete. Instance cells show status / median seconds / median audited objective. Interpret non-PASS objectives together with their status.
- SKIPPED means omitted under a predefined rule; it counts as neither a timeout nor another unsuccessful run and is excluded from timing statistics.
- Unselected instance and solver cells retain historical results for reference. Statistics above include only measurements from participating solvers in this run; historical cells are excluded from the current ranking.

## Test information and statistics

Only the run from the most recent report update is shown; historical records remain in their original run directories.

<!-- benchmark-run:pds -->
| Test information | Value |
| --- | --- |
| State / progress | Completed; 9/9 instances |
| Time (UTC) | 2026-10-25T22:21:06.080695+00:00 |
| Solvers in this run | lp_simplex, SCIP / SoPlex, HiGHS, Clp, GLPK |
| Case selection | pds full suite; pds-20, pds-30, pds-40, pds-50, pds-60, pds-70, pds-80, pds-90, pds-100 |
| System / CPU | macOS-27.0.1-arm64-arm-64bit-Mach-O / Apple M1 |
| Memory / logical CPUs | 17179869184 bytes / 8 |
| Python / NumPy / SciPy | 3.14.2 / 2.5.3 / 1.18.1 |
| lp_simplex version | 4115c34c078d6942394cc225cc33645e18006d8b; working tree modified |
| Build / compiler | Release / Apple clang version 21.0.0 (clang-2100.3.34.2) |
| SCIP / LP backend | SCIP version 10.0.3 [precision: 8 byte] [memory: block] [mode: optimized] [LP solver: SoPlex 8.0.2] [GitHash: d409edf9f6] |
| HiGHS | 1.12.0 / scipy-private-binding |
| Algorithms (lp / scip / highs) | auto / default / choose |
| Repeats / time limit / tolerance / work budget | 1 / 300 s / 1e-07 / 300000 |
| Resources | Serial; one thread per solver; no memory cap or exclusive machine access |
| Raw records | tmp/benchmarks/2026-10-26/pds (configuration, environment, SHA-256, commands, logs, solutions, and audit residuals) |
| Clp | Coin LP version 1.17.11, build Mar 11 2026 |
| GLPK | GLPSOL--GLPK LP/MIP Solver 5.0 |
| Clp / GLPK algorithms | dual simplex / primal simplex; presolve on |
| GLPK skip policy | Skip remaining GLPK runs after its first timeout; Skipped runs: 6 |

| Solver | Solved instances | Timeouts | Other unsuccessful runs | Instances solved by both this solver and lp_simplex | Geometric means on common set / s (lp_simplex / solver) | Relative time (lp_simplex / solver) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| lp_simplex | 5/9 | 4 | 0 | — | — | — |
| SCIP / SoPlex | 6/9 | 3 | 0 | 5 | 127.25 / 54.4635 | 2.33643 |
| HiGHS | 9/9 | 0 | 0 | 5 | 127.25 / 3.65749 | 34.7916 |
| Clp | 9/9 | 0 | 0 | 5 | 127.25 / 3.58432 | 35.5018 |
| GLPK | 2/9 | 1 | 0 | 2 | 66.2865 / 49.4924 | 1.33933 |

Reproduce these cases and parameters in a new run directory:

```sh
cmake -S . -B build/benchmark -DCMAKE_BUILD_TYPE=Release -DLP_SIMPLEX_BUILD_BENCHMARKS=ON

cmake --build build/benchmark --target lp_solve_mps -j 4

python3 data/benchmark.py --dataset pds --cases pds-20 pds-30 pds-40 pds-50 pds-60 pds-70 pds-80 pds-90 pds-100 --binary build/benchmark/data/lp_solve_mps --scip scip --clp clp --glpk glpsol --repeats 1 --seconds 300.0 --iterations 300000 --lp-algorithm auto --scip-algorithm default --highs-algorithm choose --solvers lp_simplex scip highs clp glpk
```
<!-- /benchmark-run -->

## Latest per-instance results

| Instance | Size (rows / columns / nonzeros) | Reference objective | lp_simplex (status / s / objective) | SCIP / SoPlex (status / s / objective) | HiGHS (status / s / objective) | Clp (status / s / objective) | GLPK (status / s / objective) |
| --- | --- | ---: | --- | --- | --- | --- | --- |
| pds-20 | 33874 / 105728 / 230200 | 23821658640 | PASS / 48.968 / 23821658640 | PASS / 8.267 / 23821658640 | PASS / 1.281 / 23821658640 | PASS / 0.895 / 23821658640 | PASS / 25.575 / 23821658640 |
| pds-30 | 49944 / 154998 / 337144 | 21385445736 | PASS / 89.730 / 21385445736 | PASS / 25.784 / 21385445736 | PASS / 2.388 / 21385445736 | PASS / 2.088 / 21385445736 | PASS / 95.776 / 21385445736 |
| pds-40 | 66844 / 212859 / 462128 | 18855198824 | PASS / 149.349 / 18855198824 | PASS / 70.382 / 18855198824 | PASS / 4.168 / 18855198824 | PASS / 4.331 / 18855198824 | TIMEOUT / — / — |
| pds-50 | 83060 / 270095 / 585114 | 16603525724 | PASS / 220.027 / 16603525724 | PASS / 152.372 / 16603525724 | PASS / 6.803 / 16603525724 | PASS / 7.856 / 16603525724 | SKIPPED / — / — |
| pds-60 | 99431 / 329643 / 712779 | 14265904407 | PASS / 231.078 / 14265904407 | PASS / 209.629 / 14265904407 | PASS / 7.544 / 14265904407 | PASS / 9.300 / 14265904407 | SKIPPED / — / — |
| pds-70 | 114944 / 382311 / 825771 | 12241162812 | TIMEOUT / — / — | PASS / 253.354 / 12241162812 | PASS / 11.645 / 12241162812 | PASS / 16.601 / 12241162812 | SKIPPED / — / — |
| pds-80 | 129181 / 426278 / 919524 | 11469077462 | TIMEOUT / — / — | TIMEOUT / — / — | PASS / 16.834 / 11469077462 | PASS / 18.394 / 11469077462 | SKIPPED / — / — |
| pds-90 | 142823 / 466671 / 1005359 | 11087561635 | TIMEOUT / — / — | TIMEOUT / — / — | PASS / 16.376 / 11087561635 | PASS / 16.391 / 11087561635 | SKIPPED / — / — |
| pds-100 | 156243 / 505360 / 1086785 | 10928229968 | TIMEOUT / — / — | TIMEOUT / — / — | PASS / 21.527 / 10928229968 | PASS / 18.424 / 10928229968 | SKIPPED / — / — |
