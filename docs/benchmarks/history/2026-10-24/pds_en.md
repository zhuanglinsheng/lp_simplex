> Historical benchmark snapshot: 2026-10-24 (Asia/Shanghai). Archived before the 2026-10-26 rerun. Original UTC timestamps are retained. The raw run directories referenced below are absent from this checkout; this snapshot preserves the published reports only. [Current results](../../README_en.md) · [History index](../README_en.md)

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

<!-- benchmark-run:pds-accepted-auto-20261024-once -->
| Test information | Value |
| --- | --- |
| State / progress | Completed; 9/9 instances |
| Time (UTC) | 2026-10-24T05:59:44.178021+00:00 |
| Solvers in this run | lp_simplex, SCIP / SoPlex, HiGHS, Clp, GLPK |
| Case selection | pds full suite; pds-20, pds-30, pds-40, pds-50, pds-60, pds-70, pds-80, pds-90, pds-100 |
| System / CPU | macOS-27.0.1-arm64-arm-64bit-Mach-O / Apple M1 |
| Memory / logical CPUs | 17179869184 bytes / 8 |
| Python / NumPy / SciPy | 3.14.2 / 2.5.3 / 1.18.1 |
| lp_simplex version | d4101511d0a7ed27c04e25406bf1707dc63a0aaf; working tree modified |
| Build / compiler | Release / Apple clang version 21.0.0 (clang-2100.3.34.2) |
| SCIP / LP backend | SCIP version 10.0.3 [precision: 8 byte] [memory: block] [mode: optimized] [LP solver: SoPlex 8.0.2] [GitHash: d409edf9f6] |
| HiGHS | 1.12.0 / scipy-private-binding |
| Algorithms (lp / scip / highs) | auto / default / choose |
| Repeats / time limit / tolerance / work budget | 1 / 300 s / 1e-07 / 300000 |
| Resources | Serial; one thread per solver; no memory cap or exclusive machine access |
| Raw records | tmp/benchmarks/pds-auto-refined-five-solvers-300s-final; tmp/benchmarks/pds-accepted-auto-20261024-once; tmp/benchmarks/pds-accepted-auto-20261024-lp-only-v2 |
| Clp | Coin LP version 1.17.11, build Mar 11 2026 |
| GLPK | GLPSOL--GLPK LP/MIP Solver 5.0 |
| Clp / GLPK algorithms | dual simplex / primal simplex; presolve on |
| GLPK skip policy | Skipped runs: 0 |

| Solver | Solved instances | Timeouts | Other unsuccessful runs | Instances solved by both this solver and lp_simplex | Geometric means on common set / s (lp_simplex / solver) | Relative time (lp_simplex / solver) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| lp_simplex | 9/9 | 0 | 0 | — | — | — |
| SCIP / SoPlex | 6/9 | 3 | 0 | 6 | 80.4532 / 68.5872 | 1.17301 |
| HiGHS | 9/9 | 0 | 0 | 9 | 111.852 / 7.20217 | 15.5303 |
| Clp | 9/9 | 0 | 0 | 9 | 111.852 / 7.38672 | 15.1423 |
| GLPK | 2/9 | 7 | 0 | 2 | 42.2868 / 53.8081 | 0.785881 |

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
| pds-20 | 33874 / 105728 / 230200 | 23821658640 | PASS / 34.717 / 23821658640 | PASS / 9.058 / 23821658640 | PASS / 1.479 / 23821658640 | PASS / 1.061 / 23821658640 | PASS / 26.418 / 23821658640 |
| pds-30 | 49944 / 154998 / 337144 | 21385445736 | PASS / 51.508 / 21385445736 | PASS / 27.050 / 21385445736 | PASS / 2.537 / 21385445736 | PASS / 2.204 / 21385445736 | PASS / 109.598 / 21385445736 |
| pds-40 | 66844 / 212859 / 462128 | 18855198824 | PASS / 92.013 / 18855198824 | PASS / 61.441 / 18855198824 | PASS / 4.596 / 18855198824 | PASS / 5.133 / 18855198824 | TIMEOUT / — / — |
| pds-50 | 83060 / 270095 / 585114 | 16603525724 | PASS / 109.296 / 16603525724 | PASS / 137.162 / 16603525724 | PASS / 5.687 / 16603525724 | PASS / 6.534 / 16603525724 | TIMEOUT / — / — |
| pds-60 | 99431 / 329643 / 712779 | 14265904407 | PASS / 102.854 / 14265904407 | PASS / 207.310 / 14265904407 | PASS / 7.527 / 14265904407 | PASS / 8.322 / 14265904407 | TIMEOUT / — / — |
| pds-70 | 114944 / 382311 / 825771 | 12241162812 | PASS / 146.615 / 12241162812 | PASS / 243.189 / 12241162812 | PASS / 8.533 / 12241162812 | PASS / 12.022 / 12241162812 | TIMEOUT / — / — |
| pds-80 | 129181 / 426278 / 919524 | 11469077462 | PASS / 190.654 / 11469077462 | TIMEOUT / — / — | PASS / 12.701 / 11469077462 | PASS / 16.419 / 11469077462 | TIMEOUT / — / — |
| pds-90 | 142823 / 466671 / 1005359 | 11087561635 | PASS / 236.439 / 11087561635 | TIMEOUT / — / — | PASS / 27.825 / 11087561635 | PASS / 18.074 / 11087561635 | TIMEOUT / — / — |
| pds-100 | 156243 / 505360 / 1086785 | 10928229968 | PASS / 224.168 / 10928229968 | TIMEOUT / — / — | PASS / 23.413 / 10928229968 | PASS / 28.121 / 10928229968 | TIMEOUT / — / — |
