> Historical benchmark snapshot: 2026-10-24 (Asia/Shanghai). Archived before the 2026-10-26 rerun. Original UTC timestamps are retained. The raw run directories referenced below are absent from this checkout; this snapshot preserves the published reports only. [Current results](../../README_en.md) · [History index](../README_en.md)

# Railway continuous LP relaxations Benchmark

[简体中文](./railway.md) | English | [Project home](../../../../README_en.md) | [Benchmark summary](./README_en.md)

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

<!-- benchmark-run:railway-auto-current-3600s-once -->
| Test information | Value |
| --- | --- |
| State / progress | Completed; 7/7 instances |
| Time (UTC) | 2026-10-23T18:29:21.672935+00:00 |
| Solvers in this run | lp_simplex, SCIP / SoPlex, HiGHS, Clp, GLPK |
| Case selection | railway full suite; rail507, rail516, rail582, rail2536, rail2586, rail4284, rail4872 |
| System / CPU | macOS-27.0.1-arm64-arm-64bit-Mach-O / Apple M1 |
| Memory / logical CPUs | 17179869184 bytes / 8 |
| Python / NumPy / SciPy | 3.14.2 / 2.5.3 / 1.18.1 |
| lp_simplex version | 00efbe80257a0530ee491cb2c44494d13b27c46c; working tree modified |
| Build / compiler | Release / Apple clang version 21.0.0 (clang-2100.3.34.2) |
| SCIP / LP backend | SCIP version 10.0.3 [precision: 8 byte] [memory: block] [mode: optimized] [LP solver: SoPlex 8.0.2] [GitHash: d409edf9f6] |
| HiGHS | 1.12.0 / scipy-private-binding |
| Algorithms (lp / scip / highs) | auto / default / choose |
| Repeats / time limit / tolerance / work budget | 1 / 3600 s / 1e-07 / 10000000 |
| Resources | Serial; one thread per solver; no memory cap or exclusive machine access |
| Raw records | tmp/benchmarks/railway-auto-current-3600s-once (configuration, environment, SHA-256, commands, logs, solutions, and audit residuals) |
| Clp | Coin LP version 1.17.11, build Mar 11 2026 |
| GLPK | GLPSOL--GLPK LP/MIP Solver 5.0 |
| Clp / GLPK algorithms | dual simplex / primal simplex; presolve on |
| GLPK skip policy | Skip remaining GLPK runs after its first timeout; Skipped runs: 3 |

| Solver | Solved instances | Timeouts | Other unsuccessful runs | Instances solved by both this solver and lp_simplex | Geometric means on common set / s (lp_simplex / solver) | Relative time (lp_simplex / solver) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| lp_simplex | 7/7 | 0 | 0 | — | — | — |
| SCIP / SoPlex | 7/7 | 0 | 0 | 7 | 27.7296 / 112.886 | 0.245642 |
| HiGHS | 7/7 | 0 | 0 | 7 | 27.7296 / 82.1388 | 0.337594 |
| Clp | 7/7 | 0 | 0 | 7 | 27.7296 / 50.2903 | 0.55139 |
| GLPK | 3/7 | 1 | 0 | 3 | 1.5798 / 3.45865 | 0.456767 |

Reproduce these cases and parameters in a new run directory:

```sh
cmake -S . -B build/benchmark -DCMAKE_BUILD_TYPE=Release -DLP_SIMPLEX_BUILD_BENCHMARKS=ON

cmake --build build/benchmark --target lp_solve_mps -j 4

python3 data/benchmark.py --dataset railway --cases rail507 rail516 rail582 rail2536 rail2586 rail4284 rail4872 --binary build/benchmark/data/lp_solve_mps --scip scip --clp clp --glpk glpsol --repeats 1 --seconds 3600.0 --iterations 10000000 --lp-algorithm auto --scip-algorithm default --highs-algorithm choose --solvers lp_simplex scip highs clp glpk
```
<!-- /benchmark-run -->

## Latest per-instance results

| Instance | Size (rows / columns / nonzeros) | Reference objective | lp_simplex (status / s / objective) | SCIP / SoPlex (status / s / objective) | HiGHS (status / s / objective) | Clp (status / s / objective) | GLPK (status / s / objective) |
| --- | --- | ---: | --- | --- | --- | --- | --- |
| rail507 | 507 / 63009 / 409349 | — | PASS / 1.346 / 172.145566677 | PASS / 5.936 / 172.145566677 | PASS / 4.206 / 172.145566677 | PASS / 3.806 / 172.145566677 | PASS / 4.002 / 172.145566677 |
| rail516 | 516 / 47311 / 314896 | — | PASS / 1.720 / 182 | PASS / 14.740 / 182 | PASS / 1.365 / 182 | PASS / 0.922 / 182 | PASS / 2.401 / 182 |
| rail582 | 582 / 55515 / 401708 | — | PASS / 1.702 / 209.712232957 | PASS / 7.835 / 209.712232957 | PASS / 4.824 / 209.712232957 | PASS / 2.372 / 209.712232957 | PASS / 4.306 / 209.712232957 |
| rail2536 | 2536 / 1081841 / 10993311 | — | PASS / 115.476 / 688.399223322 | PASS / 862.544 / 688.399223322 | PASS / 1341.432 / 688.399223322 | PASS / 391.993 / 688.399223322 | TIMEOUT / — / — |
| rail2586 | 2586 / 920683 / 8008776 | — | PASS / 153.442 / 935.921820568 | PASS / 475.296 / 935.921820568 | PASS / 589.461 / 935.921820568 | PASS / 402.475 / 935.921820568 | SKIPPED / — / — |
| rail4284 | 4284 / 1092610 / 11279748 | — | PASS / 326.133 / 1054.05449127 | PASS / 1015.431 / 1054.05449127 | PASS / 1155.374 / 1054.05449127 | PASS / 815.786 / 1054.05449127 | SKIPPED / — / — |
| rail4872 | 4872 / 968672 / 9244093 | — | PASS / 553.309 / 1509.63940162 | PASS / 818.543 / 1509.63940162 | PASS / 997.292 / 1509.63940162 | PASS / 759.785 / 1509.63940162 | SKIPPED / — / — |
