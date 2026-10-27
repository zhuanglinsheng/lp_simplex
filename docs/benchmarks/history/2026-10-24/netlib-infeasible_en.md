> Historical benchmark snapshot: 2026-10-24 (Asia/Shanghai). Archived before the 2026-10-26 rerun. Original UTC timestamps are retained. The raw run directories referenced below are absent from this checkout; this snapshot preserves the published reports only. [Current results](../../README_en.md) · [History index](../README_en.md)

# Netlib infeasible Benchmark

[简体中文](./netlib-infeasible.md) | English | [Project home](../../../../README_en.md) | [Benchmark summary](./README_en.md)

Updated by `data/benchmark.py` after each instance: generate an MPS, run the configured solvers, and clean up the MPS.

## Test and statistical criteria

- Timing is wall-clock time for model loading, solving, and solution export. It excludes MPS generation and deletion, format compatibility processing, shared adapter initialization, and independent audits. CLI timing includes process startup; HiGHS uses a Python binding, so overhead differs on small instances.
- An instance is solved only when every repetition is PASS within the time limit. Take each instance's median time, then an ordinary geometric mean, with no shift or failure penalty. Empty sets are shown as `—`.
- Each solver is compared separately with lp_simplex on their pairwise intersection of solved instances. Geometric mean times are listed as lp_simplex / solver; sets may differ between rows. Relative time is lp_simplex time divided by the solver time; above 1 means this solver is faster. The lp_simplex row omits self-comparisons.
- Successes count instances; timeouts and other unsuccessful results count individual runs. Statistics apply only to the specified run and case selection. In-progress results do not establish full-suite conclusions.
- NOT_RUN / PENDING / INCOMPLETE mean never tested / not yet tested in this run / repetitions incomplete. Instance cells show status / median seconds / median audited objective. Interpret non-PASS objectives together with their status.
- SKIPPED means omitted under a predefined rule; it counts as neither a timeout nor another unsuccessful run and is excluded from timing statistics.
- Unselected instance and solver cells retain historical results for reference. Statistics above include only measurements from participating solvers in this run; historical cells are excluded from the current ranking.

- PASS means an explicit infeasible status within the time limit, against the Netlib infeasible classification. Farkas certificates are not independently verified; PASS is a status check, not a certificate audit.
- WRONG_STATUS means optimal or unbounded; AMBIGUOUS means infeasible or unbounded; SOLVER_FAILED / ERROR means no definitive decision or an execution failure. TIMEOUT takes precedence.
- Original MPS files are passed to the solvers unchanged. greenbea contains multiple BOUNDS vectors and duplicate bounds; pang also exposes reader compatibility differences. Related failures include input compatibility effects. HiGHS read warnings are recorded in reader_status; solver logs are retained separately.

## Test information and statistics

Only the run from the most recent report update is shown; historical records remain in their original run directories.

<!-- benchmark-run:netlib-infeasible-90s-final -->
| Test information | Value |
| --- | --- |
| State / progress | Completed; 29/29 instances |
| Time (UTC) | 2026-10-24T06:55:30.641244+00:00 |
| Solvers in this run | lp_simplex, SCIP / SoPlex, HiGHS, Clp, GLPK |
| Case selection | netlib-infeasible full suite; bgdbg1, bgetam, bgindy, bgprtr, box1, ceria3d, chemcom, cplex1, cplex2, ex72a, ex73a, forest6, galenet, gosh, gran, greenbea, itest2, itest6, klein1, klein2, klein3, mondou2, pang, pilot4i, qual, reactor, refinery, vol1, woodinfe |
| System / CPU | macOS-27.0.1-arm64-arm-64bit-Mach-O / Apple M1 |
| Memory / logical CPUs | 17179869184 bytes / 8 |
| Python / NumPy / SciPy | 3.14.2 / 2.5.3 / 1.18.1 |
| lp_simplex version | d4101511d0a7ed27c04e25406bf1707dc63a0aaf; working tree modified |
| Build / compiler | Release / Apple clang version 21.0.0 (clang-2100.3.34.2) |
| SCIP / LP backend | SCIP version 10.0.3 [precision: 8 byte] [memory: block] [mode: optimized] [LP solver: SoPlex 8.0.2] [GitHash: d409edf9f6] |
| HiGHS | 1.12.0 / scipy-private-binding |
| Algorithms (lp / scip / highs) | auto / default / choose |
| Repeats / time limit / tolerance / work budget | 1 / 90 s / 1e-07 / 300000 |
| Resources | Serial; one thread per solver; no memory cap or exclusive machine access |
| Raw records | lp_simplex: tmp/benchmarks/netlib-infeasible-certificates-latest-once; SCIP / HiGHS / Clp / GLPK: tmp/benchmarks/netlib-infeasible-90s-final |
| Clp | Coin LP version 1.17.11, build Mar 11 2026 |
| GLPK | GLPSOL--GLPK LP/MIP Solver 5.0 |
| Supplemental run | tmp/benchmarks/netlib-infeasible-greenbea-serial; greenbea only; corrected HiGHS read-warning handling; original attempt retained at tmp/benchmarks/netlib-infeasible-90s-once |
| Clp / GLPK algorithms | dual simplex / primal simplex; presolve on |
| GLPK skip policy | Skip remaining GLPK runs after its first timeout; Skipped runs: 0 |

| Solver | Solved instances | Timeouts | Other unsuccessful runs | Instances solved by both this solver and lp_simplex | Geometric means on common set / s (lp_simplex / solver) | Relative time (lp_simplex / solver) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| lp_simplex | 29/29 | 0 | 0 | — | — | — |
| SCIP / SoPlex | 28/29 | 0 | 1 | 28 | 0.0281044 / 0.0542667 | 0.517893 |
| HiGHS | 29/29 | 0 | 0 | 29 | 0.0279029 / 0.00336188 | 8.2998 |
| Clp | 29/29 | 0 | 0 | 29 | 0.0279029 / 0.0295637 | 0.943822 |
| GLPK | 27/29 | 0 | 2 | 27 | 0.0277222 / 0.0161276 | 1.71892 |

| Solver | Wrong decisions | Ambiguous decisions | Execution or numerical failures | Skipped runs |
| --- | ---: | ---: | ---: | ---: |
| lp_simplex | 0 | 0 | 0 | 0 |
| SCIP / SoPlex | 0 | 1 | 0 | 0 |
| HiGHS | 0 | 0 | 0 | 0 |
| Clp | 0 | 0 | 0 | 0 |
| GLPK | 0 | 0 | 2 | 0 |

| Instance | Solver | Result | Termination information |
| --- | --- | --- | --- |
| cplex1 | SCIP / SoPlex | AMBIGUOUS | problem is solved [infeasible or unbounded] |
| pang | GLPK | SOLVER_FAILED | not optimal |
| greenbea | GLPK | SOLVER_FAILED | not optimal |

Reproduce these cases and parameters in a new run directory:

```sh
cmake -S . -B build/benchmark -DCMAKE_BUILD_TYPE=Release -DLP_SIMPLEX_BUILD_BENCHMARKS=ON

cmake --build build/benchmark --target lp_solve_mps -j 4

python3 data/benchmark.py --dataset netlib-infeasible --cases bgdbg1 bgetam bgindy bgprtr box1 ceria3d chemcom cplex1 cplex2 ex72a ex73a forest6 galenet gosh gran greenbea itest2 itest6 klein1 klein2 klein3 mondou2 pang pilot4i qual reactor refinery vol1 woodinfe --binary build/benchmark/data/lp_solve_mps --scip scip --clp clp --glpk glpsol --repeats 1 --seconds 90.0 --iterations 300000 --lp-algorithm auto --scip-algorithm default --highs-algorithm choose --solvers lp_simplex scip highs clp glpk
```
<!-- /benchmark-run -->

## Latest per-instance results

| Instance | Size (rows / columns / nonzeros) | Expected status | lp_simplex (status / s / objective) | SCIP / SoPlex (status / s / objective) | HiGHS (status / s / objective) | Clp (status / s / objective) | GLPK (status / s / objective) |
| --- | --- | ---: | --- | --- | --- | --- | --- |
| bgdbg1 | 348 / 407 / 1440 | INFEASIBLE | PASS / 0.037 / — | PASS / 0.038 / — | PASS / 0.001 / — | PASS / 0.017 / — | PASS / 0.010 / — |
| bgetam | 400 / 688 / 2409 | INFEASIBLE | PASS / 0.016 / — | PASS / 0.043 / — | PASS / 0.004 / — | PASS / 0.014 / — | PASS / 0.015 / — |
| bgindy | 2671 / 10116 / 65502 | INFEASIBLE | PASS / 0.042 / — | PASS / 0.093 / — | PASS / 0.039 / — | PASS / 0.039 / — | PASS / 0.136 / — |
| bgprtr | 20 / 34 / 64 | INFEASIBLE | PASS / 0.014 / — | PASS / 0.034 / — | PASS / 0.001 / — | PASS / 0.014 / — | PASS / 0.007 / — |
| box1 | 231 / 261 / 651 | INFEASIBLE | PASS / 0.016 / — | PASS / 0.038 / — | PASS / 0.001 / — | PASS / 0.014 / — | PASS / 0.008 / — |
| ceria3d | 3576 / 824 / 17602 | INFEASIBLE | PASS / 0.024 / — | PASS / 0.051 / — | PASS / 0.012 / — | PASS / 0.137 / — | PASS / 0.019 / — |
| chemcom | 288 / 720 / 1566 | INFEASIBLE | PASS / 0.016 / — | PASS / 0.044 / — | PASS / 0.003 / — | PASS / 0.014 / — | PASS / 0.011 / — |
| cplex1 | 3005 / 3221 / 8944 | INFEASIBLE | PASS / 0.023 / — | AMBIGUOUS / — / — | PASS / 0.089 / — | PASS / 0.096 / — | PASS / 0.044 / — |
| cplex2 | 224 / 221 / 1058 | INFEASIBLE | PASS / 0.027 / — | PASS / 0.061 / — | PASS / 0.005 / — | PASS / 0.017 / — | PASS / 0.011 / — |
| ex72a | 197 / 215 / 467 | INFEASIBLE | PASS / 0.018 / — | PASS / 0.043 / — | PASS / 0.001 / — | PASS / 0.015 / — | PASS / 0.009 / — |
| ex73a | 193 / 211 / 457 | INFEASIBLE | PASS / 0.015 / — | PASS / 0.043 / — | PASS / 0.001 / — | PASS / 0.014 / — | PASS / 0.009 / — |
| forest6 | 66 / 95 / 210 | INFEASIBLE | PASS / 0.016 / — | PASS / 0.048 / — | PASS / 0.001 / — | PASS / 0.015 / — | PASS / 0.007 / — |
| galenet | 8 / 8 / 16 | INFEASIBLE | PASS / 0.014 / — | PASS / 0.036 / — | PASS / 0.000 / — | PASS / 0.014 / — | PASS / 0.007 / — |
| gosh | 3792 / 10733 / 97231 | INFEASIBLE | PASS / 2.484 / — | PASS / 0.158 / — | PASS / 0.051 / — | PASS / 0.783 / — | PASS / 0.395 / — |
| gran | 2658 / 2520 / 20106 | INFEASIBLE | PASS / 0.028 / — | PASS / 0.052 / — | PASS / 0.012 / — | PASS / 0.102 / — | PASS / 0.035 / — |
| greenbea | 2393 / 5405 / 30883 | INFEASIBLE | PASS / 0.040 / — | PASS / 0.075 / — | PASS / 0.021 / — | PASS / 0.186 / — | SOLVER_FAILED / — / — |
| itest2 | 9 / 4 / 17 | INFEASIBLE | PASS / 0.021 / — | PASS / 0.042 / — | PASS / 0.000 / — | PASS / 0.013 / — | PASS / 0.007 / — |
| itest6 | 11 / 8 / 20 | INFEASIBLE | PASS / 0.027 / — | PASS / 0.035 / — | PASS / 0.000 / — | PASS / 0.013 / — | PASS / 0.008 / — |
| klein1 | 54 / 54 / 696 | INFEASIBLE | PASS / 0.028 / — | PASS / 0.043 / — | PASS / 0.002 / — | PASS / 0.014 / — | PASS / 0.009 / — |
| klein2 | 477 / 54 / 4585 | INFEASIBLE | PASS / 0.035 / — | PASS / 0.126 / — | PASS / 0.009 / — | PASS / 0.034 / — | PASS / 0.020 / — |
| klein3 | 994 / 88 / 12107 | INFEASIBLE | PASS / 0.092 / — | PASS / 0.113 / — | PASS / 0.030 / — | PASS / 0.104 / — | PASS / 0.076 / — |
| mondou2 | 312 / 604 / 1208 | INFEASIBLE | PASS / 0.022 / — | PASS / 0.045 / — | PASS / 0.002 / — | PASS / 0.016 / — | PASS / 0.010 / — |
| pang | 361 / 460 / 2652 | INFEASIBLE | PASS / 0.023 / — | PASS / 0.059 / — | PASS / 0.008 / — | PASS / 0.020 / — | SOLVER_FAILED / — / — |
| pilot4i | 410 / 1000 / 5141 | INFEASIBLE | PASS / 0.020 / — | PASS / 0.049 / — | PASS / 0.003 / — | PASS / 0.108 / — | PASS / 0.049 / — |
| qual | 323 / 464 / 1646 | INFEASIBLE | PASS / 0.025 / — | PASS / 0.116 / — | PASS / 0.005 / — | PASS / 0.025 / — | PASS / 0.014 / — |
| reactor | 318 / 637 / 2420 | INFEASIBLE | PASS / 0.015 / — | PASS / 0.043 / — | PASS / 0.002 / — | PASS / 0.014 / — | PASS / 0.013 / — |
| refinery | 323 / 464 / 1626 | INFEASIBLE | PASS / 0.028 / — | PASS / 0.055 / — | PASS / 0.005 / — | PASS / 0.022 / — | PASS / 0.012 / — |
| vol1 | 323 / 464 / 1646 | INFEASIBLE | PASS / 0.045 / — | PASS / 0.053 / — | PASS / 0.008 / — | PASS / 0.034 / — | PASS / 0.012 / — |
| woodinfe | 35 / 89 / 140 | INFEASIBLE | PASS / 0.014 / — | PASS / 0.042 / — | PASS / 0.000 / — | PASS / 0.015 / — | PASS / 0.007 / — |
