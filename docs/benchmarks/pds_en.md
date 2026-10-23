# PDS Benchmark

[简体中文](./pds.md) | English | [Project home](../../README_en.md) | [Benchmark summary](./README_en.md)

Updated by `data/benchmark.py` after each instance: generate an MPS, run the configured solvers, and clean up the MPS.

## Test and statistical criteria

- Timing is wall-clock time for model loading, solving, and solution export. It excludes MPS generation and deletion, format compatibility processing, shared adapter initialization, and independent audits. CLI timing includes process startup; HiGHS uses a Python binding, so overhead differs on small instances.
- PASS requires an optimal solver status and successful audits of original-model rows, bounds, and objective values. Reference objectives are also checked when available, within tolerance times `1 + abs(reference)`. Missing references are shown as `—`; their absence does not establish independent optimality certification.
- An instance is solved only when every repetition is PASS within the time limit. Take each instance's median time, then an ordinary geometric mean, with no shift or failure penalty. Empty sets are shown as `—`.
- Each solver is compared separately with lp_simplex on their pairwise intersection of solved instances. Geometric mean times are listed as lp_simplex / solver; sets may differ between rows. Relative time is lp_simplex time divided by the solver time; above 1 means this solver is faster. The lp_simplex row omits self-comparisons.
- Successes count instances; timeouts and other unsuccessful results count individual runs. Statistics apply only to the specified run and case selection. In-progress results do not establish full-suite conclusions.
- NOT_RUN / PENDING / INCOMPLETE mean never tested / not yet tested in this run / repetitions incomplete. Instance cells show status / median seconds / median audited objective. Interpret non-PASS objectives together with their status.
- SKIPPED means omitted under a predefined rule; it counts as neither a timeout nor another unsuccessful run and is excluded from timing statistics.
- Each row retains the latest record for that instance. Unselected instances are preserved. Environments and settings may differ between runs, so speed rankings are not pooled across runs.

## Test information and statistics

<!-- benchmark-run:pds-auto-refined-five-solvers-300s-final -->
| Test information | Value |
| --- | --- |
| State / progress | Completed; 9/9 instances |
| Time (UTC) | 2026-10-23T15:41:47.536059+00:00 |
| Solvers in this run | lp_simplex, SCIP / SoPlex, HiGHS, Clp, GLPK |
| Case selection | pds full suite; pds-20, pds-30, pds-40, pds-50, pds-60, pds-70, pds-80, pds-90, pds-100 |
| System / CPU | macOS-27.0.1-arm64-arm-64bit-Mach-O / Apple M1 |
| Memory / logical CPUs | 17179869184 bytes / 8 |
| Python / NumPy / SciPy | 3.14.2 / 2.5.3 / 1.18.1 |
| lp_simplex version | 1ab947ebcbbbdff01b10b99994d58bc1a7cdd3e6; working tree modified |
| Build / compiler | Release / Apple clang version 21.0.0 (clang-2100.3.34.2) |
| SCIP / LP backend | SCIP version 10.0.3 [precision: 8 byte] [memory: block] [mode: optimized] [LP solver: SoPlex 8.0.2] [GitHash: d409edf9f6] |
| HiGHS | 1.12.0 / scipy-private-binding |
| Algorithms (lp / scip / highs) | auto / default / choose |
| Repeats / time limit / tolerance / work budget | 1 / 300 s / 1e-07 / 300000 |
| Resources | Serial; one thread per solver; no memory cap or exclusive machine access |
| Raw records | tmp/benchmarks/pds-auto-refined-five-solvers-300s-final (configuration, environment, SHA-256, commands, logs, solutions, and audit residuals) |
| Clp | Coin LP version 1.17.11, build Mar 11 2026 |
| GLPK | GLPSOL--GLPK LP/MIP Solver 5.0 |
| Clp / GLPK algorithms | dual simplex / primal simplex; presolve on |

| Solver | Solved instances | Timeouts | Other unsuccessful runs | Instances solved by both this solver and lp_simplex | Geometric means on common set / s (lp_simplex / solver) | Relative time (lp_simplex / solver) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| lp_simplex | 9/9 | 0 | 0 | — | — | — |
| SCIP / SoPlex | 6/9 | 3 | 0 | 6 | 106.83 / 69.9679 | 1.52685 |
| HiGHS | 9/9 | 0 | 0 | 9 | 143.535 / 7.30131 | 19.6588 |
| Clp | 9/9 | 0 | 0 | 9 | 143.535 / 7.31296 | 19.6275 |
| GLPK | 2/9 | 7 | 0 | 2 | 64.4913 / 53.6572 | 1.20191 |

Reproduce these cases and parameters in a new run directory:

```sh
python3 data/benchmark.py --dataset pds --cases pds-100 pds-20 pds-30 pds-40 pds-50 pds-60 pds-70 pds-80 pds-90 --binary /Users/zhuanglinsheng/Documents/Github/lp_simplex/tmp/build-auto-implementation/data/netlib/test_netlib --scip /opt/homebrew/bin/scip --clp /opt/homebrew/bin/clp --glpk /opt/homebrew/bin/glpsol --repeats 1 --seconds 300.0 --iterations 300000 --lp-algorithm auto --scip-algorithm default --highs-algorithm choose --solvers lp_simplex scip highs clp glpk
```
<!-- /benchmark-run -->

## Latest per-instance results

| Instance | Size (rows / columns / nonzeros) | Reference objective | lp_simplex (status / s / objective) | SCIP / SoPlex (status / s / objective) | HiGHS (status / s / objective) | Clp (status / s / objective) | GLPK (status / s / objective) |
| --- | --- | ---: | --- | --- | --- | --- | --- |
| pds-20 | 33874 / 105728 / 230200 | 23821658640 | PASS / 59.107 / 23821658640 | PASS / 9.113 / 23821658640 | PASS / 1.440 / 23821658640 | PASS / 0.960 / 23821658640 | PASS / 27.309 / 23821658640 |
| pds-30 | 49944 / 154998 / 337144 | 21385445736 | PASS / 70.367 / 21385445736 | PASS / 30.301 / 21385445736 | PASS / 2.948 / 21385445736 | PASS / 2.226 / 21385445736 | PASS / 105.426 / 21385445736 |
| pds-40 | 66844 / 212859 / 462128 | 18855198824 | PASS / 122.208 / 18855198824 | PASS / 61.441 / 18855198824 | PASS / 4.596 / 18855198824 | PASS / 5.133 / 18855198824 | TIMEOUT / — / — |
| pds-50 | 83060 / 270095 / 585114 | 16603525724 | PASS / 130.494 / 16603525724 | PASS / 137.162 / 16603525724 | PASS / 5.687 / 16603525724 | PASS / 6.534 / 16603525724 | TIMEOUT / — / — |
| pds-60 | 99431 / 329643 / 712779 | 14265904407 | PASS / 124.617 / 14265904407 | PASS / 207.310 / 14265904407 | PASS / 7.527 / 14265904407 | PASS / 8.322 / 14265904407 | TIMEOUT / — / — |
| pds-70 | 114944 / 382311 / 825771 | 12241162812 | PASS / 179.846 / 12241162812 | PASS / 243.189 / 12241162812 | PASS / 8.533 / 12241162812 | PASS / 12.022 / 12241162812 | TIMEOUT / — / — |
| pds-80 | 129181 / 426278 / 919524 | 11469077462 | PASS / 251.145 / 11469077462 | TIMEOUT / — / — | PASS / 12.701 / 11469077462 | PASS / 16.419 / 11469077462 | TIMEOUT / — / — |
| pds-90 | 142823 / 466671 / 1005359 | 11087561635 | PASS / 289.316 / 11087561635 | TIMEOUT / — / — | PASS / 27.825 / 11087561635 | PASS / 18.074 / 11087561635 | TIMEOUT / — / — |
| pds-100 | 156243 / 505360 / 1086785 | 10928229968 | PASS / 239.419 / 10928229968 | TIMEOUT / — / — | PASS / 23.413 / 10928229968 | PASS / 28.121 / 10928229968 | TIMEOUT / — / — |
