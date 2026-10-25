# Netlib feasible Benchmark

[简体中文](./netlib.md) | English | [Project home](../../README_en.md) | [Benchmark summary](./README_en.md)

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

<!-- benchmark-run:netlib-accepted-auto-20261024-once -->
| Test information | Value |
| --- | --- |
| State / progress | Completed; 114/114 instances |
| Time (UTC) | 2026-10-24T05:33:23.851120+00:00 |
| Solvers in this run | lp_simplex, SCIP / SoPlex, HiGHS, Clp, GLPK |
| Case selection | netlib full suite; 25fv47, 80bau3b, adlittle, afiro, agg, agg2, agg3, bandm, beaconfd, blend, bnl1, bnl2, boeing1, boeing2, bore3d, brandy, capri, cre-a, cre-b, cre-c, cre-d, cycle, czprob, d2q06c, d6cube, degen2, degen3, dfl001, e226, etamacro, fffff800, finnis, fit1d, fit1p, fit2d, fit2p, forplan, ganges, gfrd-pnc, greenbea, greenbeb, grow15, grow22, grow7, israel, kb2, ken-07, ken-11, ken-13, ken-18, lotfi, maros, maros-r7, modszk1, nesm, osa-07, osa-14, osa-30, osa-60, pds-02, pds-06, pds-10, pds-20, perold, pilot, pilot.ja, pilot.we, pilot4, pilot87, pilotnov, qap12, qap15, qap8, recipe, sc105, sc205, sc50a, sc50b, scagr25, scagr7, scfxm1, scfxm2, scfxm3, scorpion, scrs8, scsd1, scsd6, scsd8, sctap1, sctap2, sctap3, seba, share1b, share2b, shell, ship04l, ship04s, ship08l, ship08s, ship12l, ship12s, sierra, stair, standata, standgub, standmps, stocfor1, stocfor2, stocfor3, truss, tuff, vtp.base, wood1p, woodw |
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
| Raw records | tmp/benchmarks/netlib-accepted-auto-20261024-once (configuration, environment, SHA-256, commands, logs, solutions, and audit residuals) |
| Clp | Coin LP version 1.17.11, build Mar 11 2026 |
| GLPK | GLPSOL--GLPK LP/MIP Solver 5.0 |
| Clp / GLPK algorithms | dual simplex / primal simplex; presolve on |
| GLPK skip policy | Skip remaining GLPK runs after its first timeout; Skipped runs: 64 |

| Solver | Solved instances | Timeouts | Other unsuccessful runs | Instances solved by both this solver and lp_simplex | Geometric means on common set / s (lp_simplex / solver) | Relative time (lp_simplex / solver) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| lp_simplex | 114/114 | 0 | 0 | — | — | — |
| SCIP / SoPlex | 109/114 | 4 | 1 | 109 | 0.0878237 / 0.15624 | 0.562108 |
| HiGHS | 112/114 | 2 | 0 | 112 | 0.0956086 / 0.0278547 | 3.43241 |
| Clp | 111/114 | 1 | 2 | 111 | 0.0998753 / 0.0493459 | 2.02399 |
| GLPK | 47/114 | 1 | 2 | 47 | 0.0944219 / 0.0659676 | 1.43134 |

Reproduce these cases and parameters in a new run directory:

```sh
cmake -S . -B build/benchmark -DCMAKE_BUILD_TYPE=Release -DLP_SIMPLEX_BUILD_BENCHMARKS=ON

cmake --build build/benchmark --target lp_solve_mps -j 4

python3 data/benchmark.py --dataset netlib --cases 25fv47 80bau3b adlittle afiro agg agg2 agg3 bandm beaconfd blend bnl1 bnl2 boeing1 boeing2 bore3d brandy capri cre-a cre-b cre-c cre-d cycle czprob d2q06c d6cube degen2 degen3 dfl001 e226 etamacro fffff800 finnis fit1d fit1p fit2d fit2p forplan ganges gfrd-pnc greenbea greenbeb grow15 grow22 grow7 israel kb2 ken-07 ken-11 ken-13 ken-18 lotfi maros-r7 maros modszk1 nesm osa-07 osa-14 osa-30 osa-60 pds-02 pds-06 pds-10 pds-20 perold pilot.ja pilot pilot.we pilot4 pilot87 pilotnov qap12 qap15 qap8 recipe sc105 sc205 sc50a sc50b scagr25 scagr7 scfxm1 scfxm2 scfxm3 scorpion scrs8 scsd1 scsd6 scsd8 sctap1 sctap2 sctap3 seba share1b share2b shell ship04l ship04s ship08l ship08s ship12l ship12s sierra stair standata standgub standmps stocfor1 stocfor2 stocfor3 truss tuff vtp.base wood1p woodw --binary build/benchmark/data/lp_solve_mps --scip scip --clp clp --glpk glpsol --repeats 1 --seconds 90.0 --iterations 300000 --lp-algorithm auto --scip-algorithm default --highs-algorithm choose --solvers lp_simplex scip highs clp glpk
```
<!-- /benchmark-run -->

## Latest per-instance results

| Instance | Size (rows / columns / nonzeros) | Reference objective | lp_simplex (status / s / objective) | SCIP / SoPlex (status / s / objective) | HiGHS (status / s / objective) | Clp (status / s / objective) | GLPK (status / s / objective) |
| --- | --- | ---: | --- | --- | --- | --- | --- |
| 25fv47 | 821 / 1571 / 10400 | 5501.845888 | PASS / 0.892 / 5501.84588829 | PASS / 0.416 / 5501.84588829 | PASS / 0.114 / 5501.84588829 | PASS / 0.122 / 5501.84588829 | PASS / 0.114 / 5501.84588829 |
| 80bau3b | 2262 / 9799 / 21002 | 987224.1924 | PASS / 0.347 / 987224.192409 | PASS / 0.457 / 987224.192409 | PASS / 0.104 / 987224.192409 | PASS / 0.171 / 987224.192409 | PASS / 0.474 / 987224.192409 |
| adlittle | 56 / 97 / 383 | 225494.9632 | PASS / 0.015 / 225494.963162 | PASS / 0.040 / 225494.963162 | PASS / 0.002 / 225494.963162 | PASS / 0.013 / 225494.963162 | PASS / 0.008 / 225494.963162 |
| afiro | 27 / 32 / 83 | -464.7531429 | PASS / 0.016 / -464.753142857 | PASS / 0.039 / -464.753142857 | PASS / 0.001 / -464.753142857 | PASS / 0.013 / -464.753142857 | PASS / 0.007 / -464.753142857 |
| agg | 488 / 163 / 2410 | -35991767.29 | PASS / 0.017 / -35991767.2866 | PASS / 0.045 / -35991767.2866 | PASS / 0.004 / -35991767.2866 | PASS / 0.016 / -35991767.2866 | PASS / 0.010 / -35991767.2866 |
| agg2 | 516 / 302 / 4284 | -20239252.36 | PASS / 0.018 / -20239252.356 | PASS / 0.047 / -20239252.356 | PASS / 0.006 / -20239252.356 | PASS / 0.015 / -20239252.356 | PASS / 0.013 / -20239252.356 |
| agg3 | 516 / 302 / 4300 | 10312115.94 | PASS / 0.019 / 10312115.9351 | PASS / 0.049 / 10312115.9351 | PASS / 0.006 / 10312115.9351 | PASS / 0.021 / 10312115.9351 | PASS / 0.013 / 10312115.9351 |
| bandm | 305 / 472 / 2494 | -158.6280185 | PASS / 0.022 / -158.62801845 | PASS / 0.054 / -158.62801845 | PASS / 0.008 / -158.62801845 | PASS / 0.017 / -158.62801845 | PASS / 0.013 / -158.62801845 |
| beaconfd | 173 / 262 / 3375 | 33592.48581 | PASS / 0.015 / 33592.4858072 | PASS / 0.039 / 33592.4858072 | PASS / 0.004 / 33592.4858072 | PASS / 0.013 / 33592.4858072 | PASS / 0.010 / 33592.4858072 |
| blend | 74 / 83 / 491 | -30.81214985 | PASS / 0.015 / -30.8121498458 | PASS / 0.041 / -30.8121498458 | PASS / 0.002 / -30.8121498458 | PASS / 0.013 / -30.8121498459 | PASS / 0.007 / -30.8121498458 |
| bnl1 | 643 / 1175 / 5121 | 1977.629562 | PASS / 0.069 / 1977.62956152 | PASS / 0.111 / 1977.62956152 | PASS / 0.021 / 1977.62956152 | PASS / 0.033 / 1977.62956152 | PASS / 0.032 / 1977.62956152 |
| bnl2 | 2324 / 3489 / 13999 | 1811.23654 | PASS / 0.103 / 1811.23654036 | PASS / 0.214 / 1811.23654036 | PASS / 0.044 / 1811.23654036 | PASS / 0.044 / 1811.23654036 | PASS / 0.210 / 1811.23654036 |
| boeing1 | 351 / 384 / 3485 | -335.2135675 | PASS / 0.026 / -335.213567507 | PASS / 0.062 / -335.213567507 | PASS / 0.074 / -335.213567507 | PASS / 0.021 / -335.213567507 | PASS / 0.017 / -335.213567507 |
| boeing2 | 166 / 143 / 1196 | -315.018728 | PASS / 0.015 / -315.018728015 | PASS / 0.042 / -315.018728015 | PASS / 0.003 / -315.018728015 | PASS / 0.016 / -315.018728015 | PASS / 0.010 / -315.018728015 |
| bore3d | 233 / 315 / 1429 | 1373.080394 | PASS / 0.016 / 1373.08039421 | PASS / 0.043 / 1373.08039421 | PASS / 0.002 / 1373.08039421 | PASS / 0.014 / 1373.08039421 | PASS / 0.008 / 1373.08039421 |
| brandy | 220 / 249 / 2148 | 1518.509896 | PASS / 0.018 / 1518.50989649 | PASS / 0.048 / 1518.50989649 | PASS / 0.005 / 1518.50989649 | PASS / 0.015 / 1518.50989649 | PASS / 0.011 / 1518.50989649 |
| capri | 271 / 353 / 1767 | 2690.012914 | PASS / 0.020 / 2690.01291377 | PASS / 0.047 / 2690.01291377 | PASS / 0.004 / 2690.01291377 | PASS / 0.018 / 2690.01291377 | PASS / 0.014 / 2690.01291377 |
| cre-a | 3516 / 4067 / 14987 | 23595407.06 | PASS / 0.195 / 23595407.061 | PASS / 0.291 / 23595407.061 | PASS / 0.053 / 23595407.061 | PASS / 0.050 / 23595407.061 | PASS / 0.240 / 23595407.061 |
| cre-b | 9648 / 72447 / 256095 | 23129639.89 | PASS / 11.693 / 23129639.8868 | PASS / 3.259 / 23129639.8868 | PASS / 0.995 / 23129639.8868 | PASS / 1.028 / 23129639.8868 | PASS / 6.065 / 23129639.8868 |
| cre-c | 3068 / 3678 / 13244 | 25275116.14 | PASS / 0.156 / 25275116.1409 | PASS / 0.189 / 25275116.1409 | PASS / 0.047 / 25275116.1409 | PASS / 0.050 / 25275116.1409 | PASS / 0.172 / 25275116.1409 |
| cre-d | 8926 / 69980 / 242646 | 24454969.76 | PASS / 3.319 / 24454969.7645 | PASS / 1.602 / 24454969.7645 | PASS / 0.667 / 24454969.7645 | PASS / 0.416 / 24454969.7645 | PASS / 3.236 / 24454969.7645 |
| cycle | 1903 / 2857 / 20720 | -5.226393025 | PASS / 0.056 / -5.22639302489 | PASS / 0.135 / -5.22639302489 | PASS / 0.140 / -5.22639302489 | PASS / 0.034 / -5.22639302489 | PASS / 0.047 / -5.22639302489 |
| czprob | 929 / 3523 / 10669 | 2185196.699 | PASS / 0.043 / 2185196.69886 | PASS / 0.099 / 2185196.69886 | PASS / 0.028 / 2185196.69886 | PASS / 0.023 / 2185196.69886 | PASS / 0.042 / 2185196.69886 |
| d2q06c | 2171 / 5167 / 32417 | 122784.2108 | PASS / 1.483 / 122784.210814 | PASS / 2.111 / 122784.210814 | PASS / 0.551 / 122784.210814 | PASS / 0.799 / 122784.210814 | PASS / 0.960 / 122784.210814 |
| d6cube | 415 / 6184 / 37704 | 315.4916667 | PASS / 3.058 / 315.491666664 | PASS / 0.287 / 315.491666667 | PASS / 0.108 / 315.491666667 | PASS / 0.064 / 315.491666667 | PASS / 1.148 / 315.491666667 |
| degen2 | 444 / 534 / 3978 | -1435.178 | PASS / 0.033 / -1435.178 | PASS / 0.070 / -1435.178 | PASS / 0.012 / -1435.178 | PASS / 0.021 / -1435.178 | PASS / 0.022 / -1435.178 |
| degen3 | 1503 / 1818 / 24646 | -987.294 | PASS / 0.337 / -987.294 | PASS / 0.658 / -987.294 | PASS / 0.135 / -987.294 | PASS / 0.342 / -987.294 | PASS / 0.244 / -987.294 |
| dfl001 | 6071 / 12230 / 35632 | 11266396.05 | PASS / 46.897 / 11266396.0467 | PASS / 10.815 / 11266396.0467 | PASS / 5.056 / 11266396.0467 | PASS / 4.190 / 11266396.0467 | PASS / 16.098 / 11266396.0467 |
| e226 | 223 / 282 / 2578 | -11.63892907 | PASS / 0.022 / -11.6389290664 | PASS / 0.055 / -11.6389290664 | PASS / 0.007 / -11.6389290664 | PASS / 0.022 / -11.6389290664 | PASS / 0.015 / -11.6389290664 |
| etamacro | 400 / 688 / 2409 | -755.7152334 | PASS / 0.028 / -755.715233352 | PASS / 0.060 / -755.715233352 | PASS / 0.009 / -755.715233301 | PASS / 0.018 / -755.715233339 | PASS / 0.016 / -755.715233139 |
| fffff800 | 524 / 854 / 6227 | 555679.5648 | PASS / 0.024 / 555679.564817 | PASS / 0.116 / 555679.564817 | PASS / 0.056 / 555679.564817 | PASS / 0.018 / 555679.564817 | PASS / 0.015 / 555679.564817 |
| finnis | 497 / 614 / 2310 | 172791.0656 | PASS / 0.018 / 172791.065596 | PASS / 0.051 / 172791.065596 | PASS / 0.005 / 172791.065596 | PASS / 0.017 / 172791.065596 | PASS / 0.014 / 172791.065596 |
| fit1d | 24 / 1026 / 13404 | -9146.378092 | PASS / 0.024 / -9146.37809242 | PASS / 0.110 / -9146.37809242 | PASS / 0.014 / -9146.37809242 | PASS / 0.018 / -9146.37809242 | PASS / 0.032 / -9146.37809242 |
| fit1p | 627 / 1677 / 9868 | 9146.378092 | PASS / 0.071 / 9146.37809242 | PASS / 0.162 / 9146.37809242 | PASS / 0.032 / 9146.37809242 | PASS / 0.034 / 9146.37809242 | PASS / 0.029 / 9146.37809242 |
| fit2d | 25 / 10500 / 129018 | -68464.29329 | PASS / 0.189 / -68464.2932938 | PASS / 0.587 / -68464.2932938 | PASS / 0.167 / -68464.2932938 | PASS / 0.093 / -68464.2932938 | PASS / 1.736 / -68464.2932938 |
| fit2p | 3000 / 13525 / 50284 | 68464.29329 | PASS / 3.075 / 68464.2932938 | PASS / 1.045 / 68464.2932938 | PASS / 0.692 / 68464.2932938 | PASS / 1.084 / 68464.2932938 | PASS / 1.050 / 68464.2932938 |
| forplan | 161 / 421 / 4563 | -664.218961272 | PASS / 0.022 / -664.218961272 | PASS / 0.057 / -664.218961272 | PASS / 0.008 / -664.218961272 | PASS / 0.016 / -664.218961272 | PASS / 0.014 / -664.218961272 |
| ganges | 1309 / 1681 / 6912 | -109585.7361 | PASS / 0.031 / -109585.736129 | PASS / 0.078 / -109585.736129 | PASS / 0.014 / -109585.736129 | PASS / 0.021 / -109585.736129 | PASS / 0.032 / -109585.736129 |
| gfrd-pnc | 616 / 1092 / 2377 | 6902236 | PASS / 0.020 / 6902235.99955 | PASS / 0.053 / 6902235.99955 | PASS / 0.008 / 6902235.99955 | PASS / 0.018 / 6902235.99955 | PASS / 0.018 / 6902235.99955 |
| greenbea | 2392 / 5405 / 30877 | -72555248.13 | PASS / 0.604 / -72555248.1298 | AUDIT_FAILED / — / -72555248.1298 | PASS / 0.198 / -72555248.1298 | PASS / 0.429 / -72555248.1299 | AUDIT_FAILED / — / -72555248.1298 |
| greenbeb | 2392 / 5405 / 30877 | -4302260.261 | PASS / 2.063 / -4302260.26121 | PASS / 1.127 / -4302260.26121 | PASS / 0.346 / -4302260.26121 | PASS / 0.268 / -4302260.26121 | PASS / 0.238 / -4302260.26121 |
| grow15 | 300 / 645 / 5620 | -106870941.3 | PASS / 0.070 / -106870941.294 | PASS / 0.097 / -106870941.294 | PASS / 0.028 / -106870941.294 | PASS / 0.022 / -106870941.294 | PASS / 0.023 / -106870941.294 |
| grow22 | 440 / 946 / 8252 | -160834336.5 | PASS / 0.163 / -160834336.483 | PASS / 0.182 / -160834336.483 | PASS / 0.069 / -160834336.483 | PASS / 0.032 / -160834336.483 | PASS / 0.037 / -160834336.483 |
| grow7 | 140 / 301 / 2612 | -47787811.81 | PASS / 0.028 / -47787811.8147 | PASS / 0.057 / -47787811.8147 | PASS / 0.008 / -47787811.8147 | PASS / 0.050 / -47787811.8147 | AUDIT_FAILED / — / -47787811.8147 |
| israel | 174 / 142 / 2269 | -896644.8219 | PASS / 0.018 / -896644.821863 | PASS / 0.051 / -896644.821863 | PASS / 0.004 / -896644.821863 | PASS / 0.016 / -896644.821863 | PASS / 0.011 / -896644.821863 |
| kb2 | 43 / 41 / 286 | -1749.90013 | PASS / 0.015 / -1749.90012991 | PASS / 0.037 / -1749.90012991 | PASS / 0.001 / -1749.90012991 | PASS / 0.013 / -1749.90012991 | PASS / 0.009 / -1749.90012991 |
| ken-07 | 2426 / 3602 / 8404 | -679520443.4 | PASS / 0.049 / -679520443.382 | PASS / 0.102 / -679520443.382 | PASS / 0.027 / -679520443.382 | PASS / 0.026 / -679520443.382 | PASS / 0.060 / -679520443.382 |
| ken-11 | 14694 / 21349 / 49058 | -6972382263 | PASS / 0.832 / -6972382262.52 | PASS / 0.870 / -6972382262.52 | PASS / 0.216 / -6972382262.52 | PASS / 0.139 / -6972382262.52 | PASS / 2.215 / -6972382262.52 |
| ken-13 | 28632 / 42659 / 97246 | -10257394789 | PASS / 3.521 / -10257394789.5 | PASS / 3.788 / -10257394789.5 | PASS / 0.638 / -10257394789.5 | PASS / 0.327 / -10257394789.5 | PASS / 16.654 / -10257394789.5 |
| ken-18 | 105127 / 154699 / 358171 | -52217025287 | PASS / 55.307 / -52217025287.4 | PASS / 42.327 / -52217025287.4 | PASS / 2.808 / -52217025287.4 | PASS / 1.767 / -52217025287.4 | TIMEOUT / — / — |
| lotfi | 153 / 308 / 1078 | -25.26470606 | PASS / 0.018 / -25.2647060619 | PASS / 0.046 / -25.2647060619 | PASS / 0.003 / -25.2647060619 | PASS / 0.015 / -25.2647060619 | SKIPPED / — / — |
| maros | 846 / 1443 / 9614 | -58063.7437 | PASS / 0.082 / -58063.7437011 | PASS / 0.103 / -58063.7437011 | PASS / 0.034 / -58063.7437011 | PASS / 0.044 / -58063.7437011 | SKIPPED / — / — |
| maros-r7 | 3136 / 9408 / 144848 | 1497185.166 | PASS / 1.340 / 1497185.16648 | PASS / 0.820 / 1497185.16648 | PASS / 0.480 / 1497185.16648 | PASS / 0.860 / 1497185.16648 | SKIPPED / — / — |
| modszk1 | 687 / 1620 / 3168 | 320.6197291 | PASS / 0.046 / 320.619729065 | PASS / 0.086 / 320.619729064 | PASS / 0.013 / 320.619729065 | PASS / 0.048 / 320.619729064 | SKIPPED / — / — |
| nesm | 662 / 2923 / 13288 | 14076036.49 | PASS / 0.143 / 14076036.4876 | PASS / 0.343 / 14076036.4876 | PASS / 0.125 / 14076036.4876 | PASS / 0.052 / 14076036.4876 | SKIPPED / — / — |
| osa-07 | 1118 / 23949 / 143694 | 535722.5173 | PASS / 0.245 / 535722.517299 | PASS / 2.237 / 535722.517299 | PASS / 0.255 / 535722.517299 | PASS / 0.080 / 535722.517299 | SKIPPED / — / — |
| osa-14 | 2337 / 52460 / 314760 | 1106462.845 | PASS / 0.561 / 1106462.84474 | PASS / 28.853 / 1106462.84474 | PASS / 0.821 / 1106462.84474 | PASS / 0.198 / 1106462.84474 | SKIPPED / — / — |
| osa-30 | 4350 / 100024 / 600138 | 2142139.873 | PASS / 1.836 / 2142139.87321 | TIMEOUT / — / — | PASS / 3.323 / 2142139.87321 | PASS / 0.373 / 2142139.87321 | SKIPPED / — / — |
| osa-60 | 10280 / 232966 / 1397793 | 4044072.503 | PASS / 8.267 / 4044072.50316 | TIMEOUT / — / — | PASS / 14.653 / 4044072.50316 | PASS / 1.068 / 4044072.50316 | SKIPPED / — / — |
| pds-02 | 2953 / 7535 / 16390 | 28857862010 | PASS / 0.101 / 28857862010 | PASS / 0.151 / 28857862010 | PASS / 0.077 / 28857862010 | PASS / 0.043 / 28857862010 | SKIPPED / — / — |
| pds-06 | 9881 / 28655 / 62524 | 27761037600 | PASS / 1.073 / 27761037600 | PASS / 0.749 / 27761037600 | PASS / 0.313 / 27761037600 | PASS / 0.131 / 27761037600 | SKIPPED / — / — |
| pds-10 | 16558 / 48763 / 106436 | 26727094976 | PASS / 5.080 / 26727094976 | PASS / 2.178 / 26727094976 | PASS / 0.482 / 26727094976 | PASS / 0.277 / 26727094976 | SKIPPED / — / — |
| pds-20 | 33874 / 105728 / 230200 | 23821658640 | PASS / 35.770 / 23821658640 | PASS / 8.644 / 23821658640 | PASS / 1.612 / 23821658640 | PASS / 0.958 / 23821658640 | SKIPPED / — / — |
| perold | 625 / 1376 / 6018 | -9380.755278 | PASS / 0.109 / -9380.75527824 | PASS / 0.274 / -9380.75527824 | PASS / 0.047 / -9380.75527824 | PASS / 0.084 / -9380.75527824 | SKIPPED / — / — |
| pilot | 1441 / 3652 / 43167 | -557.4897293 | PASS / 2.313 / -557.489729284 | PASS / 1.684 / -557.489729285 | PASS / 0.829 / -557.489729291 | PASS / 0.648 / -557.489729284 | SKIPPED / — / — |
| pilot.ja | 940 / 1988 / 14698 | -6113.136466 | PASS / 0.179 / -6113.13646558 | PASS / 0.253 / -6113.13646558 | PASS / 0.066 / -6113.13646558 | PASS / 0.158 / -6113.13646558 | SKIPPED / — / — |
| pilot.we | 722 / 2789 / 9126 | -2720107.533 | PASS / 0.195 / -2720107.53284 | PASS / 0.431 / -2720107.53284 | PASS / 0.100 / -2720107.53284 | AUDIT_FAILED / — / -2720107.53927 | SKIPPED / — / — |
| pilot4 | 410 / 1000 / 5141 | -2581.13926 | PASS / 0.043 / -2581.13925888 | PASS / 0.096 / -2581.13925888 | PASS / 0.024 / -2581.13925888 | PASS / 0.038 / -2581.13925888 | SKIPPED / — / — |
| pilot87 | 2030 / 4883 / 73152 | 301.7103473 | PASS / 19.856 / 301.710347789 | PASS / 4.816 / 301.710347333 | PASS / 2.902 / 301.710347333 | PASS / 2.657 / 301.710347333 | SKIPPED / — / — |
| pilotnov | 975 / 2172 / 13057 | -4497.276188 | PASS / 0.117 / -4497.27618822 | PASS / 0.216 / -4497.27618822 | PASS / 0.113 / -4497.27618822 | AUDIT_FAILED / — / -4497.27618822 | SKIPPED / — / — |
| qap12 | 3192 / 8856 / 38304 | 522.8943506 | PASS / 30.436 / 522.894342332 | TIMEOUT / — / — | TIMEOUT / — / — | PASS / 28.297 / 522.894350559 | SKIPPED / — / — |
| qap15 | 6330 / 22275 / 94950 | 1040.994041 | PASS / 56.071 / 1040.9940074 | TIMEOUT / — / — | TIMEOUT / — / — | TIMEOUT / — / — | SKIPPED / — / — |
| qap8 | 912 / 1632 / 7296 | 203.5 | PASS / 4.918 / 203.5 | PASS / 0.996 / 203.5 | PASS / 0.385 / 203.5 | PASS / 0.313 / 203.5 | SKIPPED / — / — |
| recipe | 91 / 180 / 663 | -266.616 | PASS / 0.014 / -266.616 | PASS / 0.038 / -266.616 | PASS / 0.002 / -266.616 | PASS / 0.014 / -266.616 | SKIPPED / — / — |
| sc105 | 105 / 103 / 280 | -52.20206121 | PASS / 0.015 / -52.2020612117 | PASS / 0.036 / -52.2020612117 | PASS / 0.001 / -52.2020612117 | PASS / 0.013 / -52.2020612117 | SKIPPED / — / — |
| sc205 | 205 / 203 / 551 | -52.20206121 | PASS / 0.016 / -52.2020612117 | PASS / 0.041 / -52.2020612117 | PASS / 0.002 / -52.2020612117 | PASS / 0.014 / -52.2020612117 | SKIPPED / — / — |
| sc50a | 50 / 48 / 130 | -64.57507706 | PASS / 0.013 / -64.5750770586 | PASS / 0.043 / -64.5750770586 | PASS / 0.001 / -64.5750770586 | PASS / 0.012 / -64.5750770586 | SKIPPED / — / — |
| sc50b | 50 / 48 / 118 | -70 | PASS / 0.016 / -70 | PASS / 0.037 / -70 | PASS / 0.001 / -70 | PASS / 0.014 / -70 | SKIPPED / — / — |
| scagr25 | 471 / 500 / 1554 | -14753433.06 | PASS / 0.021 / -14753433.0608 | PASS / 0.047 / -14753433.0608 | PASS / 0.006 / -14753433.0608 | PASS / 0.015 / -14753433.0608 | SKIPPED / — / — |
| scagr7 | 129 / 140 / 420 | -2331389.824 | PASS / 0.014 / -2331389.82433 | PASS / 0.046 / -2331389.82433 | PASS / 0.002 / -2331389.82433 | PASS / 0.015 / -2331389.82433 | SKIPPED / — / — |
| scfxm1 | 330 / 457 / 2589 | 18416.75903 | PASS / 0.022 / 18416.7590283 | PASS / 0.053 / 18416.7590283 | PASS / 0.007 / 18416.7590283 | PASS / 0.017 / 18416.7590283 | SKIPPED / — / — |
| scfxm2 | 660 / 914 / 5183 | 36660.26156 | PASS / 0.037 / 36660.261565 | PASS / 0.074 / 36660.261565 | PASS / 0.015 / 36660.261565 | PASS / 0.026 / 36660.261565 | SKIPPED / — / — |
| scfxm3 | 990 / 1371 / 7777 | 54901.25455 | PASS / 0.061 / 54901.2545498 | PASS / 0.097 / 54901.2545498 | PASS / 0.024 / 54901.2545498 | PASS / 0.034 / 54901.2545497 | SKIPPED / — / — |
| scorpion | 388 / 358 / 1426 | 1878.124823 | PASS / 0.016 / 1878.12482274 | PASS / 0.041 / 1878.12482274 | PASS / 0.003 / 1878.12482274 | PASS / 0.014 / 1878.12482274 | SKIPPED / — / — |
| scrs8 | 490 / 1169 / 3182 | 904.2969538 | PASS / 0.024 / 904.296953801 | PASS / 0.073 / 904.296953801 | PASS / 0.008 / 904.296953801 | PASS / 0.018 / 904.296953799 | SKIPPED / — / — |
| scsd1 | 77 / 760 / 2388 | 8.666666674 | PASS / 0.018 / 8.66666667433 | PASS / 0.047 / 8.66666667433 | PASS / 0.004 / 8.66666667433 | PASS / 0.015 / 8.66666667433 | SKIPPED / — / — |
| scsd6 | 147 / 1350 / 4316 | 50.50000008 | PASS / 0.025 / 50.5000000771 | PASS / 0.067 / 50.500000078 | PASS / 0.008 / 50.5000000783 | PASS / 0.019 / 50.5000000783 | SKIPPED / — / — |
| scsd8 | 397 / 2750 / 8584 | 904.9999999 | PASS / 0.073 / 904.999999925 | PASS / 0.116 / 904.999999925 | PASS / 0.038 / 904.999999925 | PASS / 0.037 / 904.999999925 | SKIPPED / — / — |
| sctap1 | 300 / 480 / 1692 | 1412.25 | PASS / 0.018 / 1412.25 | PASS / 0.052 / 1412.25 | PASS / 0.005 / 1412.25 | PASS / 0.014 / 1412.25 | SKIPPED / — / — |
| sctap2 | 1090 / 1880 / 6714 | 1724.807143 | PASS / 0.028 / 1724.80714286 | PASS / 0.075 / 1724.80714286 | PASS / 0.014 / 1724.80714286 | PASS / 0.020 / 1724.80714286 | SKIPPED / — / — |
| sctap3 | 1480 / 2480 / 8874 | 1424 | PASS / 0.039 / 1424 | PASS / 0.094 / 1424 | PASS / 0.017 / 1424 | PASS / 0.021 / 1424 | SKIPPED / — / — |
| seba | 515 / 1028 / 4352 | 15711.6 | PASS / 0.021 / 15711.6 | PASS / 0.047 / 15711.6 | PASS / 0.006 / 15711.6 | PASS / 0.014 / 15711.6 | SKIPPED / — / — |
| share1b | 117 / 225 / 1151 | -76589.31858 | PASS / 0.017 / -76589.3185792 | PASS / 0.044 / -76589.3185792 | PASS / 0.003 / -76589.3185792 | PASS / 0.014 / -76589.3185792 | SKIPPED / — / — |
| share2b | 96 / 79 / 694 | -415.7322407 | PASS / 0.016 / -415.732240741 | PASS / 0.038 / -415.732240741 | PASS / 0.002 / -415.732240741 | PASS / 0.015 / -415.732240741 | SKIPPED / — / — |
| shell | 536 / 1775 / 3556 | 1208825346 | PASS / 0.023 / 1208825346 | PASS / 0.058 / 1208825346 | PASS / 0.012 / 1208825346 | PASS / 0.020 / 1208825346 | SKIPPED / — / — |
| ship04l | 402 / 2118 / 6332 | 1793324.538 | PASS / 0.025 / 1793324.53797 | PASS / 0.067 / 1793324.53797 | PASS / 0.009 / 1793324.53797 | PASS / 0.017 / 1793324.53797 | SKIPPED / — / — |
| ship04s | 402 / 1458 / 4352 | 1798714.7 | PASS / 0.020 / 1798714.70045 | PASS / 0.054 / 1798714.70045 | PASS / 0.007 / 1798714.70045 | PASS / 0.016 / 1798714.70045 | SKIPPED / — / — |
| ship08l | 778 / 4283 / 12802 | 1909055.211 | PASS / 0.040 / 1909055.21139 | PASS / 0.089 / 1909055.21139 | PASS / 0.020 / 1909055.21139 | PASS / 0.025 / 1909055.21139 | SKIPPED / — / — |
| ship08s | 778 / 2387 / 7114 | 1920098.211 | PASS / 0.025 / 1920098.21053 | PASS / 0.064 / 1920098.21053 | PASS / 0.011 / 1920098.21053 | PASS / 0.018 / 1920098.21053 | SKIPPED / — / — |
| ship12l | 1151 / 5427 / 16170 | 1470187.919 | PASS / 0.057 / 1470187.91933 | PASS / 0.118 / 1470187.91933 | PASS / 0.029 / 1470187.91933 | PASS / 0.136 / 1470187.91933 | SKIPPED / — / — |
| ship12s | 1151 / 2763 / 8178 | 1489236.134 | PASS / 0.025 / 1489236.13441 | PASS / 0.074 / 1489236.13441 | PASS / 0.013 / 1489236.13441 | PASS / 0.024 / 1489236.13441 | SKIPPED / — / — |
| sierra | 1227 / 2036 / 7302 | 15394362.18 | PASS / 0.031 / 15394362.1836 | PASS / 0.076 / 15394362.1836 | PASS / 0.016 / 15394362.1836 | PASS / 0.020 / 15394362.1836 | SKIPPED / — / — |
| stair | 356 / 467 / 3856 | -251.2669512 | PASS / 0.033 / -251.266951193 | PASS / 0.072 / -251.266951193 | PASS / 0.013 / -251.266951193 | PASS / 0.022 / -251.266951193 | SKIPPED / — / — |
| standata | 359 / 1075 / 3031 | 1257.6995 | PASS / 0.016 / 1257.6995 | PASS / 0.048 / 1257.6995 | PASS / 0.005 / 1257.6995 | PASS / 0.016 / 1257.6995 | SKIPPED / — / — |
| standgub | 361 / 1184 / 3139 | 1257.6995 | PASS / 0.018 / 1257.6995 | PASS / 0.044 / 1257.6995 | PASS / 0.005 / 1257.6995 | PASS / 0.016 / 1257.6995 | SKIPPED / — / — |
| standmps | 467 / 1075 / 3679 | 1406.0175 | PASS / 0.017 / 1406.0175 | PASS / 0.053 / 1406.0175 | PASS / 0.006 / 1406.0175 | PASS / 0.016 / 1406.0175 | SKIPPED / — / — |
| stocfor1 | 117 / 111 / 447 | -41131.97622 | PASS / 0.027 / -41131.9762194 | PASS / 0.041 / -41131.9762194 | PASS / 0.002 / -41131.9762194 | PASS / 0.014 / -41131.9762194 | SKIPPED / — / — |
| stocfor2 | 2157 / 2031 / 8343 | -39024.40854 | PASS / 0.107 / -39024.4085379 | PASS / 0.113 / -39024.4085379 | PASS / 0.027 / -39024.4085379 | PASS / 0.030 / -39024.4085379 | SKIPPED / — / — |
| stocfor3 | 16675 / 15695 / 64875 | -39976.78394 | PASS / 4.939 / -39976.7839436 | PASS / 5.138 / -39976.7839436 | PASS / 0.382 / -39976.7839436 | PASS / 0.327 / -39976.7839436 | SKIPPED / — / — |
| truss | 1000 / 8806 / 27836 | 458815.8472 | PASS / 16.510 / 458815.847186 | PASS / 1.825 / 458815.847186 | PASS / 1.940 / 458815.847186 | PASS / 1.655 / 458815.847186 | SKIPPED / — / — |
| tuff | 333 / 587 / 4520 | 0.292147765 | PASS / 0.021 / 0.292147765094 | PASS / 0.071 / 0.292147765094 | PASS / 0.008 / 0.292147765094 | PASS / 0.018 / 0.292147765094 | SKIPPED / — / — |
| vtp.base | 198 / 203 / 908 | 129831.4625 | PASS / 0.015 / 129831.462461 | PASS / 0.037 / 129831.462461 | PASS / 0.002 / 129831.462461 | PASS / 0.013 / 129831.462461 | SKIPPED / — / — |
| wood1p | 244 / 2594 / 70215 | 1.442902412 | PASS / 0.069 / 1.44290241157 | PASS / 0.206 / 1.44290241157 | PASS / 0.086 / 1.44290241157 | PASS / 0.042 / 1.44290241157 | SKIPPED / — / — |
| woodw | 1098 / 8405 / 37474 | 1.304476333 | PASS / 0.076 / 1.30447633308 | PASS / 0.174 / 1.30447633308 | PASS / 0.066 / 1.30447633308 | PASS / 0.040 / 1.30447633308 | SKIPPED / — / — |
