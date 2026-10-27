> Historical results: 2026-10-26 (Asia/Shanghai), before the AUTO review. These runs include failed or timed-out cases and do not validate the current optimization. Original tables and UTC timestamps are retained. Only Markdown reports are archived; raw run directories mentioned below have been removed. [Current reports](../../README_en.md) · [History](../README_en.md)

# Netlib feasible Benchmark

[简体中文](./netlib.md) | English | [Project home](../../../../README_en.md) | [Benchmark summary](./README_en.md)

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

<!-- benchmark-run:netlib-glpk-completed -->
| Test information | Value |
| --- | --- |
| State / progress | Completed; 114/114 instances |
| Time (UTC) | 2026-10-25T21:49:44.718550+00:00 |
| Solvers in this run | lp_simplex, SCIP / SoPlex, HiGHS, Clp, GLPK |
| Case selection | netlib full suite; 25fv47, 80bau3b, adlittle, afiro, agg, agg2, agg3, bandm, beaconfd, blend, bnl1, bnl2, boeing1, boeing2, bore3d, brandy, capri, cre-a, cre-b, cre-c, cre-d, cycle, czprob, d2q06c, d6cube, degen2, degen3, dfl001, e226, etamacro, fffff800, finnis, fit1d, fit1p, fit2d, fit2p, forplan, ganges, gfrd-pnc, greenbea, greenbeb, grow15, grow22, grow7, israel, kb2, ken-07, ken-11, ken-13, ken-18, lotfi, maros, maros-r7, modszk1, nesm, osa-07, osa-14, osa-30, osa-60, pds-02, pds-06, pds-10, pds-20, perold, pilot, pilot.ja, pilot.we, pilot4, pilot87, pilotnov, qap12, qap15, qap8, recipe, sc105, sc205, sc50a, sc50b, scagr25, scagr7, scfxm1, scfxm2, scfxm3, scorpion, scrs8, scsd1, scsd6, scsd8, sctap1, sctap2, sctap3, seba, share1b, share2b, shell, ship04l, ship04s, ship08l, ship08s, ship12l, ship12s, sierra, stair, standata, standgub, standmps, stocfor1, stocfor2, stocfor3, truss, tuff, vtp.base, wood1p, woodw |
| System / CPU | macOS-27.0.1-arm64-arm-64bit-Mach-O / Apple M1 |
| Memory / logical CPUs | 17179869184 bytes / 8 |
| Python / NumPy / SciPy | 3.14.2 / 2.5.3 / 1.18.1 |
| lp_simplex version | 4115c34c078d6942394cc225cc33645e18006d8b; working tree modified |
| Build / compiler | Release / Apple clang version 21.0.0 (clang-2100.3.34.2) |
| SCIP / LP backend | SCIP version 10.0.3 [precision: 8 byte] [memory: block] [mode: optimized] [LP solver: SoPlex 8.0.2] [GitHash: d409edf9f6] |
| HiGHS | 1.12.0 / scipy-private-binding |
| Algorithms (lp / scip / highs) | auto / default / choose |
| Repeats / time limit / tolerance / work budget | 1 / 90 s / 1e-07 / 300000 |
| Resources | Serial; one thread per solver; no memory cap or exclusive machine access |
| Raw records | tmp/benchmarks/2026-10-26/netlib-glpk-completed (configuration, environment, SHA-256, commands, logs, solutions, and audit residuals) |
| Clp | Coin LP version 1.17.11, build Mar 11 2026 |
| GLPK | GLPSOL--GLPK LP/MIP Solver 5.0 |
| Supplemental run | Composite report: tmp/benchmarks/2026-10-26/netlib (2026-10-25T21:49:44.718550+00:00, 506 executed records retained); tmp/benchmarks/2026-10-26/netlib-glpk-supplement (2026-10-26T04:54:09.150074+00:00, 64 GLPK records replacing skips). Separate runs on the same host, with identical GLPK/HiGHS-reader binaries and budgets; only benchmark skip handling changed. Original runs are preserved; per-record source_run and component metadata retain provenance. |
| Clp / GLPK algorithms | dual simplex / primal simplex; presolve on |

| Solver | Solved instances | Timeouts | Other unsuccessful runs | Instances solved by both this solver and lp_simplex | Geometric means on common set / s (lp_simplex / solver) | Relative time (lp_simplex / solver) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| lp_simplex | 112/114 | 1 | 1 | — | — | — |
| SCIP / SoPlex | 109/114 | 4 | 1 | 107 | 0.087884 / 0.163895 | 0.536223 |
| HiGHS | 112/114 | 2 | 0 | 110 | 0.0961808 / 0.0287287 | 3.3479 |
| Clp | 111/114 | 1 | 2 | 109 | 0.10042 / 0.0517313 | 1.9412 |
| GLPK | 110/114 | 2 | 2 | 108 | 0.0953151 / 0.0797174 | 1.19566 |

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
| 25fv47 | 821 / 1571 / 10400 | 5501.845888 | PASS / 0.298 / 5501.84588829 | PASS / 0.457 / 5501.84588829 | PASS / 0.115 / 5501.84588829 | PASS / 0.124 / 5501.84588829 | PASS / 0.130 / 5501.84588829 |
| 80bau3b | 2262 / 9799 / 21002 | 987224.1924 | PASS / 0.353 / 987224.192409 | PASS / 0.433 / 987224.192409 | PASS / 0.109 / 987224.192409 | PASS / 0.173 / 987224.192409 | PASS / 0.486 / 987224.192409 |
| adlittle | 56 / 97 / 383 | 225494.9632 | PASS / 0.017 / 225494.963162 | PASS / 0.042 / 225494.963162 | PASS / 0.003 / 225494.963162 | PASS / 0.016 / 225494.963162 | PASS / 0.008 / 225494.963162 |
| afiro | 27 / 32 / 83 | -464.7531429 | PASS / 0.016 / -464.753142857 | PASS / 0.038 / -464.753142857 | PASS / 0.001 / -464.753142857 | PASS / 0.015 / -464.753142857 | PASS / 0.009 / -464.753142857 |
| agg | 488 / 163 / 2410 | -35991767.29 | PASS / 0.020 / -35991767.2866 | PASS / 0.045 / -35991767.2866 | PASS / 0.005 / -35991767.2866 | PASS / 0.016 / -35991767.2866 | PASS / 0.011 / -35991767.2866 |
| agg2 | 516 / 302 / 4284 | -20239252.36 | PASS / 0.018 / -20239252.356 | PASS / 0.046 / -20239252.356 | PASS / 0.008 / -20239252.356 | PASS / 0.018 / -20239252.356 | PASS / 0.016 / -20239252.356 |
| agg3 | 516 / 302 / 4300 | 10312115.94 | PASS / 0.019 / 10312115.9351 | PASS / 0.049 / 10312115.9351 | PASS / 0.007 / 10312115.9351 | PASS / 0.015 / 10312115.9351 | PASS / 0.016 / 10312115.9351 |
| bandm | 305 / 472 / 2494 | -158.6280185 | PASS / 0.023 / -158.62801845 | PASS / 0.055 / -158.62801845 | PASS / 0.008 / -158.62801845 | PASS / 0.019 / -158.62801845 | PASS / 0.017 / -158.62801845 |
| beaconfd | 173 / 262 / 3375 | 33592.48581 | PASS / 0.015 / 33592.4858072 | PASS / 0.039 / 33592.4858072 | PASS / 0.008 / 33592.4858072 | PASS / 0.017 / 33592.4858072 | PASS / 0.010 / 33592.4858072 |
| blend | 74 / 83 / 491 | -30.81214985 | PASS / 0.029 / -30.8121498458 | PASS / 0.055 / -30.8121498458 | PASS / 0.002 / -30.8121498458 | PASS / 0.017 / -30.8121498459 | PASS / 0.013 / -30.8121498458 |
| bnl1 | 643 / 1175 / 5121 | 1977.629562 | PASS / 0.101 / 1977.62956152 | PASS / 0.114 / 1977.62956152 | PASS / 0.022 / 1977.62956152 | PASS / 0.035 / 1977.62956152 | PASS / 0.032 / 1977.62956152 |
| bnl2 | 2324 / 3489 / 13999 | 1811.23654 | PASS / 0.099 / 1811.23654036 | PASS / 0.210 / 1811.23654036 | PASS / 0.044 / 1811.23654036 | PASS / 0.045 / 1811.23654036 | PASS / 0.211 / 1811.23654036 |
| boeing1 | 351 / 384 / 3485 | -335.2135675 | PASS / 0.026 / -335.213567507 | PASS / 0.057 / -335.213567507 | PASS / 0.012 / -335.213567507 | PASS / 0.018 / -335.213567507 | PASS / 0.015 / -335.213567507 |
| boeing2 | 166 / 143 / 1196 | -315.018728 | PASS / 0.017 / -315.018728015 | PASS / 0.189 / -315.018728015 | PASS / 0.003 / -315.018728015 | PASS / 0.014 / -315.018728015 | PASS / 0.010 / -315.018728015 |
| bore3d | 233 / 315 / 1429 | 1373.080394 | PASS / 0.019 / 1373.08039421 | PASS / 0.046 / 1373.08039421 | PASS / 0.003 / 1373.08039421 | PASS / 0.015 / 1373.08039421 | PASS / 0.011 / 1373.08039421 |
| brandy | 220 / 249 / 2148 | 1518.509896 | PASS / 0.020 / 1518.50989649 | PASS / 0.059 / 1518.50989649 | PASS / 0.005 / 1518.50989649 | PASS / 0.016 / 1518.50989649 | PASS / 0.016 / 1518.50989649 |
| capri | 271 / 353 / 1767 | 2690.012914 | PASS / 0.019 / 2690.01291377 | PASS / 0.052 / 2690.01291377 | PASS / 0.004 / 2690.01291377 | PASS / 0.018 / 2690.01291377 | PASS / 0.011 / 2690.01291377 |
| cre-a | 3516 / 4067 / 14987 | 23595407.06 | PASS / 0.266 / 23595407.061 | PASS / 0.424 / 23595407.061 | PASS / 0.064 / 23595407.061 | PASS / 0.089 / 23595407.061 | PASS / 0.356 / 23595407.061 |
| cre-b | 9648 / 72447 / 256095 | 23129639.89 | PASS / 19.914 / 23129639.8868 | PASS / 4.816 / 23129639.8868 | PASS / 1.373 / 23129639.8868 | PASS / 1.490 / 23129639.8868 | PASS / 8.909 / 23129639.8868 |
| cre-c | 3068 / 3678 / 13244 | 25275116.14 | PASS / 0.193 / 25275116.1409 | PASS / 0.256 / 25275116.1409 | PASS / 0.062 / 25275116.1409 | PASS / 0.068 / 25275116.1409 | PASS / 0.236 / 25275116.1409 |
| cre-d | 8926 / 69980 / 242646 | 24454969.76 | SOLVER_FAILED / — / 24454997.1997 | PASS / 1.681 / 24454969.7645 | PASS / 0.669 / 24454969.7645 | PASS / 0.413 / 24454969.7645 | PASS / 3.102 / 24454969.7645 |
| cycle | 1903 / 2857 / 20720 | -5.226393025 | PASS / 0.055 / -5.22639302489 | PASS / 0.112 / -5.22639302489 | PASS / 0.138 / -5.22639302489 | PASS / 0.028 / -5.22639302489 | PASS / 0.046 / -5.22639302489 |
| czprob | 929 / 3523 / 10669 | 2185196.699 | PASS / 0.041 / 2185196.69886 | PASS / 0.089 / 2185196.69886 | PASS / 0.027 / 2185196.69886 | PASS / 0.023 / 2185196.69886 | PASS / 0.042 / 2185196.69886 |
| d2q06c | 2171 / 5167 / 32417 | 122784.2108 | PASS / 1.481 / 122784.210814 | PASS / 2.117 / 122784.210814 | PASS / 0.546 / 122784.210814 | PASS / 0.769 / 122784.210814 | PASS / 0.943 / 122784.210814 |
| d6cube | 415 / 6184 / 37704 | 315.4916667 | PASS / 0.263 / 315.491666667 | PASS / 0.283 / 315.491666667 | PASS / 0.108 / 315.491666667 | PASS / 0.064 / 315.491666667 | PASS / 1.132 / 315.491666667 |
| degen2 | 444 / 534 / 3978 | -1435.178 | PASS / 0.032 / -1435.178 | PASS / 0.070 / -1435.178 | PASS / 0.013 / -1435.178 | PASS / 0.021 / -1435.178 | PASS / 0.023 / -1435.178 |
| degen3 | 1503 / 1818 / 24646 | -987.294 | PASS / 0.328 / -987.294 | PASS / 0.662 / -987.294 | PASS / 0.128 / -987.294 | PASS / 0.121 / -987.294 | PASS / 0.184 / -987.294 |
| dfl001 | 6071 / 12230 / 35632 | 11266396.05 | TIMEOUT / — / — | PASS / 14.119 / 11266396.0467 | PASS / 5.152 / 11266396.0467 | PASS / 4.147 / 11266396.0467 | PASS / 16.157 / 11266396.0467 |
| e226 | 223 / 282 / 2578 | -11.63892907 | PASS / 0.036 / -11.6389290664 | PASS / 0.080 / -11.6389290664 | PASS / 0.012 / -11.6389290664 | PASS / 0.028 / -11.6389290664 | PASS / 0.017 / -11.6389290664 |
| etamacro | 400 / 688 / 2409 | -755.7152334 | PASS / 0.036 / -755.715233352 | PASS / 0.086 / -755.715233352 | PASS / 0.008 / -755.715233301 | PASS / 0.022 / -755.715233339 | PASS / 0.019 / -755.715233139 |
| fffff800 | 524 / 854 / 6227 | 555679.5648 | PASS / 0.035 / 555679.564817 | PASS / 0.154 / 555679.564817 | PASS / 0.017 / 555679.564817 | PASS / 0.028 / 555679.564817 | PASS / 0.018 / 555679.564817 |
| finnis | 497 / 614 / 2310 | 172791.0656 | PASS / 0.023 / 172791.065596 | PASS / 0.065 / 172791.065596 | PASS / 0.006 / 172791.065596 | PASS / 0.025 / 172791.065596 | PASS / 0.024 / 172791.065596 |
| fit1d | 24 / 1026 / 13404 | -9146.378092 | PASS / 0.035 / -9146.37809242 | PASS / 0.151 / -9146.37809242 | PASS / 0.026 / -9146.37809242 | PASS / 0.029 / -9146.37809242 | PASS / 0.043 / -9146.37809242 |
| fit1p | 627 / 1677 / 9868 | 9146.378092 | PASS / 0.089 / 9146.37809242 | PASS / 0.216 / 9146.37809242 | PASS / 0.036 / 9146.37809242 | PASS / 0.051 / 9146.37809242 | PASS / 0.038 / 9146.37809242 |
| fit2d | 25 / 10500 / 129018 | -68464.29329 | PASS / 0.256 / -68464.2932938 | PASS / 0.775 / -68464.2932938 | PASS / 0.202 / -68464.2932938 | PASS / 0.120 / -68464.2932938 | PASS / 2.390 / -68464.2932938 |
| fit2p | 3000 / 13525 / 50284 | 68464.29329 | PASS / 3.967 / 68464.2932938 | PASS / 1.331 / 68464.2932938 | PASS / 0.871 / 68464.2932938 | PASS / 1.301 / 68464.2932938 | PASS / 1.288 / 68464.2932938 |
| forplan | 161 / 421 / 4563 | -664.218961272 | PASS / 0.024 / -664.218961272 | PASS / 0.064 / -664.218961272 | PASS / 0.008 / -664.218961272 | PASS / 0.041 / -664.218961272 | PASS / 0.020 / -664.218961272 |
| ganges | 1309 / 1681 / 6912 | -109585.7361 | PASS / 0.040 / -109585.736129 | PASS / 0.091 / -109585.736129 | PASS / 0.017 / -109585.736129 | PASS / 0.026 / -109585.736129 | PASS / 0.039 / -109585.736129 |
| gfrd-pnc | 616 / 1092 / 2377 | 6902236 | PASS / 0.024 / 6902235.99955 | PASS / 0.070 / 6902235.99955 | PASS / 0.011 / 6902235.99955 | PASS / 0.020 / 6902235.99955 | PASS / 0.026 / 6902235.99955 |
| greenbea | 2392 / 5405 / 30877 | -72555248.13 | PASS / 0.691 / -72555248.1298 | AUDIT_FAILED / — / -72555248.1298 | PASS / 0.175 / -72555248.1298 | PASS / 0.395 / -72555248.1299 | AUDIT_FAILED / — / -72555248.1298 |
| greenbeb | 2392 / 5405 / 30877 | -4302260.261 | PASS / 2.061 / -4302260.26121 | PASS / 1.116 / -4302260.26121 | PASS / 0.333 / -4302260.26121 | PASS / 0.264 / -4302260.26121 | PASS / 0.232 / -4302260.26121 |
| grow15 | 300 / 645 / 5620 | -106870941.3 | PASS / 0.069 / -106870941.294 | PASS / 0.100 / -106870941.294 | PASS / 0.026 / -106870941.294 | PASS / 0.026 / -106870941.294 | PASS / 0.024 / -106870941.294 |
| grow22 | 440 / 946 / 8252 | -160834336.5 | PASS / 0.157 / -160834336.483 | PASS / 0.184 / -160834336.483 | PASS / 0.069 / -160834336.483 | PASS / 0.032 / -160834336.483 | PASS / 0.037 / -160834336.483 |
| grow7 | 140 / 301 / 2612 | -47787811.81 | PASS / 0.027 / -47787811.8147 | PASS / 0.056 / -47787811.8147 | PASS / 0.008 / -47787811.8147 | PASS / 0.021 / -47787811.8147 | AUDIT_FAILED / — / -47787811.8147 |
| israel | 174 / 142 / 2269 | -896644.8219 | PASS / 0.019 / -896644.821863 | PASS / 0.043 / -896644.821863 | PASS / 0.004 / -896644.821863 | PASS / 0.015 / -896644.821863 | PASS / 0.011 / -896644.821863 |
| kb2 | 43 / 41 / 286 | -1749.90013 | PASS / 0.014 / -1749.90012991 | PASS / 0.038 / -1749.90012991 | PASS / 0.002 / -1749.90012991 | PASS / 0.016 / -1749.90012991 | PASS / 0.008 / -1749.90012991 |
| ken-07 | 2426 / 3602 / 8404 | -679520443.4 | PASS / 0.048 / -679520443.382 | PASS / 0.096 / -679520443.382 | PASS / 0.026 / -679520443.382 | PASS / 0.025 / -679520443.382 | PASS / 0.057 / -679520443.382 |
| ken-11 | 14694 / 21349 / 49058 | -6972382263 | PASS / 0.824 / -6972382262.52 | PASS / 0.820 / -6972382262.52 | PASS / 0.207 / -6972382262.52 | PASS / 0.134 / -6972382262.52 | PASS / 2.142 / -6972382262.52 |
| ken-13 | 28632 / 42659 / 97246 | -10257394789 | PASS / 3.449 / -10257394789.5 | PASS / 3.543 / -10257394789.5 | PASS / 0.614 / -10257394789.5 | PASS / 0.313 / -10257394789.5 | PASS / 16.288 / -10257394789.5 |
| ken-18 | 105127 / 154699 / 358171 | -52217025287 | PASS / 50.033 / -52217025287.4 | PASS / 41.176 / -52217025287.4 | PASS / 2.706 / -52217025287.4 | PASS / 1.700 / -52217025287.4 | TIMEOUT / — / — |
| lotfi | 153 / 308 / 1078 | -25.26470606 | PASS / 0.018 / -25.2647060619 | PASS / 0.042 / -25.2647060619 | PASS / 0.003 / -25.2647060619 | PASS / 0.015 / -25.2647060619 | PASS / 0.016 / -25.2647060619 |
| maros | 846 / 1443 / 9614 | -58063.7437 | PASS / 0.075 / -58063.7437011 | PASS / 0.100 / -58063.7437011 | PASS / 0.032 / -58063.7437011 | PASS / 0.036 / -58063.7437011 | PASS / 0.062 / -58063.7437011 |
| maros-r7 | 3136 / 9408 / 144848 | 1497185.166 | PASS / 1.270 / 1497185.16648 | PASS / 0.796 / 1497185.16648 | PASS / 0.465 / 1497185.16648 | PASS / 0.838 / 1497185.16648 | PASS / 2.049 / 1497185.16648 |
| modszk1 | 687 / 1620 / 3168 | 320.6197291 | PASS / 0.030 / 320.619729065 | PASS / 0.077 / 320.619729064 | PASS / 0.013 / 320.619729065 | PASS / 0.020 / 320.619729064 | PASS / 0.037 / 320.619729064 |
| nesm | 662 / 2923 / 13288 | 14076036.49 | PASS / 0.137 / 14076036.4876 | PASS / 0.324 / 14076036.4876 | PASS / 0.120 / 14076036.4876 | PASS / 0.048 / 14076036.4876 | PASS / 0.233 / 14076036.4876 |
| osa-07 | 1118 / 23949 / 143694 | 535722.5173 | PASS / 0.227 / 535722.517299 | PASS / 1.997 / 535722.517299 | PASS / 0.247 / 535722.517299 | PASS / 0.078 / 535722.517299 | PASS / 0.395 / 535722.517299 |
| osa-14 | 2337 / 52460 / 314760 | 1106462.845 | PASS / 0.591 / 1106462.84474 | PASS / 32.420 / 1106462.84474 | PASS / 0.894 / 1106462.84474 | PASS / 0.208 / 1106462.84474 | PASS / 2.577 / 1106462.84474 |
| osa-30 | 4350 / 100024 / 600138 | 2142139.873 | PASS / 1.819 / 2142139.87321 | TIMEOUT / — / — | PASS / 5.323 / 2142139.87321 | PASS / 0.407 / 2142139.87321 | PASS / 6.213 / 2142139.87321 |
| osa-60 | 10280 / 232966 / 1397793 | 4044072.503 | PASS / 11.012 / 4044072.50316 | TIMEOUT / — / — | PASS / 16.993 / 4044072.50316 | PASS / 1.288 / 4044072.50316 | PASS / 18.052 / 4044072.50316 |
| pds-02 | 2953 / 7535 / 16390 | 28857862010 | PASS / 0.098 / 28857862010 | PASS / 0.178 / 28857862010 | PASS / 0.084 / 28857862010 | PASS / 0.043 / 28857862010 | PASS / 0.123 / 28857862010 |
| pds-06 | 9881 / 28655 / 62524 | 27761037600 | PASS / 1.203 / 27761037600 | PASS / 0.865 / 27761037600 | PASS / 0.300 / 27761037600 | PASS / 0.164 / 27761037600 | PASS / 1.194 / 27761037600 |
| pds-10 | 16558 / 48763 / 106436 | 26727094976 | PASS / 6.024 / 26727094976 | PASS / 2.489 / 26727094976 | PASS / 0.549 / 26727094976 | PASS / 0.335 / 26727094976 | PASS / 4.529 / 26727094976 |
| pds-20 | 33874 / 105728 / 230200 | 23821658640 | PASS / 62.991 / 23821658640 | PASS / 10.227 / 23821658640 | PASS / 1.830 / 23821658640 | PASS / 1.100 / 23821658640 | PASS / 41.190 / 23821658640 |
| perold | 625 / 1376 / 6018 | -9380.755278 | PASS / 0.143 / -9380.75527824 | PASS / 0.361 / -9380.75527824 | PASS / 0.055 / -9380.75527824 | PASS / 0.102 / -9380.75527824 | PASS / 0.080 / -9380.75527824 |
| pilot | 1441 / 3652 / 43167 | -557.4897293 | PASS / 2.845 / -557.489729284 | PASS / 2.022 / -557.489729285 | PASS / 1.009 / -557.489729291 | PASS / 0.791 / -557.489729284 | PASS / 1.691 / -557.489706157 |
| pilot.ja | 940 / 1988 / 14698 | -6113.136466 | PASS / 0.202 / -6113.13646558 | PASS / 0.312 / -6113.13646558 | PASS / 0.074 / -6113.13646558 | PASS / 0.193 / -6113.13646558 | PASS / 0.158 / -6113.13646558 |
| pilot.we | 722 / 2789 / 9126 | -2720107.533 | PASS / 0.235 / -2720107.53284 | PASS / 0.531 / -2720107.53284 | PASS / 0.129 / -2720107.53284 | AUDIT_FAILED / — / -2720107.53927 | PASS / 0.098 / -2720107.53284 |
| pilot4 | 410 / 1000 / 5141 | -2581.13926 | PASS / 0.054 / -2581.13925888 | PASS / 0.111 / -2581.13925888 | PASS / 0.031 / -2581.13925888 | PASS / 0.049 / -2581.13925888 | PASS / 0.036 / -2581.13925888 |
| pilot87 | 2030 / 4883 / 73152 | 301.7103473 | PASS / 22.582 / 301.710347789 | PASS / 5.332 / 301.710347333 | PASS / 3.706 / 301.710347333 | PASS / 3.286 / 301.710347333 | PASS / 4.782 / 301.710347333 |
| pilotnov | 975 / 2172 / 13057 | -4497.276188 | PASS / 0.143 / -4497.27618822 | PASS / 0.258 / -4497.27618822 | PASS / 0.139 / -4497.27618822 | AUDIT_FAILED / — / -4497.27618822 | PASS / 0.082 / -4497.27618822 |
| qap12 | 3192 / 8856 / 38304 | 522.8943506 | PASS / 38.479 / 522.894342332 | TIMEOUT / — / — | TIMEOUT / — / — | PASS / 33.494 / 522.894350559 | PASS / 26.808 / 522.894350559 |
| qap15 | 6330 / 22275 / 94950 | 1040.994041 | PASS / 66.745 / 1040.9940074 | TIMEOUT / — / — | TIMEOUT / — / — | TIMEOUT / — / — | TIMEOUT / — / — |
| qap8 | 912 / 1632 / 7296 | 203.5 | PASS / 8.282 / 203.5 | PASS / 1.162 / 203.5 | PASS / 0.443 / 203.5 | PASS / 0.375 / 203.5 | PASS / 0.384 / 203.5 |
| recipe | 91 / 180 / 663 | -266.616 | PASS / 0.015 / -266.616 | PASS / 0.049 / -266.616 | PASS / 0.002 / -266.616 | PASS / 0.018 / -266.616 | PASS / 0.014 / -266.616 |
| sc105 | 105 / 103 / 280 | -52.20206121 | PASS / 0.019 / -52.2020612117 | PASS / 0.041 / -52.2020612117 | PASS / 0.001 / -52.2020612117 | PASS / 0.013 / -52.2020612117 | PASS / 0.015 / -52.2020612117 |
| sc205 | 205 / 203 / 551 | -52.20206121 | PASS / 0.023 / -52.2020612117 | PASS / 0.041 / -52.2020612117 | PASS / 0.004 / -52.2020612117 | PASS / 0.019 / -52.2020612117 | PASS / 0.014 / -52.2020612117 |
| sc50a | 50 / 48 / 130 | -64.57507706 | PASS / 0.015 / -64.5750770586 | PASS / 0.047 / -64.5750770586 | PASS / 0.002 / -64.5750770586 | PASS / 0.019 / -64.5750770586 | PASS / 0.010 / -64.5750770586 |
| sc50b | 50 / 48 / 118 | -70 | PASS / 0.016 / -70 | PASS / 0.050 / -70 | PASS / 0.001 / -70 | PASS / 0.020 / -70 | PASS / 0.012 / -70 |
| scagr25 | 471 / 500 / 1554 | -14753433.06 | PASS / 0.023 / -14753433.0608 | PASS / 0.060 / -14753433.0608 | PASS / 0.006 / -14753433.0608 | PASS / 0.016 / -14753433.0608 | PASS / 0.019 / -14753433.0608 |
| scagr7 | 129 / 140 / 420 | -2331389.824 | PASS / 0.015 / -2331389.82433 | PASS / 0.040 / -2331389.82433 | PASS / 0.002 / -2331389.82433 | PASS / 0.015 / -2331389.82433 | PASS / 0.015 / -2331389.82433 |
| scfxm1 | 330 / 457 / 2589 | 18416.75903 | PASS / 0.030 / 18416.7590283 | PASS / 0.069 / 18416.7590283 | PASS / 0.007 / 18416.7590283 | PASS / 0.026 / 18416.7590283 | PASS / 0.039 / 18416.7590283 |
| scfxm2 | 660 / 914 / 5183 | 36660.26156 | PASS / 0.050 / 36660.261565 | PASS / 0.091 / 36660.261565 | PASS / 0.016 / 36660.261565 | PASS / 0.033 / 36660.261565 | PASS / 0.034 / 36660.261565 |
| scfxm3 | 990 / 1371 / 7777 | 54901.25455 | PASS / 0.078 / 54901.2545498 | PASS / 0.114 / 54901.2545498 | PASS / 0.027 / 54901.2545498 | PASS / 0.037 / 54901.2545497 | PASS / 0.048 / 54901.2545498 |
| scorpion | 388 / 358 / 1426 | 1878.124823 | PASS / 0.018 / 1878.12482274 | PASS / 0.055 / 1878.12482274 | PASS / 0.004 / 1878.12482274 | PASS / 0.021 / 1878.12482274 | PASS / 0.023 / 1878.12482274 |
| scrs8 | 490 / 1169 / 3182 | 904.2969538 | PASS / 0.026 / 904.296953801 | PASS / 0.074 / 904.296953801 | PASS / 0.011 / 904.296953801 | PASS / 0.022 / 904.296953799 | PASS / 0.036 / 904.296953801 |
| scsd1 | 77 / 760 / 2388 | 8.666666674 | PASS / 0.020 / 8.66666667433 | PASS / 0.051 / 8.66666667433 | PASS / 0.004 / 8.66666667433 | PASS / 0.023 / 8.66666667433 | PASS / 0.012 / 8.66666667433 |
| scsd6 | 147 / 1350 / 4316 | 50.50000008 | PASS / 0.024 / 50.5000000771 | PASS / 0.076 / 50.500000078 | PASS / 0.016 / 50.5000000783 | PASS / 0.020 / 50.5000000783 | PASS / 0.025 / 50.5000000783 |
| scsd8 | 397 / 2750 / 8584 | 904.9999999 | PASS / 0.081 / 904.999999925 | PASS / 0.142 / 904.999999925 | PASS / 0.043 / 904.999999925 | PASS / 0.051 / 904.999999925 | PASS / 0.054 / 904.999999925 |
| sctap1 | 300 / 480 / 1692 | 1412.25 | PASS / 0.025 / 1412.25 | PASS / 0.078 / 1412.25 | PASS / 0.006 / 1412.25 | PASS / 0.025 / 1412.25 | PASS / 0.022 / 1412.25 |
| sctap2 | 1090 / 1880 / 6714 | 1724.807143 | PASS / 0.036 / 1724.80714286 | PASS / 0.090 / 1724.80714286 | PASS / 0.013 / 1724.80714286 | PASS / 0.022 / 1724.80714286 | PASS / 0.048 / 1724.80714286 |
| sctap3 | 1480 / 2480 / 8874 | 1424 | PASS / 0.040 / 1424 | PASS / 0.104 / 1424 | PASS / 0.018 / 1424 | PASS / 0.035 / 1424 | PASS / 0.063 / 1424 |
| seba | 515 / 1028 / 4352 | 15711.6 | PASS / 0.020 / 15711.6 | PASS / 0.061 / 15711.6 | PASS / 0.006 / 15711.6 | PASS / 0.016 / 15711.6 | PASS / 0.029 / 15711.6 |
| share1b | 117 / 225 / 1151 | -76589.31858 | PASS / 0.018 / -76589.3185792 | PASS / 0.049 / -76589.3185792 | PASS / 0.003 / -76589.3185792 | PASS / 0.015 / -76589.3185792 | PASS / 0.017 / -76589.3185792 |
| share2b | 96 / 79 / 694 | -415.7322407 | PASS / 0.017 / -415.732240741 | PASS / 0.044 / -415.732240741 | PASS / 0.002 / -415.732240741 | PASS / 0.022 / -415.732240741 | PASS / 0.017 / -415.732240741 |
| shell | 536 / 1775 / 3556 | 1208825346 | PASS / 0.028 / 1208825346 | PASS / 0.080 / 1208825346 | PASS / 0.019 / 1208825346 | PASS / 0.023 / 1208825346 | PASS / 0.028 / 1208825346 |
| ship04l | 402 / 2118 / 6332 | 1793324.538 | PASS / 0.034 / 1793324.53797 | PASS / 0.075 / 1793324.53797 | PASS / 0.010 / 1793324.53797 | PASS / 0.023 / 1793324.53797 | PASS / 0.032 / 1793324.53797 |
| ship04s | 402 / 1458 / 4352 | 1798714.7 | PASS / 0.028 / 1798714.70045 | PASS / 0.062 / 1798714.70045 | PASS / 0.007 / 1798714.70045 | PASS / 0.020 / 1798714.70045 | PASS / 0.027 / 1798714.70045 |
| ship08l | 778 / 4283 / 12802 | 1909055.211 | PASS / 0.049 / 1909055.21139 | PASS / 0.126 / 1909055.21139 | PASS / 0.021 / 1909055.21139 | PASS / 0.027 / 1909055.21139 | PASS / 0.054 / 1909055.21139 |
| ship08s | 778 / 2387 / 7114 | 1920098.211 | PASS / 0.026 / 1920098.21053 | PASS / 0.077 / 1920098.21053 | PASS / 0.012 / 1920098.21053 | PASS / 0.025 / 1920098.21053 | PASS / 0.042 / 1920098.21053 |
| ship12l | 1151 / 5427 / 16170 | 1470187.919 | PASS / 0.068 / 1470187.91933 | PASS / 0.137 / 1470187.91933 | PASS / 0.038 / 1470187.91933 | PASS / 0.031 / 1470187.91933 | PASS / 0.082 / 1470187.91933 |
| ship12s | 1151 / 2763 / 8178 | 1489236.134 | PASS / 0.043 / 1489236.13441 | PASS / 0.080 / 1489236.13441 | PASS / 0.020 / 1489236.13441 | PASS / 0.022 / 1489236.13441 | PASS / 0.036 / 1489236.13441 |
| sierra | 1227 / 2036 / 7302 | 15394362.18 | PASS / 0.050 / 15394362.1836 | PASS / 0.083 / 15394362.1836 | PASS / 0.023 / 15394362.1836 | PASS / 0.028 / 15394362.1836 | PASS / 0.058 / 15394362.1836 |
| stair | 356 / 467 / 3856 | -251.2669512 | PASS / 0.052 / -251.266951193 | PASS / 0.089 / -251.266951193 | PASS / 0.024 / -251.266951193 | PASS / 0.039 / -251.266951193 | PASS / 0.034 / -251.266951193 |
| standata | 359 / 1075 / 3031 | 1257.6995 | PASS / 0.030 / 1257.6995 | PASS / 0.057 / 1257.6995 | PASS / 0.006 / 1257.6995 | PASS / 0.026 / 1257.6995 | PASS / 0.019 / 1257.6995 |
| standgub | 361 / 1184 / 3139 | 1257.6995 | PASS / 0.019 / 1257.6995 | PASS / 0.059 / 1257.6995 | PASS / 0.005 / 1257.6995 | PASS / 0.020 / 1257.6995 | PASS / 0.016 / 1257.6995 |
| standmps | 467 / 1075 / 3679 | 1406.0175 | PASS / 0.020 / 1406.0175 | PASS / 0.066 / 1406.0175 | PASS / 0.006 / 1406.0175 | PASS / 0.030 / 1406.0175 | PASS / 0.025 / 1406.0175 |
| stocfor1 | 117 / 111 / 447 | -41131.97622 | PASS / 0.017 / -41131.9762194 | PASS / 0.045 / -41131.9762194 | PASS / 0.002 / -41131.9762194 | PASS / 0.014 / -41131.9762194 | PASS / 0.019 / -41131.9762194 |
| stocfor2 | 2157 / 2031 / 8343 | -39024.40854 | PASS / 0.121 / -39024.4085379 | PASS / 0.131 / -39024.4085379 | PASS / 0.033 / -39024.4085379 | PASS / 0.032 / -39024.4085379 | PASS / 0.112 / -39024.4085379 |
| stocfor3 | 16675 / 15695 / 64875 | -39976.78394 | PASS / 5.434 / -39976.7839436 | PASS / 1.536 / -39976.7839436 | PASS / 0.388 / -39976.7839436 | PASS / 0.382 / -39976.7839436 | PASS / 6.917 / -39976.7839436 |
| truss | 1000 / 8806 / 27836 | 458815.8472 | PASS / 16.906 / 458815.847186 | PASS / 1.710 / 458815.847186 | PASS / 1.911 / 458815.847186 | PASS / 1.606 / 458815.847186 | PASS / 0.875 / 458815.847186 |
| tuff | 333 / 587 / 4520 | 0.292147765 | PASS / 0.019 / 0.292147765094 | PASS / 0.064 / 0.292147765094 | PASS / 0.007 / 0.292147765094 | PASS / 0.015 / 0.292147765094 | PASS / 0.027 / 0.292147765094 |
| vtp.base | 198 / 203 / 908 | 129831.4625 | PASS / 0.013 / 129831.462461 | PASS / 0.037 / 129831.462461 | PASS / 0.002 / 129831.462461 | PASS / 0.013 / 129831.462461 | PASS / 0.015 / 129831.462461 |
| wood1p | 244 / 2594 / 70215 | 1.442902412 | PASS / 0.069 / 1.44290241157 | PASS / 0.194 / 1.44290241157 | PASS / 0.083 / 1.44290241157 | PASS / 0.040 / 1.44290241157 | PASS / 0.122 / 1.44290241157 |
| woodw | 1098 / 8405 / 37474 | 1.304476333 | PASS / 0.073 / 1.30447633308 | PASS / 0.158 / 1.30447633308 | PASS / 0.063 / 1.30447633308 | PASS / 0.039 / 1.30447633308 | PASS / 0.110 / 1.30447633308 |
