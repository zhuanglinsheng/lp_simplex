# 基准测试汇总

简体中文 | [English](./README_en.md) | [项目主页](../../README.md) | [文档](../README.md) | [历史性能](./history/README.md)

由 `data/benchmark.py` 自动生成；仅包含测试信息、对比统计及复现说明。

## 测试信息与统计口径

逐案例生成临时 MPS，本轮所有求解器完成全部重复后清理，再同步更新数据集报告和本汇总。

成功案例须全部重复 PASS 且每次均在限时内。先取重复耗时中位数，再计算普通几何平均，不移位、不对失败加惩罚；空集合为 `—`。各求解器分别与 lp_simplex 使用双方共同成功案例比较；几何平均按“lp_simplex / 本求解器”列出，各行集合可能不同；lp_simplex 行不列自比较统计。相对耗时为 lp_simplex 耗时除以本求解器耗时，大于 1 表示本求解器更快。成功数按案例计，超时与其他未通过数按重复运行计。

计时含模型读取、求解和解导出，不含 MPS 生成/删除、格式兼容处理、公共适配器初始化和独立审核。串行、单线程，不保证机器独占；CLI 与 Python binding 开销有差异。

可行问题的 PASS 要求求解器报告最优且原模型可行性、目标复算通过；有参考目标时额外核对。缺少参考值不代表独立最优性认证。进行中或子集结果不构成全集结论，同一机器、求解器、输入和参数相同的测量直接复用，并按共同成功案例重算比较。

Netlib infeasible 单独统计：PASS 只核对明确的不可行状态，不代表独立不可行证书验证。

## Netlib feasible

[逐案例结果与完整报告](netlib.md)

展示当前版本的逐例结果及对比统计。

<!-- benchmark-run:accepted-netlib -->
| 测试信息 | 值 |
| --- | --- |
| 状态 / 进度 | 已完成；114/114 个案例 |
| 时间（UTC） | 2026-10-28T12:58:45.277298+00:00 |
| 对比求解器 | lp_simplex, SCIP / SoPlex, HiGHS, Clp, GLPK |
| 案例范围 | netlib 全集；25fv47, 80bau3b, adlittle, afiro, agg, agg2, agg3, bandm, beaconfd, blend, bnl1, bnl2, boeing1, boeing2, bore3d, brandy, capri, cre-a, cre-b, cre-c, cre-d, cycle, czprob, d2q06c, d6cube, degen2, degen3, dfl001, e226, etamacro, fffff800, finnis, fit1d, fit1p, fit2d, fit2p, forplan, ganges, gfrd-pnc, greenbea, greenbeb, grow15, grow22, grow7, israel, kb2, ken-07, ken-11, ken-13, ken-18, lotfi, maros, maros-r7, modszk1, nesm, osa-07, osa-14, osa-30, osa-60, pds-02, pds-06, pds-10, pds-20, perold, pilot, pilot.ja, pilot.we, pilot4, pilot87, pilotnov, qap12, qap15, qap8, recipe, sc105, sc205, sc50a, sc50b, scagr25, scagr7, scfxm1, scfxm2, scfxm3, scorpion, scrs8, scsd1, scsd6, scsd8, sctap1, sctap2, sctap3, seba, share1b, share2b, shell, ship04l, ship04s, ship08l, ship08s, ship12l, ship12s, sierra, stair, standata, standgub, standmps, stocfor1, stocfor2, stocfor3, truss, tuff, vtp.base, wood1p, woodw |
| 系统 / CPU | macOS-27.0.1-arm64-arm-64bit-Mach-O / Apple M1 |
| 内存 / 逻辑 CPU | 17179869184 bytes / 8 |
| Python / NumPy / SciPy | 3.14.2 / 2.5.3 / 1.18.1 |
| lp_simplex 版本 | b6f657753fe886c7de801ccff83be0ab4129a3ad；工作区有改动 |
| 构建 / 编译器 | Release / Apple clang version 21.0.0 (clang-2100.3.34.2) |
| SCIP / LP 后端 | SCIP version 10.0.3 [precision: 8 byte] [memory: block] [mode: optimized] [LP solver: SoPlex 8.0.2] [GitHash: d409edf9f6] |
| HiGHS | 1.12.0 / scipy-private-binding |
| 算法（lp / scip / highs） | auto / default / choose |
| 重复 / 单次限时 / 容差 / 工作预算 | 1 / 90 秒 / 1e-07 / 300000 |
| 资源 | 串行；每个求解器单线程；未限制内存，不保证机器独占 |
| 原始记录 | 逐例测量与指纹保存在本报告中；临时文件已清理 |
| Clp | Coin LP version 1.17.11, build Mar 11 2026 |
| GLPK | GLPSOL--GLPK LP/MIP Solver 5.0 |
| Clp / GLPK 算法 | dual simplex / primal simplex；presolve on |

| 求解器 | 成功案例 | 超时次数 | 其他未通过次数 | 与 lp_simplex 共同成功数 | 共同案例几何平均/秒（lp_simplex / 本求解器） | 相对耗时（lp_simplex / 各求解器） |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| lp_simplex | 114/114 | 0 | 0 | — | — | — |
| SCIP / SoPlex | 109/114 | 4 | 1 | 109 | 0.0822446 / 0.174508 | 0.471295 |
| HiGHS | 112/114 | 2 | 0 | 112 | 0.0908197 / 0.0310268 | 2.92714 |
| Clp | 111/114 | 1 | 2 | 111 | 0.0948068 / 0.0549174 | 1.72635 |
| GLPK | 110/114 | 2 | 2 | 110 | 0.0916918 / 0.0865553 | 1.05934 |

复现本轮案例与参数（新建运行目录）：

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DLP_SIMPLEX_BUILD_BENCHMARKS=ON

cmake --build build --target lp_solve_mps -j 4

python3 data/benchmark.py --dataset netlib --cases 25fv47 80bau3b adlittle afiro agg agg2 agg3 bandm beaconfd blend bnl1 bnl2 boeing1 boeing2 bore3d brandy capri cre-a cre-b cre-c cre-d cycle czprob d2q06c d6cube degen2 degen3 dfl001 e226 etamacro fffff800 finnis fit1d fit1p fit2d fit2p forplan ganges gfrd-pnc greenbea greenbeb grow15 grow22 grow7 israel kb2 ken-07 ken-11 ken-13 ken-18 lotfi maros-r7 maros modszk1 nesm osa-07 osa-14 osa-30 osa-60 pds-02 pds-06 pds-10 pds-20 perold pilot.ja pilot pilot.we pilot4 pilot87 pilotnov qap12 qap15 qap8 recipe sc105 sc205 sc50a sc50b scagr25 scagr7 scfxm1 scfxm2 scfxm3 scorpion scrs8 scsd1 scsd6 scsd8 sctap1 sctap2 sctap3 seba share1b share2b shell ship04l ship04s ship08l ship08s ship12l ship12s sierra stair standata standgub standmps stocfor1 stocfor2 stocfor3 truss tuff vtp.base wood1p woodw --binary build/data/lp_solve_mps --scip scip --clp clp --glpk glpsol --repeats 1 --seconds 90.0 --iterations 300000 --lp-algorithm auto --scip-algorithm default --highs-algorithm choose --solvers lp_simplex scip highs clp glpk
```
<!-- /benchmark-run -->

## Netlib infeasible

[逐案例结果与完整报告](netlib-infeasible.md)

展示当前版本的逐例结果及对比统计。

<!-- benchmark-run:accepted-infeasible -->
| 测试信息 | 值 |
| --- | --- |
| 状态 / 进度 | 已完成；29/29 个案例 |
| 时间（UTC） | 2026-10-28T13:06:50.478741+00:00 |
| 对比求解器 | lp_simplex, SCIP / SoPlex, HiGHS, Clp, GLPK |
| 案例范围 | netlib-infeasible 全集；bgdbg1, bgetam, bgindy, bgprtr, box1, ceria3d, chemcom, cplex1, cplex2, ex72a, ex73a, forest6, galenet, gosh, gran, greenbea, itest2, itest6, klein1, klein2, klein3, mondou2, pang, pilot4i, qual, reactor, refinery, vol1, woodinfe |
| 系统 / CPU | macOS-27.0.1-arm64-arm-64bit-Mach-O / Apple M1 |
| 内存 / 逻辑 CPU | 17179869184 bytes / 8 |
| Python / NumPy / SciPy | 3.14.2 / 2.5.3 / 1.18.1 |
| lp_simplex 版本 | b6f657753fe886c7de801ccff83be0ab4129a3ad；工作区有改动 |
| 构建 / 编译器 | Release / Apple clang version 21.0.0 (clang-2100.3.34.2) |
| SCIP / LP 后端 | SCIP version 10.0.3 [precision: 8 byte] [memory: block] [mode: optimized] [LP solver: SoPlex 8.0.2] [GitHash: d409edf9f6] |
| HiGHS | 1.12.0 / scipy-private-binding |
| 算法（lp / scip / highs） | auto / default / choose |
| 重复 / 单次限时 / 容差 / 工作预算 | 1 / 90 秒 / 1e-07 / 300000 |
| 资源 | 串行；每个求解器单线程；未限制内存，不保证机器独占 |
| 原始记录 | 逐例测量与指纹保存在本报告中；临时文件已清理 |
| Clp | Coin LP version 1.17.11, build Mar 11 2026 |
| GLPK | GLPSOL--GLPK LP/MIP Solver 5.0 |
| Clp / GLPK 算法 | dual simplex / primal simplex；presolve on |

| 求解器 | 成功案例 | 超时次数 | 其他未通过次数 | 与 lp_simplex 共同成功数 | 共同案例几何平均/秒（lp_simplex / 本求解器） | 相对耗时（lp_simplex / 各求解器） |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| lp_simplex | 29/29 | 0 | 0 | — | — | — |
| SCIP / SoPlex | 28/29 | 0 | 1 | 28 | 0.0279978 / 0.0424539 | 0.659487 |
| HiGHS | 29/29 | 0 | 0 | 29 | 0.0279817 / 0 | — |
| Clp | 29/29 | 0 | 0 | 29 | 0.0279817 / 0.0241099 | 1.16059 |
| GLPK | 27/29 | 0 | 2 | 27 | 0.0268443 / 0.0136329 | 1.96908 |

| 求解器 | 错误判定次数 | 判定不明确次数 | 运行或数值失败次数 | 跳过次数 |
| --- | ---: | ---: | ---: | ---: |
| lp_simplex | 0 | 0 | 0 | 0 |
| SCIP / SoPlex | 0 | 1 | 0 | 0 |
| HiGHS | 0 | 0 | 0 | 0 |
| Clp | 0 | 0 | 0 | 0 |
| GLPK | 0 | 0 | 2 | 0 |

| 实例 | 求解器 | 结果 | 终止信息 |
| --- | --- | --- | --- |
| cplex1 | SCIP / SoPlex | AMBIGUOUS | AMBIGUOUS |
| greenbea | GLPK | SOLVER_FAILED | SOLVER_FAILED |
| pang | GLPK | SOLVER_FAILED | SOLVER_FAILED |

复现本轮案例与参数（新建运行目录）：

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DLP_SIMPLEX_BUILD_BENCHMARKS=ON

cmake --build build --target lp_solve_mps -j 4

python3 data/benchmark.py --dataset netlib-infeasible --cases bgdbg1 bgetam bgindy bgprtr box1 ceria3d chemcom cplex1 cplex2 ex72a ex73a forest6 galenet gosh gran greenbea itest2 itest6 klein1 klein2 klein3 mondou2 pang pilot4i qual reactor refinery vol1 woodinfe --binary build/data/lp_solve_mps --scip scip --clp clp --glpk glpsol --repeats 1 --seconds 90.0 --iterations 300000 --lp-algorithm auto --scip-algorithm default --highs-algorithm choose --solvers lp_simplex scip highs clp glpk
```
<!-- /benchmark-run -->

## PDS

[逐案例结果与完整报告](pds.md)

展示当前版本的逐例结果及对比统计。

<!-- benchmark-run:accepted-pds -->
| 测试信息 | 值 |
| --- | --- |
| 状态 / 进度 | 已完成；9/9 个案例 |
| 时间（UTC） | 2026-10-28T12:53:58.313619+00:00 |
| 对比求解器 | lp_simplex, SCIP / SoPlex, HiGHS, Clp, GLPK |
| 案例范围 | pds 全集；pds-20, pds-30, pds-40, pds-50, pds-60, pds-70, pds-80, pds-90, pds-100 |
| 系统 / CPU | macOS-27.0.1-arm64-arm-64bit-Mach-O / Apple M1 |
| 内存 / 逻辑 CPU | 17179869184 bytes / 8 |
| Python / NumPy / SciPy | 3.14.2 / 2.5.3 / 1.18.1 |
| lp_simplex 版本 | b6f657753fe886c7de801ccff83be0ab4129a3ad；工作区有改动 |
| 构建 / 编译器 | Release / Apple clang version 21.0.0 (clang-2100.3.34.2) |
| SCIP / LP 后端 | SCIP version 10.0.3 [precision: 8 byte] [memory: block] [mode: optimized] [LP solver: SoPlex 8.0.2] [GitHash: d409edf9f6] |
| HiGHS | 1.12.0 / scipy-private-binding |
| 算法（lp / scip / highs） | auto / default / choose |
| 重复 / 单次限时 / 容差 / 工作预算 | 1 / 300 秒 / 1e-07 / 300000 |
| 资源 | 串行；每个求解器单线程；未限制内存，不保证机器独占 |
| 原始记录 | 逐例测量与指纹保存在本报告中；临时文件已清理 |
| Clp | Coin LP version 1.17.11, build Mar 11 2026 |
| GLPK | GLPSOL--GLPK LP/MIP Solver 5.0 |
| Clp / GLPK 算法 | dual simplex / primal simplex；presolve on |
| GLPK 跳过策略 | 任一案例超时后跳过后续测试；已跳过 6 次 |

| 求解器 | 成功案例 | 超时次数 | 其他未通过次数 | 与 lp_simplex 共同成功数 | 共同案例几何平均/秒（lp_simplex / 本求解器） | 相对耗时（lp_simplex / 各求解器） |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| lp_simplex | 9/9 | 0 | 0 | — | — | — |
| SCIP / SoPlex | 6/9 | 3 | 0 | 6 | 16.6697 / 70.3678 | 0.236894 |
| HiGHS | 9/9 | 0 | 0 | 9 | 26.4508 / 7.08912 | 3.73118 |
| Clp | 9/9 | 0 | 0 | 9 | 26.4508 / 7.23811 | 3.65438 |
| GLPK | 2/9 | 1 | 0 | 2 | 5.73654 / 49.4921 | 0.115908 |

复现本轮案例与参数（新建运行目录）：

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DLP_SIMPLEX_BUILD_BENCHMARKS=ON

cmake --build build --target lp_solve_mps -j 4

python3 data/benchmark.py --dataset pds --cases pds-20 pds-30 pds-40 pds-50 pds-60 pds-70 pds-80 pds-90 pds-100 --binary build/data/lp_solve_mps --scip scip --clp clp --glpk glpsol --repeats 1 --seconds 300.0 --iterations 300000 --lp-algorithm auto --scip-algorithm default --highs-algorithm choose --solvers lp_simplex scip highs clp glpk
```
<!-- /benchmark-run -->

## Railway 连续 LP 松弛

[逐案例结果与完整报告](railway.md)

展示当前版本的逐例结果及对比统计。

<!-- benchmark-run:accepted-railway -->
| 测试信息 | 值 |
| --- | --- |
| 状态 / 进度 | 已完成；7/7 个案例 |
| 时间（UTC） | 2026-10-28T13:07:22.428932+00:00 |
| 对比求解器 | lp_simplex, SCIP / SoPlex, HiGHS, Clp, GLPK |
| 案例范围 | railway 全集；rail507, rail516, rail582, rail2536, rail2586, rail4284, rail4872 |
| 系统 / CPU | macOS-27.0.1-arm64-arm-64bit-Mach-O / Apple M1 |
| 内存 / 逻辑 CPU | 17179869184 bytes / 8 |
| Python / NumPy / SciPy | 3.14.2 / 2.5.3 / 1.18.1 |
| lp_simplex 版本 | b6f657753fe886c7de801ccff83be0ab4129a3ad；工作区有改动 |
| 构建 / 编译器 | Release / Apple clang version 21.0.0 (clang-2100.3.34.2) |
| SCIP / LP 后端 | SCIP version 10.0.3 [precision: 8 byte] [memory: block] [mode: optimized] [LP solver: SoPlex 8.0.2] [GitHash: d409edf9f6] |
| HiGHS | 1.12.0 / scipy-private-binding |
| 算法（lp / scip / highs） | auto / default / choose |
| 重复 / 单次限时 / 容差 / 工作预算 | 1 / 3600 秒 / 1e-07 / 10000000 |
| 资源 | 串行；每个求解器单线程；未限制内存，不保证机器独占 |
| 原始记录 | 逐例测量与指纹保存在本报告中；临时文件已清理 |
| Clp | Coin LP version 1.17.11, build Mar 11 2026 |
| GLPK | GLPSOL--GLPK LP/MIP Solver 5.0 |
| Clp / GLPK 算法 | dual simplex / primal simplex；presolve on |
| GLPK 跳过策略 | 任一案例超时后跳过后续测试；已跳过 3 次 |

| 求解器 | 成功案例 | 超时次数 | 其他未通过次数 | 与 lp_simplex 共同成功数 | 共同案例几何平均/秒（lp_simplex / 本求解器） | 相对耗时（lp_simplex / 各求解器） |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| lp_simplex | 7/7 | 0 | 0 | — | — | — |
| SCIP / SoPlex | 7/7 | 0 | 0 | 7 | 34.2404 / 128.945 | 0.265543 |
| HiGHS | 7/7 | 0 | 0 | 7 | 34.2404 / 90.0054 | 0.380426 |
| Clp | 7/7 | 0 | 0 | 7 | 34.2404 / 53.5764 | 0.639095 |
| GLPK | 3/7 | 1 | 0 | 3 | 1.42014 / 3.19814 | 0.444052 |

复现本轮案例与参数（新建运行目录）：

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DLP_SIMPLEX_BUILD_BENCHMARKS=ON

cmake --build build --target lp_solve_mps -j 4

python3 data/benchmark.py --dataset railway --cases rail507 rail516 rail582 rail2536 rail2586 rail4284 rail4872 --binary build/data/lp_solve_mps --scip scip --clp clp --glpk glpsol --repeats 1 --seconds 3600.0 --iterations 10000000 --lp-algorithm auto --scip-algorithm default --highs-algorithm choose --solvers lp_simplex scip highs clp glpk
```
<!-- /benchmark-run -->

## 复现说明

在仓库根目录运行。需 CMake、C 编译器、BLAS/LAPACK，并安装 SCIP、Clp、GLPK、NumPy、SciPy，以及 highspy 或兼容的 SciPy HiGHS 绑定。PDS 转换器由脚本临时编译。

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DLP_SIMPLEX_BUILD_BENCHMARKS=ON

cmake --build build --target lp_solve_mps -j 4

python3 data/benchmark.py --dataset netlib --binary build/data/lp_solve_mps

python3 data/benchmark.py --dataset pds --binary build/data/lp_solve_mps

python3 data/benchmark.py --dataset railway --binary build/data/lp_solve_mps
```

省略 `--cases` 测试所选数据集全集；短程验证可加 `--cases <实例名> --repeats 1 --seconds 10`。simple_examples 仅用于功能测试，不参与基准。上方各轮的复现命令保留其案例与具体参数，二进制路径可换为新构建结果；复核版本及 SHA-256 后再比较耗时。

默认每轮新建 `tmp/benchmarks/<UTC 时间>/`，保存配置、环境、命令、结果、日志与审核信息。正常结束、异常、Ctrl+C 或 SIGTERM 清理临时 MPS；SIGKILL/断电残留可随运行目录删除。

直接从已保存的 Markdown 更新报告和汇总，不需要临时日志或再次求解：

```sh
python3 data/benchmark.py --refresh-reports
```

如仍保留原始运行目录，也可从中重建报告：

```sh
python3 data/benchmark.py --report-only tmp/benchmarks/<运行目录>
```

`--docs-dir` 可指定报告目录；重建汇总会保留其他数据集已保存的测试信息和统计。更多参数与审核说明见 [数据与基准说明](../../data/README.md)。
