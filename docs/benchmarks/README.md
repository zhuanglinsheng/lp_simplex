# 基准测试汇总

简体中文 | [English](./README_en.md) | [项目主页](../../README.md) | [文档](../README.md)

由 `data/benchmark.py` 自动生成；仅包含测试信息、对比统计及复现说明。

## 测试信息与统计口径

逐案例生成临时 MPS，本轮所有求解器完成全部重复后清理，再同步更新数据集报告和本汇总。

成功案例须全部重复 PASS 且每次均在限时内。先取重复耗时中位数，再计算普通几何平均，不移位、不对失败加惩罚；空集合为 `—`。各求解器分别与 lp_simplex 使用双方共同成功案例比较；几何平均按“lp_simplex / 本求解器”列出，各行集合可能不同；lp_simplex 行不列自比较统计。相对耗时为 lp_simplex 耗时除以本求解器耗时，大于 1 表示本求解器更快。成功数按案例计，超时与其他未通过数按重复运行计。

计时含模型读取、求解和解导出，不含 MPS 生成/删除、格式兼容处理、公共适配器初始化和独立审核。串行、单线程，不保证机器独占；CLI 与 Python binding 开销有差异。

PASS 要求求解器报告最优且原模型可行性、目标复算通过；有参考目标时额外核对。缺少参考值不代表独立最优性认证。进行中或子集结果不构成全集结论，不同运行分别统计，不混合环境、配置或案例集合。

## Netlib feasible

[逐案例结果与完整报告](netlib.md)

<!-- benchmark-run:netlib-auto-current-90s-once -->
| 测试信息 | 值 |
| --- | --- |
| 状态 / 进度 | 已完成；114/114 个案例 |
| 时间（UTC） | 2026-10-23T17:35:04.295362+00:00 |
| 参与求解器 | lp_simplex, SCIP / SoPlex, HiGHS, Clp, GLPK |
| 案例范围 | netlib 全集；25fv47, 80bau3b, adlittle, afiro, agg, agg2, agg3, bandm, beaconfd, blend, bnl1, bnl2, boeing1, boeing2, bore3d, brandy, capri, cre-a, cre-b, cre-c, cre-d, cycle, czprob, d2q06c, d6cube, degen2, degen3, dfl001, e226, etamacro, fffff800, finnis, fit1d, fit1p, fit2d, fit2p, forplan, ganges, gfrd-pnc, greenbea, greenbeb, grow15, grow22, grow7, israel, kb2, ken-07, ken-11, ken-13, ken-18, lotfi, maros, maros-r7, modszk1, nesm, osa-07, osa-14, osa-30, osa-60, pds-02, pds-06, pds-10, pds-20, perold, pilot, pilot.ja, pilot.we, pilot4, pilot87, pilotnov, qap12, qap15, qap8, recipe, sc105, sc205, sc50a, sc50b, scagr25, scagr7, scfxm1, scfxm2, scfxm3, scorpion, scrs8, scsd1, scsd6, scsd8, sctap1, sctap2, sctap3, seba, share1b, share2b, shell, ship04l, ship04s, ship08l, ship08s, ship12l, ship12s, sierra, stair, standata, standgub, standmps, stocfor1, stocfor2, stocfor3, truss, tuff, vtp.base, wood1p, woodw |
| 系统 / CPU | macOS-27.0.1-arm64-arm-64bit-Mach-O / Apple M1 |
| 内存 / 逻辑 CPU | 17179869184 bytes / 8 |
| Python / NumPy / SciPy | 3.14.2 / 2.5.3 / 1.18.1 |
| lp_simplex 版本 | 8ec4f0dc66f059824b0b38b36fac7632e062114f；工作区有改动 |
| 构建 / 编译器 | Release / Apple clang version 21.0.0 (clang-2100.3.34.2) |
| SCIP / LP 后端 | SCIP version 10.0.3 [precision: 8 byte] [memory: block] [mode: optimized] [LP solver: SoPlex 8.0.2] [GitHash: d409edf9f6] |
| HiGHS | 1.12.0 / scipy-private-binding |
| 算法（lp / scip / highs） | auto / default / choose |
| 重复 / 单次限时 / 容差 / 工作预算 | 1 / 90 秒 / 1e-07 / 300000 |
| 资源 | 串行；每个求解器单线程；未限制内存，不保证机器独占 |
| 原始记录 | tmp/benchmarks/netlib-auto-current-90s-once（配置、环境、SHA-256、命令、日志、解和审核残差） |
| Clp | Coin LP version 1.17.11, build Mar 11 2026 |
| GLPK | GLPSOL--GLPK LP/MIP Solver 5.0 |
| Clp / GLPK 算法 | dual simplex / primal simplex；presolve on |

| 求解器 | 成功案例 | 超时次数 | 其他未通过次数 | 与 lp_simplex 共同成功数 | 共同案例几何平均/秒（lp_simplex / 本求解器） | 相对耗时（lp_simplex / 各求解器） |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| lp_simplex | 114/114 | 0 | 0 | — | — | — |
| SCIP / SoPlex | 109/114 | 4 | 1 | 109 | 0.113925 / 0.205502 | 0.554373 |
| HiGHS | 112/114 | 2 | 0 | 112 | 0.123357 / 0.0346396 | 3.56117 |
| Clp | 111/114 | 1 | 2 | 111 | 0.129037 / 0.0618844 | 2.08513 |
| GLPK | 110/114 | 2 | 2 | 110 | 0.122373 / 0.0797 | 1.53543 |

复现本轮案例与参数（新建运行目录）：

```sh
python3 data/benchmark.py --dataset netlib --cases 25fv47 80bau3b adlittle afiro agg agg2 agg3 bandm beaconfd blend bnl1 bnl2 boeing1 boeing2 bore3d brandy capri cre-a cre-b cre-c cre-d cycle czprob d2q06c d6cube degen2 degen3 dfl001 e226 etamacro fffff800 finnis fit1d fit1p fit2d fit2p forplan ganges gfrd-pnc greenbea greenbeb grow15 grow22 grow7 israel kb2 ken-07 ken-11 ken-13 ken-18 lotfi maros-r7 maros modszk1 nesm osa-07 osa-14 osa-30 osa-60 pds-02 pds-06 pds-10 pds-20 perold pilot.ja pilot pilot.we pilot4 pilot87 pilotnov qap12 qap15 qap8 recipe sc105 sc205 sc50a sc50b scagr25 scagr7 scfxm1 scfxm2 scfxm3 scorpion scrs8 scsd1 scsd6 scsd8 sctap1 sctap2 sctap3 seba share1b share2b shell ship04l ship04s ship08l ship08s ship12l ship12s sierra stair standata standgub standmps stocfor1 stocfor2 stocfor3 truss tuff vtp.base wood1p woodw --binary /Users/zhuanglinsheng/Documents/Github/lp_simplex/tmp/netlib-90s-snapshot/build/data/netlib/test_netlib --scip /opt/homebrew/bin/scip --clp /opt/homebrew/bin/clp --glpk /opt/homebrew/bin/glpsol --repeats 1 --seconds 90.0 --iterations 300000 --lp-algorithm auto --scip-algorithm default --highs-algorithm choose --solvers lp_simplex scip highs clp glpk
```
<!-- /benchmark-run -->

## PDS

[逐案例结果与完整报告](pds.md)

仅展示最近一次报告更新对应的运行；历史记录保存在各轮原始运行目录中。

<!-- benchmark-run:pds-auto-refined-five-solvers-300s-final -->
| 测试信息 | 值 |
| --- | --- |
| 状态 / 进度 | 已完成；9/9 个案例 |
| 时间（UTC） | 2026-10-23T15:41:47.536059+00:00 |
| 参与求解器 | lp_simplex, SCIP / SoPlex, HiGHS, Clp, GLPK |
| 案例范围 | pds 全集；pds-20, pds-30, pds-40, pds-50, pds-60, pds-70, pds-80, pds-90, pds-100 |
| 系统 / CPU | macOS-27.0.1-arm64-arm-64bit-Mach-O / Apple M1 |
| 内存 / 逻辑 CPU | 17179869184 bytes / 8 |
| Python / NumPy / SciPy | 3.14.2 / 2.5.3 / 1.18.1 |
| lp_simplex 版本 | 1ab947ebcbbbdff01b10b99994d58bc1a7cdd3e6；工作区有改动 |
| 构建 / 编译器 | Release / Apple clang version 21.0.0 (clang-2100.3.34.2) |
| SCIP / LP 后端 | SCIP version 10.0.3 [precision: 8 byte] [memory: block] [mode: optimized] [LP solver: SoPlex 8.0.2] [GitHash: d409edf9f6] |
| HiGHS | 1.12.0 / scipy-private-binding |
| 算法（lp / scip / highs） | auto / default / choose |
| 重复 / 单次限时 / 容差 / 工作预算 | 1 / 300 秒 / 1e-07 / 300000 |
| 资源 | 串行；每个求解器单线程；未限制内存，不保证机器独占 |
| 原始记录 | tmp/benchmarks/pds-auto-refined-five-solvers-300s-final（配置、环境、SHA-256、命令、日志、解和审核残差） |
| Clp | Coin LP version 1.17.11, build Mar 11 2026 |
| GLPK | GLPSOL--GLPK LP/MIP Solver 5.0 |
| Clp / GLPK 算法 | dual simplex / primal simplex；presolve on |

| 求解器 | 成功案例 | 超时次数 | 其他未通过次数 | 与 lp_simplex 共同成功数 | 共同案例几何平均/秒（lp_simplex / 本求解器） | 相对耗时（lp_simplex / 各求解器） |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| lp_simplex | 9/9 | 0 | 0 | — | — | — |
| SCIP / SoPlex | 6/9 | 3 | 0 | 6 | 106.83 / 69.9679 | 1.52685 |
| HiGHS | 9/9 | 0 | 0 | 9 | 143.535 / 7.30131 | 19.6588 |
| Clp | 9/9 | 0 | 0 | 9 | 143.535 / 7.31296 | 19.6275 |
| GLPK | 2/9 | 7 | 0 | 2 | 64.4913 / 53.6572 | 1.20191 |

复现本轮案例与参数（新建运行目录）：

```sh
python3 data/benchmark.py --dataset pds --cases pds-100 pds-20 pds-30 pds-40 pds-50 pds-60 pds-70 pds-80 pds-90 --binary /Users/zhuanglinsheng/Documents/Github/lp_simplex/tmp/build-auto-implementation/data/netlib/test_netlib --scip /opt/homebrew/bin/scip --clp /opt/homebrew/bin/clp --glpk /opt/homebrew/bin/glpsol --repeats 1 --seconds 300.0 --iterations 300000 --lp-algorithm auto --scip-algorithm default --highs-algorithm choose --solvers lp_simplex scip highs clp glpk
```
<!-- /benchmark-run -->

## Railway 连续 LP 松弛

[逐案例结果与完整报告](railway.md)

仅展示最近一次报告更新对应的运行；历史记录保存在各轮原始运行目录中。

尚未运行；环境、版本、算法、限时和重复次数尚未记录。

| 求解器 | 成功案例 | 超时次数 | 其他未通过次数 | 与 lp_simplex 共同成功数 | 共同案例几何平均/秒（lp_simplex / 本求解器） | 相对耗时（lp_simplex / 各求解器） |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| lp_simplex | 0/7 | 0 | 0 | — | — | — |
| SCIP / SoPlex | 0/7 | 0 | 0 | 0 | — / — | — |
| HiGHS | 0/7 | 0 | 0 | 0 | — / — | — |
| Clp | 0/7 | 0 | 0 | 0 | — / — | — |
| GLPK | 0/7 | 0 | 0 | 0 | — / — | — |

## 复现说明

在仓库根目录运行。需 CMake、C 编译器、BLAS/LAPACK，并安装 SCIP、Clp、GLPK、NumPy、SciPy，以及 highspy 或兼容的 SciPy HiGHS 绑定。PDS 转换器由脚本临时编译。

```sh
cmake -S . -B build/benchmark -DCMAKE_BUILD_TYPE=Release -DLP_SIMPLEX_BUILD_BENCHMARKS=ON
cmake --build build/benchmark -j 4
python3 data/benchmark.py --dataset netlib --binary build/benchmark/data/netlib/test_netlib
python3 data/benchmark.py --dataset pds --binary build/benchmark/data/netlib/test_netlib
python3 data/benchmark.py --dataset railway --binary build/benchmark/data/netlib/test_netlib
```

省略 `--cases` 测试所选数据集全集；短程验证可加 `--cases <实例名> --repeats 1 --seconds 10`。simple_examples 仅用于功能测试，不参与基准。上方各轮的复现命令保留其案例与具体参数，二进制路径可换为新构建结果；复核版本及 SHA-256 后再比较耗时。

默认每轮新建 `tmp/benchmarks/<UTC 时间>/`，保存配置、环境、命令、结果、日志与审核信息。正常结束、异常、Ctrl+C 或 SIGTERM 清理临时 MPS；SIGKILL/断电残留可随运行目录删除。

仅重建报告和本汇总，不需要 MPS 或再次求解：

```sh
python3 data/benchmark.py --report-only tmp/benchmarks/<运行目录>
```

`--docs-dir` 可指定报告目录；重建汇总会保留其他数据集已保存的测试信息和统计。更多参数与审核说明见 [数据与基准说明](../../data/README.md)。
