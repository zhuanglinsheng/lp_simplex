> 历史测试记录：2026-10-26（北京时间），AUTO 复查前的运行，包含未通过或超时案例，不代表当前优化已通过验收。保留原数值表和 UTC 时间；仅归档 Markdown 报告，下文引用的原始运行目录已清理。[当前报告](../../README.md) · [历史索引](../README.md)

# Netlib infeasible 基准测试

简体中文 | [English](./netlib-infeasible_en.md) | [项目主页](../../../../README.md) | [基准测试汇总](./README.md)

由 `data/benchmark.py` 逐案例生成 MPS、执行求解器对比测试、清理 MPS 后更新。

## 测试与统计口径

- 计时为读取、求解和解导出的墙钟时间；不含 MPS 生成/删除、格式兼容处理、适配器公共初始化和独立审核。CLI 含进程启动，HiGHS 使用 Python binding，小案例开销有差异。
- 成功案例须全部重复为 PASS 且每次均在限时内；先取各案例重复耗时中位数，再计算普通几何平均（不移位、不对失败加惩罚）。空集合显示 `—`。
- 每个求解器分别与 lp_simplex 比较，仅使用双方共同成功的案例。几何平均耗时按“lp_simplex / 本求解器”列出；各行集合可能不同。相对耗时为 lp_simplex 耗时除以本求解器耗时，大于 1 表示本求解器更快；lp_simplex 行不列自比较统计。
- 成功数按案例计，超时与其他未通过数按重复运行计；统计仅属于对应运行的明确案例范围，进行中结果不构成全集结论。
- NOT_RUN / PENDING / INCOMPLETE 表示从未测试 / 本轮尚未测试 / 未完成全部重复。逐例单元格为“状态 / 秒数中位数 / 复算目标中位数”；非 PASS 的目标值仍须结合状态阅读。
- SKIPPED 表示按预设规则跳过，不计为超时或其他未通过，也不参与耗时统计。
- 未选中的案例及求解器单元格保留历史结果，仅作参考；上方统计只包含本次参与求解器的实测记录，不将历史单元格计入本次排名。

- PASS 表示明确报告不可行且在限时内结束；预期状态来自 Netlib infeasible 分类。未独立验证 Farkas 证书，PASS 不代表证书审核通过。
- WRONG_STATUS 表示报告最优或无界；AMBIGUOUS 表示不可行或无界；SOLVER_FAILED / ERROR 表示未得到明确判定或运行失败。TIMEOUT 优先于其他状态。
- 原始 MPS 按原样交给求解器。greenbea 含多个 BOUNDS 向量及重复界，pang 存在读取器兼容性差异；相关失败包含输入兼容性影响。HiGHS 的读取警告保存在 reader_status，实际求解日志另行保留。

## 测试信息与统计

仅展示最近一次报告更新对应的运行；历史记录保存在各轮原始运行目录中。

<!-- benchmark-run:netlib-infeasible -->
| 测试信息 | 值 |
| --- | --- |
| 状态 / 进度 | 已完成；29/29 个案例 |
| 时间（UTC） | 2026-10-25T22:20:02.419269+00:00 |
| 参与求解器 | lp_simplex, SCIP / SoPlex, HiGHS, Clp, GLPK |
| 案例范围 | netlib-infeasible 全集；bgdbg1, bgetam, bgindy, bgprtr, box1, ceria3d, chemcom, cplex1, cplex2, ex72a, ex73a, forest6, galenet, gosh, gran, greenbea, itest2, itest6, klein1, klein2, klein3, mondou2, pang, pilot4i, qual, reactor, refinery, vol1, woodinfe |
| 系统 / CPU | macOS-27.0.1-arm64-arm-64bit-Mach-O / Apple M1 |
| 内存 / 逻辑 CPU | 17179869184 bytes / 8 |
| Python / NumPy / SciPy | 3.14.2 / 2.5.3 / 1.18.1 |
| lp_simplex 版本 | 4115c34c078d6942394cc225cc33645e18006d8b；工作区有改动 |
| 构建 / 编译器 | Release / Apple clang version 21.0.0 (clang-2100.3.34.2) |
| SCIP / LP 后端 | SCIP version 10.0.3 [precision: 8 byte] [memory: block] [mode: optimized] [LP solver: SoPlex 8.0.2] [GitHash: d409edf9f6] |
| HiGHS | 1.12.0 / scipy-private-binding |
| 算法（lp / scip / highs） | auto / default / choose |
| 重复 / 单次限时 / 容差 / 工作预算 | 1 / 90 秒 / 1e-07 / 300000 |
| 资源 | 串行；每个求解器单线程；未限制内存，不保证机器独占 |
| 原始记录 | tmp/benchmarks/2026-10-26/netlib-infeasible（配置、环境、SHA-256、命令、日志、解和审核残差） |
| Clp | Coin LP version 1.17.11, build Mar 11 2026 |
| GLPK | GLPSOL--GLPK LP/MIP Solver 5.0 |
| Clp / GLPK 算法 | dual simplex / primal simplex；presolve on |
| GLPK 跳过策略 | 任一案例超时后跳过后续测试；已跳过 0 次 |

| 求解器 | 成功案例 | 超时次数 | 其他未通过次数 | 与 lp_simplex 共同成功数 | 共同案例几何平均/秒（lp_simplex / 本求解器） | 相对耗时（lp_simplex / 各求解器） |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| lp_simplex | 29/29 | 0 | 0 | — | — | — |
| SCIP / SoPlex | 28/29 | 0 | 1 | 28 | 0.0211895 / 0.0425177 | 0.498368 |
| HiGHS | 29/29 | 0 | 0 | 29 | 0.0210961 / 0.00311983 | 6.76195 |
| Clp | 29/29 | 0 | 0 | 29 | 0.0210961 / 0.0242736 | 0.869098 |
| GLPK | 27/29 | 0 | 2 | 27 | 0.0208329 / 0.0135933 | 1.53259 |

| 求解器 | 错误判定次数 | 判定不明确次数 | 运行或数值失败次数 | 跳过次数 |
| --- | ---: | ---: | ---: | ---: |
| lp_simplex | 0 | 0 | 0 | 0 |
| SCIP / SoPlex | 0 | 1 | 0 | 0 |
| HiGHS | 0 | 0 | 0 | 0 |
| Clp | 0 | 0 | 0 | 0 |
| GLPK | 0 | 0 | 2 | 0 |

| 实例 | 求解器 | 结果 | 终止信息 |
| --- | --- | --- | --- |
| cplex1 | SCIP / SoPlex | AMBIGUOUS | problem is solved [infeasible or unbounded] |
| greenbea | GLPK | SOLVER_FAILED | not optimal |
| pang | GLPK | SOLVER_FAILED | not optimal |

复现本轮案例与参数（新建运行目录）：

```sh
cmake -S . -B build/benchmark -DCMAKE_BUILD_TYPE=Release -DLP_SIMPLEX_BUILD_BENCHMARKS=ON

cmake --build build/benchmark --target lp_solve_mps -j 4

python3 data/benchmark.py --dataset netlib-infeasible --cases bgdbg1 bgetam bgindy bgprtr box1 ceria3d chemcom cplex1 cplex2 ex72a ex73a forest6 galenet gosh gran greenbea itest2 itest6 klein1 klein2 klein3 mondou2 pang pilot4i qual reactor refinery vol1 woodinfe --binary build/benchmark/data/lp_solve_mps --scip scip --clp clp --glpk glpsol --repeats 1 --seconds 90.0 --iterations 300000 --lp-algorithm auto --scip-algorithm default --highs-algorithm choose --solvers lp_simplex scip highs clp glpk
```
<!-- /benchmark-run -->

## 逐例最新结果

| 实例 | 规模（行 / 列 / 非零元） | 预期状态 | lp_simplex（状态 / 秒 / 目标） | SCIP / SoPlex（状态 / 秒 / 目标） | HiGHS（状态 / 秒 / 目标） | Clp（状态 / 秒 / 目标） | GLPK（状态 / 秒 / 目标） |
| --- | --- | ---: | --- | --- | --- | --- | --- |
| bgdbg1 | 348 / 407 / 1440 | INFEASIBLE | PASS / 0.014 / — | PASS / 0.034 / — | PASS / 0.001 / — | PASS / 0.012 / — | PASS / 0.008 / — |
| bgetam | 400 / 688 / 2409 | INFEASIBLE | PASS / 0.014 / — | PASS / 0.036 / — | PASS / 0.003 / — | PASS / 0.013 / — | PASS / 0.013 / — |
| bgindy | 2671 / 10116 / 65502 | INFEASIBLE | PASS / 0.037 / — | PASS / 0.076 / — | PASS / 0.038 / — | PASS / 0.037 / — | PASS / 0.126 / — |
| bgprtr | 20 / 34 / 64 | INFEASIBLE | PASS / 0.013 / — | PASS / 0.033 / — | PASS / 0.001 / — | PASS / 0.012 / — | PASS / 0.007 / — |
| box1 | 231 / 261 / 651 | INFEASIBLE | PASS / 0.014 / — | PASS / 0.035 / — | PASS / 0.001 / — | PASS / 0.012 / — | PASS / 0.007 / — |
| ceria3d | 3576 / 824 / 17602 | INFEASIBLE | PASS / 0.019 / — | PASS / 0.046 / — | PASS / 0.011 / — | PASS / 0.133 / — | PASS / 0.018 / — |
| chemcom | 288 / 720 / 1566 | INFEASIBLE | PASS / 0.014 / — | PASS / 0.035 / — | PASS / 0.003 / — | PASS / 0.013 / — | PASS / 0.010 / — |
| cplex1 | 3005 / 3221 / 8944 | INFEASIBLE | PASS / 0.019 / — | AMBIGUOUS / — / — | PASS / 0.087 / — | PASS / 0.090 / — | PASS / 0.044 / — |
| cplex2 | 224 / 221 / 1058 | INFEASIBLE | PASS / 0.022 / — | PASS / 0.050 / — | PASS / 0.005 / — | PASS / 0.016 / — | PASS / 0.011 / — |
| ex72a | 197 / 215 / 467 | INFEASIBLE | PASS / 0.013 / — | PASS / 0.034 / — | PASS / 0.001 / — | PASS / 0.012 / — | PASS / 0.008 / — |
| ex73a | 193 / 211 / 457 | INFEASIBLE | PASS / 0.013 / — | PASS / 0.035 / — | PASS / 0.001 / — | PASS / 0.012 / — | PASS / 0.007 / — |
| forest6 | 66 / 95 / 210 | INFEASIBLE | PASS / 0.013 / — | PASS / 0.035 / — | PASS / 0.001 / — | PASS / 0.013 / — | PASS / 0.007 / — |
| galenet | 8 / 8 / 16 | INFEASIBLE | PASS / 0.013 / — | PASS / 0.032 / — | PASS / 0.000 / — | PASS / 0.011 / — | PASS / 0.006 / — |
| gosh | 3792 / 10733 / 97231 | INFEASIBLE | PASS / 1.453 / — | PASS / 0.097 / — | PASS / 0.050 / — | PASS / 0.721 / — | PASS / 0.327 / — |
| gran | 2658 / 2520 / 20106 | INFEASIBLE | PASS / 0.020 / — | PASS / 0.045 / — | PASS / 0.011 / — | PASS / 0.089 / — | PASS / 0.031 / — |
| greenbea | 2393 / 5405 / 30883 | INFEASIBLE | PASS / 0.029 / — | PASS / 0.058 / — | PASS / 0.018 / — | PASS / 0.168 / — | SOLVER_FAILED / — / — |
| itest2 | 9 / 4 / 17 | INFEASIBLE | PASS / 0.012 / — | PASS / 0.033 / — | PASS / 0.000 / — | PASS / 0.011 / — | PASS / 0.007 / — |
| itest6 | 11 / 8 / 20 | INFEASIBLE | PASS / 0.012 / — | PASS / 0.032 / — | PASS / 0.000 / — | PASS / 0.011 / — | PASS / 0.006 / — |
| klein1 | 54 / 54 / 696 | INFEASIBLE | PASS / 0.020 / — | PASS / 0.036 / — | PASS / 0.002 / — | PASS / 0.012 / — | PASS / 0.008 / — |
| klein2 | 477 / 54 / 4585 | INFEASIBLE | PASS / 0.029 / — | PASS / 0.056 / — | PASS / 0.007 / — | PASS / 0.026 / — | PASS / 0.021 / — |
| klein3 | 994 / 88 / 12107 | INFEASIBLE | PASS / 0.086 / — | PASS / 0.099 / — | PASS / 0.026 / — | PASS / 0.077 / — | PASS / 0.060 / — |
| mondou2 | 312 / 604 / 1208 | INFEASIBLE | PASS / 0.014 / — | PASS / 0.035 / — | PASS / 0.002 / — | PASS / 0.014 / — | PASS / 0.010 / — |
| pang | 361 / 460 / 2652 | INFEASIBLE | PASS / 0.021 / — | PASS / 0.048 / — | PASS / 0.008 / — | PASS / 0.019 / — | SOLVER_FAILED / — / — |
| pilot4i | 410 / 1000 / 5141 | INFEASIBLE | PASS / 0.015 / — | PASS / 0.036 / — | PASS / 0.003 / — | PASS / 0.026 / — | PASS / 0.010 / — |
| qual | 323 / 464 / 1646 | INFEASIBLE | PASS / 0.019 / — | PASS / 0.042 / — | PASS / 0.005 / — | PASS / 0.018 / — | PASS / 0.010 / — |
| reactor | 318 / 637 / 2420 | INFEASIBLE | PASS / 0.014 / — | PASS / 0.035 / — | PASS / 0.002 / — | PASS / 0.015 / — | PASS / 0.009 / — |
| refinery | 323 / 464 / 1626 | INFEASIBLE | PASS / 0.023 / — | PASS / 0.045 / — | PASS / 0.004 / — | PASS / 0.017 / — | PASS / 0.011 / — |
| vol1 | 323 / 464 / 1646 | INFEASIBLE | PASS / 0.028 / — | PASS / 0.047 / — | PASS / 0.008 / — | PASS / 0.022 / — | PASS / 0.011 / — |
| woodinfe | 35 / 89 / 140 | INFEASIBLE | PASS / 0.015 / — | PASS / 0.033 / — | PASS / 0.000 / — | PASS / 0.013 / — | PASS / 0.007 / — |
