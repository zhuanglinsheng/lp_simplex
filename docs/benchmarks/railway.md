# Railway 连续 LP 松弛 基准测试

简体中文 | [English](./railway_en.md) | [项目主页](../../README.md) | [基准测试汇总](./README.md)

由 `data/benchmark.py` 逐案例生成 MPS、执行求解器对比测试、清理 MPS 后更新。

## 测试与统计口径

- 计时为读取、求解和解导出的墙钟时间；不含 MPS 生成/删除、格式兼容处理、适配器公共初始化和独立审核。CLI 含进程启动，HiGHS 使用 Python binding，小案例开销有差异。
- PASS 要求求解器报告最优，原模型行约束、变量界和目标复算通过；有参考目标时再核对参考值（误差限为容差乘以 `1 + abs(reference)`）。无参考值显示 `—`，不代表独立最优性认证。
- 成功案例须全部重复为 PASS 且每次均在限时内；先取各案例重复耗时中位数，再计算普通几何平均（不移位、不对失败加惩罚）。空集合显示 `—`。
- 每个求解器分别与 lp_simplex 比较，仅使用双方共同成功的案例。几何平均耗时按“lp_simplex / 本求解器”列出；各行集合可能不同。相对耗时为 lp_simplex 耗时除以本求解器耗时，大于 1 表示本求解器更快；lp_simplex 行不列自比较统计。
- 成功数按案例计，超时与其他未通过数按重复运行计；统计仅属于对应运行的明确案例范围，进行中结果不构成全集结论。
- NOT_RUN / PENDING / INCOMPLETE 表示从未测试 / 本轮尚未测试 / 未完成全部重复。逐例单元格为“状态 / 秒数中位数 / 复算目标中位数”；非 PASS 的目标值仍须结合状态阅读。
- SKIPPED 表示按预设规则跳过，不计为超时或其他未通过，也不参与耗时统计。
- 未选中的案例及求解器单元格保留历史结果，仅作参考；上方统计只包含本次参与求解器的实测记录，不将历史单元格计入本次排名。

## 测试信息与统计

仅展示最近一次报告更新对应的运行；历史记录保存在各轮原始运行目录中。

<!-- benchmark-run:railway-auto-current-3600s-once -->
| 测试信息 | 值 |
| --- | --- |
| 状态 / 进度 | 已完成；7/7 个案例 |
| 时间（UTC） | 2026-10-23T18:29:21.672935+00:00 |
| 参与求解器 | lp_simplex, SCIP / SoPlex, HiGHS, Clp, GLPK |
| 案例范围 | railway 全集；rail507, rail516, rail582, rail2536, rail2586, rail4284, rail4872 |
| 系统 / CPU | macOS-27.0.1-arm64-arm-64bit-Mach-O / Apple M1 |
| 内存 / 逻辑 CPU | 17179869184 bytes / 8 |
| Python / NumPy / SciPy | 3.14.2 / 2.5.3 / 1.18.1 |
| lp_simplex 版本 | 00efbe80257a0530ee491cb2c44494d13b27c46c；工作区有改动 |
| 构建 / 编译器 | Release / Apple clang version 21.0.0 (clang-2100.3.34.2) |
| SCIP / LP 后端 | SCIP version 10.0.3 [precision: 8 byte] [memory: block] [mode: optimized] [LP solver: SoPlex 8.0.2] [GitHash: d409edf9f6] |
| HiGHS | 1.12.0 / scipy-private-binding |
| 算法（lp / scip / highs） | auto / default / choose |
| 重复 / 单次限时 / 容差 / 工作预算 | 1 / 3600 秒 / 1e-07 / 10000000 |
| 资源 | 串行；每个求解器单线程；未限制内存，不保证机器独占 |
| 原始记录 | tmp/benchmarks/railway-auto-current-3600s-once（配置、环境、SHA-256、命令、日志、解和审核残差） |
| Clp | Coin LP version 1.17.11, build Mar 11 2026 |
| GLPK | GLPSOL--GLPK LP/MIP Solver 5.0 |
| Clp / GLPK 算法 | dual simplex / primal simplex；presolve on |
| GLPK 跳过策略 | 任一案例超时后跳过后续测试；已跳过 3 次 |

| 求解器 | 成功案例 | 超时次数 | 其他未通过次数 | 与 lp_simplex 共同成功数 | 共同案例几何平均/秒（lp_simplex / 本求解器） | 相对耗时（lp_simplex / 各求解器） |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| lp_simplex | 7/7 | 0 | 0 | — | — | — |
| SCIP / SoPlex | 7/7 | 0 | 0 | 7 | 27.7296 / 112.886 | 0.245642 |
| HiGHS | 7/7 | 0 | 0 | 7 | 27.7296 / 82.1388 | 0.337594 |
| Clp | 7/7 | 0 | 0 | 7 | 27.7296 / 50.2903 | 0.55139 |
| GLPK | 3/7 | 1 | 0 | 3 | 1.5798 / 3.45865 | 0.456767 |

复现本轮案例与参数（新建运行目录）：

```sh
cmake -S . -B build/benchmark -DCMAKE_BUILD_TYPE=Release -DLP_SIMPLEX_BUILD_BENCHMARKS=ON

cmake --build build/benchmark --target lp_solve_mps -j 4

python3 data/benchmark.py --dataset railway --cases rail507 rail516 rail582 rail2536 rail2586 rail4284 rail4872 --binary build/benchmark/data/lp_solve_mps --scip scip --clp clp --glpk glpsol --repeats 1 --seconds 3600.0 --iterations 10000000 --lp-algorithm auto --scip-algorithm default --highs-algorithm choose --solvers lp_simplex scip highs clp glpk
```
<!-- /benchmark-run -->

## 逐例最新结果

| 实例 | 规模（行 / 列 / 非零元） | 参考目标 | lp_simplex（状态 / 秒 / 目标） | SCIP / SoPlex（状态 / 秒 / 目标） | HiGHS（状态 / 秒 / 目标） | Clp（状态 / 秒 / 目标） | GLPK（状态 / 秒 / 目标） |
| --- | --- | ---: | --- | --- | --- | --- | --- |
| rail507 | 507 / 63009 / 409349 | — | PASS / 1.346 / 172.145566677 | PASS / 5.936 / 172.145566677 | PASS / 4.206 / 172.145566677 | PASS / 3.806 / 172.145566677 | PASS / 4.002 / 172.145566677 |
| rail516 | 516 / 47311 / 314896 | — | PASS / 1.720 / 182 | PASS / 14.740 / 182 | PASS / 1.365 / 182 | PASS / 0.922 / 182 | PASS / 2.401 / 182 |
| rail582 | 582 / 55515 / 401708 | — | PASS / 1.702 / 209.712232957 | PASS / 7.835 / 209.712232957 | PASS / 4.824 / 209.712232957 | PASS / 2.372 / 209.712232957 | PASS / 4.306 / 209.712232957 |
| rail2536 | 2536 / 1081841 / 10993311 | — | PASS / 115.476 / 688.399223322 | PASS / 862.544 / 688.399223322 | PASS / 1341.432 / 688.399223322 | PASS / 391.993 / 688.399223322 | TIMEOUT / — / — |
| rail2586 | 2586 / 920683 / 8008776 | — | PASS / 153.442 / 935.921820568 | PASS / 475.296 / 935.921820568 | PASS / 589.461 / 935.921820568 | PASS / 402.475 / 935.921820568 | SKIPPED / — / — |
| rail4284 | 4284 / 1092610 / 11279748 | — | PASS / 326.133 / 1054.05449127 | PASS / 1015.431 / 1054.05449127 | PASS / 1155.374 / 1054.05449127 | PASS / 815.786 / 1054.05449127 | SKIPPED / — / — |
| rail4872 | 4872 / 968672 / 9244093 | — | PASS / 553.309 / 1509.63940162 | PASS / 818.543 / 1509.63940162 | PASS / 997.292 / 1509.63940162 | PASS / 759.785 / 1509.63940162 | SKIPPED / — / — |
