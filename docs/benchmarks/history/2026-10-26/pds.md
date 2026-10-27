> 历史测试记录：2026-10-26（北京时间），AUTO 复查前的运行，包含未通过或超时案例，不代表当前优化已通过验收。保留原数值表和 UTC 时间；仅归档 Markdown 报告，下文引用的原始运行目录已清理。[当前报告](../../README.md) · [历史索引](../README.md)

# PDS 基准测试

简体中文 | [English](./pds_en.md) | [项目主页](../../../../README.md) | [基准测试汇总](./README.md)

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

<!-- benchmark-run:pds -->
| 测试信息 | 值 |
| --- | --- |
| 状态 / 进度 | 已完成；9/9 个案例 |
| 时间（UTC） | 2026-10-25T22:21:06.080695+00:00 |
| 参与求解器 | lp_simplex, SCIP / SoPlex, HiGHS, Clp, GLPK |
| 案例范围 | pds 全集；pds-20, pds-30, pds-40, pds-50, pds-60, pds-70, pds-80, pds-90, pds-100 |
| 系统 / CPU | macOS-27.0.1-arm64-arm-64bit-Mach-O / Apple M1 |
| 内存 / 逻辑 CPU | 17179869184 bytes / 8 |
| Python / NumPy / SciPy | 3.14.2 / 2.5.3 / 1.18.1 |
| lp_simplex 版本 | 4115c34c078d6942394cc225cc33645e18006d8b；工作区有改动 |
| 构建 / 编译器 | Release / Apple clang version 21.0.0 (clang-2100.3.34.2) |
| SCIP / LP 后端 | SCIP version 10.0.3 [precision: 8 byte] [memory: block] [mode: optimized] [LP solver: SoPlex 8.0.2] [GitHash: d409edf9f6] |
| HiGHS | 1.12.0 / scipy-private-binding |
| 算法（lp / scip / highs） | auto / default / choose |
| 重复 / 单次限时 / 容差 / 工作预算 | 1 / 300 秒 / 1e-07 / 300000 |
| 资源 | 串行；每个求解器单线程；未限制内存，不保证机器独占 |
| 原始记录 | tmp/benchmarks/2026-10-26/pds（配置、环境、SHA-256、命令、日志、解和审核残差） |
| Clp | Coin LP version 1.17.11, build Mar 11 2026 |
| GLPK | GLPSOL--GLPK LP/MIP Solver 5.0 |
| Clp / GLPK 算法 | dual simplex / primal simplex；presolve on |
| GLPK 跳过策略 | 任一案例超时后跳过后续测试；已跳过 6 次 |

| 求解器 | 成功案例 | 超时次数 | 其他未通过次数 | 与 lp_simplex 共同成功数 | 共同案例几何平均/秒（lp_simplex / 本求解器） | 相对耗时（lp_simplex / 各求解器） |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| lp_simplex | 5/9 | 4 | 0 | — | — | — |
| SCIP / SoPlex | 6/9 | 3 | 0 | 5 | 127.25 / 54.4635 | 2.33643 |
| HiGHS | 9/9 | 0 | 0 | 5 | 127.25 / 3.65749 | 34.7916 |
| Clp | 9/9 | 0 | 0 | 5 | 127.25 / 3.58432 | 35.5018 |
| GLPK | 2/9 | 1 | 0 | 2 | 66.2865 / 49.4924 | 1.33933 |

复现本轮案例与参数（新建运行目录）：

```sh
cmake -S . -B build/benchmark -DCMAKE_BUILD_TYPE=Release -DLP_SIMPLEX_BUILD_BENCHMARKS=ON

cmake --build build/benchmark --target lp_solve_mps -j 4

python3 data/benchmark.py --dataset pds --cases pds-20 pds-30 pds-40 pds-50 pds-60 pds-70 pds-80 pds-90 pds-100 --binary build/benchmark/data/lp_solve_mps --scip scip --clp clp --glpk glpsol --repeats 1 --seconds 300.0 --iterations 300000 --lp-algorithm auto --scip-algorithm default --highs-algorithm choose --solvers lp_simplex scip highs clp glpk
```
<!-- /benchmark-run -->

## 逐例最新结果

| 实例 | 规模（行 / 列 / 非零元） | 参考目标 | lp_simplex（状态 / 秒 / 目标） | SCIP / SoPlex（状态 / 秒 / 目标） | HiGHS（状态 / 秒 / 目标） | Clp（状态 / 秒 / 目标） | GLPK（状态 / 秒 / 目标） |
| --- | --- | ---: | --- | --- | --- | --- | --- |
| pds-20 | 33874 / 105728 / 230200 | 23821658640 | PASS / 48.968 / 23821658640 | PASS / 8.267 / 23821658640 | PASS / 1.281 / 23821658640 | PASS / 0.895 / 23821658640 | PASS / 25.575 / 23821658640 |
| pds-30 | 49944 / 154998 / 337144 | 21385445736 | PASS / 89.730 / 21385445736 | PASS / 25.784 / 21385445736 | PASS / 2.388 / 21385445736 | PASS / 2.088 / 21385445736 | PASS / 95.776 / 21385445736 |
| pds-40 | 66844 / 212859 / 462128 | 18855198824 | PASS / 149.349 / 18855198824 | PASS / 70.382 / 18855198824 | PASS / 4.168 / 18855198824 | PASS / 4.331 / 18855198824 | TIMEOUT / — / — |
| pds-50 | 83060 / 270095 / 585114 | 16603525724 | PASS / 220.027 / 16603525724 | PASS / 152.372 / 16603525724 | PASS / 6.803 / 16603525724 | PASS / 7.856 / 16603525724 | SKIPPED / — / — |
| pds-60 | 99431 / 329643 / 712779 | 14265904407 | PASS / 231.078 / 14265904407 | PASS / 209.629 / 14265904407 | PASS / 7.544 / 14265904407 | PASS / 9.300 / 14265904407 | SKIPPED / — / — |
| pds-70 | 114944 / 382311 / 825771 | 12241162812 | TIMEOUT / — / — | PASS / 253.354 / 12241162812 | PASS / 11.645 / 12241162812 | PASS / 16.601 / 12241162812 | SKIPPED / — / — |
| pds-80 | 129181 / 426278 / 919524 | 11469077462 | TIMEOUT / — / — | TIMEOUT / — / — | PASS / 16.834 / 11469077462 | PASS / 18.394 / 11469077462 | SKIPPED / — / — |
| pds-90 | 142823 / 466671 / 1005359 | 11087561635 | TIMEOUT / — / — | TIMEOUT / — / — | PASS / 16.376 / 11087561635 | PASS / 16.391 / 11087561635 | SKIPPED / — / — |
| pds-100 | 156243 / 505360 / 1086785 | 10928229968 | TIMEOUT / — / — | TIMEOUT / — / — | PASS / 21.527 / 10928229968 | PASS / 18.424 / 10928229968 | SKIPPED / — / — |
