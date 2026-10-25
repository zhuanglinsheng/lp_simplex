# 数据与基准测试

简体中文 | [English](./README_en.md) | [项目主页](../README.md) | [文档](../docs/README.md)

本目录维护测试数据、按需生成 MPS 的脚本，以及 lp_simplex、SCIP/SoPlex、HiGHS、Clp 和 GLPK 的性能对比入口。以下命令均在仓库根目录执行。

## 数据集

| 数据集 | 内容 | 性能报告 |
|---|---|---|
| [Netlib](netlib/README.md) | 114 个可行 LP，另有不可行模型和预处理版本；此组仅测试可行原模型 | [逐例结果](../docs/benchmarks/netlib.md) |
| Netlib infeasible | 29 个不可行 LP，独立统计不可行性判定 | [逐例结果](../docs/benchmarks/netlib-infeasible.md) |
| [PDS](pds/README.md) | 9 个大规模 LP，约 10.6万—50.5万变量 | [逐例结果](../docs/benchmarks/pds.md) |
| [Railway](railway/README.md) | 7 个铁路集合覆盖问题，测试其连续 LP 松弛，约 4.7万—109.3万变量 | [逐例结果](../docs/benchmarks/railway.md) |
| [小型示例](simple_examples/README.md) | 9 个功能与回归测试模型 | 不纳入性能基准 |

Netlib、PDS 和 Railway 只保存 raw 数据及转换所需文件，各目录的 `prepare.py` 按需生成 MPS。`simple_examples` 直接保存小型 MPS 和配套 C 示例；如何建模和调用求解器见[使用指南](../docs/usage.md)。

命令行默认采用正式数值配置，不受数值 `LP_SIMPLEX_*` 实验开关影响。复现环境实验时须显式添加 `--experimental-environment`；该标志只用于单模型驱动，不作为默认基准配置。详见[正式配置与实验环境边界](../docs/runtime-configuration.md)。

## 计时口径

`benchmark.py` 保留原有 `wall_seconds`：CLI 的启动、读取、求解、点导出，以及适配器读取日志/点文件；不含适配器初始化或独立原模型审核。CLI 运行另记录 `cli_wall_seconds`（启动至退出，不含之后的解析）、`cli_user_seconds`、`cli_system_seconds` 及二者之和 `cli_cpu_seconds`。CPU 来自隔离 worker 中该 CLI 及其子进程的 `getrusage` 差值；平台不支持或计数无效时为 `null`，不能按零处理。HiGHS Python 绑定没有 CLI 子进程，不提供这些字段。

`compare_performance.py` 的原汇总仍使用 `wall_seconds`，新增 `cli_timing` 分别汇总上述四项；只对双方重复次数完整、全部通过且度量有效的案例计算比例。数值 PASS 与计时可用性分开。驱动器日志的 `time:` 是求解加一次原模型审核的进程 CPU，范围不同，不代替完整 CLI 计时。硬截止与原有成功判据保持不变；子进程被杀且没有完整测量时，CPU 保持未知。原始日志、点审核和每次数据均应保留。

## 准备运行环境

使用 CMake、C 编译器和 BLAS/LAPACK 构建单模型测试程序：

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DLP_SIMPLEX_BUILD_BENCHMARKS=ON

cmake --build build --target lp_solve_mps -j 4
```

构建产物为 `build/data/lp_solve_mps`。单独运行它只需准备好 MPS；从 raw 生成模型还需要 Python 3.9+，PDS 转换器由 Python 调用 C 编译器临时构建。

仅当选中 SCIP、Clp 或 GLPK 时，才需安装对应可执行文件，并确保它位于 PATH 中（也可用 `--scip`、`--clp` 或 `--glpk` 指定路径）。Python 环境需包含 NumPy、SciPy，以及 highspy 或兼容的 SciPy HiGHS 绑定。例如：

macOS 可用 `brew install clp glpk` 安装 Clp 和 GLPK。

```sh
python3 -m pip install numpy scipy highspy
```

## 运行完整基准

PDS 和 Railway 按案例编号从小到大执行，作为由小规模到大规模的测试顺序。GLPK 一旦超时，后续 GLPK 测试均标记为 `SKIPPED`；其他求解器继续运行。跳过不计为超时或失败，不参与耗时统计。

分别运行三个数据集的全部基准案例：

```sh
python3 data/benchmark.py --dataset netlib
python3 data/benchmark.py --dataset pds
python3 data/benchmark.py --dataset railway
```

默认只运行 lp_simplex，每个案例测试 1 次；Netlib 每次限时 60 秒，PDS 和 Railway 每次限时 120 秒。Netlib 仅测试 114 个可行原模型；`simple_examples` 不参与性能基准。自定义构建目录时，用 `--binary PATH` 指定 `lp_solve_mps`。

脚本逐案例执行以下流程：

1. 从 raw 数据生成本案例的临时 MPS。
2. 串行运行本轮选中的求解器及全部重复，记录耗时、状态并审核原模型上的解。
3. 删除本案例的临时 MPS。
4. 保存运行记录，同步更新数据集报告和中英文汇总，再处理下一案例。

基准不需要预先生成整个数据集，临时 MPS 位于本轮运行目录内。手动生成的 MPS 不参与这套清理流程。各数据集的报告维护逐案例最新结果，测试信息和统计按轮次保留。

调试单个案例时，可以减少重复次数和限时：

```sh
python3 data/benchmark.py --dataset netlib --cases afiro --repeats 1 --seconds 10
```

`--cases` 接收案例名，例如 Netlib 的 `afiro`、PDS 的 `pds-30`、Railway 的 `rail507`，不带目录和文件后缀。该选项仅限制本轮测试范围。

中断后继续测试时，保留旧运行目录，用 `--cases` 选择尚未完成的案例并启动新一轮运行；沿用原算法、求解器、预算和限时，且 `--run-dir` 必须使用新目录。已完成且未选中的案例行会保留，各轮统计分别记录。脚本没有原目录断点续跑功能，也不会自动跳过案例内已经完成的求解器或重复；中断的案例会从头运行。`--report-only` 只重建报告，不继续求解。

## 单独检查一个模型

[`lp_solve_mps.c`](lp_solve_mps.c) 用于单模型求解和预期结果校验。例如，先生成 Afiro，再求解并打印变量：

```sh
python3 data/netlib/prepare.py feasible/afiro

./build/data/lp_solve_mps --expected -464.7531429 --print-solution data/netlib/feasible/afiro.mps
```

这里 `prepare.py` 需要类别前缀 `feasible/`，以区分 Netlib 中不同类别的同名模型。生成的 MPS 会保留在数据目录中；其他数据集的生成方式见各自 README。

该驱动通用于 Netlib、PDS、Railway 及其他 MPS 模型，不按模型名或目录推断预期结果。用 `--expected VALUE` 指定参考目标，或用 `--infeasible` 预期不可行；未指定时只要求求解成功。数据准备、参考目标选择和独立审核由 `benchmark.py` 管理。

退出码 0 表示匹配预期结果，1 表示读取、求解或结果校验失败，2 表示用法错误。因此，预期不可行且验证通过时也返回 0。

## 查看结果与重建报告

[基准测试汇总](../docs/benchmarks/README.md)包含测试环境、对比统计、统计口径和复现命令，并链接三个数据集的逐例报告。PASS 要求求解器报告最优，且原模型可行性、目标复算及可用参考目标的核对通过；它不等同于完整的独立最优性证明。

每轮配置、环境、结果、日志和解向量保存在 `tmp/benchmarks/<运行目录>/`，由 Git 忽略。需要长期留存时请另行备份。`--run-dir` 可指定新的运行目录，`--docs-dir` 可指定报告目录。

从已有运行记录重建报告时，无需重新求解。将下面路径替换为实际运行目录：

```sh
python3 data/benchmark.py --report-only tmp/benchmarks/<运行目录>
```

脚本同步更新对应数据集的中英文报告和汇总；不完整的运行会显示尚未完成的结果。

## 命令行参数参考

### 求解器基准：benchmark.py

| 参数 | 默认值 | 说明 |
| --- | --- | --- |
| `--dataset` | `netlib` | `netlib`、`pds` 或 `railway`；不包含 simple_examples |
| `--cases` | 所选数据集全集 | 一个或多个 raw 案例名 |
| `--binary` | `build/data/lp_solve_mps` | lp_simplex 单模型驱动路径 |
| `--scip` | `scip` | SCIP 可执行文件名称或路径 |
| `--clp` | `clp` | Clp 可执行文件名称或路径 |
| `--glpk` | `glpsol` | GLPK 可执行文件名称或路径 |
| `--glpk-gate-case` | 不启用 | 先测试指定案例；若 GLPK 未通过，则跳过其余案例的 GLPK 测试。跳过单独记录，不计为超时或失败 |
| `--solvers` | `lp_simplex` | 选择本轮执行的一个或多个求解器，可单独测试其他求解器；仅检查选中求解器的可执行文件 |
| `--repeats` | 1 | 每个案例、每个求解器的重复次数 |
| `--seconds` | Netlib: 60；PDS/Railway: 120 | 每次运行的墙钟秒数上限 |
| `--iterations` | 300000 | lp_simplex 工作预算 |
| `--lp-algorithm` | `auto` | `auto`、`dual-revised`、`pan-bda`、`tableau` |
| `--highs-algorithm` | `choose` | `choose`、`simplex`、`ipm` |
| `--scip-algorithm` | `default` | `default` 或 `dual` |
| `--run-dir` | `tmp/benchmarks/<UTC 时间>` | 必须是新目录，保存本轮原始记录 |
| `--docs-dir` | `docs/benchmarks` | 报告输出目录 |
| `--report-only DIR` | 不启用 | 从已有记录重建报告，不重新求解 |

日常优化默认只测 lp_simplex，每例一次。需要刷新对照时通过 `--solvers lp_simplex highs clp` 等参数选择；也可用 `--solvers scip` 单独刷新其他求解器。未选中的报告单元格保留历史参考值，本次统计不包含这些历史值。HiGHS binding、NumPy 和 SciPy 仍用于读取模型及独立审核，即使不测试 HiGHS 也需要安装。

### 单模型检查：lp_solve_mps

```sh
./build/data/lp_solve_mps [选项] MODEL.mps
```

| 参数 | 默认值 | 说明 |
| --- | --- | --- |
| `--algorithm` | `auto` | `auto`、`dual-revised`、`pan-bda`、`tableau` |
| `--criteria` | `dantzig` | 手动 tableau 使用 `bland/dantzig`，Pan 使用 `dantzig/normalized`；AUTO/dual 自动定价 |
| `--iterations` | 300000 | 正整数，覆盖所选算法的工作预算 |
| `--no-presolve` | 不启用该开关 | 传入时关闭预处理 |
| `--tolerance` | `1e-7` | 仅控制参考目标的相对比较，不改变求解器 primal/dual tolerance |
| `--expected VALUE` | 不核对参考目标 | 指定参考目标 |
| `--infeasible` | 关闭 | 显式预期模型不可行 |
| `--print-solution` | 关闭 | 打印原始变量 |
| `--progress` | 关闭 | 将原始模型上下界和 gap 及时输出到 stderr |
| `--progress-interval S` | 5 个 CPU 秒 | 正数，设置报告间隔并启用进度 |
| `--save-point PATH` | 不导出 | 求解过程中原子保存新的可行目标改善；不保证最优 |
| `--help` | — | 查看当前可用参数 |

以上为测试工具的参数。C API 的求解选项及结果字段见[参数与结果说明](../docs/parameters.md)。当前工具支持的参数也可通过 `python3 data/benchmark.py --help` 和 `./build/data/lp_solve_mps --help` 查看。

运行不可行性基准：

```sh
python3 data/benchmark.py --dataset netlib-infeasible --binary build/data/lp_solve_mps --solvers lp_simplex scip highs clp glpk --repeats 1 --seconds 90
```
