# Netlib 数据集

简体中文 | [English](./README_en.md) | [项目主页](../../README.md) | [数据总览](../README.md)

收录 114 个可行 LP、29 个不可行 LP，以及 Gurobi 10 预处理版本。可行模型与不可行模型分别运行、分别统计。参考目标值见 [feasible_gurobi_1e-8.csv](feasible_gurobi_1e-8.csv)。

## 准备模型

需要 Python 3.9+。在仓库根目录执行：

```sh
python3 data/netlib/prepare.py feasible/afiro feasible/degen2
python3 data/netlib/prepare.py --all
```

名称需带类别，因为不同类别可能有同名模型。默认恢复到原路径，例如 `data/netlib/feasible/afiro.mps`；用 `--output-dir DIR` 可指定其他输出目录。

公共驱动 [`../lp_solve_mps.c`](../lp_solve_mps.c) 用于单案例求解和参考结果校验。构建后可执行：

```sh
./build/data/lp_solve_mps --expected -464.7531429 data/netlib/feasible/afiro.mps
```

五个求解器的批量对比由 [`../benchmark.py`](../benchmark.py) 统一调度。

## 数据保存

- `raw/{feasible,infeasible,netlib_grbp}/*.mps.gz`：仓库既有 MPS 的无损 gzip 存档。
- `SHA256SUMS`：压缩输入的校验值。
- 生成的 MPS 由 Git 忽略，成功生成后才替换已有输出。

原始集合见 [Netlib](https://www.netlib.org/lp/)。性能结果见 [Netlib 基准报告](../../docs/benchmarks/netlib.md)。

不可行性判定基准见 [Netlib infeasible](../../docs/benchmarks/netlib-infeasible.md)。每例一次、限时 90 秒；PASS 只核对不可行状态，不代表独立证书验证。
