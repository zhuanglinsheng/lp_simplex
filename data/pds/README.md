# PDS 数据集

简体中文 | [English](./README_en.md) | [项目主页](../../README.md) | [数据总览](../README.md)

收录 Mittelmann PDS 目录的 9 个 LP 实例。PDS 指 Patient Distribution System，即伤员分配与疏散问题。`pds-20` 与 Netlib 同名模型相同，较小的 `pds-02/06/10` 位于 Netlib 目录。

## 实例规模

| 实例 | 约束数 | 变量数 | 约束矩阵非零元 |
|---|---:|---:|---:|
| pds-20 | 33,874 | 105,728 | 230,200 |
| pds-30 | 49,944 | 154,998 | 337,144 |
| pds-40 | 66,844 | 212,859 | 462,128 |
| pds-50 | 83,060 | 270,095 | 585,114 |
| pds-60 | 99,431 | 329,643 | 712,779 |
| pds-70 | 114,944 | 382,311 | 825,771 |
| pds-80 | 129,181 | 426,278 | 919,524 |
| pds-90 | 142,823 | 466,671 | 1,005,359 |
| pds-100 | 156,243 | 505,360 | 1,086,785 |

非零元仅统计约束矩阵，不包含目标系数和显式零值。完整来源、大小与校验信息见 `manifest.csv` 和 `SHA256SUMS`。

## 准备模型

需要 Python 3.9+ 与 C 编译器。在仓库根目录执行：

```sh
python3 data/pds/prepare.py pds-30
python3 data/pds/prepare.py --all
```

原始格式为 bzip2 压缩的 Netlib MPC。脚本离线解压，通过 `CC`（默认 `cc`）临时编译 `raw/emps.c`，再转换为 MPS。默认输出到 `data/pds/mps/`，可用 `--output-dir DIR` 修改。中间文件自动清理，失败不会覆盖已有输出。

## 来源与结果

- [Mittelmann PDS 数据目录](https://plato.asu.edu/ftp/lptestset/pds/)及[官方格式说明](https://plato.asu.edu/ftp/lptestset/pds/00README)。
- [Netlib MPC → MPS 转换器](https://www.netlib.org/lp/data/emps.c)，源码保存在 `raw/emps.c`。
- [PDS 基准报告](../../docs/benchmarks/pds.md)及[统一测试入口](../README.md)。
