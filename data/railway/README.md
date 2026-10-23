# Railway 数据集

简体中文 | [English](./README_en.md) | [项目主页](../../README.md) | [数据总览](../README.md)

收录 OR-Library 全部 7 个铁路集合覆盖实例，源于意大利铁路的实际应用。本项目求解其连续 LP 松弛，变量界为 `[0, 1]`。

基准测试和报告按名称中的数字递增排列：`rail507`、`rail516`、`rail582`、`rail2536`、`rail2586`、`rail4284`、`rail4872`。

## 实例规模

| 实例 | 约束数 | 变量数 | 矩阵非零元 |
|---|---:|---:|---:|
| rail507 | 507 | 63,009 | 409,349 |
| rail516 | 516 | 47,311 | 314,896 |
| rail582 | 582 | 55,515 | 401,708 |
| rail2536 | 2,536 | 1,081,841 | 10,993,311 |
| rail2586 | 2,586 | 920,683 | 8,008,776 |
| rail4284 | 4,284 | 1,092,610 | 11,279,748 |
| rail4872 | 4,872 | 968,672 | 9,244,093 |

百万规模指变量数量；实例名称中的数字对应约束数。非零元不包括目标系数。来源、文件大小和校验值见 `manifest.csv` 与 `SHA256SUMS`。

## 准备模型

需要 Python 3.9+，无需第三方包。在仓库根目录执行：

```sh
python3 data/railway/prepare.py rail507
python3 data/railway/prepare.py --all
```

`raw/rail*.gz` 是官方压缩数据，其内容为 Railway 专用列格式，并非 MPS。脚本读取列成本及覆盖行，生成 fixed-column MPS 连续模型。默认输出到 `data/railway/mps/`，可用 `--output-dir DIR` 修改。

大实例会生成较大的 MPS，建议按需准备。压缩输入保留原始字节，生成文件由 Git 忽略；转换失败不会覆盖已有输出。

## 来源与结果

- [OR-Library 数据说明](https://people.brunel.ac.uk/~mastjjb/jeb/orlib/scpinfo.html)及[文件目录](https://people.brunel.ac.uk/~mastjjb/jeb/orlib/files/)。
- `raw/scpinfo.html` 保留官方格式说明，`manifest.csv` 记录每个实例的下载来源和解析结果。
- [Railway 基准报告](../../docs/benchmarks/railway.md)及[统一测试入口](../README.md)。五个求解器使用同一连续 LP 松弛模型。
