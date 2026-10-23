# Railway 铁路集合覆盖测试集

下载日期：2026-10-23。收录 OR-Library 全部 7 个 `rail` 实例。
这些实例源于意大利铁路的实际集合覆盖应用，由 Paolo Nobili 提供。

## 来源与官方说明

- [OR-Library 官方说明](https://people.brunel.ac.uk/~mastjjb/jeb/orlib/scpinfo.html)
- [官方文件目录](https://people.brunel.ac.uk/~mastjjb/jeb/orlib/files/)
- `raw/scpinfo.html`：下载的官方说明原文，未修改；这是所有 set-covering 数据的说明，其中包含 railway 专节。
- 每个实例的直接下载 URL、日期、文件大小和校验值见 `manifest.csv`。

## 实例规模

以下统计来自下载后的完整文本解析，非零元不包括目标系数。

| 实例 | 约束数 | 变量数 | 矩阵非零元 |
|---|---:|---:|---:|
| rail507 | 507 | 63,009 | 409,349 |
| rail516 | 516 | 47,311 | 314,896 |
| rail582 | 582 | 55,515 | 401,708 |
| rail2536 | 2,536 | 1,081,841 | 10,993,311 |
| rail2586 | 2,586 | 920,683 | 8,008,776 |
| rail4284 | 4,284 | 1,092,610 | 11,279,748 |
| rail4872 | 4,872 | 968,672 | 9,244,093 |

rail2536 与 rail4284 超过百万变量；rail2586 与 rail4872 接近百万变量。
实例名称中的数字对应约束数。百万规模指变量数量，矩阵非零元最高约 1,128 万。
仓库仅保留原始压缩包；MPS 按需生成，不纳入版本控制。

## 目录

- `raw/rail*.gz`：官方 gzip 压缩文件，保留原始字节。
- `raw/scpinfo.html`：官方来源、规模、结构特点与文件格式说明。
- `prepare.py`：直接流式读取 gzip，生成 fixed-column MPS 连续 LP 松弛。
- `mps/`：脚本生成的 MPS，由 Git 忽略。
- `manifest.csv`：原始数据来源、下载日期、实测规模、大小、SHA-256 和解析检查结果；
  `text_bytes` / `text_sha256` 用于校验按需解压的文本。
- `SHA256SUMS`：本目录数据及说明的本地校验值，不是发布者提供的签名。

## 文件格式

这些文件采用 OR-Library 的 railway 专用列格式，**不是 MPS**。
项目现有 MPS 读取器不能直接读取解压后的文本。在仓库根目录执行
（Python 3.9+，无需第三方包或 C 编译器）：

```sh
python3 data/railway/prepare.py rail507
python3 data/railway/prepare.py --all
```

默认输出到 `mps/`；`--output-dir DIR` 可指定其他目录。脚本直接读取压缩包，
不会保存解压文本。所有实例均生成变量界为 `[0, 1]` 的连续 LP 松弛。
每列的成本、覆盖行和完整性均在转换时检查；失败不会覆盖已有输出。
大实例的 MPS 很大，建议显式指定所需实例。仅改扩展名不构成格式转换。

按空白分隔读取整数，不应假定每条列记录只占一行：

```text
m n
c_1 k_1 row_1 row_2 ... row_k1
c_2 k_2 row_1 row_2 ... row_k2
...
c_n k_n row_1 row_2 ... row_kn
```

`m` 为行数，`n` 为列数；`c_j` 为列成本，`k_j` 为该列覆盖的行数；
后续行号使用从 1 开始的编号，对应矩阵系数 1，未列出的系数为 0。
这与官方同页其他 set-covering 实例的行格式不同。

## 后续求解口径

原始集合覆盖模型为：

```text
minimize  sum_j c_j * x_j
subject to sum_{j: column j covers row i} x_j >= 1,  for each row i
           x_j in {0, 1}
```

本项目目前求解连续 LP，因此计划使用明确保留上下界的 LP 松弛：
将最后一行替换为 `0 <= x_j <= 1`。
三求解器对比时，lp_simplex、SCIP/SoPlex、HiGHS 必须读取同一个连续模型，
不能将 SCIP 的整数求解结果与我们的 LP 松弛结果混比。
rail507 的 MPS 已经由独立读取器逐列对照原始文本，验证全部成本、行索引、
系数、变量界和 RHS。此处未填入未经验证的参考最优值。

## 完整性检查

本次已逐个读取 gzip 到文件末尾以检查 CRC，并完整解析每个实例：

- 文件头的行列数与官方说明一致；恰好包含声明数量的列记录，无多余数据。
- 成本均为 1 或 2；每列覆盖 1 至 12 行；行号均在 `[1, m]` 内。
- 每列没有重复行号，每行至少被一列覆盖。
- 按全部列记录累计非零元，结果见上表和 `manifest.csv`。
- 三个小实例的 gzip 解压结果与另外下载的官方 `.txt` 字节完全一致。

从仓库根目录复核本地文件：

```sh
(cd data/railway && shasum -a 256 -c SHA256SUMS)
```

从已保存的压缩包重新解压单个实例，例如：

```sh
gzip -dc data/railway/raw/rail2536.gz > /tmp/rail2536.txt
```
