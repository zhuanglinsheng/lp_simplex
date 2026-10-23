# PDS LP 测试集

下载日期：2026-10-22。收录 Mittelmann PDS 目录的全部 9 个实例：
`pds-20`、`pds-30`、`pds-40`、`pds-50`、`pds-60`、`pds-70`、
`pds-80`、`pds-90` 和 `pds-100`。

PDS 是 Patient Distribution System（伤员分配/疏散）问题。
`pds-20` 与 `../netlib/feasible/pds-20.mps` 内容相同；其余 8 个实例
用于扩展现有 Netlib 集合。较小的 `pds-02/06/10` 仍位于 Netlib 目录。

## 实例规模

下表为当前项目读取器的统计，约束矩阵非零元不包括目标系数和显式零值。
仓库仅保留约 24 MiB 的压缩原文件；全部生成的 MPS 约 295 MiB，不纳入版本控制。

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

## 目录

- `raw/`：9 个原始 `.bz2` 文件及官方 `00README`。
- `raw/emps.c`：Netlib 官方 MPC → MPS 转换器源码，未修改。
- `prepare.py`：按需编译转换器并生成 MPS，仅依赖 Python 3.9+ 与 C 编译器。
- `mps/`：脚本生成的 fixed-column MPS，由 Git 忽略。
- `manifest.csv`：来源、文件大小及模型规模。
- `SHA256SUMS`：保留数据文件的本地 SHA-256 校验值。

原始文件有两层编码：**bzip2 → Netlib MPC → MPS**。在仓库根目录运行：

```sh
python3 data/pds/prepare.py pds-30
python3 data/pds/prepare.py --all
```

脚本完全离线工作，使用 Python 标准库解压，并通过 `CC`（默认 `cc`）编译
本地 `raw/emps.c`。转换器及中间 MPC 位于临时目录，结束后自动清理。
`--output-dir DIR` 可指定其他输出目录。MPS 成功生成后才替换目标文件；
压缩数据损坏或编译失败不会覆盖已有结果。单独解压或改扩展名不能得到 MPS。

```sh
(cd data/pds && shasum -a 256 -c SHA256SUMS)
```

## 使用

先按上方命令生成所需实例，再在仓库根目录执行：

```sh
./build/tools/test_netlib --algorithm dual-revised --iterations 300000 \
    data/pds/mps/pds-30.mps
```

最初整理只验证压缩包、格式转换和当时源码 MPS 读取器的解析，未进行
全量求解基准测试，也未生成新的最优目标值。`test_netlib` 对没有参考值的
实例只要求求解成功；要验证目标值，应提供独立核实的 `--expected VALUE`。

## 来源

- [Mittelmann PDS 下载目录](https://plato.asu.edu/ftp/lptestset/pds/)
- [官方来源与格式说明](https://plato.asu.edu/ftp/lptestset/pds/00README)
- [Netlib emps.c](https://www.netlib.org/lp/data/emps.c)
- W. Carolan, J. Hill, J. Kennington, S. Niemi, S. Wichmann,
  “An Empirical Evaluation of the KORBX Algorithms for Military Airlift
  Applications”, Operations Research 38 (1990), 240–248。
