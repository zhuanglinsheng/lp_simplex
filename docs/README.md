# lp_simplex documentation

This directory contains the design notes, file-format documentation, and user
guides for `lp_simplex`.

## Tutorials

- [基于 CSC 的稀疏 dual revised simplex](./sparse-dual-revised-simplex.md)：
  介绍 CSC 模型、隐式逻辑列、稀疏基分解、Harris BFRT 以及 dual revised
  simplex 的完整迭代过程。

## Benchmarks

- [CSC sparse dual revised simplex 性能测试](./dual-revised-simplex-benchmark.md)：
  记录 NETLIB 长耗时样本的重复计时、Gurobi 仓库参考数据、数值残差和扩展
  正确性回归。
