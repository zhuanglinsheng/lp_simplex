# lp_simplex documentation

This directory contains the design notes, file-format documentation, and user
guides for `lp_simplex`.

## Tutorials

- [基于 CSC 的稀疏 dual revised simplex](./sparse-dual-revised-simplex.md)：
  介绍 CSC 模型、隐式逻辑列、稀疏基分解、Harris BFRT 以及 dual revised
  simplex 的完整迭代过程。
- [Pan 退化控制与主求解流程融合](./pan-degeneracy-control.md)：
  介绍 deficient-face 对偶方向、退化触发状态机、符号扰动、循环保护及其
  在 sparse dual revised simplex 主循环中的实现。
- [Hypersparse dual simplex 架构](./hypersparse-solver-architecture.md)：
  定义 presolve/postsolve、packed sparse vector、pricing engine、稀疏基更新和
  跨组件成本策略的边界。

## Benchmarks

- [CSC sparse dual revised simplex 性能测试](./dual-revised-simplex-benchmark.md)：
	记录 30 个代表性 NETLIB 模型的当前求解时间、迭代数、正确性结果以及
	Gurobi 仓库参考数据。
