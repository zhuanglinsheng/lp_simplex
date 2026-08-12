# lp_simplex 文档

- [整体算法说明](./algorithm.md)：当前实现的权威说明，覆盖公开入口、
  presolve、dual revised simplex、basis factor、退化控制、tableau 与
  postsolve，并给出端到端和单次 pivot 流程图。
- [Pan/BDA 动态亏基算法](./pan-bda.md)：独立 Pan 求解路径的数学模型、
  Phase I/II、不变量、动态亏基更新、稀疏正交分解、终止认证与实现映射。
- [NETLIB 基准快照](./dual-revised-simplex-benchmark.md)：带测试日期的历史结果，
  用于性能回归与问题定位，不作为当前算法行为的定义。

Presolve、hypersparse、dual/basis 和 CSC 的共同流程统一收录在
`algorithm.md`；Pan/BDA 因具有独立的标准形、可行性不变量和动态亏基线性代数，
单列专文。`algorithm.md` 第 9 节只说明它与 dual-revised 中 Pan-inspired
anti-stalling 的边界，并链接到完整定义。这样既避免重复，又不把两个名称相近、
数学机制不同的算法混为一谈。
