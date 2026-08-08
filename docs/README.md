# lp_simplex 文档

- [整体算法说明](./algorithm.md)：当前实现的权威说明，覆盖公开入口、
  presolve、dual revised simplex、basis factor、退化控制、tableau 与
  postsolve，并给出端到端和单次 pivot 流程图。
- [NETLIB 基准快照](./dual-revised-simplex-benchmark.md)：带测试日期的历史结果，
  用于性能回归与问题定位，不作为当前算法行为的定义。

过去分散的 presolve、hypersparse、Pan、dual/basis 和 CSC 教程已经合并进
`algorithm.md`。这样“当前实现”和“未来设计”不会继续分散在多份相互重复的
文档中。
