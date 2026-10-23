# 项目文档

简体中文 | [English](./README_en.md) | [项目主页](../README.md) | [算法总览](./algorithm/README.md)

本文档目录提供算法原理、求解结果解释和性能测试的阅读入口。构建、安装及基本使用见[项目主页](../README.md)。

## 使用求解器

[使用指南](./usage.md)以 `data/simple_examples/production_planning.mps` 为例，依次展示数学模型、C API 手动建模、MPS 读取和独立程序的编译运行，并链接其他完整示例。

## 算法原理

[算法文集](./algorithm/README.md)按数学模型、假设、推导与最优性条件组织，覆盖预处理、两阶段原始单纯型、对偶修正单纯型、抗退化、NNLS、Pan/BDA、切锥恢复、PDHG 及基矩阵运算。每篇均提供中英文版本，可在标题下切换语言。

建议先阅读总览中的统一记号与原始—对偶条件，再按关注的方法阅读专题。

## 求解结果

[统一结果契约](./parameters.md)解释成功与失败状态、目标值、原始可行性、对偶诊断所属模型空间及缺失值，供调用求解接口和判断返回结果时查阅。

## 性能测试

[基准测试汇总](./benchmarks/README.md)提供 lp_simplex、SCIP/SoPlex 与 HiGHS 的测试信息、对比统计和复现说明。各数据集的逐案例结果见：

- [Netlib feasible](./benchmarks/netlib.md)
- [PDS](./benchmarks/pds.md)
- [Railway 连续 LP 松弛](./benchmarks/railway.md)

数据准备、运行参数和审核口径见[数据与基准测试说明](../data/README.md)。逐次运行的配置、原始结果、日志和解向量保存在本地 `tmp/benchmarks/`，由 Git 忽略。
