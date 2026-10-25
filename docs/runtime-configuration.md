# 正式配置与实验环境边界

简体中文 | [English](./runtime-configuration_en.md) | [文档](./README.md) | [参数与结果](./parameters.md)

`lp_simplex_solve` 是正式求解入口。四种公开算法均使用公开 Options 和确定的数值默认值，忽略库的数值 `LP_SIMPLEX_*` 环境开关。默认 AUTO 不需要设置环境变量，也不会因上一次实验留下的开关改变预解、缩放、扰动、定价或后端偏好。`LP_SIMPLEX_PROFILE` 继续只控制诊断；第三方 BLAS 等运行库自己的线程环境不属于这套开关。

## 调用与迁移

```c
struct lp_simplex_Options options;
struct lp_simplex_Result result;
int state;
lp_simplex_default_options(&options, lp_simplex_ALGORITHM_AUTO);
/* Set supported options such as presolve, budget, tolerances and callback. */
state = lp_simplex_solve(model, &options, point, &result);
```

需要复现实验环境时，明确改用 `lp_simplex_solve_experimental(model, &options, point, &result)`。该入口在调用开始读取一次数值开关，之后使用求解局部快照；求解器不修改进程环境，回调中改变环境也不改变当前数值配置。调用方负责进程环境并发修改的同步；数值快照不会使环境写入变得线程安全。实验开关不是正式的稳定参数或新增公开算法。

命令行默认也使用正式入口。旧的环境实验须增加显式选项，例如：

```sh
LP_SIMPLEX_DISABLE_SCALING=1 build/data/lp_solve_mps \
  --experimental-environment --algorithm dual-revised --no-presolve MODEL.mps
```

驱动输出 `configuration: production` 或 `configuration: experimental-environment`；PROFILE 根事件也标明模式。没有增加 Options、Progress 或 Result 字段，既有结构布局不变。以前通过普通入口设置数值环境开关的程序，应迁移到实验入口；普通默认求解仍保留历史“环境为空”时的数值行为。

## 配置由谁拥有

| 边界 | 责任 |
| --- | --- |
| 公开求解根入口 | 创建 `simplex_SolveContext`，解析 `simplex_RunConfig`，固定本次配置 |
| AUTO 主控 | 私有独立调用缺少配置时建立正式默认配置；已有根配置沿模型作用域继承 |
| 预解、缩放及子模型 | 继承同一根快照，不重新解释进程环境 |
| Dual | 将解析后的 `simplex_DualControls` 复制到迭代状态；热循环读取字段 |
| LU | configured 创建接口只接收具体 KLU 偏好，不接收 AUTO/实验模式标签 |
| BDA | 手动入口使用根快照；AUTO 完整路线仍由适配器选择确定的配置 |
| 原模型 NNLS 初始化 | 接收已解析数值控制并明确使用标准形；不再从裸基构造器间接读实验环境 |
| 失败比较 | 键使用实际消费的规范化 Dual 控制；PROFILE 和不被 Dual 消费的 PAN 设置不能改变比较资格 |

AUTO 的显式实验入口可启用预解/缩放与 Dual 实验开关，但不把 PAN 环境重新混入其完整 BDA 路线工厂。原模型 NNLS 初始化使用根快照中的数值设置，模型表示仍由初始化器明确选择。

低层未绑定配置的私有兼容入口仍支持历史实验调用，供内部测试与研究使用；它们不代表正式公开接口。组合新流程时须创建并配置共享求解上下文，或使用明确的 configured 数值接口。有效的基、点或模型身份不等同于允许重新读取环境。

## 历史实验语义

存在型开关的字符串 `0` 仍表示启用；这只适用于显式实验入口。Pan 的整数预算、epoch 和 block 参数沿用有界解析：错误/负数禁用，溢出饱和，block 再限制在支持范围；权重要求有限且合法。标准入口不解析这些数值环境控制。所有入口仍使用相同的输入验证、累计预算、原模型审核和结果契约，实验设置不能降低认证要求。

切换数值算法、模型表示、oracle 或后端时仍应分别解释对应层次。环境快照只是配置来源，不是求解状态、候选点、履历或新的算法路线。实现为 `simplex_run_config.c/.h`，验证见 `test_environment`；逐次证据和限制在[开发记录](./development/README.md)。
