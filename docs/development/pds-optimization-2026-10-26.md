# PDS 优化：原则复核与验收（2026-10-26）

2026-10-27 当前状态：本文以下分支／候选哈希及验收描述属于 10 月 26 日的历史快照；相关候选现已汇集为 main 的未提交工作区改动，不再有独立 PDS 分支。最新离线归因及暂停状态见[更新后的深度分析](pds-performance-analysis-2026-10-26.md#7-2026-10-27隐含上界中的搜索容差污染)。中间原始文件不作为仓库交付。

旧轮次曾保留原始本机证据；按用户后续要求，中间测试文件清理，历史性能目录只保留基准测试 Markdown。下文的原始证据路径与保留描述属于当时快照，不表示当前仍保留这些文件。

本轮从 `c8298eb` 的文档状态、`559100b` 的数值基线出发。初始分层候选为 `be6a6a0`；随后加入同次回滚后列索引复用，当前待验收候选冻结为 `375bbb6`，位于 `codex/pds-kernel-optimization`；在性能验收完成前不成为 main 默认策略。[先前归因](pds-performance-analysis-2026-10-26.md)与[恢复候选反例](pds-auto-review-2026-10-26.md)保留原样。

## 已理解并逐项复核的原则

| 原则与来源 | 本次实现如何遵守 | 验证依据 |
| --- | --- | --- |
| [AUTO 分层与单一状态 owner](../auto-policy.md) | 改动位于 BDA 的 PDHG、恢复内核与服务适配；AUTO 路线及主控状态不改。runtime 仅追加诊断名称。 | 架构检查；AUTO 契约、状态、能力及组合测试 |
| [累计预算与产出权限](../pan-services.md) | 修复使用父 PDHG 的剩余授权和累计计数；子工作已包含于父汇总。仅可行性停滞时，仿射与锥修复各一次外层校正，内部 NNLS 保留原有有界多步授权。 | services、bounded、数值政策及恢复触发测试 |
| [独立认证](../numerical-policy.md) | 停滞只产生机会；修复仅产出点。经审核的可行改进点可以交接，最终仍需原生 KKT 与原模型审核。 | 已知最优、不可行、无界、结果契约及原模型恢复回归 |
| [数值尺度与保护](../numerical-policy.md) | 容差、补偿求和、非有限值检查、目标上限及修复前残差减半条件不改；不读取模型名称、参考答案或墙钟期限。 | kernels、数值政策与跨数据集原模型审核 |
| [配置快照](../runtime-configuration.md) | 正式入口仍使用确定默认配置；不增加环境开关，不在回调后重读数值环境。 | environment、配置快照与 AUTO 入口测试 |
| [身份及映射生命周期](../pan-services.md) | 复用原有绑定表示和服务 begin/finish；不借用过期事件，不将同维度视为同模型，也不改 presolve journal。 | services 身份、配置、绑定、重置和回调测试 |
| [诊断不控制数值行为](../auto-policy.md) | 新触发种类来自实际进展观察和显式服务参数，不从履历槽位推导资格；履历只记录。 | 诊断关闭、容量拒绝、清空和同址重置回归 |
| [资源失败与事务](../numerical-policy.md) | 不新增按模型尺寸分配的必需数组；单轮失败沿原路径回滚，资源失败不伪造不可行或奇异。 | storage、故障注入、服务拒绝及恢复回滚测试 |
| [C90 与接口边界](../parameters.md) | 不更改公开 Options、Result、Progress 或算法种类；变化只在私有接口及追加诊断枚举。 | Release/Debug-ASAN、完整/最小依赖配置；格式及架构检查 |

“原则成立”指上述控制边界、认证与资源语义保持，并有对应回归覆盖；有限样本测试不能证明所有 LP 都更快或任何平台均无缺陷。

## 选择的数值优化

加权综合指标和固定工作模型的未加权可行性残差分别保留进展历史。改变预条件权重仅重置综合指标的历史；物理可行性历史不因此丢失。沿用 8192 步窗口与 1% 累计进展标准，未增加按数据集调节的门槛。

综合指标停滞仍使用原有有限恢复搜索。仅可行性停滞时，授权一次仿射校正和一次锥校正；单次外层尝试不等于一个 NNLS 工作单位。新增两行三变量例验证：内部 NNLS 必须可以使用多个工作单位；另一个带上界的例验证单轮失败完整回滚、完整搜索仍可成功。

恢复后只有独立审核确认的可行改进点才允许提前交接；交接不宣称最优。原生算法继续优化并认证。目标上限、相对上次实际修复至少减半的残差、原始/逐行尺度检查、父剩余预算及最终认证均保留。

这些是针对可观测数值状态的通用规则，不是保证每例加速的定理。之前 qap15 的额外成本仍必须纳入验收。

## 修复路径中的重复工作优化

初始候选 `be6a6a0` 的两对 qap15 控制中，基线两次通过（78.000 / 87.951 秒），候选一次通过（89.272 秒）、一次在 90 秒超时；流程因此停止，没有据此合入默认策略。

复查发现：单次仿射修复失败并将点完整回滚后，锥修复会对完全相同的模型、同一点重复选择自由面独立列。`375bbb6` 仅在第一次、唯一一次仿射选列成功的条件下，保留列索引提示并让紧随的锥修复借用；未知选择和多轮仿射不复用。锥修复仍为实际锥表示重新加列、分解和审核，未跨模型继承活因子、点、乘子或认证。它复用的是同次调用的索引工作区，替代原有同尺寸分配，不增加必需工作区峰值。

新增 `test_bounded_restore_storage` 将根工作区分配逐一注入失败，覆盖空／非空自由面、锥修复成功和一轮失败回滚；检查共享预算、原始残差或精确回滚、泄漏和错误释放。Release 定向八项及更新后的四配置全部通过，共 704 项；完整性能验收仍未完成。复用不会调整任何数值政策门槛，也不根据数据集选择路径。

## 没有采用的算术候选

先测了两个保持算术结果的实现候选，避免直接把新策略与局部性能变化混在一起：

| 候选 | 正确性检查 | 性能证据 | 处置 |
| --- | --- | --- | --- |
| 将列变量更新与残差 scatter 合并 | 179 项功能测试及逐位差分通过 | 实际 PDS 矩阵八对内核微测，候选/基线 CPU 比中位数约 1.247（20）和 1.170（70） | 撤销，保留补丁与原始记录 |
| 将逐行残差校正与对偶更新合并 | 179 项功能测试及补充差分通过 | pds-20 两对完整求解均通过；CPU 中位数约快 1%，墙钟中位数约慢 2% | 未证明稳定收益，撤销 |

两个算术候选均不在当前冻结候选中。不能用一次较快结果覆盖反序配对的较慢结果，也不能把内核微测的比例当作完整求解加速比。日期归档保留各自实际范围、原始记录与撤销补丁。

## 验收状态

功能与内存检查使用同一冻结源码，按 Release/Debug-ASAN × 完整/最小依赖四配置检查。性能阶段使用独立构建与输入 SHA-256，先检查 qap15 的两对反序控制，再按原参数覆盖全部 159 个 LP 案例：Netlib 114、不可行集 29、PDS 9、Railway 7。正式限时、容差和工作预算不延长。

四配置结果为 Release 完整依赖 180 项、Debug-ASAN 完整依赖 180 项、Release 最小依赖 172 项、Debug-ASAN 最小依赖 172 项，合计 704 项全部通过。归档保留各配置完整构建／测试日志。

当前候选的 qap15 两对反序控制全部通过，但墙钟中位数由基线约 53.430 秒增至候选约 62.754 秒，仍慢约 17.45%；这项代价没有消失，也不能用不同轮次的绝对计时宣称索引复用已经带来完整求解加速。

随后 PDS 批次仅产生七个完整结果：pds-20／30／40／50／60／70 通过（39.556／115.841／101.809／136.096／143.895／191.910 秒），pds-80 在原 300 秒限时超时；pds-90 被中断而无完整结果，pds-100 未运行。Netlib、不可行集和 Railway 尚未开始。这些记录不能换算成九例全集通过率，单次 pds-70 通过也不是稳定加速证明。

机器负载反复升高，影响大案例计时。用户于 2026-10-26 明确要求稍后安排性能测试；测试和监控进程已停止，没有后台续跑或自动调度。候选保留在独立分支，main、默认 build 及原五求解器完整报告保持原基线。下一轮需要在可比较环境补齐原参数全集与跨数据集控制，尤其检查 qap15 的额外成本；当前不能宣称普遍加速或性能无回退。

历史证据目录保存源码与构建指纹、补丁、704 项日志、两轮 qap15 控制、PDS 超时和中断记录、机器负载以及被撤销的算术候选；SHA-256 清单覆盖归档文件。

## English summary

The candidate preserves algorithm ownership, shared work accounting, model identity, configuration snapshots, output-only diagnostics, conservative resource/rollback handling and independent final certification. Independent physical-feasibility history permits a bounded one-affine/one-cone opportunity without suppressing the original bounded multi-work NNLS. Existing merit-stall search and all final tolerances remain unchanged. No dataset name, reference answer or wall-time threshold enters numerical policy.

Two arithmetic-fusion attempts were tested and rejected rather than presented as improvements: primal/column fusion was slower in eight-pair matrix-kernel measurements; row/dual fusion did not demonstrate stable wall-time benefit in two full pds-20 pairs. Their raw results and patches are retained. All 704 functionality/memory tests passed across four configurations. Full cross-dataset performance acceptance remains incomplete: both alternating qap15 pairs passed, but candidate median wall time was 17.45% slower; PDS has only six passes and one full timeout, with two cases unfinished. The user deferred further performance testing, and all related processes stopped. Main, the default build and the complete five-solver reports remain at baseline. Raw successes, failures, interruptions and provenance are kept locally outside the repository.
