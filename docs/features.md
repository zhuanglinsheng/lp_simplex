# 算法与功能索引

简体中文 | [English](./features_en.md) | [文档](./README.md) | [算法](./algorithm/README.md)

按调用条件、输入产出和实现入口查找能力。只有公开算法可通过 `Options.algorithm` 选择；内部角色、模型表示和后端不构成新的公开算法。测试列为覆盖入口，不是性能保证。

## 公开求解方式

| 名称与用途 | 调用方与启用条件 | 输入 → 产出 | 限制 | 实现入口 | 对应测试 |
| --- | --- | --- | --- | --- | --- |
| TABLEAU：两阶段原始单纯型 | 公开选择；AUTO 小模型恢复 | LP → 点与求解状态 | 稠密化；不导出对偶乘子 | [tableau_solver](../src/tableau/simplex_tableau_solver.c) | [result_contract](../test/test_result_contract.c) |
| DUAL_REVISED：稀疏对偶修正单纯型 | 公开选择；AUTO 主路线及最终恢复 | LP、可选内部基提示 → 点、状态及诊断 | 基提示需重建审核；公开入口不接受热启动参数 | [dual](../src/dual/simplex_dual.c) | [regressions](../test/test_regressions.c) |
| PAN_BDA：动态亏基优化 | 公开选择；内部可从可行种子启动 | LP、可选内部种子 → 点、状态及诊断 | 默认标准形；原生界等变体由内部路线或手动设置选择 | [pan](../src/bda/simplex_pan.c) | [pan](../test/test_pan.c)、[warmstart](../test/test_pan_warmstart.c) |
| AUTO：跨算法调度 | 公开选择 | LP → 累计预算内的求解结果及保留点 | 启发式路线选择；不保证最快 | [auto](../src/core/simplex_auto.c) | [auto_policy](../test/test_auto_policy.c)、[auto_flow](../test/test_auto_flow.c) |

## 内部路线与子过程

| 名称与用途 | 调用方与启用条件 | 输入 → 产出 | 限制 | 实现入口 | 对应测试 |
| --- | --- | --- | --- | --- | --- |
| 原生上下界混合路线 | AUTO hybrid；手动 `PAN_NATIVE_BOUNDS` 设置 | 盒约束 LP、可选可行点 → 候选或认证结果 | AUTO 检查后端及工作区；失败保留有效点 | [bounded](../src/bda/simplex_pan_bounded.c) | [pan_bounded](../test/test_pan_bounded.c) |
| 显式／隐式切锥与标准形变体 | AUTO 路线；PAN 手动设置 | 当前表示 → 方向、候选或求解结果 | AUTO 从隐式切到显式需明确的 oracle 实现故障证据；转换受资源限制 | [route_solve](../src/core/simplex_auto_routes.c)、[cone](../src/bda/simplex_pan_cone.c)、[implicit](../src/bda/simplex_pan_cone_implicit.c) | [auto_policy](../test/test_auto_policy.c)、[cone_oracle](../test/test_pan_cone_oracle.c) |
| 配对单列专用路线 | AUTO singleton，结构及稠密资源检查通过 | 特殊等式 LP → 输入模型上经审核的点与分类 | 零下界、非负变量；逐行正负正成本单列配对；辅助对偶不可行须结合结构证明才能映射为输入无界 | [singleton_dual](../src/dual/simplex_singleton_dual.c) | [singleton_classification](../test/test_singleton_classification.c)、[auto_policy](../test/test_auto_policy.c) |
| 活动列覆盖搜索 | AUTO；宽模型、预算、结构和工作区满足条件 | 非负覆盖 LP → 候选或全量认证结果 | 仅支持规定的 GE、成本和变量界；子问题失败不证明原 LP 不可行 | [cover](../src/dual/simplex_cover.c) | [cover](../test/test_cover.c)、[stage_report](../test/test_stage_report.c) |
| 原模型可行性初始化 | 预求解外层提供回调；AUTO Dual 让出且无可行点时尝试 | 原模型 → 原模型及约化模型审核后的种子 | 有展开规模门槛与独立上限；不能宣称终止状态 | [original_seed](../src/core/simplex.c) | [auto_flow](../test/test_auto_flow.c)、[execution](../test/test_execution.c) |
| NNLS Phase I 与混合 NNLS | PAN 初始化、切锥及可行性恢复 | 最小二乘／带界辅助问题 → 可行性候选和活动集 | 辅助问题结果需按 LP 审核 | [phase1](../src/bda/simplex_pan_phase1.c)、[cone](../src/bda/simplex_pan_cone.c) | [pan_phase1](../test/test_pan_phase1.c)、[mixed_nnls](../test/test_pan_mixed_nnls.c) |
| PDHG 一阶推进 | AUTO 原生 hybrid 默认配置包含；手动 `PAN_NATIVE_PDHG` | 稀疏 LP → 近似点／乘子及交接信号 | 内部有限预算子过程；近似点不等于认证成功 | [pdhg](../src/bda/simplex_pan_pdhg.c) | [pdhg_kernels](../test/test_pan_pdhg_kernels.c)、[pan_bounded](../test/test_pan_bounded.c) |
| 支撑压缩、可行性修复与对偶修整 | 原生 PAN 按状态调用；AUTO 配置含恢复与对偶修整 | 候选点／乘子 → 经检查的改进候选 | 不能仅凭辅助残差宣称最优 | [recovery](../src/bda/simplex_pan_recovery.c)、[first_order](../src/bda/simplex_pan_first_order.c) | [pan_recovery](../test/test_pan_recovery.c)、[pan_bounded](../test/test_pan_bounded.c) |
| 其他手动 PAN 变体 | 显式实验入口捕获环境设置，选择预测、块 NNLS、Newton QR、中心化等 | 当前点／基 → 可选搜索或修复 | 非正式稳定配置；部分操作依赖 SPQR，AUTO 配置隔离相关手动实验 | [config](../src/bda/simplex_pan_config.c)、[配置边界](./runtime-configuration.md) | [environment](../test/test_environment.c)、[pan_trial](../test/test_pan_trial.c) |

## 模型处理与数值机制

| 名称与用途 | 调用方与启用条件 | 输入 → 产出 | 限制 | 实现入口 | 对应测试 |
| --- | --- | --- | --- | --- | --- |
| BDA 配置契约 | AUTO 适配层、手动兼容入口 | `PAN_BDA` 身份及独立的模型/方向/oracle 控制 → 只读配置 | 私有接口；组合需校验；QR 后端独立 | [config](../src/bda/simplex_pan_config.h) | [bda_config](../test/test_bda_config.c) |
| 预求解与后求解 | 求解入口，`presolve=1` 默认开启 | LP → 约化 LP、还原记录 → 原模型点 | 含空行列、固定变量、单列、活动区间、双元代换及支配规则；部分规则有工作上限 | [presolve](../src/presolve/simplex_presolve.c)、[postsolve](../src/presolve/simplex_presolve_postsolve.c) | [presolve](../test/test_presolve.c)、[dominance](../test/test_presolve_dominance.c) |
| 缩放与标准形转换 | Dual 缩放默认启用；PAN／tableau 按路线转换 | 输入 LP → 工作表示及逆映射 | 残差空间和单位不能混用；界展开可能增大模型 | [scaling](../src/runtime/simplex_scaling.c)、[standard](../src/bda/simplex_pan_standard.c)、[transform](../src/tableau/simplex_transform.c) | [scaling](../test/test_scaling.c)、[pan_bounded](../test/test_pan_bounded.c) |
| 定价、比值检验与翻界 | tableau、Dual、PAN 内核 | 当前基／活动集 → 主元或方向选择 | 公开算法与 pricing 组合受验证；Devex 属 Dual 内部原始清理机制 | [dual_pricing](../src/dual/simplex_dual_pricing.c)、[primal_pricing](../src/dual/simplex_dual_primal_pricing.c)、[pivot](../src/tableau/simplex_pivot.c) | [dual_pricing](../test/test_dual_pricing.c)、[primal_pricing](../test/test_primal_pricing.c) |
| 退化、防循环与有限扰动 | Dual 内核按停滞／重复状态触发；tableau 可选 Bland | 进展与基状态 → 面处理、字典序、扰动或恢复 | 工作目标／界变化不提供原模型终止证明 | [degeneracy](../src/dual/simplex_degeneracy.c)、[recovery](../src/dual/simplex_recovery.c) | [degeneracy](../test/test_degeneracy.c)、[recovery](../test/test_recovery.c) |
| 进展观察与活跃状态续跑 | AUTO 中 Dual 探索 | 可比进展、重复状态 → 续跑或让出 | 续跑发生在同一活跃调用内；检查点重启另行重建 | [explore](../src/dual/simplex_degeneracy.c)、[kernel](../src/dual/simplex_dual_kernel.c) | [auto_policy](../test/test_auto_policy.c)、[stage_report](../test/test_stage_report.c) |
| 基检查点与导入审核 | AUTO Dual 探索／恢复；内部调用 | 模型绑定的基提示 → 重建、分类或拒绝 | 地址、版本及执行身份匹配；不保存活跃分解 | [checkpoint](../src/dual/simplex_dual_checkpoint.c)、[entry](../src/dual/simplex_dual.c)、[kernel](../src/dual/simplex_dual_kernel.c) | [checkpoint](../test/test_checkpoint.c)、[dual_seed](../test/test_dual_seed.c) |
| 后端恢复与重复失败抑制 | KLU 数值失败时尝试内部 LU；共享账本限制已知重复冷启动 | 失败及尝试条件 → 恢复或跳过 | 合计消费原授权；未知条件不冒充已失败；账本有限容量 | [backend](../src/dual/simplex_dual_backend.c)、[dual_retry](../src/dual/simplex_dual_retry.c)、[retry](../src/runtime/simplex_retry.c) | [retry](../test/test_retry.c)、[execution](../test/test_execution.c) |
| 原模型重试与候选归档 | 外层后求解审核按条件触发；AUTO 保留可行点 | 还原点及失败 → 原模型重算或恢复已存点 | 剩余预算内执行；保留点不意味着最优 | [solve](../src/core/simplex.c)、[progress](../src/runtime/simplex_progress.c) | [execution](../test/test_execution.c)、[progress](../test/test_progress.c) |

## 线性代数基础

| 名称与用途 | 调用方与启用条件 | 输入 → 产出 | 限制 | 实现入口 | 对应测试 |
| --- | --- | --- | --- | --- | --- |
| 最小范数乘子直接 QR 恢复 | 标准形 BDA 重建后仍未通过乘子审核 | 当前原始基方程 → 独立 QR 候选及原始残差 | 32 MiB 新增工作区；不是 LP 认证；失败不提交乘子 | [basis_qr](../src/bda/simplex_pan_basis_qr.c) | [minimum_norm](../test/test_pan_minimum_norm.c) |
| 稀疏 LU、基更新与精化 | Dual；KLU 可用时使用，否则内部 LU | 方基 → 分解、FTRAN／BTRAN | KLU 是可选构建依赖；更新需残差及重建保护 | [basis](../src/basis/simplex_basis.c)、[sparse_lu](../src/basis/simplex_sparse_lu.c) | [basis_update](../test/test_basis_update.c) |
| 动态亏基 QR、增删列与判秩 | PAN；SPQR 可用时稀疏后端，否则稠密后备 | 独立活动列 → 投影、最小二乘、乘子 | 后端能力与资源不同；部分高级操作无稠密实现 | [basis](../src/bda/simplex_pan_basis.c)、[spqr](../src/bda/simplex_pan_basis_spqr.c)、[linalg](../src/bda/linalg/) | [dynamic_qr](../test/test_pan_dynamic_sparse_qr.c)、[multifrontal_qr](../test/test_pan_multifrontal_qr.c)、[rank_repair](../test/test_pan_rank_repair.c) |
| CSC／CSR 与稀疏工作向量 | 稀疏内核 | 矩阵／索引 → 稀疏乘积及工作集 | 内部基础设施，不是独立求解算法 | [csc](../src/matrix/simplex_csc.c)、[sparse_vector](../src/matrix/simplex_sparse_vector.c) | [sparse_vector](../test/test_sparse_vector.c) |

## 认证、接口与诊断

| 名称与用途 | 调用方与启用条件 | 输入 → 产出 | 限制 | 实现入口 | 对应测试 |
| --- | --- | --- | --- | --- | --- |
| 可行性、最优性与结果审核 | 内核、AUTO、最终返回 | 点／乘子／状态 → 残差及接纳结论 | 认证归内核；公开对偶残差不是独立最优性证明 | [problem](../src/runtime/simplex_problem.c)、[auto](../src/core/simplex_auto.c)、[solve](../src/core/simplex.c) | [result_contract](../test/test_result_contract.c)、[auto_flow](../test/test_auto_flow.c) |
| 浮点不可行证据验证 | Dual 候选不可行出口 | 行组合与变量界 → 矛盾验证 | 验证失败表示未确认，不表示可行 | [certificate](../src/runtime/simplex_certificate.c) | [infeasibility_certificate](../test/test_infeasibility_certificate.c) |
| 精确基行验证 | 浮点无法确认且构建启用 GMP | 当前基 → 精确有理数矛盾验证 | 结构基规模／位数受限；缺少 GMP 不降低证据要求 | [basis_certificate](../src/runtime/simplex_basis_certificate.c) | [infeasibility_certificate](../test/test_infeasibility_certificate.c) |
| 观测、进度回调与因果履历 | 阶段报告；公开可选回调；内部履历 | 点、界及执行事件 → 进度／来源记录 | 缺失值保留未知；显示履历容量有限；不凭小 gap 终止 | [stage](../src/runtime/simplex_stage.c)、[progress](../src/runtime/simplex_progress.c)、[execution](../src/runtime/simplex_execution.c) | [stage_report](../test/test_stage_report.c)、[progress](../test/test_progress.c)、[execution](../test/test_execution.c) |
| C 模型接口与 MPS 读取 | 公开 API 显式调用 | 行列、变量界、目标或 MPS → 连续 LP 模型 | 仅连续 LP；MPS 支持范围以解析器为准，不执行整数优化 | [model](../src/core/model.c)、[mps](../src/core/mps.c)、[API](../include/lp_simplex/) | [simple_examples](../test/test_simple_examples.c)、[regressions](../test/test_regressions.c)、[model_storage](../test/test_model_storage.c) |

参数默认值、环境设置和资源上限见[数值政策](./numerical-policy.md)；AUTO 顺序和恢复契约见[AUTO 调度说明](./auto-policy.md)。新增或移除能力、修改默认启用条件或公开范围时，同步更新本清单及英文版。

架构、阶段契约、数值政策与验证记录分别见 [AUTO](./auto-policy.md)、[BDA](./pan-services.md)、[数值政策](./numerical-policy.md)和[开发记录](./development/README.md)。
