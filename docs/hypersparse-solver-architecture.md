# Hypersparse dual simplex 架构

本文定义 sparse dual revised simplex 的模块边界。目标是让 presolve、pricing、
基更新和退化控制共享稀疏状态，而不是在一个求解循环中维护彼此独立的稠密数组。

## 模型生命周期

求解流程分为 original、presolved、solver 和 postsolve 四个空间。

presolve 不修改用户模型；每项 reduction 都保存原行列映射和恢复数据。solver 只读取
presolved CSC/CSR，求解结束后由 postsolve 恢复原变量并在原模型上计算目标值。

当前可逆规则包括：

- 固定列消去
- 空列按目标方向选界
- 可行空行删除
- 矛盾空行检测
- singleton 行界收紧
- 一般行 implied-bound propagation
- 活动区间不可行
  与冗余判定
- 同支撑平行行的支配和冲突检测。

行列删除时增量维护 degree；
bound propagation 使用向外舍入的安全界决定 reduction，随后 crash 在 reduced
matrix 上计算精确工作界。

presolved model 直接构造 CSC/CSR，不分配 `m × n`
稠密矩阵。后续 aggregation 和 substitution 继续复用同一 reduction record 与
postsolve 栈。

## 稀疏向量协议

`simplex_SparseVector` 同时维护 packed `index/value`、dense scatter、generation
marker 和索引到 packed slot 的常数时间映射。`clear` 只访问上一代活动位置。

当前 BTRAN 结果 `rho` 已打包后交给 CSR pricing；后续 LU 内核应直接产生活动索引，消除从稠密结果重新扫描的过渡成本。

## Pricing

pricing engine 接收 packed `rho`、非基状态集合和 reduced costs，负责：

1. 在 hypersparse CSR scatter 与 dense CSC pricing 间按实际访问量切换；
2. 维护 partial-pricing 分区和周期性全局认证；
3. 一次生成多个稳定候选；
4. 将候选交给 Harris/BFRT；
5. 把 pivot 强度、列稀疏度和预测 fill 反馈给 pivot 选择。

当前 pricing 已有独立模块和生命周期，统一拥有 nonbasic index/slot、candidate
arrays、packed `rho/alpha` 与 ratio 状态，并已完成第 1 项。partial block 不能
破坏全局对偶可行性，未扫描分区必须由 reduced-cost bound 或周期性完整
pricing 认证。

## 基分解与更新

基模块对上层只暴露 factorize、FTRAN、BTRAN、update 和统计信息。更新策略将
从固定长度 product-form eta 演进为稀疏 Forrest–Tomlin 类更新。reinversion
同时考虑 L/U 与 update 的 fill、pivot 相对稳定性、近期 FTRAN/BTRAN 成本、
iterative refinement 次数，以及 compact Pan core 的阶数变化。

## 全流程联动

presolve 结构用于 crash 和 pricing 分区；crash 初始基用于 symbolic ordering；
Pan、BFRT 和 perturbation 共享退化状态；pricing 候选携带稳定性与预测 fill；
basis policy 根据更新成本决定接受 update 或 reinvert。组件通过显式统计结构
通信，不读取其他模块的私有数组。

## 正确性约束

- 所有 reduction 必须可 postsolve；
- partial pricing 必须定期进行全局对偶可行性认证；
- perturbation 必须在返回前移除并重新计算 reduced costs；
- BFRT 的 bound flip 必须同步 primal values、状态和候选集合；
- 稀疏 LU update 必须接受残差认证，失败时触发 reinversion；
- 总体性能按代表性集合评价，但错误最优、错误不可行和错误无界均是正确性失败。
