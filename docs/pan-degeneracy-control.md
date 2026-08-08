# Pan 退化控制与 sparse dual revised simplex 融合

## 目标

高度退化的 LP 会在同一个顶点附近产生大量基交换。对偶单纯形仍然保持对偶可行，但对偶目标的实际改变量接近零：

\[
\Delta z_D = \theta\,v_p \approx 0,
\]

其中 \(v_p\) 是离基变量相对目标边界的违反量，\(\theta\) 是对偶比率
检验得到的步长。此时继续只依赖普通的 pivot magnitude tie-break，可能在
多个表示同一顶点的基之间停滞。

本求解器把 Pan 的 basis-deficiency 思路融合为一种自动激活的
**bounded-deficiency face mode**。正常迭代使用 CSC、稀疏 LU、eta update、
dual steepest edge 和 Harris ratio test；观测到持续退化后，在当前离基列所
定义的亏缺面上启用符号扰动，并允许基求解器使用亏缺结构对应的紧凑分解。

## 亏缺面的计算

设当前基为

\[
B=[b_1,\ldots,b_m].
\]

选择位置 \(p\) 离基后，暂时保留的列为

\[
C=B\setminus\{b_p\}\in\mathbb R^{m\times(m-1)}.
\]

对偶 revised simplex 计算

\[
B^T\rho=e_p.
\]

因此

\[
b_i^T\rho=0\quad(i\ne p),
\]

即 \(\rho\in\operatorname{null}(C^T)\)。这正是从 rank-\((m-1)\) deficient
basis 所关联的面离开的对偶方向。实现不需要显式构造矩形矩阵或稠密
\(Q\)：现有 BTRAN 已经给出该亏缺面的正交方向。

沿 \(\rho\) 计算每个非基列的 tableau 系数：

\[
\alpha_j=a_j^T\rho.
\]

随后仍执行有界变量对偶比率检验。进入列通过 FTRAN 验证：

\[
d=B^{-1}a_q,
\qquad d_p=\alpha_q.
\]

所以 Pan 模式不会绕过现有的相对 pivot 检查、bound flipping 或残差检查。

## 退化观测量

每次成功 pivot 后记录：

\[
g_k=|\theta_k v_{p,k}|.
\]

当

\[
g_k\le 10\epsilon_D
\]

时，把该 pivot 记为退化 pivot，其中 \(\epsilon_D\) 是 dual feasibility
tolerance。控制器同时维护：

- 连续退化 pivot 数；
- 最近 16 次 pivot 的退化位图；
- Pan 模式进入次数；
- 亏缺面符号扰动次数；
- 最近 64 个求解状态的指纹。

以下任一条件成立时激活 Pan 模式：

1. 连续 12 次 pivot 的实际对偶目标改变量接近零；
2. 最近 16 次 pivot 中至少 15 次退化。

这是对求解过程行为的判断，不依赖模型名称、矩阵规模或 NETLIB 分类。

## 符号扰动

在普通模式下，Harris ratio window 内选择 pivot magnitude 最大的进入列，
以改善数值稳定性。

Pan 模式只在最小比率接近零时改变同率候选的顺序。中小型基按变量的稳定
全局编号选择第一个候选，相当于对零 reduced cost 施加不实际写入浮点数组
的字典序微扰：

\[
\bar c_j(\varepsilon)=\bar c_j+
\varepsilon^{j+1}.
\]

实际使用的 \(\theta\) 仍是未扰动问题的比率，因此不会为了 anti-stalling
而引入真实的 dual infeasibility。非零比率仍使用 Harris window 和最大
pivot 规则。

对于至少 4096 行的大型基，Pan 模式同时使用 deficient-core 结构评分：

1. 在 DSE 得分达到当前最大值 80% 的离基行中，优先选择结构基列离基；
2. 对同率进入列，首先最小化结构核心阶数的变化
   \(\Delta r=[q\in A]-[p\in A]\)；
3. 若 \(\Delta r\) 相同，则最大化
   \(|\alpha_q|/(1+\operatorname{nnz}(a_q))\)，兼顾 pivot 强度和预测 fill-in；
4. 最后才用变量编号作确定性 tie-break。

因此大型问题上的 Pan pivot 不只是防循环，而会主动让结构列退出、逻辑列
进入，使后续 reinversion 的 deficient core 变小。DSE 的 80% 窗口保留了
离基行的全局几何质量，避免为了降一阶而选择严重违反 steepest-edge 方向的行。

取得一次非退化进展后，控制器立即退出 Pan 模式，并进入 32 次 pivot 的
cooldown。这样符号扰动只负责离开退化面，不长期替代 DSE 的全局几何
选择。

一次激活最多连续执行 64 个退化 pivot。若求解过程中已经成功构造过紧凑核，
预算耗尽后退出 Pan 模式并冷却 32 次 pivot，随后若停滞仍然存在则允许重新
进入；若始终未形成紧凑核，预算耗尽说明本轮 face 探索没有产生线性代数降维
收益，后续改回确定性的 DSE/Harris 路径。检测到真实状态重复时也会在本次
求解中停止 Pan 扰动。这个策略允许有实际降维价值的 Pan 路径重复工作，同时
限制仅改变 pivot 轨迹、却不能缩小分解维度的长期探索。

## 紧凑基分解

模型变换后的逻辑列为 \(-e_i\)。在一次 reinversion 时，把基列分成结构列
集合 \(S\) 和逻辑列集合 \(L\)。令 \(D\) 为逻辑列对应的行，\(R\) 为其补集。
若结构基列数为 \(r\)，则 \(|R|=|S|=r\)。按行列分块后，基矩阵等价于

\[
B=\begin{bmatrix}
A_{R,S}&0\\
A_{D,S}&-I
\end{bmatrix}.
\]

因此只需对 \(r\times r\) 的核心

\[
K=A_{R,S}
\]

做稀疏 LU。这个消元是精确的块分解，不改变原 LP，也不删除约束。

对 FTRAN \(Bx=b\)：

\[
Kx_S=b_R,\qquad x_L=A_{D,S}x_S-b_D.
\]

对 BTRAN \(B^Ty=g\)：

\[
y_D=-g_L,\qquad
K^Ty_R=g_S-A_{D,S}^Ty_D.
\]

核心矩阵直接从原 CSC 的结构列抽取。`row_to_core` 把原约束行映射到核心行，
`core_basis` 和 `core_position` 保存结构列及其基位置；恢复步骤仍按原基位置
输出，所以现有 eta 更新不需要改变。

紧凑求解后以原始完整基计算残差。相对残差超过 \(10^{-11}\) 时最多执行两次
迭代精化，修正方程仍由同一个核心 LU 求解。

## 激活与成本模型

紧凑分解不会在 dual crash 中启用。退化控制器首次确认 stalling 后，给基求解器一次持续授权；此后每次正常 reinversion 都重新检查当前结构，而不是只在 64 次符号扰动预算内降维。

当前同时满足以下条件才采用紧凑核心：

- 已经由实际目标进展观测触发 Pan 模式；
- 完整基至少有 4096 行；
- 核心阶数不超过完整基阶数的三分之一。

第一个条件避免在 crash 阶段改变建基轨迹；规模条件把持续降维用于恢复成本能够被大规模稀疏分解收益覆盖的问题。核心进入紧凑模式后不再设置 512 阶
下限，可以继续缩小到真正的低维 deficient-core 表示；旧的下限会在 Pan 最成功
地降低核心后反而恢复完整基。若核心增长超过三分之一阈值，当前精确 factor
与 eta 链仍可继续使用，到下一次自然 reinversion 再构造完整基。

## 循环保护

首次进入 Pan face 时，对以下状态计算双指纹：

- 有序 basis index 数组；
- 所有变量的 basic/lower/upper/fixed/free 状态。

此后 basis exchange 和 bound flip 只对变化的槽位执行 XOR 增量更新，不再扫描完整 basis/status。控制器保存最近 64 个状态。如果同一状态再次出现，本次求解永久停止 Pan
扰动，继续使用确定性的 DSE/Harris 路径。指纹只用于拒绝探索，不用于宣告最优、不可行或无界，因此哈希碰撞最多导致少用一次 Pan 模式，不影响求解结论。

## 与主求解流程的关系

主循环顺序如下：

1. DSE 选择 primal-infeasible 的离基位置；
2. BTRAN 得到亏缺面的方向 \(\rho\)，并精确校准该行 DSE 权重；
3. CSC pricing 计算 \(\alpha=A^T\rho\)；
4. 根据控制器状态执行普通 Harris 选择或 Pan 符号扰动；
5. FTRAN 验证进入列和相对 pivot；
6. 执行 bound flipping 或 basis exchange；
7. 使用精确 DSE 递推更新全部权重；
8. 更新 eta factor、退化统计和求解状态指纹；
9. 必要时 reinversion，并重新计算 primal values 与 reduced costs。

Pan 模式不在 dual crash 中激活。Crash 的任务是建立初始 dual-feasible
basis，此时还没有足够的主迭代历史判断 stalling。reinversion、精度恢复和
dual feasibility 恢复也优先于退化控制。

## 数值不变量

融合实现始终保持以下条件：

\[
Bx_B=-N x_N,
\]

\[
B^T\pi=c_B,
\]

以及非基变量相对于其 lower/upper/free 状态的 reduced-cost 符号条件。

Pan 模式不会近似 FTRAN、BTRAN 或 DSE 权重。紧凑核心来自逻辑列的精确块消元；相对 pivot tolerance、残差重算和 reinversion 条件与普通模式相同。

## 代码结构

- `simplex_degeneracy.c/.h`：退化窗口、激活/退出、cooldown、状态指纹；
- `simplex_dual.c`：目标步观测、亏缺面 ratio tie perturbation、统计接入；
- `simplex_basis.c`：核心行/列映射、块 FTRAN/BTRAN、迭代精化和 eta update；
- `simplex_sparse_lu.c`：从 CSC 抽取核心子矩阵并完成稀疏 LU；
- `simplex_csc.c`：继续承担列点积和进入列装配。

设置 `LP_SIMPLEX_PROFILE=1` 时，求解器额外报告：

- `pan_degenerate`：观测到的退化 pivot 数；
- `pan_activations`：Pan 模式激活次数；
- `pan_probes`：符号扰动选择次数；
- `factor_core`：当前实际 LU 阶数与完整基阶数；
- `compact`：紧凑分解次数以及本次求解出现过的最小、最大核心阶数。

## 代表性结果

rank/fill-aware 联合 pivot 和持续低维核心启用后，再结合紧凑基双右端 FTRAN
与 CSC/CSR 混合稀疏定价，固定 30 题总时间从 46.532680 秒降到
40.909731 秒，30/30 保持正确。大型退化模型的 pivot 变化如下：

| 模型 | 调整前迭代数 | 当前迭代数 | 当前时间 |
|---|---:|---:|---:|
| `cre-b` | 18559 | 15847 | 8.213842 s |
| `cre-d` | 10850 | 10316 | 4.306680 s |
| `stocfor3` | 17107 | 16912 | 12.330118 s |

三题迭代总数和 30 题总时间均下降。`qap8`、`pilot.we` 等中小型基保持原确定性路径，分别以 5501 和 2265 次迭代通过。

测试对照可设置 `LP_SIMPLEX_DISABLE_PAN=1`，它只关闭退化控制，不改变其他
dual revised simplex 组件。

eta update 按方向结构自适应存储。Pan 紧凑基使用稀疏索引和值；大型完整基仅在方向非零数不超过基阶数八分之一时使用稀疏形式，其余情况保留 BLAS
稠密更新。该选择不截断数值，只省略严格等于零的项。

## 参考文献

Ping-Qi Pan, *A Basis-Deficiency-Allowing Variation of the Simplex Method for
Linear Programming*, Computers & Mathematics with Applications, 36(3),
33–53, 1998, received November 1997 and accepted December 1997,
DOI: `10.1016/S0898-1221(98)00127-8`.

Ping-Qi Pan, *A Primal Deficient-Basis Simplex Algorithm for Linear
Programming*, Applied Mathematics and Computation, 196(2), 898–912, 2008,
DOI: `10.1016/j.amc.2007.07.030`.



