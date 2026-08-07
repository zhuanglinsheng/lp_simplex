# 基于 CSC 的稀疏 dual revised simplex

本文介绍有界线性规划的稀疏 dual revised simplex，包括 CSC 矩阵、基矩阵
分解、dual steepest-edge 定价、Harris bound-flipping ratio test、换基更新和
数值校正。算法始终保持对偶可行，通过消除基本变量的原始可行性违反得到最优
基解。

## 1. 有界线性规划

考虑最小化问题：

\[
\begin{aligned}
\min\quad & c^T x \\
\text{s.t.}\quad & l^r \le Ax \le u^r, \\
& l \le x \le u,
\end{aligned}
\]

其中 \(A\in\mathbb R^{m\times n}\)。等式和单侧不等式都通过行界表示：

```text
a_i x  = b_i    row_lower[i] = b_i,  row_upper[i] = b_i
a_i x <= b_i    row_lower[i] = -INF, row_upper[i] = b_i
a_i x >= b_i    row_lower[i] = b_i,  row_upper[i] = INF
```

引入行变量 \(y\)：

\[
Ax-y=0,\qquad l^r\le y\le u^r.
\]

增广矩阵和成本向量为

\[
\bar A=[A\; -I],\qquad \bar c=(c,0).
\]

统一变量编号如下：

```text
0 ... n-1       结构变量，对应 A 的列
n ... n+m-1     逻辑行变量，对应 -I 的列
```

逻辑列 `n+i` 只有一个元素 `(row=i, value=-1)`。这些列由列访问接口即时生成，
不写入 CSC。

## 2. CSC 矩阵

CSC（Compressed Sparse Column）将同一列的非零元素连续存储：

```c
struct simplex_CscMatrix {
	int rows;
	int columns;
	int nonzeros;
	int *column_start; /* columns + 1 */
	int *row_index;    /* nonzeros */
	double *value;     /* nonzeros */
};
```

第 `j` 列位于：

```c
start = A->column_start[j];
end   = A->column_start[j + 1];
```

CSC 必须保持以下不变量：

- `column_start[0] == 0`；
- `column_start[columns] == nonzeros`；
- 同一列的 `row_index` 严格递增；
- 重复坐标已合并；
- 不保存精确零和低于导入阈值的系数；
- 求解期间不修改矩阵结构和数值。

### 2.1 Triplet 转换为 CSC

输入解析器先产生 triplet：

```text
(row[k], column[k], value[k]), k = 0 ... nz-1
```

转换算法为：

1. 统计每列 triplet 数量；
2. 对列计数做前缀和，得到 `column_start`；
3. 复制 `column_start` 为写指针 `next`；
4. 将 triplet 写入 `next[column]++` 指定的位置；
5. 在每列内部按 `row_index` 排序；
6. 合并相同行号，对重复值求和；
7. 删除合并后绝对值不超过导入阈值的元素；
8. 压紧 `row_index` 和 `value`，重新计算 `column_start`。

合并重复项必须发生在删除小元素之前。例如两个同坐标系数 `1` 和 `-1` 应被
合并并删除，而不能作为两个独立非零元参与点积。

### 2.2 列访问

CSC 的基本列操作为：

```c
double csc_column_dot(
		const struct simplex_CscMatrix *A,
		int column, const double *vector);

void csc_scatter_column(
		const struct simplex_CscMatrix *A, int column,
		double *dense, int *index, int *count);
```

列点积只访问该列非零元：

```c
sum = 0.;
for (k = A->column_start[j]; k < A->column_start[j + 1]; k++)
	sum += vector[A->row_index[k]] * A->value[k];
```

求解器使用统一的增广列接口：

```text
get_column(j):
    if j < n:
        return CSC column A[:, j]
    else:
        return one entry (row=j-n, value=-1)
```

FTRAN、定价和基矩阵提取都通过该接口读取列。

## 3. 基和变量状态

基由 \(m\) 个增广列组成：

\[
B=[\bar A_{b_0},\bar A_{b_1},\ldots,\bar A_{b_{m-1}}].
\]

核心状态包括：

```c
int *basis;                  /* basis position -> variable */
int *basis_position;         /* variable -> position, 非基本变量为 -1 */
unsigned char *status;       /* BASIC, LOWER, UPPER, FIXED, FREE */
double *primal_value;
double *lower;
double *upper;
double *reduced_cost;
```

`basis` 与 `basis_position` 必须互为逆映射。固定变量不参与定价。自由变量不能
普通地停在某个边界：它作为非基本变量时必须满足 reduced cost 为零，否则应在
dual Phase I 中移入基或进行等价拆分。

给定非基本变量值 \(x_N\)，基本变量满足：

\[
x_B=B^{-1}(-\bar A_Nx_N).
\]

对偶乘子和 reduced cost 为：

\[
B^T\pi=c_B,
\qquad
r_j=c_j-a_j^T\pi.
\]

对于最小化问题，非基本变量的对偶可行条件为：

\[
\begin{array}{ll}
x_j=l_j:& r_j\ge -\varepsilon_d,\\
x_j=u_j:& r_j\le \varepsilon_d,\\
l_j=u_j:& r_j\text{ 不受符号限制},\\
x_j\text{ 自由}:& |r_j|\le\varepsilon_d.
\end{array}
\]

对非固定、非自由的非基本变量定义：

\[
s_j=\begin{cases}
+1,&x_j=l_j,\\
-1,&x_j=u_j,
\end{cases}
\qquad
z_j=s_jr_j.
\]

于是对偶可行条件统一为 \(z_j\ge -\varepsilon_d\)。

## 4. 稀疏基矩阵分解

算法反复求解：

\[
Bd=a_q \quad\text{(FTRAN)},
\qquad
B^T\rho=e_p \quad\text{(BTRAN)}.
\]

求解器不显式计算 \(B^{-1}\)，而是维护带行列置换的稀疏 LU：

\[
PBQ=LU.
\]

符号分析使用近似最小度或 Markowitz 类排序抑制 fill-in；数值分解使用阈值
主元策略，在稀疏度和稳定性之间取舍。重新分解时，根据 `basis[p]` 提取 CSC
基列：结构变量复制 \(A\) 的对应列，逻辑变量生成一个 `-1`。

基模块提供统一接口：

```c
int basis_factorize(struct simplex_Basis *basis);
int basis_ftran(struct simplex_Basis *basis, struct simplex_Vector *rhs);
int basis_btran(struct simplex_Basis *basis, struct simplex_Vector *rhs);
int basis_update(struct simplex_Basis *basis,
		int leaving_position, const struct simplex_Vector *direction);
```

换基使用 product-form eta 更新维护基逆的作用。每个 eta 保存
\(d=B^{-1}a_q\) 与离基位置；FTRAN 按正序应用 \(E^{-1}\)，BTRAN 按逆序应用
\(E^{-T}\)。更新模块同时记录：

- 更新次数；
- 更新向量非零数；
- 最小主元；
- 因子增长；
- FTRAN/BTRAN 残差；
- LU 的非零元和 fill ratio。

达到更新次数上限、主元过小、因子显著增长或残差超限时，立即从当前
`basis` 重新构造 \(B\) 并数值分解。重新分解不是固定周期动作；它由更新质量
和线性求解残差共同触发。

## 5. 初始基和 dual Phase I

全逻辑基使用所有 \(-I\) 列，因此 \(B=-I\)。逻辑变量成本为零，此时
\(\pi=0\)，结构变量 reduced cost 为 \(r_j=c_j\)。结构变量的初始非基状态按
对偶可行符号选择：

```text
c[j] >= 0 and lower[j] finite  -> LOWER
c[j] <  0 and upper[j] finite  -> UPPER
lower[j] == upper[j]           -> FIXED
```

随后由 `y=A*x` 计算基本逻辑变量。行变量对行界的违反构成 dual simplex 需要
消除的原始不可行性。

当全逻辑基不能给出对偶可行状态时，dual Phase I 使用 cost shifting：

1. 由当前基求 \(\pi\) 和真实 reduced cost；
2. 对符号违反的非基本变量建立临时成本 shift，使 shifted reduced cost 满足
   对偶可行条件；
3. 以 shifted cost 运行同一个 dual simplex 内核；
4. 每次换基后收缩可移除的 shift，并更新 Phase I 对偶不可行度；
5. 所有 shift 回到零且真实 reduced cost 对偶可行时结束 Phase I；
6. 如果必要 shift 无法通过允许的 pivot 消除，则返回 dual Phase I 失败；只有
   同时取得原始可行解或无界射线时才能报告原问题无界，不能仅凭对偶不可行
   推断最终状态。

Crash procedure 在数值分解前从结构列和逻辑列中选择近似三角、主元较大且
fill-in 较小的列。Crash basis 必须经过完整的 primal/dual residual 计算，不能
仅根据组合结构宣布可行。

## 6. 离基选择：dual steepest-edge

Dual simplex 保持对偶可行，并选择违反自身上下界的基本变量离基。基本位置
\(i\) 的原始违反为：

\[
v_i=\begin{cases}
l_{b_i}-x_{b_i},&x_{b_i}<l_{b_i},\\
x_{b_i}-u_{b_i},&x_{b_i}>u_{b_i},\\
0,&\text{otherwise}.
\end{cases}
\]

Dual steepest-edge 使用权重

\[
\gamma_i=\|e_i^TB^{-1}\|_2^2
\]

并按归一化违反量选择离基位置：

\[
p=\arg\max_i \frac{v_i^2}{\gamma_i}.
\]

权重随换基递推更新。每次选定离基行并完成
\(B^T\rho=e_p\) 后，直接用

\[
\gamma_p\leftarrow\rho^T\rho
\]

精确校准本轮实际使用的权重，再以该值执行 DSE 递推。这样把权重修正集中在
访问到的行上，不需要在每次重新分解后额外执行 \(m\) 次 BTRAN。检测到对偶
状态漂移并重新 crash 时，可把权重重置为 Devex 初值，随后继续惰性校准。

若不存在超过 primal tolerance 的基本变量违反，同时当前基保持对偶可行，
则当前基最优。

设离基变量的目标边界为 \(\bar x_p\)，定义修复方向：

\[
\kappa=\begin{cases}
+1,&x_{b_p}<l_{b_p},\\
-1,&x_{b_p}>u_{b_p}.
\end{cases}
\]

## 7. BTRAN 和 CSC 定价

求解：

\[
B^T\rho=e_p.
\]

对每个非基本结构列，用 CSC 点积计算 pivot row 系数：

\[
\alpha_j=\rho^Ta_j.
\]

非基本逻辑列 `n+i` 的系数为：

\[
\alpha_{n+i}=-\rho_i.
\]

定义：

\[
g_j=s_j\alpha_j,
\qquad
z_j=s_jr_j.
\]

只有满足

\[
\kappa g_j<-\varepsilon_{pivot}
\]

的非基本变量能够修复离基行。它的对偶 breakpoint 为：

\[
\theta_j=\frac{z_j}{-\kappa g_j}.
\]

完整列定价只遍历 CSC 非零元和非基本逻辑列，主要成本为
\(O(\operatorname{nnz}(A)+m)\)。超宽模型可使用分块 partial pricing，但每个
分块必须保存上次扫描状态，并定期执行完整扫描验证没有遗漏合格变量。

## 8. Harris bound-flipping ratio test

Harris 两遍比值检验将数值稳定性和步长选择分开：

1. 第一遍使用放宽的 dual tolerance，确定允许的 breakpoint 区间；
2. 第二遍只考察该区间内的候选，优先选择绝对值较大的 \(|\alpha_j|\)；
3. 对具有有限两侧边界的候选执行 bound flipping；
4. 在累计翻转不能继续吸收离基行违反时，确定真正的入基变量 \(q\)。

变量从一个边界翻到另一个边界时：

\[
\Delta x_j=\begin{cases}
u_j-l_j,&LOWER\rightarrow UPPER,\\
l_j-u_j,&UPPER\rightarrow LOWER.
\end{cases}
\]

所有 flip 对基本变量的影响可合并为：

\[
h=\sum_{j\in F}a_j\Delta x_j,
\qquad
x_B\leftarrow x_B-B^{-1}h.
\]

`h` 通过 CSC 列稀疏累加，只执行一次 FTRAN。状态更新在 FTRAN 成功后统一
提交；若线性求解失败，可以丢弃 flip 列表并在重新分解后重试，不留下部分修改
的状态。

如果存在原始违反却没有合格候选，问题原始不可行。BTRAN 结果 \(\rho\)、离基
方向和变量界共同构成 Farkas certificate，返回前应使用未缩放模型验证证书。

## 9. FTRAN、换基和解更新

读取入基列并求解：

\[
Bd=a_q.
\]

理论上

\[
d_p=\rho^Ta_q=\alpha_q.
\]

若两者的绝对或相对差异超过更新容差，说明基分解误差过大。此时重新分解
当前基，重新执行 BTRAN、ratio test 和 FTRAN，不能继续使用旧候选结果。

完成 BFRT 的批量 flip 更新后，令：

\[
\Delta=\frac{x_{b_p}-\bar x_p}{d_p}.
\]

基本解更新为：

\[
x_{b_i}\leftarrow x_{b_i}-d_i\Delta,\quad i\ne p.
\]

离基变量设置为 \(\bar x_p\)，入基变量的新值为 \(x_q+\Delta\)。随后原子地
更新映射：

```text
basis_position[leaving] = -1
basis[p] = entering
basis_position[entering] = p
status[leaving] = target bound status
status[entering] = BASIC
```

最后用 \(d\) 和位置 \(p\) 追加 eta 更新，并递推 steepest-edge 权重。

## 10. 对偶更新

若入基变量的 breakpoint 为 \(\theta\)，对偶步长为：

\[
\tau=-\kappa\theta.
\]

更新对偶乘子和 reduced cost：

\[
\pi\leftarrow\pi+\tau\rho,
\qquad
r_j\leftarrow r_j-\tau\alpha_j.
\]

入基变量满足 \(r_q=0\)。离基变量的新 reduced cost 必须符合其目标边界的
符号条件。增量更新后，把数值落在 dual tolerance 内的符号违反压回零，但不能
把明显违反静默截断。

完整 reduced cost 校正通过以下步骤完成：

\[
B^T\pi=c_B,
\qquad
r_j=c_j-a_j^T\pi.
\]

其中 \(A^T\pi\) 逐 CSC 列计算。校正发生在重新分解后、检测到 drift 时以及
宣布最优之前。

普通 pivot 之间使用上述公式增量维护 \(x_B\)、\(\pi\) 和 \(r\)，不重复执行
完整 FTRAN、BTRAN 与全列定价。重新分解、检测到 dual drift、关键枢轴一致性
失败以及终止验收会触发全量校正；稳定的更新链不承担固定周期校正的成本。

## 11. 稀疏工作向量

FTRAN、BTRAN、列组合和更新向量使用“稠密值 + 稀疏索引”表示：

```c
struct simplex_Vector {
	double *value; /* length m */
	int *index;    /* touched indices */
	int count;
	int *mark;
	int generation;
};
```

首次写入位置 `i` 时将其追加到 `index`。清空向量只遍历已触碰位置，无需对
长度 \(m\) 的数组执行 `memset`。达到一定密度后切换为稠密循环，避免索引间接
访问成本。

稀疏向量删除元素时使用独立的 drop tolerance。该阈值必须显著小于 primal、
dual 和 pivot tolerance，并结合行列 scaling 解释；否则过早丢弃元素会改变
ratio test 候选集合。

## 12. 数值稳定性

不同判断使用不同容差：

```text
primal_feasibility_tolerance
dual_feasibility_tolerance
pivot_tolerance
factor_pivot_threshold
update_residual_tolerance
drop_tolerance
```

数值控制包括：

- 行列 equilibration scaling；
- 稀疏 LU 阈值主元；
- Harris 两遍 BFRT；
- 小主元拒绝和候选重选；
- cost perturbation 与最终去扰动；
- 基本解、对偶乘子和 reduced cost 的周期性重算；
- FTRAN/BTRAN 迭代改进；
- 基更新质量驱动的重新分解；
- 最优、不可行和无界证书的未缩放验证。

最优性验证至少计算：

\[
\begin{aligned}
e_p &= \max_i\{l_i-x_i,\ x_i-u_i,\ 0\},\\
e_d &= \max_{j\in N}\{-s_jr_j,\ 0\},\\
e_r &= \|Ax-y\|_\infty.
\end{aligned}
\]

只有三个误差都在各自容差内，才能返回最优状态。

## 13. 完整主循环

```text
build immutable CSC A and bound arrays
build crash basis and nonbasic bound statuses
factorize B
run dual Phase I until true reduced costs are dual feasible
compute and correct primal values, pi and reduced costs

while iteration < iteration_limit:
    p = dual_steepest_edge_choose_primal_violation()

    if p == NONE:
        recompute primal values, pi, reduced costs and residuals
        if primal feasible and dual feasible and residual feasible:
            unscale solution and verify
            return OPTIMAL
        refactorize B and continue

    target, kappa = violated_bound_and_direction(p)
    rho = BTRAN(unit_vector(p))

    scan nonbasic CSC and logical columns
    compute alpha, signed reduced costs and Harris breakpoints
    flips, q, theta = harris_bound_flipping_ratio_test()

    if q == NONE:
        build and verify Farkas certificate
        return INFEASIBLE

    apply all flips with one accumulated FTRAN
    d = FTRAN(column(q))

    if d[p] disagrees with alpha[q] or pivot is too small:
        rollback flips
        refactorize B
        continue

    update primal basic values
    update pi and reduced costs
    atomically replace basis column p by q
    append product-form eta update
    update dual steepest-edge weights

    if factor/update diagnostics require reinversion:
        factorize current B
        recompute primal values, pi and reduced costs

return ITERATION_LIMIT
```

实现 flip 和换基状态时应使用事务式工作区：所有候选修改先写入临时列表，只有
FTRAN、主元检查和 basis update 全部成功后才提交。这样，重新分解和候选重选
不会破坏 `basis`、变量状态与数值向量之间的一致性。

## 14. 复杂度和诊断指标

一次完整定价扫描约访问 \(\operatorname{nnz}(A)+m\) 个系数。BTRAN、FTRAN
和 basis update 的成本由 LU 因子及更新向量的非零数决定，无法仅用 \(m,n\)
表示；排序和主元选择导致的 fill-in 是主要性能变量。

求解日志至少记录：

```text
rows / columns / input nnz
simplex and dual Phase I iterations
pricing, BTRAN, FTRAN and update time
number of bound flips
number of refactorizations
L/U nonzeros and fill ratio
rejected pivots and numerical retries
maximum primal, dual and residual error
```

CSC 的作用是让原始矩阵列访问与完整定价保持稀疏；稀疏 LU、eta 更新、
dual steepest-edge 和 Harris BFRT 则共同决定整个 dual revised simplex 的效率
与可靠性。
