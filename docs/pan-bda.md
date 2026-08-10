# Pan/BDA 动态亏基算法

本文定义仓库中 `lp_simplex_ALGORITHM_PAN_BDA` 的实际实现。它对应 Pan 的
basis-deficiency-allowing generalized simplex，而不是 `src/dual/` 中仅用于
Harris 同率排序的 Pan-inspired anti-stalling 策略。

## 1. 求解形式与模型变换

公开模型先稀疏转换为

\[
  \min d^Tz+d_0,\qquad Az=c,\qquad z\ge0.
\]

转换覆盖自由变量、单下界、单上界、双侧界、固定变量以及等式、`<=`、`>=`
行：自由变量拆成正负两列；有限下界作平移；单上界用 \(x=u-z\)；双侧界除
平移外增加一条宽度等式和 slack；不等式增加带正确符号的 slack。每个变换列
记录原变量、比例和位移，求解后直接恢复原变量。

按照 Pan 的对偶记号，算法等价地处理

\[
  \max b^Ty,\qquad Ay=c,\qquad y\ge0,
\]

其中 \(b=-d\)。维护的活动列集 \(B\) 始终线性无关，但列数
\(k=|B|\) 只要求 \(k\le m\)，不要求是传统的 \(m\) 阶方基。

## 2. Phase I：可靠可行亏基

Phase I 求解非负最小二乘

\[
  \min_{y\ge0}\frac12\|Ay-c\|_2^2.
\]

实现采用 Lawson--Hanson active-set NNLS：以残差 \(r=c-Ay\) 为梯度，选择
归一化相关性 \(a_j^Tr/\|a_j\|_2\) 最大的列进入 passive set；只有通过
Cholesky 秩检验的列才能进入。无约束最小二乘解若含非正分量，就沿当前可行解
到该解的线段前进，所有先触及零的列同时删除。于是 Phase I 结束时直接得到

\[
  Ay=c,\quad y\ge0,\quad B\text{ 满列秩},\quad |B|\le m.
\]

若 NNLS 的 KKT entering 条件已经满足而残差仍超出
`primal_tolerance * (1 + ||c||)`，返回 `Infeasibility`；Phase I 和 Phase II
共享同一个迭代预算。

## 3. Phase II：Pan 动态亏基主循环

给定可行的 \((B,y_B)\)，首先求满足

\[
  B^Tx=b_B
\]

的最小二范数乘子 \(x=B(B^TB)^{-1}b_B\)。对非活动列计算违反量

\[
  v_j=b_j-a_j^Tx=-d_j-a_j^Tx.
\]

默认定价选择原始违反量 `v_j` 最大的列（Dantzig violation pricing）；它在
大型退化实例上明显减少了活动基交换次数。设置
`options.pricing = lp_simplex_PRICING_PAN_NORMALIZED` 可切换回
`v_j / ||a_j||` 归一化定价。若不存在超过 `dual_tolerance` 的违反，则当前
可行解最优。

对进入列 \(a_p\) 解最小二乘投影 \(B\delta\approx a_p\)：

1. 若投影残差大于秩阈值，则 \(a_p\notin\mathcal R(B)\)。将它以零值加入
   \(B\)，基维数从 \(k\) 增至 \(k+1\)，原始可行点不变。
2. 若 \(a_p=B\delta\)，可行方向为
   \(y_p=\theta,\ y_B(\theta)=y_B-\theta\delta\)。若所有
   \(\delta_i\le0\)，目标可无限改善，返回 `Unboundedness`。
3. 否则执行最小比值检验
   \[
     \theta^*=\min_{\delta_i>0}\frac{y_i}{\delta_i}.
   \]
   所有达到同一最小比值的集合 \(Q\) 同时离基，再加入 \(p\)。新基维数为
   \(k-|Q|+1\)，因此退化时可以下降而不是用任意零值列补成方基。

扩基、单列交换和多列同时删除都是普通路径，而不是停滞达到阈值后才启用的
补丁。这正是 Pan 动态亏基的核心。

## 4. Multifrontal Householder 锚点与动态正交更新

稀疏后端将归一化活动基 \(\widehat B\) 保持为 CSC，并维护真正的

\[
  \widehat B=QR,
  \qquad Q^TQ=I,
  \qquad \widehat b_i=b_i/\|b_i\|_2,
\]

其中 \(Q\) 由两层正交变换组成：SuiteSparseQR multifrontal sparse-Householder
锚点，以及锚点之后的局部 Givens/尾部 Householder 变换日志。扩基直接计算
\(Q^Ta\)，对未消去尾部构造稳定 Householder，并把反射向量保留在 Q 中；不再用
\(Q=\widehat B R^{-1}\)，因此最小二乘不退化为 \(R^{-1}R^{-T}\) 半正规方程。

减基从 \(R\) 删除目标列，用稳定 Givens 沿受影响的上 Hessenberg bulge 做局部
retriangularization；同一旋转以常数大小记录到 Q 日志，不扫描原始行空间。
`R` 的非零元采用稳定整数句柄节点池，同时维护排序行链和列链。行旋转先把两个
模式线性归并到预分配的连续 scratch，再原位复用原节点；只有真正产生 fill 时
才从 free-list 取得节点。相邻列交换只访问两列非零行的并集，单边非零节点直接
改变列归属，不再删除、搜索并重插。节点还将活动态 `row` 与空闲态 `free_next`
联合存储，并把数值放在首字段，使节点由 40 字节压缩到 32 字节。由 `R` 模式
增量维护消去树父链、前沿大小和
dirty ancestor path，因此结构修改只标记修改列到根的路径。

链式双视图适合结构修改，却不适合重复三角求解。后端因此另维护按结构代次失效
的 packed CSR 数值视图：一次活动基修改后的首次求解线性打包，此后的最小二乘、
最小范数和迭代精化都在连续 `column/value` 数组上运行，并缓存对角元。容量按
几何增长且在后端生命周期内复用；除容量增长外，求解与 Givens 热路径不做临时
分配。这是
“结构所有权”和“数值遍历”两种数据结构的职责分离，不是密度阈值调参。

Q/Qᵀ应用按正确方向组合 multifrontal Householder 锚点和局部变换。SPQR 锚点
内部使用超节点前沿、块 Householder/WY 与 BLAS-3。后端累计实际局部变换应用
工作量，并与上一次 SPQR 报告的真实 factor flop count 比较；只有可消除的日志
成本达到 multifrontal 重构成本时才重建锚点。该决策不依赖更新次数或用户阈值。
三角方程后向误差若超过维数缩放的机器精度界，也会请求稳定锚点重构。

监测数据包括 Householder 正交恒等式误差、三角方程后向误差、`R` fill、最大
前沿、局部 retriangularization、锚点重构和 SPQR 实测 flop count。所有算法
残差仍由原始稀疏列独立重算并精化；最终 SPQR certification 保持不变。

最终成功返回前，算法在原标准形上独立重算 \(c-Ay\) 和非负性。若维护值未
通过认证，则用 SuiteSparseQR 对当前亏基做 rank-revealing 稀疏正交
恢复；只有 QR 解满列秩、非负且重新计算的残差通过容差时才接受，否则返回
`PrecisionError`。QR 不会添加虚构基列，也不会把亏基补成方阵。

SuiteSparse 不可用时，构建系统保留原来的归一化稠密 Gram/LAPACK fallback；
`-DLP_SIMPLEX_USE_SPQR=OFF` 可显式测试该路径。动态 QR 普通更新不调用 CHOLMOD；
SuiteSparseQR 仍作为最终独立认证与恢复后端。

## 5. 代码结构与使用

| 文件 | 职责 |
|---|---|
| `simplex_pan_standard.c` | 原模型到非负等式标准形及解恢复 |
| `simplex_pan_phase1.c` | active-set NNLS Phase I |
| `linalg/pan_dynamic_sparse_qr.c` | 稳定节点池、行/列双视图、packed 求解缓存、Givens 事件和稳定 downdate |
| `linalg/pan_multifrontal_qr.c` | sparse-Householder 锚点、局部正交日志、消去树/前沿和成本模型 |
| `simplex_pan_basis_spqr.c` | 活动基适配、真正 Q/Qᵀ求解、后向误差和 SPQR 认证 |
| `simplex_pan_basis.c` | SuiteSparse 不可用时的稠密兼容后端 |
| `simplex_pan.c` | Phase II 定价、扩基、最小比值、多列离基、终止状态 |

公开选择方式：

```c
struct lp_simplex_Options options;
lp_simplex_default_options(&options, lp_simplex_ALGORITHM_PAN_BDA);
options.iteration_limit = 100000;
lp_simplex_solve(model, &options, x, &result);
```

默认使用 `lp_simplex_PRICING_DANTZIG`。需要复现原归一化轨迹时：

```c
options.pricing = lp_simplex_PRICING_PAN_NORMALIZED;
```

命令行差分诊断：

```sh
./build/tools/test_netlib --algorithm pan-bda data/netlib/feasible/afiro.mps
./build/tools/test_netlib --algorithm pan-bda --criteria normalized \
    data/netlib/feasible/degen2.mps
./build/tools/test_netlib --algorithm pan-bda --no-presolve \
    data/netlib/feasible/degen2.mps
```

设置 `LP_SIMPLEX_PROFILE=1` 会打印总时间、Phase I+II 迭代数、最终
`rank/rows`、后端名称、总分解次数、Phase I 的迭代/分解拆分、bordered extension、
downdate、Givens 旋转、增量符号更新和迭代精化次数。

## 6. 当前版本的验证与效率快照

下表是在 Apple Silicon、Release 构建、默认 \(10^{-7}\) 容差下的单次诊断
快照，用来说明算法行为而非跨机器性能承诺：

| 模型 | Presolve | Pan 迭代 | 最终动态基 | Pan 时间 | Dual revised 时间 | 原模型最大 primal residual |
|---|---:|---:|---:|---:|---:|---:|
| `afiro` | 开 | 42 | 41 / 50 | 0.00145 s | 0.00034 s | `1.42e-14` |
| `adlittle` | 开 | 266 | 146 / 147 | 0.01126 s | 未测 | `3.05e-6` |
| `degen2` | 开 | 1392 | 814 / 855 | 约 1.28--1.55 s（热态） | 0.01885 s | `2.00e-7` |

`degen2` 的 814/855 直接表明活动基确实允许亏秩 41。Presolve 回代和
outward-rounded 界会把原模型残差放大到约 `2e-7`；此前关闭 presolve 的差分
测试可达到 `2.29e-14`。类似地，`adlittle` 开启 presolve 时本次回代残差为
`3.05e-6`，关闭后为 `2.27e-13`，所以该差异不应归因于 Pan 的动态基更新。

`degen2` 的当前 Householder 路径执行 1392 次 bordered extension、578 次
downdate、116336 次 Givens 和 174 次成本/误差驱动锚点重构。稳定节点原位复用
加 packed 三角求解把 Release 热态从节点全量重建版本的约 50 秒降到
`1.28--1.55` 秒。正交恒等式误差约 `4.05e-15`，最终三角后向误差约
`1.75e-16`；优化前后迭代数、fill、锚点数和误差监测值一致。

`degen3` 的 4000 次截断画像达到 rank 2328：4000 次 bordered extension、1672
次 downdate、约 99.4 万次 Givens 和 398 次 multifrontal 锚点重构。稳定节点版
当前约 40.25 秒（节点压缩前约 43.1 秒），正交误差约 `1.11e-14`，后向误差约
`2.1e-16`。一次把节点池改为
行优先物理排列的实验反而恶化到约 55.8 秒，原因是列交换失去局部性，已经撤销。
这说明同一节点排列无法同时优化行、列两个方向；大型前沿的下一层优化应是独立
的分段 packed 数值行/超节点存储，在一次活动基更新边界与稳定符号节点同步。
因此当前实现虽已消除全量行重建和重复链式求解，但仍不应宣称“极致完成”。

专用单元测试覆盖退化最优、等式、不可行、零列无界以及自由/上下界变换；
`test_simple_examples` 对三种求解路径运行全部小型 MPS，CTest 另含 Pan/BDA
的 `afiro` smoke test。

## 7. 参考

- P.-Q. Pan, “A basis-deficiency-allowing variation of the simplex method,”
  *Computers & Mathematics with Applications*, 1998。
- 现代 BDA 论文中的 Algorithm 1 明确说明其与 Pan primal BDA 等价；本文用它
  交叉核对了动态扩基、最小比值和多列同时离基的流程。
