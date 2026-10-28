# qap15：等式依赖与 presolve 分析（2026-10-28）

[开发记录](README.md) | [Netlib 正式结果](../benchmarks/netlib.md) | [数值政策](../numerical-policy.md)

初始诊断发现可精确验证的 632 条独立冗余等式。第 1—6 节记录离线分析，第 7—9 节记录实现、回归与准入决定。依赖约化因完整求解退步改为显式实验选项，默认关闭；第 9 节记录当天较早的默认 AUTO Netlib **114/114**，qap15 用时 **53.38 秒**，见[该轮历史报告](../benchmarks/history/2026-10-28/netlib.md)。后续 AUTO 动态特征修订及新一轮验证见 [PDS 调度复查](pds-performance-analysis-2026-10-26.md#15-auto-动态特征复查2026-10-28)。独立 Dual 的 qap15 当时仍在 90 秒超时；混合 AUTO 的成功不表示独立 Dual 已解决该例。

## 1. 输入及实验边界

Apple M1、Release、KLU/SPQR/GMP 开启，当前 main 工作区。共享库 SHA-256 为 `4a179cb99f77274e18fc1c977c3260a23d348e7a5d613720239e03ee290b4369`。有限步实验使用生产 Dual、10000 工作、默认容差、单线程，开启输出型 `LP_SIMPLEX_PROFILE`。两轮都用同一个共享库，没有更改数值环境控制。

HiGHS 使用现有 SciPy binding 的 1.12.0，单线程、presolve on、serial dual simplex、原始/对偶容差 `1e-7`、随机种子 0、90 秒保护。只进行了本案例的诊断，没有重跑其他求解器或全套基准。该次 HiGHS 在 55252 次迭代后超时，日志目标约 888.133；它仍有原始不可行量，不是最优解或可行目标。90.09 秒求解时间除以迭代数约为 1.63 ms，包含启动阶段，不是独立内核微基准。

以下时间均是单次诊断。机器未独占，且 PROFILE 有观测开销；不同输入的基路径、退化和 FTRAN 次数也不同，不能把每次主元的平均值视为同一个计算的指令级加速。

## 2. 数据结构

| 属性 | 实测值 |
| --- | ---: |
| 行 / 列 / 非零元 | 6330 / 22275 / 94950 |
| 等式行 | 6330 |
| 每行非零元 | 15 |
| 变量非负 / 原始显式上界 | 全部 / 全部无穷 |
| 度为 30 的列 | 225 |
| 度为 4 的列 | 22050 |
| 正 / 负单位系数 | 88650 / 6300 |
| RHS 为 1 / 0 的行 | 30 / 6300 |
| 零成本列 | 6525 |

可将 225 个高连接变量记为 `x(i,k)`，其余变量记为两个不同位置的配对 `y(ik,jl)`；这只是对输入矩阵的解释，诊断检测使用系数和邻接，不使用行列名称。每个配对变量在两个端点的边际等式中合计出现四次。

30 条分配等式分别约束每组 x 之和为 1。每个 x 另连接 28 条边际等式，均为 `sum(14 个 y) - x = 0`。这 28 条可以分成两组 14 条：两种方向的边际和描述同一组配对变量。这同时解释了局部等式依赖、跨组依赖及零成本/退化结构；不能仅由这些计数断言某种反循环策略一定更快。

## 3. 两边现有 presolve 做了什么

HiGHS 单独 presolve 约 0.114 秒，输出仍是 6330 / 22275 / 94950，成本、行界、变量界均未改变。本次默认设置下没有获得可用于解释求解差距的大规模约化。[HiGHS 1.12.0 的 HPresolve 实现](https://github.com/ERGO-Code/HiGHS/blob/v1.12.0/highs/presolve/HPresolve.cpp)包含依赖等式搜索及资源保护；源码具备规则与这次实际删行是两回事。

我们约 0.052 秒，未删行列，但收紧了 22275 个界。分配等式及非负性推出 `x <= 1`，局部边际和进一步推出 `y <= x <= 1`。因此隐含常数上界传播已有覆盖，不能把它重复算作待实现收益。当前缺口是跨行等式依赖识别。

所有行宽 15，变量度为 4 或 30，现有单行、重复行、doubleton 和无填充代换规则无法直接识别这些长的守恒关系。尤其不能以“变量度只有 4”为理由随意消元：宽行代换仍会产生填充，还要保留被消变量的界。

## 4. 可验证的 632 条独立依赖

构造了 661 个稀疏整数见证，共 19800 个非零系数：

1. **225 个局部边际恒等式。** 每个 x 的第一组 14 条等式之和，减去第二组 14 条之和，所有 y 与 x 都精确相消。
2. **435 个跨组恒等式。** 对 30 个分配组的每一对 a、b，取“端点属于 a、目标属于 b”的边际行和，减去反方向行和，再加分配行 a、减分配行 b。配对变量、共享端点和 RHS 全部相消。
3. **1 个总量恒等式。** 两类分配行的总和相减。

将见证组成 W，用 int64 稀疏乘法逐条验证 `W A = 0`、`W b = 0`。本例中计数及系数很小，整数计算没有溢出。再在素数 `1000000007` 下消元，找到 632 条独立见证及对应的 632 个主元行位置。

这里模运算不是近似可行性判据：见证与输入的乘积已经用整数精确验证；选中子式模素数非零，证明该整数子式在有理数上可逆。因此删除这些主元位置对应的原方程后，被删除方程可由保留方程严格推出，变量和目标均未改变。只删除有证据的集合，不宣称模秩发现了所有可能依赖。

| 输入 | 已验证可删除行 | 删行后的行 / 列 / 非零元 |
| --- | ---: | --- |
| qap8 | 170 | 742 / 1632 / 5936 |
| qap12 | 398 | 2794 / 8856 / 33528 |
| qap15 | 632 | 5698 / 22275 / 85470 |

相同邻接检测和见证算法在三例上通过，代码没有写入案例名、15 或预期删除数。它证明该结构族上的适用性，尚不能替代更广泛的通用 presolve 回归。

独立浮点 SPQR 对 qap15 的转置给出秩估计 5698，与上述结果一致；但耗时约 7.34 CPU 秒，峰值 CHOLMOD 记账内存约 671 MB、R 约 889 万非零元。它仅作诊断交叉检查，不应把全模型 QR 或浮点秩阈值直接放进默认删行规则。

## 5. 填充与后续求解成本

只删除上述行：行数减少约 9.98%，非零元减少约 9.98%，没有矩阵填充。随后再用我们的现有 presolve，另外删除 3 行、543 列，剩余 **5695 行 / 21732 列**；其中 540 列走已有 singleton projection，说明消除依赖也能暴露后续规则的机会。HiGHS 对同一删行模型的额外约化是 3 行、3 列，剩余 5695 / 22272 / 85425。

相反，若直接用每个共享 x 的一条边际式消掉 x，并向其余关联行代入，本例变成 **6105 / 22050 / 164640**。变量虽少 225 个，非零元反而增加约 73.4%。这是明确的填充反例，不能用“删除列数更多”来评价聚合器。

| 10000 次 Dual 主元的诊断 | 原输入 | 632 条依赖删行后，再执行现有 presolve |
| --- | ---: | ---: |
| Dual 内核 CPU / 秒 | 30.316746 | 21.152248 |
| 内核平均 ms / 主元 | 3.0317 | 2.1152 |
| FTRAN 秒 / 调用数 | 23.356441 / 33548 | 14.041780 / 23049 |
| BTRAN 秒 / 调用数 | 2.857243 / 10022 | 3.967909 / 10022 |
| 分解秒 / 次数 | 0.108599 / 22 | 0.142394 / 22 |
| 定价及 ratio / 秒 | 1.001838 | 1.499732 |
| flat / zero-theta 主元 | 2652 / 114 | 2370 / 590 |

两个实验都在 10000 工作处正常返回迭代上限，均未解到最优。内核 CPU 比约 1.43，但一部分来自 FTRAN 调用次数及轨迹改变。不能据此声称完整 qap15 已解决，也不能把退化总量改善与所有退化指标混为一谈：删行后 zero-theta 反而更多。

原输入 FTRAN 占内核约 77%，定价约 3.3%，重分解约 0.36%。直接优化定价或继续减少分解次数的可得收益有限；后续应检查 FTRAN 的 eta 链访问、稀疏到稠密切换、DSE 附加 RHS 和 bound-flip RHS 次数。[HiGHS 的 Dual 实现](https://github.com/ERGO-Code/HiGHS/blob/v1.12.0/highs/simplex/HEkkDual.cpp)也使用 DSE 更新及按密度选择的稀疏定价，不能把已有机制简单列为我们完全缺失的能力。

## 6. 应落地的通用方向

优先补充**有稀疏整数见证的等式依赖约化**，放在现有局部规则稳定之后，成功删行后重新进入原有队列。检测基于精确 signed-unit 系数、共享列和边际邻接；不按数据集、名称、目标答案或墙钟选择策略。先支持便宜的局部关系，再以明确的内存、候选数、消元工作和系数增长预算搜索跨组关系。上述 Python 模运算实验不是生产实现成本，不能承诺当前 3.61 秒消元成本直接可接受。

删除一个集合前须验证完整的系数关系及 RHS，并验证对应主元子式独立。近似平行、RHS 不匹配、整数溢出、预算或分配拒绝时保留原模型；不得据此返回不可行。支持一般浮点系数应单独设计误差界，不能把本例的整数论证替换成固定 epsilon。映射与工作区先准备再发布，后处理仍审核完整原模型。

删行保持所有原变量及目标；任何 reduced primal 仍需过全部原方程，reduced dual 可在删除行处补零形成原模型乘子，但跨阶段基交接还需要正确处理原模型的逻辑列/基维度，不能把变量映射正确误当成基映射完成。后续暴露的 singleton/projection 继续使用现有 journal 及逆序恢复。

验收需包括 RHS 扰动使依赖失效、近似单位系数、行列重排和重命名、一般无该结构的模型、内存/工作拒绝，以及完整 Netlib/不可行集和原模型 primal/dual 审核。有限步收益仅决定优先级；完成实现后，仍须按原 90 秒及 300000 工作验证 qap15 完整求解。

初始诊断仅保留这份分析，模型、见证数组、原始日志、临时脚本和诊断工具均已清理。

## 7. 生产候选与完整求解验证

新增 [simplex_presolve_dependencies.c](../../src/presolve/simplex_presolve_dependencies.c)，从精确 signed-unit 行的共享列和邻接关系识别守恒见证，不读取名称、参考目标或模型规模分段。整数 RHS 的准入上限是 `1073741823`，加法检查 `long` 溢出。完整整数乘积验证先于模素数 `65521` 下的独立性选择；两项职责没有合并。

可选自有存储总上限为 32 MiB，搜索工作为 `512*(nnz+m+n)`，每次预解只搜索一次。未满足整数结构、工作/内存不足或算术拒绝时保留尚未证明可删的行；此前完整证明的独立集合可保留。搜索结束统一发布删行，更新列度，重新进入已有固定点队列和代换流程。

删行选择优先保留全同号的求和等式，因为它们能直接表达总量和隐含界；依赖验证标准不变。该顺序也经过行列重排和等式换号测试。最终 qap15 约化为 **5698 行、21672 列**：删除 632 条依赖等式，随后现有 singleton projection 删除 603 列，预解约 0.128 CPU 秒。不同删行集合会改变后续代换和数值路径，不能只比较最终行数。

同时修正 [simplex_basis.c](../../src/basis/simplex_basis.c) 的重分解成本门槛：原来的大行数保护为内部 LU 的高成本而设，却也禁用了 KLU 的已有确定性成本规则。现在 KLU 在累计 eta 应用工作达到因子工作估计时请求重分解；大基内部 LU 的策略、弱主元保护、更新长度和最终残差要求不变。在保留求和等式的排序调整之前，10000 次主元诊断由约 21.15 降至 11.30 内核 CPU 秒，FTRAN 从约 14.04 降至 5.49 秒，分解则从 22 增至 838 次、约 2.06 秒。这说明旧分解次数少并不意味着总成本低；该诊断属于中间候选，不冒充最终排序的完整性能数据。

最终候选共享库 SHA-256 为 `1ab0268fdea354d574ab5716026f2e1ddb3b3953e42961e1a3dd388fd7dc0f1e`。使用现成 benchmark 的原模型审核和 90 秒/300000 工作保护，分别测试 AUTO 与纯 Dual，**两者仍均超时**。每轮完成记录均确认源码、二进制、共享库、模型和对比求解器在运行期间不变。没有扩大限时、放宽容差或把可行点作为最优点接受。

最终排序的 10000 次主元诊断为 **8.754066 内核 CPU 秒**，其中 FTRAN 4.162210 秒 / 24301 次，BTRAN 0.802868 秒 / 10913 次，分解 1.769996 秒 / 913 次，ratio 0.522036 秒。相对初始未约化的 30.316746 秒约为 3.46 倍有限步内核提速；两次均未到最优，轨迹及每步工作不同，不能外推完整求解倍数。

验证范围：最终 Release 182/182 检查通过；最小依赖 ASAN/UBSAN 的 7 项相关检查在最后一次排序调整后通过。另有首版的完整最小依赖 ASAN/UBSAN 155/155 通过，该完整轮次不覆盖随后所有改动。新增测试独立比较删除前后的仿射行空间，覆盖两个规模、行列重排、换号、微小 RHS/系数扰动、等于模素数的 RHS 错误、部分搜索及内存/工作拒绝。基测试覆盖长内部 eta 链与大 KLU 基的成本触发重分解，并以独立矩阵残差验证数值。

因此已经落地的是安全的结构约化及一项有成本依据的后端政策修正；qap15 的完整求解问题仍未解决。正式基准文档没有被这些单案例候选报告覆盖。运行日志、临时模型、实验构建及工具均清理；没有暂存、提交或推送。

## 8. DSE 与共享定价存储复查

继续诊断发现两个独立问题，均按通用数值与存储不变量修正，没有增加案例名称判断或调整 AUTO 调度。

**DSE 的离基权重需要重新核对。** 前 10000 步的 fresh/stored 平方范数比最大约 `1.000143674`，没有严重低估；延长观察后出现 7 次低估超过四倍，最大比例 `1.6e13`。其中递推值触及 `1e-12` 底限，当前 BTRAN 平方范数约为 16。因此仅检查开头的有限步不足以排除后续递推消去误差。[HiGHS 1.12.0 的离基选择](https://github.com/ERGO-Code/HiGHS/blob/v1.12.0/highs/simplex/HEkkDual.cpp#L1386-L1448)也在 BTRAN 后核对平方范数，并拒绝低于 fresh 范数四分之一的旧权重。

现在复用已经求出的 BTRAN，写回当前范数；严重低估时更新可行性堆并重新选择，不推进主元工作或改动点、成本、界。同一基内限定最多所有行各修正一次，再做最终选择；非有限或非正范数返回精度失败，不能据此认定不可行。新增上下界两种离基方向的测试，覆盖堆选择改变、点和 merit 不变、四分之一边界及非有限/零范数拒绝。

**稠密与稀疏定价共用数组时，稠密写入没有登记坐标。** 稀疏 clear 只清登记坐标，导致稠密到稀疏切换后未触及位置残留旧值。旧测试只检查当次稀疏支撑，未检查所有未触及位置；补上完整数组与独立参考向量的比较后，修复前测试失败。修复将稠密 CSC 写入也通过已有 sparse-vector setter 登记，保持清理契约。修复后该回归通过；20000 步生产 Dual 诊断的 115 次稀疏定价完整校验全部差异为零、修复坐标数为零。此前多次较大校验差异的存储原因得到消除，但不能由此推断所有数值困难均已解决。

该 20000 步诊断仍返回迭代上限，CPU 约 33.04 秒，原模型原始残差约 `1.92e4`；权重路径改变后记录 14 次严重低估，最大比例仍为 `1.6e13`。这不是最优解，也不能与旧的不同迭代路径直接比较速度。两项修正后的 AUTO 单案例按原 90 秒保护重跑，仍超时。

另做了原生 Simplex 可行基诊断：读取此前求解器保留的 0/1 可行点（目标 1536），映射并检查到约化模型，再用逻辑基交给已有 primal cleanup，生产 context 保持默认。这是复用外部保存点的诊断，不能作为端到端性能。KLU 清理在 1108 步中出现 900 个零步、532 个 Bland 步，最终因基数值失败拒绝；内部后端重试在 124 步中出现 94 个零步、60 个 Bland 步及 3 次重复状态后拒绝。对偶残差并不要求单调下降，其升高本身不是错误证明；这里的失败证据是实际主元/基拒绝及重复状态。

仅改变同强度主元的固定基变量优先级仍遇到循环，试验已撤回。打开已有界扰动的独立诊断也未得到可接受结果，没有因此打开生产开关。后续重点应是高度退化可行顶点上的稳定基构造与原目标 primal 枢轴，而非继续压缩已经约 0.13 秒的 presolve，或通过延长限时宣称成功。

当前两项修正的 Release 182/182、最小依赖 ASAN/UBSAN 7/7 相关检查通过；格式与架构检查通过。共享库 SHA-256 为 `fb0ab3f1756c76be6dda75256bb4c0117417144d0a59d12c02c1e1f3c4cf0246`。完整 Netlib 回归结果见本节后续记录；正式报告不会用单案例或未完成的诊断覆盖。

## 9. 完整回归与默认准入

两项定价修正、KLU 政策及依赖约化全部启用的完整 Netlib 114 例回归得到 **112 PASS、2 TIMEOUT**：qap12 和 qap15 均触及原 90 秒限时。全部完成记录确认源码、二进制、共享库、模型及竞争求解器在运行期间未改变。其余 112 例通过原模型约束、变量界、复算目标和参考目标审核；这里没有重跑其他求解器。

随后对 qap12 保持同限时、同生产配置逐项隔离：

| 候选 | qap12 |
| --- | --- |
| 依赖约化与两项定价修正均启用 | TIMEOUT，90 秒 |
| 仅停用 DSE 范数修正，保留依赖约化与存储修复 | TIMEOUT，90 秒 |
| 保留 DSE 与存储修正，仅停用依赖约化 | PASS，63.16 秒 |

这些对照说明依赖约化后的路径是需要解决的性能问题；不能用更少行列或前 10000 步更快来证明默认收益。它不否定精确删行的代数证明，也没有证明所有机器上都会有相同耗时，但足以拒绝本轮默认启用。

因此将依赖约化收回为**显式实验选项**，默认对所有模型关闭，不按 qap12/qap15 名称或规模特判。新增内部 `simplex_presolve_run_configured` 接收已经解析的布尔值，普通私有预解和正式公开入口都关闭；只有实验入口读取 `LP_SIMPLEX_PRESOLVE_DEPENDENCIES`。数值内核不读取环境，也不接收实验模式标签。独立见证/仿射行空间测试保留，并增加默认关闭、显式启用实际产生约化的集成检查；环境隔离测试验证正式入口忽略残留开关。

重现实验须清空其他数值实验变量，并同时保留生产 Dual 的初始成本及单 epoch 预设，例如：

```sh
LP_SIMPLEX_PRESOLVE_DEPENDENCIES=1 LP_SIMPLEX_INITIAL_COST_PERTURB=1 \
LP_SIMPLEX_SINGLE_PERTURB_EPOCH=1 build/data/lp_solve_mps \
  --experimental-environment --algorithm dual-revised MODEL.mps
```

最后准入候选的共享库 SHA-256 为 `0dfe7cc74826e7dec8e3e6b60cb8e8dae9b33133bff1d56a8b2efa4326249da0`。Release 182/182 检查通过，新增集成检查另行通过；最小依赖 ASAN/UBSAN 8/8 相关检查通过，覆盖环境隔离。最终默认配置重跑 Netlib 全部 114 例，**114/114 PASS**，原 90 秒、300000 工作及 `1e-7` 标准不变。源码、二进制、共享库、模型及竞争求解器的运行期间变更标志全部为否。qap12 为 59.24 秒，qap15 为 53.38 秒，dfl001 为 50.06 秒。qap15 的独立原模型审核残差为 `4.6337291e-8`，复算目标 `1040.9940074016308`，与参考目标误差 `3.3598369e-5`，通过原有相对误差标准；没有放宽精度。

已用完整 AUTO 轮次更新中英文 Netlib 报告及汇总，保留同机其他求解器的已有测量。后续纯 Dual 单案例诊断没有覆盖正式 AUTO 报告：同样 90 秒保护下仍超时。独立 PROFILE 复现约 54 秒成功，轨迹显示 Dual 探索消耗 16384 工作，随后可行性初始化消耗 330 工作，最终由已有的 `native-hybrid`／PDHG 阶段完成，总工作 67342。PROFILE 是诊断复现，正式耗时仍取未开启 PROFILE 的完整基准轮次。没有新增或修改 AUTO 调度来绕开纯 Dual 的困难；当前成功属于既有混合流程。独立 Dual 的完整求解仍需研究退化与基构造，不能用 AUTO 的通过数证明该路线已经解决。用户后续明确允许混合方法，当前验收以 AUTO 的通用性和完整基准成功为准。没有暂存、提交或推送。

## English summary

qap15 contains 6330 equalities, 22275 nonnegative variables and 94950 signed-unit entries. Exact integer witnesses certify 632 independent row dependencies without fill. The bounded optional rule and subsequent projections reduce it to 5698 rows and 21672 columns. However, a full trial with this rule enabled passes only 112/114 Netlib cases: qap12 and qap15 time out. qap12 still times out with only the DSE guard disabled, but passes in 63.16 seconds when only dependency elimination is disabled. The reduction therefore requires explicit experimental opt-in and is disabled by default for every model.

Retained changes restore cost-driven KLU refactorization on large bases, check optimistic DSE norms using the existing BTRAN, and register dense pricing writes in shared sparse storage. A full-array regression fails before the storage fix and passes afterward; all 115 pricing audits in a 20000-pivot diagnostic then agree exactly. Final Release passes 182 checks, the added opt-in integration check passes separately, and eight relevant minimal-backend ASAN/UBSAN checks pass.

The earlier same-day default AUTO build passes **114/114** Netlib cases under the unchanged 90-second, 300000-work and `1e-7` limits. qap15 passes in **53.38 seconds**, qap12 in **59.24 seconds** and dfl001 in **50.06 seconds**. The complete AUTO results update the bilingual reports; competitor measurements are reused. A separate pure-Dual qap15 trial still times out at 90 seconds. A PROFILE replay confirms that AUTO uses 16384 Dual work, 330 feasibility-initialization work, and finishes through the existing native-hybrid/PDHG route, with 67342 total work. These measurements belong to the earlier build archived above. Later AUTO feature work is recorded in the linked PDS analysis. AUTO success does not establish standalone Dual success; the current acceptance direction permits hybrid methods and requires general AUTO behavior across the complete suites. No changes have been staged, committed or pushed.
