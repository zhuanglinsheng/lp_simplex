# 开发记录

本页保留影响当前设计的决策、已修复缺陷和验证范围。使用接口请读[文档首页](../README.md)；调度、阶段和数值规则分别由 [AUTO](../auto-policy.md)、[BDA](../pan-services.md)和[数值政策](../numerical-policy.md)维护。

## 2026-10-24—26：求解系统重构

| 主题 | 最终决策 |
| --- | --- |
| 算法分层 | Dual、Tableau、BDA 独立；PDHG/NNLS/切锥/恢复/对偶修整属于内部角色，表示、oracle 和后端分别配置。 |
| AUTO | 主控独占运行状态；策略只返回提案，routes 处理准入与适配。阶段共用完成协议，八个组合入口归并为一个。 |
| 交接与因果 | 数学身份独立于履历；点、乘子、基分别审核，跨回调引用绑定 token。原因、父调用和产出来源不混同。 |
| 失败消费 | 有限路线游标与精确失败请求账本分工；PDHG 修复残差基线、NNLS blocked 集和面请求按各自状态决定重复。 |
| 配置 | 正式入口不读取数值环境；显式实验入口在根捕获配置。回调改变环境不改变当前数值行为。 |
| 目录与接口 | BDA 位于 `src/bda`，共享运行设施位于 `src/runtime`；原模型审核恢复与 presolve 代数映射分开命名。删除无独立语义的私有包装。 |
| 资源与提交 | 模型、稀疏缓存、journal、LU/QR 的存储检查与提交各归所有者；资源拒绝不能冒充不可行或奇异。 |
| 外部交付 | 可重定位 CMake 安装包，独立 C90/C++ 消费者验证四种公开选择及进度回调。 |

主要缺陷修复包括：近似平行行误判；非有限 Tableau Phase I 目标误分类；Singleton 快捷路径丢失根配置/作用域；不可逆缩放；基导入与分解资源类别丢失；恢复记录部分发布；履历槽复用导致错误来源；后端选点工作不是累计值；Dual 初始化尝试上限整数溢出。回归分别归在现有 owner 的测试中，见 [test/CMakeLists.txt](../../test/CMakeLists.txt)。

最终重构版本在 macOS arm64 的 Release/Debug-ASAN、可选 KLU/SPQR/GMP 开/关四配置通过 636 项测试（163、155、163、155）。风格、架构及安装包移位检查通过，构建无警告。该结果只对应当时的源码与配置，不自动证明后续改动或其他平台。

## 性能记录

性能以各次保存的版本、配置和数据为准；[基准报告](../benchmarks/README.md)保留逐例结果，[数据说明](../../data/README.md)定义审核与计时方法。本页不将不同候选结果交叉相除。

- AUTO 接受性检查：55 例各一次配对，110 次运行通过原模型审核，工作数相同；基线/候选总耗时 845.502863/868.379459 秒，候选慢 2.7057%，在当次预设 5% 范围内。不是重复全集或稳定加速证明。
- 原生面请求变化路径的早期样本保留了约 6.42% 的组件 CPU 增加及 11—26% 的墙钟增加。后续真实 BDA 对照中，CPU 成对比值中位数约 0.988263，墙钟中位数 1.1158285/1.5194234985 秒，候选慢约 36.17%。这些慢样本不能由较快案例或接近的 CPU 结果抹去，也不冒称最终结构版本的测量。
- 2026-10-26 用户明确取消收尾时的全集性能复测，接受历史版本结果作为参考；没有因此发布新的速度提升或当前版本无回退结论。此前约两倍加速的研究目标不属于架构重构的完成证明。

## 2026-10-26：文档与仓库整理

现行规则按主题收拢，功能索引去除重复验收表，开发记录合并为本页。保留算法推导、公开接口参考、生成的基准表和原始数据；清理旧构建、临时实验、Python 缓存及系统元数据。

整理前的完整源码/文档快照、逐轮日志和本地构建已另存于仓库旁的 `lp_simplex-archive-20261026`。这些机器相关证据不进入源码版本库；需要复现时用对应版本及原始配置，不能将本页汇总替代原始计时。

本次整理未改变生产代码、C 回归测试或 CMake 构建输入；复用上述四配置证据，并重跑安装包、架构、基准工具三项检查。新增文档链接检查接入 CI，验证本地文件及 Markdown 锚点。

## 2026-10-26：PDS AUTO 复查

[复查说明与反例](pds-auto-review-2026-10-26.md)记录加权综合指标对可行性停滞的观察盲点，以及增加修复所带来的跨数据集回退。候选未合入 main；179 项功能测试通过，完整性能验收因机器争用未完成。分层候选的小样本 qap15 配对仍有约 13.5% 墙钟中位数开销。[当天完整基线与实验记录](../benchmarks/history/2026-10-26/README.md)分别归档。

## 2026-10-26：PDS 分阶段性能归因

[深度分析](pds-performance-analysis-2026-10-26.md)区分预处理自身成本、削减能力及后续收敛变化。九例 presolve CPU 均不足 0.5 秒；完整 pds-70 诊断的主要成本在 PDHG。开关预处理对 pds-20／70 的影响不同；重复自身预处理只多删少量重复行，追加 HiGHS 预处理仍能显著缩减结构。未因此调整默认策略或替换正式基准结果。

## 2026-10-26：PDS 优化候选

[原则复核与验收](pds-optimization-2026-10-26.md)记录可行性停滞观察、有界修复及同次精确回滚后的列索引复用。四配置 704 项测试通过；qap15 两对控制均通过但墙钟中位数仍慢约 17.45%，PDS 仅有七个完整结果。用户要求稍后安排性能测试，候选尚未合入 main；日期归档保留成功、超时、中断与撤销方案。

## 2026-10-27：PDS 实现复核与暂停测试

[新增离线分析](pds-performance-analysis-2026-10-26.md#7-2026-10-27隐含上界中的搜索容差污染)确认隐含界的用户容差填充、全矩阵扫描与活动面结构，以及前置 Dual 的成本限制。当前候选在 main 未提交工作区，旧候选分支已删除；PDS／Railway 求解及队列已停止，20 倍尚未实现或验证。随后[矩阵自由修整候选](pds-performance-analysis-2026-10-26.md#10-2026-10-27矩阵自由乘子修整候选)让 qap15 通过控制，完整 Release／最小依赖 ASAN 共 352 项功能测试通过；dfl001 的额外成本仍须保留。没有发布新的完整性能报告。

## 2026-10-28：Simplex 聚合预处理与退化定价

[规则消融与候选验证](pds-performance-analysis-2026-10-26.md#11-2026-10-28沿-simplex-定位聚合缺口与退化定价)将后续优化主线转到 Simplex。HiGHS 消融确认主要预处理差距在聚合器；候选补充跨行隐含界、signed-unit 邻域的无填充消元及补偿 postsolve。pds-20／70／100 的结构均缩减，最初 pilot.we 回退已通过保守准入修复，最终 Release／最小依赖 ASAN 共 352 项功能检查通过。现成 Devex 实验避免 pds-20 原目标清理后的两次 Phase I 重启，但默认定价未改；20 倍仍待验证，PDS 全集／Railway 未运行，候选未提交。

## 2026-10-28：qap15 等式依赖分析

[结构分析及生产候选](qap15-presolve-analysis-2026-10-28.md)补充了有整数见证和有限资源的依赖删行、KLU 重分解成本规则及 DSE／定价存储复查。依赖删行的固定工作成本改善没有转化为完整求解优势，因 qap12 回退而默认关闭。完整默认 AUTO Netlib 验证为 114/114；不能把早期候选超时或纯 Dual 结果当成最终 AUTO 结论。后续调度特征修改需另做冻结回归。

## 2026-10-28：AUTO 阶段相关动态特征

[复查与反例](pds-performance-analysis-2026-10-26.md#15-auto-动态特征复查2026-10-28)区分结构准入、不可变输入审核和算法内动态观察。Phase-I 原始残差反弹不能单独证明 Dual 无进展；固定输入目标势须同时越过抵消和成本漂移尺度。扩大输入进展续跑权限曾引入 qap15 超时，修订为按阶段目标使用证据，未降低阈值或扩大预算。恢复次序和观察窗口仍是有明确边界的启发式；见 [AUTO 架构](../auto-policy.md#特征提取与在线重分类)。最新冻结候选在 AUTO 下通过 Netlib 114/114、不可行集 29/29、PDS 9/9、Railway 7/7；逐例耗时有快有慢，不以全集成功证明逐例速度无回退。未提交。

## English summary

The 2026-10-24–26 refactor separated algorithm identity, model representation, oracle implementation and numerical backends. AUTO owns cross-algorithm state; policy proposes decisions and routes adapts calls. Typed handoffs, exact failed requests, token-bound provenance, deterministic production snapshots and owner-specific resource/transaction boundaries replace implicit conventions. BDA and shared runtime now have separate source directories.

The final refactor passed 636 tests across Release/ASAN and optional-backend on/off configurations on macOS arm64, including relocated C90/C++ consumers, plus style/architecture gates. These are version-specific results, not evidence for later revisions or untested platforms.

Historical performance includes a 55-pair acceptance sample (candidate +2.7057% total wall time) and retained slow native-face/BDA samples, including +36.17% wall in one later comparison. The user waived a new full performance run; no new speedup or current-version regression claim follows. This cleanup consolidates live documentation and archives machine-specific history/builds outside the repository while preserving algorithms, reference tables and source data.
