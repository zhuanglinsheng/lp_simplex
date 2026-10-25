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

## English summary

The 2026-10-24–26 refactor separated algorithm identity, model representation, oracle implementation and numerical backends. AUTO owns cross-algorithm state; policy proposes decisions and routes adapts calls. Typed handoffs, exact failed requests, token-bound provenance, deterministic production snapshots and owner-specific resource/transaction boundaries replace implicit conventions. BDA and shared runtime now have separate source directories.

The final refactor passed 636 tests across Release/ASAN and optional-backend on/off configurations on macOS arm64, including relocated C90/C++ consumers, plus style/architecture gates. These are version-specific results, not evidence for later revisions or untested platforms.

Historical performance includes a 55-pair acceptance sample (candidate +2.7057% total wall time) and retained slow native-face/BDA samples, including +36.17% wall in one later comparison. The user waived a new full performance run; no new speedup or current-version regression claim follows. This cleanup consolidates live documentation and archives machine-specific history/builds outside the repository while preserving algorithms, reference tables and source data.
