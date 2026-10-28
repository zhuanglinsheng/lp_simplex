# PDS AUTO 复查（2026-10-26）

2026-10-27 当前状态：本文以下分支／候选哈希及验收描述属于 10 月 26 日的历史快照；相关候选现已汇集为 main 的未提交工作区改动，不再有独立 PDS 分支。最新离线归因及暂停状态见[更新后的深度分析](pds-performance-analysis-2026-10-26.md#7-2026-10-27隐含上界中的搜索容差污染)。中间原始文件不作为仓库交付。

旧轮次曾保留原始本机证据；按用户后续要求，中间测试文件清理，历史性能目录只保留基准测试 Markdown。下文的原始证据路径与保留描述属于当时快照，不表示当前仍保留这些文件。

分层候选保存在实验分支 `codex/pds-auto-feasibility`，数值实现提交为 `8730c9c`。功能回归已完成，计时验收尚待空闲机器条件；本页不宣称它已经成为 main 的默认策略。

本次复查使用同一台 macOS arm64 机器、Release 构建、KLU/SPQR/GMP 开启及单线程配置。基线为 `559100b` 的数值代码；独立构建分别冻结可执行文件和动态库，配对测试交替运行，保存输入、源码和二进制 SHA-256。原模型审核、参考目标、容差及工作预算均未放宽。

[复查前的完整报告](../benchmarks/history/2026-10-26/README.md) · 被否决的完整修复候选 · [当前报告](../benchmarks/README.md)

## 发现的调度盲点

原完整测试中，PDS 的四个失败都是 300 秒超时。输入生成和校验正常，AUTO 也确实在 Dual 让出后进入原生界 BDA；并非漏测案例或漏掉整条恢复路线。

在 pds-70 的基线诊断中，Dual 在 16,384 次工作后让出。随后 PDHG 的固定工作模型可行性残差在约 65k—82k 步保持约 `5.14e-5`，明显高于 `1e-7` 的要求，但加权综合指标继续下降，权重变化又重置了综合指标的历史。因此旧条件没有给出恢复机会，直到约 105k 步才获得更优可行点。

这揭示了两个不同的问题：综合指标改善不能推出可行性分量改善；另一方面，可行性停滞也不能直接授权昂贵的完整修复搜索，因为对偶分量仍可能取得有用进展。观察信号与干预成本必须分别控制。

## 候选及反例

第一版加入独立可行性观察后，能更早恢复可行点，但继续等待 PDHG 自身认证，pds-70 诊断仍在 330 秒超时。第二版将经审核的可行改进交回原生算法优化和认证。两版都没有修改最终审核标准。

第二版在 PDS 有效，却给 Netlib 的 qap15 引入回退。诊断显示，在联合指标仍有进展时，一次新触发的恢复连续进行了 16 次原生校正，全部失败并回滚。这个候选没有作为默认策略直接采用。

| 第二版配对案例 | 基线 | 候选 | 口径 |
| --- | --- | --- | --- |
| pds-20 | 70.812 / 75.532 秒，均通过 | 59.642 / 55.983 秒，均通过 | 各两次；300 秒 / 300000 工作 |
| pds-70 | 两次均在 300 秒超时 | 210.628 / 288.158 秒，均通过 | 各两次；300000 工作 |
| qap15 | 54.993 / 65.002 秒，均通过 | 两次均在 90 秒超时 | 各两次；300000 工作 |

诊断额外启用输出专用 PROFILE，使用 120 或 330 秒观察阶段；这些诊断不参与正式的 PASS 或速度统计。不同案例和不同轮次未合并为同一加速比。尤其不能用 PDS 的改善抵消 qap15 的新增失败。

一个后续试验把整个恢复授权限制为一个工作单位，虽然 qap15 两对都通过，却截断了 PDS 所需的锥修复内部 NNLS。PDS 新增超时后已中止该候选；未完成记录单独保留，不作为全集统计。日志中的锥修复 `attempts=0` 表示第一轮成功，不代表零工作。

## 分层修订的通用规则

- 分别记录加权综合指标和固定工作模型的未加权可行性残差；权重变化仅重置前者。沿用原来的 8192 步停滞窗口、1% 累计进展标准及尝试前残差严格减半的条件。
- 仅判定可行性停滞时，仿射修复和锥修复各最多试一次校正；保留锥修复内部 NNLS 原有的有限共享预算，失败后完整回滚并继续 PDHG。联合指标停滞时保留原有的有限修复搜索。
- 恢复当前点或平均点后，只有经审核的更优可行点可以交接；最终最优性仍由独立的原生 KKT 和原模型审核确定。所有子过程共享原累计预算。

规则不读取数据集名称、尺寸区间、参考目标或墙钟期限。限制针对各修复方法的外层校正次数，避免把“可行性机会”扩大成反复换面搜索；内部 NNLS 的合法多步求解仍有预算。没有为某个案例调整残差门槛。

## 冻结验证

| 数据集 | 原完整基线 | 第二版候选 | 分层修订版 |
| --- | ---: | ---: | ---: |
| PDS | 5/9 | 8/9 | 待空闲验证 |
| Netlib | 112/114 | 113/114 | 待空闲验证 |
| Netlib infeasible | 29/29 | 29/29 | 待空闲验证 |
| Railway | 7/7 | 7/7 | 待空闲验证 |

第二版的 PDS 数字合并了两个开发案例的各两次配对与其余七例的各一次冻结验证；它不是统一重复次数的九例运行。修订版的完整验收要求按原预算单次重跑全部四个数据集，并对 qap15 另作两次配对；当前尚未取得满足条件的完整结果。每个实验保留失败记录和原始范围，不把较快的重复挑入最终报告。

第二版和分层修订版均通过 179 项 CTest；格式与架构检查通过。新增两行三变量回归：当前自由面无法修复另一行残差，必须允许锥修复内部 NNLS 使用多个工作单位；另一个带上界的案例验证单次锥校正失败回滚、完整搜索仍可成功。这些功能检查不替代上述性能验证。原始本机记录位于 `tmp/auto-pds-review`，归档保留机器结果、构建来源、候选补丁及诊断日志。

## 计时环境与验收状态

后续控制组中，同一个未修改基线的 qap15 从此前约 55—65 秒变为两次均在 90 秒超时，分层候选也两次超时。该时段八逻辑核机器的一分钟系统负载约 114；未发现本仓库遗留求解器进程，测试脚本会终止整个 worker 进程组。此轮高负载记录保留原样，不能用它区分策略影响与机器争用。

续测要求运行前负载回落，并连续记录测量期的系统负载。先检验冻结基线与候选的两对 qap15 和 pds-70，再以原有参数收集完整四个数据集。不能用延长正式限时、删去失败或选择较快重复来取得通过。

## 分层候选的实测结果与处置

机器负载回落后完成了第二轮控制验证：qap15 基线为 58.090 / 62.069 秒，候选为 67.237 / 69.132 秒，两对均通过。中位数分别为 60.080 / 68.185 秒，候选约慢 13.5%；各次累计工作分别为 107211 / 107214。工作计数接近并不能证明计算成本接近，两对小样本也不足以证明性能等价。pds-70 在原 300 秒 / 300000 工作预算内以 190.631 秒通过原模型审核。

随后 PDS 全集测量期间，系统负载再次超过 110；部分先前成功案例出现超时，已主动中断。该批次标记为 `completed=false`，原始失败和负载轨迹均保留，不能拼入完整通过率或从中挑选较快结果。分层候选记录包含控制配对、单例验证及中断全集的独立范围。

本轮复查确认了观察条件的盲点，也确认了无条件增加完整修复会伤害其他问题。分层候选仍有计时代价且未完成全集验收，因此保留在实验分支，不合入 main，不替换当前完整基准报告。默认 Release 构建保持原基线；本轮测量及空闲监测均已停止。后续采用候选前，需要在稳定负载下完成四个数据集的验证，并明确评估 qap15 的额外成本。

## English summary

AUTO did invoke native BDA after Dual yielded; the four original PDS failures were wall-time limits. A weighted joint merit can improve while feasibility stalls, and reweighting resets its history. An independent feasibility observer exposes opportunities but does not justify unrestricted repair: the second candidate improved PDS from 5/9 to 8/9, yet introduced a qap15 timeout. Two alternating pairs confirmed that regression, so the candidate was not adopted directly.

An intermediate one-work lease wrongly suppressed cone NNLS and was stopped after PDS regressions; its partial records remain separate. The layered revision permits one affine and one cone correction for a feasibility-only opportunity, retaining the bounded inner NNLS budget and rolling back on failure. Joint stagnation retains the existing bounded search. Audited feasible improvements can hand off to native optimization and certification within the shared budget. No dataset names, dimensions, reference answers or wall-time thresholds enter the rule. The original tolerance and work limits remain unchanged. All 179 functional tests pass. Timing acceptance is pending a quiet-host rerun: both variants timed out in a later control while host load reached114 on eight logical CPUs. That censored comparison is retained and is not evidence of policy-specific regression. The final table will retain the actual scope and failures.

The later qap15 pairs all passed: baseline 58.090/62.069 seconds versus candidate 67.237/69.132 seconds, a 13.5% median wall-time overhead in two pairs. The pds-70 pilot passed in 190.631 seconds under the original 300-second budget. Full PDS collection was interrupted when host load again exceeded 110; all partial failures and load samples are retained separately. The review is documented, but general performance acceptance remains incomplete. The candidate stays on the experimental branch; main and the current full benchmark reports retain the original policy. No benchmark or idle monitor remains running.
