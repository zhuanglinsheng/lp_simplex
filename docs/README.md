# 项目文档

简体中文 | [English](./README_en.md) | [项目主页](../README.md)

## 使用与查阅

| 文档 | 解决的问题 |
| --- | --- |
| [使用指南](./usage.md) | 建模、读取 MPS、编译并完成第一次求解 |
| [参数与结果](./parameters.md) | Options、变量界、状态、残差空间、进度和构建选项 |
| [配置边界](./runtime-configuration.md) | 正式默认值与显式环境实验 |
| [功能索引](./features.md) | 能力、适用条件、限制及实现入口 |
| [算法文集](./algorithm/README.md) | 数学记号、推导、假设及最优性条件 |

## 架构与数值

| 文档 | 唯一维护的主题 |
| --- | --- |
| [AUTO](./auto-policy.md) | 目录职责、调度、预算、候选和恢复链 |
| [BDA 阶段](./pan-services.md) | 子过程角色、身份、授权及重复准入 |
| [列模型交接](./column-handoff.md) | cover 的列投影、点/乘子/基权限 |
| [数值政策](./numerical-policy.md) | 内部阈值的单位、尺度和资源失败边界 |

## 数据与测试

[数据说明](../data/README.md)维护准备、命令行、审核和计时方法；[基准汇总](./benchmarks/README.md)及逐数据集报告保存版本对应的测量。[开发记录](./development/README.md)仅保留设计决策和验证历史。

## 维护约定

自有 C/C++ 使用 Tab（8 列）、Linux 大括号、80 列排版目标、最多一个空行，并以 Doxygen `@file`/`@brief` 说明文件职责。CMake/Python 使用 4 空格，YAML/JSON 使用 2 空格；统一 LF 和末尾换行。格式版本固定在仓库根目录。

```sh
python3 -m pip install clang-format==22.1.8 cmakelang==0.6.13
python3 tools/check_style.py
python3 tools/check_architecture.py
python3 tools/check_docs.py
ctest --test-dir build --output-on-failure
```

格式检查覆盖 C/C++、CMake 和脚本基础空白/语法；不声称完整格式化 Python/YAML。架构检查覆盖私有包含边界；文档检查覆盖本地文件和 Markdown 锚点；它们均不能替代数值回归。生成的基准报告通过 `data/benchmark.py --report-only` 更新，不手工修改测量数据。
