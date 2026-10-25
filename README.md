# lp_simplex

简体中文 | [English](README_en.md) | [文档](docs/README.md)

`lp_simplex` 是一个用 C 编写的连续线性规划求解器，提供 C API 和 MPS 文件输入，可作为独立库嵌入优化与科学计算程序。

与 [HiGHS](https://highs.dev/)、[SoPlex](https://soplex.zib.de/)、[Clp](https://github.com/coin-or/Clp) 和 [GLPK](https://www.gnu.org/software/glpk/) 等常用的开源 LP 求解器相比，本项目的特色是采用允许基亏缺的广义单纯型法（BDA）。它以线性无关的活动列表示可行点，允许活动基维数随求解过程增减。退化时无需补齐方基，可以直接利用问题的低维结构。

我们在此基础上使用混合求解流程：用原始—对偶混合梯度法（PDHG）推进近似解，通过最小二乘恢复可行性，再由 BDA 按需继续优化并检验最优性。这条路线将一阶迭代与动态亏基优化连接起来，并通过原生上下界处理减少模型展开。

## 性能对比

在当前基准测试中，lp_simplex 与 SCIP/SoPlex 的整体表现处于同一水平，速度落后于 HiGHS 和 Clp，但在大规模、高度退化的问题上具有较强的稳健性。完整结果见[综合性能对比报告](docs/benchmarks/README.md)。

- [Netlib 数据集](data/netlib/README.md)包含 32—232,966 个变量、24—105,127 条约束。逐例测试结果见 [Netlib 问题对比](docs/benchmarks/netlib.md)。
- 大规模问题 [PDS 数据集](data/pds/README.md)包含约 10.6万—50.5万变量、3.4万—15.6万约束。逐例测试结果见 [PDS 问题对比](docs/benchmarks/pds.md)。
- 大规模 [Railway 数据集](data/railway/README.md)问题采用连续 LP 松弛，包含约 4.7万—109.3万变量、507—4,872 条约束。逐例测试结果见 [Railway 问题对比](docs/benchmarks/railway.md)。

## 文档

[文档首页](docs/README.md)提供以下阅读入口，中英文页面可通过标题下的链接切换：

- [使用指南](docs/usage.md)：从数学模型出发，分别手动建模或读取 MPS，再编译运行独立示例。
- [算法总览](docs/algorithm/README.md)：统一数学记号、主要算法的推导与最优性条件。
- [参数与结果说明](docs/parameters.md)：模型参数、求解选项、状态与结果契约、构建安装参考。
- [数据与基准测试](data/README.md)：数据准备、性能测试与结果复现。

公开接口见[模型头文件](include/lp_simplex/model.h)和[求解头文件](include/lp_simplex/solve.h)。使用指南、参数与结果说明、算法文集、性能报告和数据说明均提供中英文版本。

## 安装与依赖

### 预编译版本

可从 [Releases](https://github.com/zhuanglinsheng/lp_simplex/releases) 下载预编译版本，也可以按以下步骤从源码安装。平台与数值库要求见[安装与依赖说明](INSTALL.md)。

### 从源码安装

源码构建需要 C 编译器、CMake 3.12+、BLAS 和 LAPACK。SuiteSparse 为可选依赖，提供 KLU 和 SPQR 稀疏运算后端；未安装时使用内置后备路径。GMP 也是可选依赖：仅在浮点不可行证据无法确认时，用于精确基行验证；可用 `-DLP_SIMPLEX_USE_GMP=OFF` 关闭。缺少 GMP 时不会降低证据要求，无法确认的情况继续数值恢复或返回失败。

在项目根目录执行，将 `/path/to/install` 替换为安装目录：

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/path/to/install
cmake --build build --target install
```

### 安装内容

- 库目录（通常 `lib/`）：共享库 `liblp_simplex.so`（Linux）或 `liblp_simplex.dylib`（macOS）。目录由 `CMAKE_INSTALL_LIBDIR` 控制。
- `include/lp_simplex/`：公开头文件 `lp_simplex.h`、`model.h`、`mps.h`、`solve.h`、`status.h` 和 `version.h`。
- 库目录下的 `cmake/lp_simplex/`：CMake 包配置与 `lp_simplex::lp_simplex` 导入目标。
- `share/licenses/lp_simplex/`：LGPL 与所引用的 GPL 条款文本。

外部工程可使用已安装的包，C 与 C++ 共用同一个目标：

```cmake
cmake_minimum_required(VERSION 3.12)
project(my_lp_client LANGUAGES C)
find_package(lp_simplex CONFIG REQUIRED)
add_executable(my_lp_client main.c)
target_link_libraries(my_lp_client PRIVATE lp_simplex::lp_simplex)
```

配置外部工程时传入 `-DCMAKE_PREFIX_PATH=/path/to/install`。使用相对的 GNUInstallDirs 安装目录时，安装前缀可整体移动，移动后以新前缀配置客户端。自定义库目录时也可设置 `lp_simplex_DIR=<前缀>/<库目录>/cmake/lp_simplex`。BLAS/LAPACK 及构建时启用的其他动态库由部署环境提供。包配置中的 `lp_simplex_WITH_KLU`、`lp_simplex_WITH_SPQR`、`lp_simplex_WITH_GMP` 反映实际探测并启用的能力。

Release 用于交付。非 MSVC 的 Debug 构建启用 AddressSanitizer，外部 Debug 程序也需要以 `-fsanitize=address` 链接；安装回归会按对应配置设置 sanitizer。

`ctest --test-dir build --output-on-failure` 包含安装移动后的独立 C/C++ 调用检查。检测到 Python 3 时加入数据准备回归，另检测到 NumPy 时加入基准工具回归；CI 安装这些测试依赖。交叉编译时不运行宿主机安装消费者测试。

## 参与与许可证

欢迎通过 [GitHub Issues](https://github.com/zhuanglinsheng/lp_simplex/issues) 报告问题或讨论改进，也欢迎提交 Pull Request。报告求解或性能问题时，请附模型、构建环境、运行参数及复现步骤。

本项目采用 [GNU Lesser General Public License v3.0（LGPL-3.0）](LICENSE)。其引用的 GNU GPL v3.0 条款文本保存在 [LICENSE.GPL-3.0](LICENSE.GPL-3.0)。

## 交付与兼容性边界

公开接口位于 `include/lp_simplex/`，`src/` 属于私有实现。调用方应使用匹配的头文件重新构建，不承诺跨版本二进制兼容。macOS arm64 已验证 Release/ASAN、可选后端开/关及安装包移位；其他平台以实际 CI 结果为准。工作预算不等于 CPU 超时或总体 RSS 限额。配置和结果含义见[配置边界](./docs/runtime-configuration.md)与[结果参考](./docs/parameters.md)，测量范围见[验证历史](./docs/development/README.md)。
