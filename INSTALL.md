# Installing lp_simplex 1.0

## 简体中文

正式包提供共享库、C90/C++ 可用的公开头文件和 CMake 导入目标。
Linux x64、macOS x64 和 macOS arm64 分别构建；请按本机架构选择。
1.0 正式二进制包启用 KLU、SPQR 和 GMP。

数值动态库由运行环境提供，压缩包不包含第三方数值库：

- Ubuntu/Debian：安装 `libopenblas-dev liblapack-dev libsuitesparse-dev libgmp-dev`。
- macOS：运行 `brew install suite-sparse gmp`，BLAS/LAPACK 使用系统实现。

不同 Linux 发行版、SuiteSparse 主版本和 BLAS 实现不保证二进制兼容。
Linux 正式包在 Ubuntu 24.04 上构建；环境不兼容时请从冻结源码构建。
包内 `BUILD-INFO.json` 保存构建提交、工具链、实际启用能力和构建时依赖路径；
这些路径用于追溯，不是要求用户复制的安装前缀。

解压后，可将包目录整体移动。外部工程使用：

```cmake
find_package(lp_simplex 1.0.0 EXACT CONFIG REQUIRED)
target_link_libraries(my_program PRIVATE lp_simplex::lp_simplex)
```

配置时传入 `-DCMAKE_PREFIX_PATH=/path/to/extracted/package`。
头文件版本为 `LP_SIMPLEX_VERSION`，加载库版本可用 `lp_simplex_version()` 查询。
升级时使用匹配的头文件重新编译，不承诺跨版本二进制兼容。

从源码构建时先安装依赖，再运行构建与测试：

```sh
sh tools/install_dependencies.sh
python3 data/netlib/prepare.py --all
sh tools/build.sh build-release
cmake --install build-release --prefix /path/to/install
```

需要 C/C++ 编译器、CMake 3.12+ 和 Make 或 Ninja。
数据准备与部分回归需要 Python 3；NumPy 可启用基准工具回归。
`tools/build.sh` 默认要求三个可选后端全部可用，不会静默降级。
也可以直接使用 CMake 关闭可选后端，以构建只依赖 BLAS/LAPACK 的版本：

```sh
cmake -S . -B build-minimal -DCMAKE_BUILD_TYPE=Release \
  -DLP_SIMPLEX_USE_KLU=OFF -DLP_SIMPLEX_USE_SPQR=OFF -DLP_SIMPLEX_USE_GMP=OFF
cmake --build build-minimal
ctest --test-dir build-minimal --output-on-failure
```

安装脚本使用系统包管理器，需要联网；Linux 系统安装可能要求 sudo。
CMake 配置本身不安装软件、不下载依赖源码。

## English

Choose the native Linux x64, macOS x64 or macOS arm64 archive. Release binaries
include the lp_simplex shared library, public C/C++ headers and CMake package.
KLU, SPQR and GMP support is enabled; third-party numerical shared libraries
must be provided by the host. On Ubuntu/Debian install `libopenblas-dev`,
`liblapack-dev`, `libsuitesparse-dev` and `libgmp-dev`; on macOS run
`brew install suite-sparse gmp` (system BLAS/LAPACK is used).

Linux archives are built on Ubuntu 24.04. Binary compatibility across Linux
distributions, SuiteSparse major versions or BLAS providers is not guaranteed.
Build from the frozen source when the host libraries are incompatible.
`BUILD-INFO.json` records the source revision, toolchain and resolved backends.

Extract and move the package anywhere, then configure a consumer with
`-DCMAKE_PREFIX_PATH=/path/to/package` and link `lp_simplex::lp_simplex`.
The CMake example and source-build commands above apply to both languages.
Rebuild clients against matching headers when upgrading; cross-version ABI
compatibility is not promised. `LP_SIMPLEX_VERSION` identifies the headers and
`lp_simplex_version()` identifies the loaded library.

The dependency installer uses apt or Homebrew and may require network access
and elevated package-manager permissions. CMake itself never installs system
packages. `tools/build.sh` requires all three optional backends by default.
Source builds require CMake 3.12+, a C compiler and Make or Ninja; C++ is used for
consumer validation. Python enables data preparation and additional tests;
NumPy enables benchmark-tool regression tests. These are not C API runtime
requirements. For full benchmark comparisons see [data documentation](data/README.md).
