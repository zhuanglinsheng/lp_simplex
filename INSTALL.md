# Installing lp_simplex 1.0

## 简体中文

正式包提供共享库、C90/C++ 可用的公开头文件和 CMake 导入目标。
Linux x64、macOS x64/arm64 和 Windows x64 分别构建；请按本机架构选择。
Linux/macOS 包启用 KLU、SPQR 和 GMP。Windows x64 使用已验证的最小依赖版本，
默认关闭三个可选后端，使用 BLAS/LAPACK 和内置算法。最终发布以各平台验收结果为准。

Linux/macOS 的数值动态库由运行环境提供，压缩包不包含第三方数值库：

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
数据准备与部分回归需要 Python 3；NumPy 和 SciPy 可启用基准工具回归。
`tools/build.sh` 默认要求三个可选后端全部可用，不会静默降级。
也可以直接使用 CMake 关闭可选后端，以构建只依赖 BLAS/LAPACK 的版本：

```sh
cmake -S . -B build-minimal -DCMAKE_BUILD_TYPE=Release \
  -DLP_SIMPLEX_USE_KLU=OFF -DLP_SIMPLEX_USE_SPQR=OFF -DLP_SIMPLEX_USE_GMP=OFF
cmake --build build-minimal
ctest --test-dir build-minimal --output-on-failure
```

Windows 源码构建使用 Visual Studio 2022 x64 开发者环境、PowerShell 7.3+、
CMake 和 Ninja。依赖版本由 `vcpkg.json` 固定，默认仅安装 BLAS/LAPACK：

```powershell
./tools/install_dependencies.ps1
python data/netlib/prepare.py --all
./tools/build.ps1
./tools/build.ps1 -BuildDir build-windows-debug -Configuration Debug
```

Windows 打包脚本收集实际链接的数值 DLL，并验证解压移位后的 C/C++ 调用。
运行 Windows 二进制包需要可用的 Microsoft Visual C++ x64 运行库；
压缩包不捆绑 Windows 系统 DLL 或 MSVC 运行库。
可选的 `-OptionalBackends` 必须同时用于依赖安装和构建；这是单独验证的扩展配置，
不属于当前 Windows 默认验收范围。关闭可选后端时，公开接口保留，能力标志为 OFF，
求解使用相应回退实现；不承诺不同后端具有相同性能或数值轨迹。

安装脚本使用系统包管理器，需要联网；Linux 系统安装可能要求 sudo。
CMake 配置本身不安装软件、不下载依赖源码。

## English

Choose the native Linux x64, macOS x64/arm64 or Windows x64 archive. Release binaries
include the lp_simplex shared library, public C/C++ headers and CMake package.
Linux/macOS enable KLU, SPQR and GMP; third-party numerical shared libraries
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
NumPy and SciPy enable benchmark-tool regression tests. These are not C API runtime
requirements. For full benchmark comparisons see [data documentation](data/README.md).

Windows x64 validation uses native MSVC 2022, PowerShell 7.3+, CMake and Ninja.
Run the PowerShell commands above to install pinned BLAS/LAPACK dependencies
and test Release/Debug. Optional KLU/SPQR/GMP backends are OFF by default;
public interfaces remain available through fallback implementations and report
those capabilities as OFF. Backend performance and numerical paths can differ.
Use `-OptionalBackends` for both PowerShell scripts only when validating that
additional configuration separately. Windows archives bundle the numerical DLLs
actually required by the library and are tested after extraction and relocation.
The host must provide the Microsoft Visual C++ x64 runtime; Windows system DLLs
and the MSVC runtime are not bundled.
The Windows baseline has passed native Release/Debug and relocated archive validation.
