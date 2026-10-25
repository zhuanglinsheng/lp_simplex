# lp_simplex

[简体中文](README.md) | English | [Documentation](docs/README_en.md)

`lp_simplex` is a continuous linear programming solver written in C. It provides a C API and MPS file input, and can be embedded as a standalone library in optimization and scientific computing applications.

Compared with common open-source LP solvers such as [HiGHS](https://highs.dev/), [SoPlex](https://soplex.zib.de/), [Clp](https://github.com/coin-or/Clp), and [GLPK](https://www.gnu.org/software/glpk/), this project emphasizes a basis-deficiency-allowing generalized simplex method (BDA). It represents feasible points using linearly independent active columns and allows the active basis dimension to grow or shrink during optimization. At degenerate points, it can directly exploit the problem’s lower-dimensional structure without completing a square basis.

Building on this, we use a hybrid solution flow: primal–dual hybrid gradient (PDHG) advances approximate solutions, least squares restores feasibility, and BDA then continues optimization and checks optimality as needed. This approach connects first-order iterations with dynamic deficient-basis optimization and reduces model expansion through native bound handling.

## Performance comparison

In the current benchmarks, lp_simplex and SCIP/SoPlex perform at a similar overall level. Although lp_simplex is slower than HiGHS and Clp, it demonstrates strong robustness on large-scale, highly degenerate problems. See the [benchmark report](docs/benchmarks/README_en.md) for full results.

- The small- and medium-scale problems in the [Netlib dataset](data/netlib/README_en.md) contain 32–232,966 variables and 24–105,127 constraints. See the [Netlib per-instance results](docs/benchmarks/netlib_en.md).
- The large-scale [PDS dataset](data/pds/README_en.md) contains approximately 106,000–505,000 variables and 34,000–156,000 constraints. See the [PDS per-instance results](docs/benchmarks/pds_en.md).
- The large-scale problems in the [Railway dataset](data/railway/README_en.md) use continuous LP relaxations, with approximately 47,000–1,093,000 variables and 507–4,872 constraints. See the [Railway per-instance results](docs/benchmarks/railway_en.md).

## Documentation

The [documentation home](docs/README_en.md) provides the following entry points. Chinese and English pages include language links below their titles:

- [Usage guide](docs/usage_en.md): follow a mathematical model through manual C modeling, MPS input, compilation, and standalone examples.
- [Algorithm overview](docs/algorithm/README_en.md): shared mathematical notation, algorithm derivations, and optimality conditions.
- [Parameters and results reference](docs/parameters_en.md): model fields, solver options, statuses and result contract, build/install reference.
- [Data and benchmarks](data/README_en.md): data preparation, performance tests, and reproduction.

The public interfaces are defined in the [model header](include/lp_simplex/model.h) and [solve header](include/lp_simplex/solve.h). The algorithm collection, result contract, benchmark summary, and data instructions are available in Chinese and English.

## Installation and Dependencies

### Precompiled releases

Download precompiled packages from [Releases](https://github.com/zhuanglinsheng/lp_simplex/releases), or build from source using the steps below. See [installation and dependencies](INSTALL.md) for platform and numerical-library requirements.

### Build from source

Building requires a C compiler, CMake 3.12+, BLAS, and LAPACK. SuiteSparse is optional and provides the KLU and SPQR sparse backends. Built-in fallbacks are used when it is unavailable. GMP is also optional: it verifies exact terminal basis-row certificates when floating verification is inconclusive. Disable it with `-DLP_SIMPLEX_USE_GMP=OFF`. Without GMP, evidence requirements remain unchanged; unverified cases continue numerical recovery or return failure.

Run from the project root, replacing `/path/to/install` with your installation directory:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/path/to/install
cmake --build build --target install
```

### Installed files

- Library directory (usually `lib/`): `liblp_simplex.so` (Linux), `liblp_simplex.dylib` (macOS), or the Windows import library `lp_simplex.lib`, controlled by `CMAKE_INSTALL_LIBDIR`.
- The Windows DLL is `bin/lp_simplex.dll`. Native MSVC builds use only BLAS/LAPACK by default; see the [installation guide](./INSTALL.md).
- `include/lp_simplex/`: public headers `lp_simplex.h`, `model.h`, `mps.h`, `solve.h`, `status.h`, and `version.h`.
- `cmake/lp_simplex/` under the library directory: package configuration and the `lp_simplex::lp_simplex` imported target.
- `share/licenses/lp_simplex/`: LGPL and referenced GPL license texts.

An external C or C++ project can use the installed target:

```cmake
cmake_minimum_required(VERSION 3.12)
project(my_lp_client LANGUAGES C)
find_package(lp_simplex CONFIG REQUIRED)
add_executable(my_lp_client main.c)
target_link_libraries(my_lp_client PRIVATE lp_simplex::lp_simplex)
```

Configure the client with `-DCMAKE_PREFIX_PATH=/path/to/install`. With relative GNUInstallDirs destinations, the installation prefix may be moved as a unit; configure clients against the new prefix after moving it. For a custom library directory, `lp_simplex_DIR=<prefix>/<library-dir>/cmake/lp_simplex` selects the package directly. The deployment environment provides BLAS/LAPACK and other enabled dynamic dependencies. `lp_simplex_WITH_KLU`, `lp_simplex_WITH_SPQR`, and `lp_simplex_WITH_GMP` report resolved build capabilities.

Use Release for delivery. Non-MSVC Debug builds enable AddressSanitizer; external Debug programs must also link with `-fsanitize=address`. The package regression sets this for the matching test configuration.

`ctest --test-dir build --output-on-failure` includes independent C/C++ consumers after installation and relocation. Python 3 enables data-preparation regression; NumPy and SciPy additionally enable benchmark-tool regression. CI installs these test dependencies. Host package execution tests are omitted for cross-compilation.

## Contributing and License

Report issues or discuss improvements through [GitHub Issues](https://github.com/zhuanglinsheng/lp_simplex/issues), or submit a pull request. For solver or performance issues, include the model, build environment, runtime parameters, and reproduction steps.

This project is licensed under the [GNU Lesser General Public License v3.0 (LGPL-3.0)](LICENSE). The GNU GPL v3.0 terms incorporated by that license are provided in [LICENSE.GPL-3.0](LICENSE.GPL-3.0).

## Delivery and compatibility boundary

The public interface is `include/lp_simplex/`; `src/` is private. Rebuild consumers against matching headers; cross-version binary compatibility is not promised. Release/ASAN and optional-backend on/off configurations, including relocated consumers, were validated on macOS arm64. Other platforms require their actual CI results. Work limits do not impose CPU timeouts or total RSS caps. See [configuration](./docs/runtime-configuration_en.md), [result interpretation](./docs/parameters_en.md) and [validation history](./docs/development/README.md).
