# lp_simplex

[简体中文](README.md) | English | [Documentation](docs/README_en.md)

`lp_simplex` is a continuous linear programming solver written in C. It provides a C API and MPS file input, and can be embedded as a standalone library in optimization and scientific computing applications.

Compared with common open-source LP solvers such as [HiGHS](https://highs.dev/), [SoPlex](https://soplex.zib.de/), [Clp](https://github.com/coin-or/Clp), and [GLPK](https://www.gnu.org/software/glpk/), this project emphasizes a basis-deficiency-allowing generalized simplex method (BDA). It represents feasible points using linearly independent active columns and allows the active basis dimension to grow or shrink during optimization. At degenerate points, it can directly exploit the problem’s lower-dimensional structure without completing a square basis.

Building on this, we explore a novel hybrid solution approach: primal–dual hybrid gradient (PDHG) advances approximate solutions, least squares restores feasibility, and BDA then continues optimization and checks optimality as needed. This approach connects first-order iterations with dynamic deficient-basis optimization and reduces model expansion through native bound handling.

## Performance comparison

In the current benchmarks, lp_simplex and SCIP/SoPlex perform at a similar overall level. Although lp_simplex is slower than HiGHS and Clp, it demonstrates strong robustness on large-scale, highly degenerate problems. See the [benchmark report](docs/benchmarks/README_en.md) for full results.

- The [Netlib dataset](data/netlib/README_en.md) contains 32–232,966 variables and 24–105,127 constraints. See the [Netlib per-instance results](docs/benchmarks/netlib_en.md).
- The large-scale [PDS dataset](data/pds/README_en.md) contains approximately 106,000–505,000 variables and 34,000–156,000 constraints. See the [PDS per-instance results](docs/benchmarks/pds_en.md).
- The large-scale [Railway dataset](data/railway/README_en.md) uses continuous LP relaxations, with approximately 47,000–1,093,000 variables and 507–4,872 constraints. See the [Railway per-instance results](docs/benchmarks/railway_en.md).

## Documentation

The [documentation home](docs/README_en.md) provides the following entry points. Chinese and English pages include language links below their titles:

- [Usage guide](docs/usage_en.md): follow a mathematical model through manual C modeling, MPS input, compilation, and standalone examples.
- [Algorithm overview](docs/algorithm/README_en.md): shared mathematical notation, algorithm derivations, and optimality conditions.
- [Parameters and results reference](docs/parameters_en.md): model fields, solver options, statuses and result contract, build/install reference.
- [Data and benchmarks](data/README_en.md): data preparation, performance tests, and reproduction.

The public interfaces are defined in the [model header](include/lp_simplex/model.h) and [solve header](include/lp_simplex/solve.h). The algorithm collection, result contract, benchmark summary, and data instructions are available in Chinese and English.

## Installation and Dependencies

### Precompiled releases

Download precompiled packages from [Releases](https://github.com/zhuanglinsheng/lp_simplex/releases), or build from source using the steps below.

### Build from source

Building requires a C compiler, CMake 3.12+, BLAS, and LAPACK. SuiteSparse is optional and provides the KLU and SPQR sparse backends. Built-in fallbacks are used when it is unavailable.

Run from the project root, replacing `/path/to/install` with your installation directory:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/path/to/install
cmake --build build --target install
```

### Installed files

- `lib/`: the shared library `liblp_simplex.so` (Linux) or `liblp_simplex.dylib` (macOS).
- `include/lp_simplex/`: public headers `lp_simplex.h`, `model.h`, `mps.h`, `solve.h`, and `status.h`.

## Contributing and License

Report issues or discuss improvements through [GitHub Issues](https://github.com/zhuanglinsheng/lp_simplex/issues), or submit a pull request. For solver or performance issues, include the model, build environment, runtime parameters, and reproduction steps.

This project is licensed under the [GNU Lesser General Public License v3.0 (LGPL-3.0)](LICENSE). The GNU GPL v3.0 terms incorporated by that license are provided in [LICENSE.GPL-3.0](LICENSE.GPL-3.0).
