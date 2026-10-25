# Project documentation

[简体中文](./README.md) | English | [Project home](../README_en.md)

## Usage and reference

| Document | Question answered |
| --- | --- |
| [Usage](./usage_en.md) | Model, read MPS, compile and solve a first problem |
| [Parameters/results](./parameters_en.md) | Options, bounds, statuses, residual spaces, progress and build options |
| [Configuration](./runtime-configuration_en.md) | Production defaults versus explicit environment experiments |
| [Capability index](./features_en.md) | Features, applicability, limits and implementation entries |
| [Algorithms](./algorithm/README_en.md) | Mathematical notation, derivations, assumptions and optimality |

## Architecture and numerics

| Document | Owned topic |
| --- | --- |
| [AUTO](./auto-policy_en.md) | Source ownership, scheduling, grants, candidates and recovery chain |
| [BDA phases](./pan-services_en.md) | Subprocess roles, identity, grants and repeat admission |
| [Column handoff](./column-handoff_en.md) | Cover projections and point/multiplier/basis permissions |
| [Numerical policy](./numerical-policy_en.md) | Internal thresholds, units, scale and resource failures |

## Data and tests

[Data instructions](../data/README_en.md) own preparation, CLI, audits and timing. [Benchmark summaries](./benchmarks/README_en.md) and dataset reports retain version-specific measurements. [Development history](./development/README.md) records design decisions and validation history only.

## Maintenance conventions

First-party C/C++ uses tabs (8 columns), Linux braces, an 80-column target, at most one blank line and Doxygen `@file`/`@brief` descriptions. CMake/Python use four spaces; YAML/JSON use two. Files use LF and a final newline. Formatter versions are pinned at the repository root.

```sh
python3 -m pip install clang-format==22.1.8 cmakelang==0.6.13
python3 tools/check_style.py
python3 tools/check_architecture.py
python3 tools/check_docs.py
ctest --test-dir build --output-on-failure
```

Style checks cover C/C++, CMake and basic script whitespace/syntax, not complete Python/YAML formatting. Architecture checks cover private include boundaries; documentation checks cover local files and Markdown anchors. None replaces numerical regression. Update generated benchmark reports with `data/benchmark.py --report-only`, without manually editing measurements.
