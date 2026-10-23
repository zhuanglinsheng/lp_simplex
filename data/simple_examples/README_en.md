# Small Examples

[简体中文](./README.md) | English | [Usage guide](../../docs/usage_en.md) | [Data overview](../README_en.md)

Each model has three matching files: `xxx_1.c` constructs it manually, `xxx_2.c`
reads it from a file, and `xxx.mps` stores the model. All 18 C files are standalone
programs using the public API, without a benchmark driver or data-preparation step.

Build and run both versions from the repository root:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DLP_SIMPLEX_BUILD_EXAMPLES=ON
cmake --build build -j 4
./build/data/simple_examples/example_production_planning_manual
./build/data/simple_examples/example_production_planning_mps data/simple_examples/production_planning.mps
ctest --test-dir build -R '^example_' --output-on-failure
```

| Model | Manual construction | MPS input | Data |
| --- | --- | --- | --- |
| `production_planning` | [production_planning_1.c](production_planning_1.c) | [production_planning_2.c](production_planning_2.c) | [MPS](production_planning.mps) |
| `mixed_inequalities` | [mixed_inequalities_1.c](mixed_inequalities_1.c) | [mixed_inequalities_2.c](mixed_inequalities_2.c) | [MPS](mixed_inequalities.mps) |
| `equality_system` | [equality_system_1.c](equality_system_1.c) | [equality_system_2.c](equality_system_2.c) | [MPS](equality_system.mps) |
| `free_variables` | [free_variables_1.c](free_variables_1.c) | [free_variables_2.c](free_variables_2.c) | [MPS](free_variables.mps) |
| `free_and_lower_bounds` | [free_and_lower_bounds_1.c](free_and_lower_bounds_1.c) | [free_and_lower_bounds_2.c](free_and_lower_bounds_2.c) | [MPS](free_and_lower_bounds.mps) |
| `variable_bounds` | [variable_bounds_1.c](variable_bounds_1.c) | [variable_bounds_2.c](variable_bounds_2.c) | [MPS](variable_bounds.mps) |
| `ranged_rows` | [ranged_rows_1.c](ranged_rows_1.c) | [ranged_rows_2.c](ranged_rows_2.c) | [MPS](ranged_rows.mps) |
| `degenerate_pivot` | [degenerate_pivot_1.c](degenerate_pivot_1.c) | [degenerate_pivot_2.c](degenerate_pivot_2.c) | [MPS](degenerate_pivot.mps) |
| `energy_system` | [energy_system_1.c](energy_system_1.c) | [energy_system_2.c](energy_system_2.c) | [MPS](energy_system.mps) |

Set `LP_SIMPLEX_BUILD_EXAMPLES=OFF` to skip the examples. With BUILD_TESTING enabled,
18 tests are registered. These models are usage examples and functional tests,
excluded from Netlib, PDS, and Railway performance benchmarks.
