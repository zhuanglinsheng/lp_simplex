# 小型示例

简体中文 | [English](./README_en.md) | [使用指南](../../docs/usage.md) | [数据总览](../README.md)

每个模型有三个同名文件：`xxx_1.c` 手动建模，`xxx_2.c` 从文件读取，`xxx.mps` 保存模型。
18 个 C 文件都是独立程序，直接使用公开 API，不需要基准驱动或数据准备脚本。

在仓库根目录构建并运行两个版本：

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DLP_SIMPLEX_BUILD_EXAMPLES=ON
cmake --build build -j 4
./build/data/simple_examples/example_production_planning_manual
./build/data/simple_examples/example_production_planning_mps data/simple_examples/production_planning.mps
ctest --test-dir build -R '^example_' --output-on-failure
```

| 模型 | 手动建模 | MPS 读取 | 数据 |
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

关闭 `LP_SIMPLEX_BUILD_EXAMPLES` 可跳过示例。开启 `BUILD_TESTING` 时注册 18 项测试。
这些模型只用于使用演示与功能回归，不参与 Netlib、PDS、Railway 性能基准。
