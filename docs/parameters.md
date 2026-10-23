# 参数与结果说明

简体中文 | [English](./parameters_en.md) | [使用指南](./usage.md) | [文档首页](./README.md)

本文集中说明建模参数、求解选项、返回状态和结果诊断，并提供构建选项与安装排查。
第一次求解模型请先阅读[使用指南](./usage.md)；本文按接口字段组织，供配置和查阅使用。
公开定义以 [model.h](../include/lp_simplex/model.h)、[solve.h](../include/lp_simplex/solve.h)
和 [status.h](../include/lp_simplex/status.h) 为准。

## 1. 模型参数

`lp_model_create(m, n)` 创建 `m` 条约束、`n` 个变量；两者均须为正整数。
默认连续变量、下界 0，目标为最小化。稠密矩阵填写或修改后调用 `lp_model_build_sparse`。

| 字段 | 含义 |
| --- | --- |
| `objective[j]` | 内部最小化形式中第 `j` 个变量的目标系数 |
| `objective_offset` | 同一目标方向下的常数项，默认 0 |
| `objective_sense` | `+1` 为最小化；最大化时按下面约定存储负目标并设为 `-1` |
| `constraints[i].coef[j]` | 第 `i` 条约束中第 `j` 个变量的系数；下标从 0 开始 |
| `constraints[i].rhs` | 约束右端常数 |
| `constraints[i].type` | `optm_CONS_T_LE` / `optm_CONS_T_EQ` / `optm_CONS_T_GE` 分别为 `<=` / `=` / `>=` |
| `bounds[j].b_type` | 有效变量界的种类，见下表 |
| `bounds[j].lb` / `.ub` | 下界 / 上界，只有相应类型启用的界参与求解 |
| `bounds[j].v_type` | 默认 `optm_VAR_T_REAL`；当前求解器拒绝整数或二进制变量 |

| 变量界类型 | 生效的界 | 示例 |
| --- | --- | --- |
| `optm_BOUND_T_FR` | 无 | 自由变量 |
| `optm_BOUND_T_LO` | `lb` | 默认 `x >= 0`，也可指定负下界 |
| `optm_BOUND_T_UP` | `ub` | 仅上界，允许向负方向无界 |
| `optm_BOUND_T_BS` | `lb` 与 `ub` | 双侧界；二者相等表示固定变量 |

下面片段在 `lp_model_create` 后、求解前使用；所有数值必须有限、界必须一致。

```c
/* X1 is free; its returned value is x[0]. */
model->bounds[0].b_type = optm_BOUND_T_FR;

/* 1 <= X2 <= 5; its returned value is x[1]. */
model->bounds[1].b_type = optm_BOUND_T_BS;
model->bounds[1].lb = 1.;
model->bounds[1].ub = 5.;

/* Equality and lower-bound row types are also supported. */
model->constraints[0].type = optm_CONS_T_EQ;
model->constraints[1].type = optm_CONS_T_GE;
```

固定变量使用 `optm_BOUND_T_BS`，并令 `lb == ub`。仅上界使用 `optm_BOUND_T_UP`，
仅下界使用 `optm_BOUND_T_LO`；类型决定哪些界有效，不能仅改 `lb/ub` 而保留错误类型。

程序化构造真正的最大化模型时，按公开约定将目标系数与目标常数取负，
并设置 `objective_sense = -1`；`result.objective` 将恢复成最大化方向。
[生产计划示例](../data/simple_examples/production_planning_1.c)已经存储负利润系数，因此若想让 API 按最大化利润返回，只需在求解前设置：

```c
/* Coefficients remain -3 and -2; do not negate them again. */
model->objective_sense = -1.;
```

此时报告值为 `+22.5`，最优变量仍为 `X1 = X2 = 4.5`。这与前面的默认最小化示例是两种目标方向约定，
比较结果时要保持一致。
读取带 `OBJSENSE MAXIMIZE` 的 MPS 时，读取器会自动处理，不要再次取负。

## 2. 求解选项：lp_simplex_Options

先调用默认配置函数，再覆盖需要调整的字段：

```c
lp_simplex_default_options(&options, lp_simplex_ALGORITHM_AUTO);
options.iteration_limit = 300000;
options.primal_tolerance = 1e-7;
options.dual_tolerance = 1e-7;
options.pivot_tolerance = 1e-9;
options.presolve = 1;
```

| 字段 | 默认值 | 作用与约束 |
| --- | --- | --- |
| `algorithm` | 由 `lp_simplex_default_options` 的第二个参数选择 | 可选算法见下表 |
| `pricing` | 随算法选择 | 只接受与算法匹配的枚举，见下表 |
| `iteration_limit` | AUTO 为 300000；其他算法为 100000 | 必须大于 0；是累计工作预算，不是墙钟秒数，AUTO 各阶段及重试共用预算 |
| `primal_tolerance` | `1e-7` | 原始可行性容差；最终原模型最大绝对约束/界违反量须不大于此值 |
| `dual_tolerance` | `1e-7` | 内核搜索和最优性判断所用的对偶容差；不能直接当作所有公开对偶诊断的绝对阈值 |
| `pivot_tolerance` | `1e-9` | 数值换基/主元相关容差 |
| `presolve` | 1 | 只能为 0 或 1：0 关闭预处理，1 开启 |

三种容差都必须有限且严格大于 0。选项全清零不能得到有效默认配置。
公开 C API 没有墙钟超时字段，单次分解也会耗时；三方基准脚本另外在进程外实施限时。

| 算法 | `algorithm` 枚举 | 允许的 `pricing`；默认项加粗 |
| --- | --- | --- |
| 自动调度 | `lp_simplex_ALGORITHM_AUTO` | **`lp_simplex_PRICING_DUAL_STEEPEST_EDGE`** |
| 对偶修正单纯型 | `lp_simplex_ALGORITHM_DUAL_REVISED` | **`lp_simplex_PRICING_DUAL_STEEPEST_EDGE`** |
| Pan/BDA | `lp_simplex_ALGORITHM_PAN_BDA` | **`lp_simplex_PRICING_DANTZIG`**、`lp_simplex_PRICING_PAN_NORMALIZED` |
| 两阶段单纯型表 | `lp_simplex_ALGORITHM_TABLEAU` | **`lp_simplex_PRICING_BLAND`**、`lp_simplex_PRICING_DANTZIG` |

切换算法时重新调用默认配置函数，再设置自定义预算与容差，以免沿用不匹配的 pricing。
内部实验环境开关不属于这套公开参数接口。

## 3. API 返回值与状态码

```c
state = lp_simplex_solve(model, &options, x, &result);
if (state == lp_simplex_EXIT_SUCCESS && result.status == lp_simplex_Success) {
    /* Use the solved point and objective. */
}
```

API 返回码 `lp_simplex_EXIT_SUCCESS` 为 0，`lp_simplex_EXIT_FAILURE` 为 -1。
必须同时检查返回码与 `result.status`。失败时的有限候选点不能当作最优解。

| `result.status` | 值 | 含义与处理 |
| --- | ---: | --- |
| `lp_simplex_Success` | 0 | 成功，通过统一的最终验收 |
| `lp_simplex_MemoryAllocError` | 1 | 分配失败；检查模型规模及稠密存储 |
| `lp_simplex_CondUnsatisfied` | 2 | 参数、输入或前提条件不满足 |
| `lp_simplex_ExceedIterLimit` | 3 | 工作预算耗尽，不能据此推断不可行或无界 |
| `lp_simplex_Singularity` | 4 | 数值分解/基奇异失败，检查尺度并尝试其他算法 |
| `lp_simplex_OverDetermination` | 5 | 后端约束系统处理失败；不能只凭此状态认定原模型不可行 |
| `lp_simplex_Unboundedness` | 6 | 求解器判定目标无界 |
| `lp_simplex_Infeasibility` | 7 | 求解器判定模型不可行 |
| `lp_simplex_Degeneracy` | 8 | 退化相关失败状态，不是最优状态 |
| `lp_simplex_PrecisionError` | 9 | 数值精度相关失败，保留模型和参数以便复现 |

## 4. 结果字段与统一验收契约

所有 `lp_simplex_solve` 算法入口（包括 AUTO、dual、Pan 和 tableau）共用
`simplex_finalize_result`，不因算法不同放宽最终原模型精度。

| 字段 | 统一含义 |
|---|---|
| status / 返回码 | Success 必须对应 EXIT_SUCCESS；候选可行点不能代替最优状态 |
| objective | 对返回 x 在原模型重算，包括目标常数与最大化方向；无有限点时 NaN |
| primal_infeasibility | 原模型行约束和变量界的最大绝对违反量；无有限点时为无穷 |
| dual_infeasibility | 明确模型空间中的 KKT 符号诊断，必须连同 dual_residual_space 解读；未提供时 NaN |
| dual_residual_space | ORIGINAL / PRESOLVED / STANDARD / SCALED / AUXILIARY / UNAVAILABLE |
| iterations | 本次调用累计工作量，包括后端与原模型重试；重试不重置预算 |

最终成功必须经过内核的原目标最优性判定，且原模型目标有限、原始残差不大于
`options.primal_tolerance`；已提供的对偶诊断也必须有限。旧入口的最终
`10 * primal_tolerance` 放宽已移除。所有容差参数必须有限且严格大于零。
失败时若仍返回有限候选点，objective 和 primal 也描述该点，状态仍明确为失败。
输入无效时结果初始化为“未提供”，不会把调用者已有的 x 当作求解结果。

双单纯型在恢复行乘子缩放后，按输入问题原单位重算对偶诊断。原生 BDA 使用
同一个检查函数；标准形 Pan 在最终点重算标准形诊断。符号残差采用归一成最小化
后的成本：对活动下界检查 reduced cost 非负，对活动上界检查非正，对自由/内点
变量检查为零；固定变量无需限制 reduced cost。活动距离使用 primal tolerance。
非等式行还检查乘子方向及非活动行乘子为零。标准形对应非负变量和等式行。

内部 dual 使用 `c-A^T*pi`，原生 BDA 使用 `c+A^T*y`；检查函数显式接收符号，
数值意义一致。内核使用自身的缩放和成本尺度容差进行搜索及最优性判定；公开
绝对符号残差不是相对间隙，不应仅用 `dual_infeasibility <= dual_tolerance`
替代内核证书检查。AUTO 也不再用这个有量纲的绝对诊断额外改变内核验收。

预处理后未重建原模型行乘子的诊断标为 PRESOLVED，不能冒称 ORIGINAL。
标准型变量/宽度行未恢复时标为 STANDARD。预处理直接求出解、tableau 未导出
乘子时标为 UNAVAILABLE；这表示缺少该公开诊断，不表示数值为零，也不自动
撤销内核的最优性结论。完整原模型对偶与 Farkas 恢复尚未统一实现。

`primal_infeasibility` 始终在原模型上计算；`dual_residual_space` 只描述
`dual_infeasibility` 所属的模型空间。

`lp_simplex_Result` 新增来源字段，C 调用方须使用更新的头文件重新编译；不保证
与旧结构体布局二进制兼容。测试覆盖四入口、原生/标准形、启用/禁用缩放、
预处理、最大化与目标常数、迭代上限候选点、非有限输入以及严格原始精度。

### 对偶诊断的模型空间

`dual_residual_space` 的完整枚举名称如下；不能根据成功状态推断模型空间。

| 枚举（前缀均为 `lp_simplex_RESIDUAL_`） | 值 | 对偶诊断所属空间 |
| --- | ---: | --- |
| `UNAVAILABLE` | 0 | 未提供；不代表残差为 0 |
| `ORIGINAL` | 1 | 原模型 |
| `PRESOLVED` | 2 | 预处理后模型 |
| `STANDARD` | 3 | 标准形 |
| `SCALED` | 4 | 缩放模型 |
| `AUXILIARY` | 5 | 辅助问题 |

MPS 诊断与基准脚本的命令行参数见[数据与基准测试说明](../data/README.md)。

## 5. 构建、安装和输入排查

| CMake 参数 | 默认值 | 作用 |
| --- | --- | --- |
| `CMAKE_BUILD_TYPE` | 由构建配置决定 | 使用 `Release` 测性能，`Debug` 做调试 |
| `LP_SIMPLEX_BUILD_EXAMPLES` | ON | 构建 18 个独立使用示例 |
| `LP_SIMPLEX_BUILD_BENCHMARKS` | ON | 构建单模型诊断/基准驱动 |
| `BUILD_TESTING` | ON | 注册回归与示例测试 |
| `LP_SIMPLEX_USE_KLU` | ON | 尝试检测并使用 SuiteSparse KLU |
| `LP_SIMPLEX_USE_SPQR` | ON | 尝试检测并使用 CHOLMOD/SPQR |
| `CMAKE_INSTALL_PREFIX` | CMake 默认安装前缀 | 在配置阶段指定共享库与头文件安装位置 |

MPS 支持 `ROWS`、`COLUMNS`、`RHS`、`RANGES`、`BOUNDS` 和 `OBJSENSE`。
名称采用 fixed-column 布局，数值支持完整的长小数字段；不要假定任意 free-format MPS 都能读取。
`INTORG/INTEND`、`BV/LI/UI` 等整数声明会导致求解入口拒绝模型，需先明确构造连续松弛。
SOS、二次目标和多个 RHS/RANGES/BOUNDS 向量不受支持。

```sh
cmake -S . -B build -DCMAKE_INSTALL_PREFIX="$PWD/install"
cmake --build build -j 4
cmake --install build
cc -std=c90 -I install/include data/simple_examples/production_planning_1.c \
  -L install/lib -Wl,-rpath,"$PWD/install/lib" -llp_simplex -lm -o solve_installed
./solve_installed
```

在配置阶段指定安装路径，使 macOS 的动态库标识与安装位置保持一致；不要只在最后安装时更换前缀。
安装的是共享库和 `include/lp_simplex/` 公开头文件；示例程序从构建目录运行。
调用程序链接时指定安装目录的 `include`、`lib`，并为动态加载器设置 rpath 或平台对应的库搜索路径。

- 找不到示例程序：确认开启 `LP_SIMPLEX_BUILD_EXAMPLES`，并使用实际的 `build/data/simple_examples/` 路径。
- 找不到共享库：核对编译时 `-L` 和运行时 rpath；更换构建目录后同步修改路径。
- MPS 读取失败：核对固定列格式、整数/二次声明，以及文件是否真正完成格式转换。
- 希望测试其他小模型：直接使用 `data/simple_examples/*.mps`；大数据准备与三方测试见[数据说明](../data/README.md)。
- 准备运行 Netlib 回归测试：在 CMake 配置前先运行 `python3 data/netlib/prepare.py --all`，
  再配置、构建并执行 `ctest --test-dir build --output-on-failure`，以便注册依赖数据文件的测试。
