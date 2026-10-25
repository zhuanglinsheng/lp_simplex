# 使用指南

简体中文 | [English](./usage_en.md) | [项目主页](../README.md) | [文档首页](./README.md)

下面以生产计划为例，依次介绍数学模型、C API 建模、MPS 读取，以及编译运行。
完整源文件和 MPS 都在 [`data/simple_examples/`](../data/simple_examples/) 中，
可以直接编译，无需从文档复制代码。

## 1. 数学模型

生产两种产品，数量分别为 `x1`、`x2`，单位利润分别为 3、2。
资源容量及生产上限给出以下连续线性规划：

```text
maximize    3 x1 + 2 x2
subject to  x1 + x2 <= 9
            3 x1 + x2 <= 18
            x1 <= 7
            x2 <= 6
            x1, x2 >= 0
```

示例采用等价的最小化形式 `min -3 x1 - 2 x2`。
下文的手动建模与 MPS 文件都使用这个最小化形式，便于直接比较结果。

最优解为 `x1 = x2 = 4.5`，求解器目标值为 `-22.5`，对应最大利润 `22.5`。


| 数学变量 | 示例名称 / MPS 名称 | C 解向量 | C 变量界 |
| --- | --- | --- | --- |
| `x1` | `X1` | `x[0]` | `model->bounds[0]` |
| `x2` | `X2` | `x[1]` | `model->bounds[1]` |

## 2. 用 C API 建立模型

手动建模的完整程序是 [`production_planning_1.c`](../data/simple_examples/production_planning_1.c)。
其中 `create_model` 按照上面的数学模型填写目标系数、约束矩阵和右端项：

该函数使用 `<lp_simplex/lp_simplex.h>` 和 `<string.h>`。

```c
static struct lp_Model *create_model(void)
{
	static const double a[4][2] = {{1., 1.}, {3., 1.}, {1., 0.}, {0., 1.}};
	static const double c[2] = {-3., -2.};
	static const double rhs[4] = {9., 18., 7., 6.};
	static const int type[4] = {
		optm_CONS_T_LE, optm_CONS_T_LE, optm_CONS_T_LE, optm_CONS_T_LE};
	struct lp_Model *model = lp_model_create(4, 2);
	int i, j;

	if (model == NULL)
		return NULL;
	/* Column j represents the same variable in C and in the MPS file. */
	memcpy(model->bounds[0].name, "X1", sizeof("X1"));
	memcpy(model->bounds[1].name, "X2", sizeof("X2"));
	for (j = 0; j < 2; j++)
		model->objective[j] = c[j];
	for (i = 0; i < 4; i++) {
		model->constraints[i].type = type[i];
		model->constraints[i].rhs = rhs[i];
		for (j = 0; j < 2; j++)
			model->constraints[i].coef[j] = a[i][j];
	}

	if (lp_model_build_sparse(model) != lp_simplex_EXIT_SUCCESS) {
		lp_model_free(model);
		return NULL;
	}
	return model;
}
```

`lp_model_create(4, 2)` 为模型分配 4 条约束、2 个变量的存储空间，变量默认连续且非负。
接下来用 `a`、`type` 和 `rhs` 分别指定每条约束的系数、比较方向和右端常数。
C 数组下标从 0 开始：`a[i][0]` 是第 `i + 1` 条约束中 `x1` 的系数，
`a[i][1]` 是同一条约束中 `x2` 的系数。

代码中的四行对应如下：

| 数组下标 `i` | 系数行 `a[i]` | 方向 `type[i]` | 右端常数 `rhs[i]` | 对应约束 |
| --- | --- | --- | ---: | --- |
| 0 | `{1, 1}` | `optm_CONS_T_LE` | 9 | `x1 + x2 <= 9` |
| 1 | `{3, 1}` | `optm_CONS_T_LE` | 18 | `3 x1 + x2 <= 18` |
| 2 | `{1, 0}` | `optm_CONS_T_LE` | 7 | `x1 <= 7` |
| 3 | `{0, 1}` | `optm_CONS_T_LE` | 6 | `x2 <= 6` |

`optm_CONS_T_LE` 表示“小于或等于”。循环将上述数组的内容写入模型：
`a[i][j]` 写入 `model->constraints[i].coef[j]`，方向和右端常数分别写入
`model->constraints[i].type` 和 `.rhs`。

目标函数单独设置：`c[0] = -3`、`c[1] = -2` 分别是 `x1`、`x2` 的目标系数，
写入 `model->objective` 后表示 `min -3 x1 - 2 x2`，不计入上面的四条约束。

最后调用 `lp_model_build_sparse(model)`，将已填写的约束矩阵转换为求解器使用的
CSC／CSR 稀疏存储。它保留非零系数及其位置，不改变模型的约束或目标。

转换先构造新缓存，成功后才替换旧缓存；分配或计数越界失败时，原稀疏数组和 `nnz` 保持不变。已经修改的稠密系数仍须成功重建后再求解。重建期间同时保留旧、新缓存，峰值内存包含两者。

模型建好后，程序分配 `model->n` 个 double 存放解，初始化选项，再调用求解器：

```c
lp_simplex_default_options(&options, lp_simplex_ALGORITHM_AUTO);
state = lp_simplex_solve(model, &options, x, &result);
```

只有返回值为 `lp_simplex_EXIT_SUCCESS` 且 `result.status == lp_simplex_Success`，才表示求解成功。

`_1.c` 是独立的手动建模程序，不读取文件。它包含完整的内存检查、求解状态检查、
目标验证和释放流程，具体代码可直接打开源文件查看。

求解入口要求矩阵、目标、目标常数、RHS 和生效的界端点为有限数，目标方向字段为 `+1` 或 `-1`。用 FR/LO/UP 类型表示缺失的界，未使用的端点不参与检查。非法输入在启动进度回调和写入解向量前返回 `CondUnsatisfied`、零工作；密集、缓存稀疏和纯稀疏输入均适用。多个矩阵表示的一致性仍由调用方维护，求解器不自动协调它们。

## 3. 已有模型文件时，直接读取 MPS

同一个模型也保存在 [`production_planning.mps`](../data/simple_examples/production_planning.mps) 中。
读取文件的完整程序是 [`production_planning_2.c`](../data/simple_examples/production_planning_2.c)。
它用下面一行替代手动填矩阵的过程，之后同样初始化选项、调用 `lp_simplex_solve`：

```c
model = lp_read_mps(argv[1]);
```

`argv[1]` 是运行时提供的 MPS 路径。返回 `NULL` 表示读取失败；读取成功后，目标、
约束、变量界和稀疏矩阵已建立，不需要再调用 `lp_model_build_sparse`。
`_2.c` 只演示文件读取，不包含手动建模分支。

本例的完整 MPS 内容如下。复制或编辑时要保留 fixed-column 名称与字段的空格位置：

```text
NAME          PRODPLAN
* Production planning profit maximization
* EXPECTED STATUS: OPTIMAL
* EXPECTED OBJECTIVE: -22.5
ROWS
 N  OBJECTIV
 L  C1
 L  C2
 L  C3
 L  C4
COLUMNS
    X1        OBJECTIV            -3
    X1        C1                   1
    X1        C2                   3
    X1        C3                   1
    X2        OBJECTIV            -2
    X2        C1                   1
    X2        C2                   1
    X2        C4                   1
RHS
    RHS1      C1                   9
    RHS1      C2                  18
    RHS1      C3                   7
    RHS1      C4                   6
ENDATA
```

| 数学模型中的内容 | MPS 中的内容 |
| --- | --- |
| 目标 `-3 x1 - 2 x2` | `N OBJECTIV` 与 `COLUMNS` 中的目标系数 |
| 四条 `<=` 约束 | `L C1` 至 `L C4` |
| 每个变量参与哪些约束 | `COLUMNS` 中按变量列出的非零系数 |
| 右端项 `9, 18, 7, 6` | `RHS` |
| 非负变量 | 未声明 `BOUNDS` 时的默认下界 0 |

手动模型与文件模型由调用者拥有，最后都要用 `lp_model_free` 释放。
解向量 `x` 单独分配和释放；稠密矩阵行指针 `constraints[i].coef` 由模型拥有，不能单独释放。

## 4. 编译并运行两个程序

需要 C 编译器、CMake 3.12+、BLAS 和 LAPACK；SuiteSparse KLU/SPQR 是可选稀疏后端。
以下命令在仓库根目录执行：

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DLP_SIMPLEX_BUILD_EXAMPLES=ON
cmake --build build -j 4

# 手动建模版本，不需要输入文件。
./build/data/simple_examples/example_production_planning_manual

# MPS 版本，显式传入同一模型的文件。
./build/data/simple_examples/example_production_planning_mps data/simple_examples/production_planning.mps
```

两个程序都应得到如下结果：

```text
objective=-22.5 primal=0 iterations=2
X1=4.5
X2=4.5
```

迭代数与残差最后几位可能随构建和算法变化。两个程序各自检查求解状态、原始可行性
和本例参考目标，全部满足才返回退出码 0。

也可以只构建这两个目标：

```sh
cmake --build build --target example_production_planning_manual example_production_planning_mps
```

若想直接用 C 编译器编译源文件，可链接已经构建的库：

```sh
cc -std=c90 -I include data/simple_examples/production_planning_1.c \
  -L build -Wl,-rpath,"$PWD/build" -llp_simplex -lm -o production_manual
cc -std=c90 -I include data/simple_examples/production_planning_2.c \
  -L build -Wl,-rpath,"$PWD/build" -llp_simplex -lm -o production_mps
./production_manual
./production_mps data/simple_examples/production_planning.mps
```

## 5. 查看其他完整示例

每个模型固定提供三个同名文件：`xxx.mps` 是数据，`xxx_1.c` 手动建模，
`xxx_2.c` 读取 MPS。所有 C 文件均为独立程序，只依赖公开 C API。
本目录仅用于使用演示与功能回归，不参与 Netlib、PDS、Railway 性能基准。

| 示例 | 内容 | 手动建模 `_1.c` | MPS 读取 `_2.c` | 数据 |
| --- | --- | --- | --- | --- |
| `production_planning` | 生产计划与负利润最小化 | [源码](../data/simple_examples/production_planning_1.c) | [源码](../data/simple_examples/production_planning_2.c) | [MPS](../data/simple_examples/production_planning.mps) |
| `mixed_inequalities` | 混合方向不等式 | [源码](../data/simple_examples/mixed_inequalities_1.c) | [源码](../data/simple_examples/mixed_inequalities_2.c) | [MPS](../data/simple_examples/mixed_inequalities.mps) |
| `equality_system` | 等式约束 | [源码](../data/simple_examples/equality_system_1.c) | [源码](../data/simple_examples/equality_system_2.c) | [MPS](../data/simple_examples/equality_system.mps) |
| `free_variables` | 自由变量 | [源码](../data/simple_examples/free_variables_1.c) | [源码](../data/simple_examples/free_variables_2.c) | [MPS](../data/simple_examples/free_variables.mps) |
| `free_and_lower_bounds` | 自由变量与负下界 | [源码](../data/simple_examples/free_and_lower_bounds_1.c) | [源码](../data/simple_examples/free_and_lower_bounds_2.c) | [MPS](../data/simple_examples/free_and_lower_bounds.mps) |
| `variable_bounds` | 单侧界与双侧界 | [源码](../data/simple_examples/variable_bounds_1.c) | [源码](../data/simple_examples/variable_bounds_2.c) | [MPS](../data/simple_examples/variable_bounds.mps) |
| `ranged_rows` | 范围约束的双不等式表示 | [源码](../data/simple_examples/ranged_rows_1.c) | [源码](../data/simple_examples/ranged_rows_2.c) | [MPS](../data/simple_examples/ranged_rows.mps) |
| `degenerate_pivot` | 退化换基 | [源码](../data/simple_examples/degenerate_pivot_1.c) | [源码](../data/simple_examples/degenerate_pivot_2.c) | [MPS](../data/simple_examples/degenerate_pivot.mps) |
| `energy_system` | 较完整的能源系统模型 | [源码](../data/simple_examples/energy_system_1.c) | [源码](../data/simple_examples/energy_system_2.c) | [MPS](../data/simple_examples/energy_system.mps) |

其他程序沿用 `example_<名称>_manual` 与 `example_<名称>_mps` 的命名方式。
例如测试自由变量：

```sh
./build/data/simple_examples/example_free_variables_manual
./build/data/simple_examples/example_free_variables_mps data/simple_examples/free_variables.mps
```

`ranged_rows_1.c` 将四条区间约束分别拆成 `>=` 与 `<=` 两条普通约束；
`ranged_rows_2.c` 读取原始 `RANGES` 段，两种表示描述同一个可行域。
`energy_system_1.c` 用非零条目列表填写矩阵，便于查阅较大示例。
开启 `BUILD_TESTING` 时，CMake 为 9 组程序注册 18 项直接建模 / MPS 读取测试：

```sh
ctest --test-dir build -R '^example_' --output-on-failure
```

模型字段、算法与容差、返回状态、结果诊断和安装排查见[参数与结果说明](./parameters.md)。
