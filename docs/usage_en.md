# Usage Guide

[简体中文](./usage.md) | English | [Project home](../README_en.md) | [Documentation](./README_en.md)

This guide follows one production-planning problem from its mathematical model to
C API construction, MPS input, and finally compilation and execution. All complete
sources and MPS files are in [`data/simple_examples/`](../data/simple_examples/),
so no source code needs to be copied out of the tutorial.

## 1. Define the mathematical model

Let `x1` and `x2` be production quantities with unit profits 3 and 2:

```text
maximize    3 x1 + 2 x2
subject to  x1 + x2 <= 9
            3 x1 + x2 <= 18
            x1 <= 7
            x2 <= 6
            x1, x2 >= 0
```

The examples use the equivalent minimization model `min -3 x1 - 2 x2`.
Its optimum is `x1 = x2 = 4.5`, with objective `-22.5` and business profit `22.5`.
Both the manual C model and the MPS file use this same minimization form.

| Mathematical variable | Example / MPS name | C solution vector | C variable bound |
| --- | --- | --- | --- |
| `x1` | `X1` | `x[0]` | `model->bounds[0]` |
| `x2` | `X2` | `x[1]` | `model->bounds[1]` |

## 2. Construct the model with the C API

The complete manual-modeling program is [`production_planning_1.c`](../data/simple_examples/production_planning_1.c).
Its `create_model` function fills the objective, constraint matrix, and right-hand sides:

The function uses `<lp_simplex/lp_simplex.h>` and `<string.h>`.

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

`lp_model_create(4, 2)` allocates storage for four constraints and two variables,
which default to continuous and nonnegative. The arrays `a`, `type`, and `rhs`
specify each constraint's coefficients, comparison direction, and right-hand side.
C array indices start at zero: `a[i][0]` is the coefficient of `x1` in constraint
`i + 1`, and `a[i][1]` is the coefficient of `x2` in that constraint:

| Array index `i` | Coefficient row `a[i]` | Direction `type[i]` | RHS `rhs[i]` | Constraint |
| --- | --- | --- | ---: | --- |
| 0 | `{1, 1}` | `optm_CONS_T_LE` | 9 | `x1 + x2 <= 9` |
| 1 | `{3, 1}` | `optm_CONS_T_LE` | 18 | `3 x1 + x2 <= 18` |
| 2 | `{1, 0}` | `optm_CONS_T_LE` | 7 | `x1 <= 7` |
| 3 | `{0, 1}` | `optm_CONS_T_LE` | 6 | `x2 <= 6` |

`optm_CONS_T_LE` means “less than or equal to.” The loop copies `a[i][j]` into
`model->constraints[i].coef[j]`, and copies the direction and right-hand side into
`model->constraints[i].type` and `.rhs`.

The objective is set separately: `c[0] = -3` and `c[1] = -2` are the objective
coefficients of `x1` and `x2`. Copying them into `model->objective` defines
`min -3 x1 - 2 x2`; the objective is not one of the four constraints.

Finally, `lp_model_build_sparse(model)` converts the populated constraint matrix
into the solver's sparse storage, retaining nonzero coefficients and their
positions without changing the constraints or objective.

Conversion constructs replacement CSC/CSR buffers before committing them. Allocation or count-overflow failure preserves the previous sparse arrays and `nnz`. Pending dense edits still require a successful rebuild before solving. Rebuilding temporarily retains both old and new caches; peak memory includes both.

The program allocates `model->n` doubles for the solution, initializes options,
and solves the model:

```c
lp_simplex_default_options(&options, lp_simplex_ALGORITHM_AUTO);
state = lp_simplex_solve(model, &options, x, &result);
```

Success requires both `state == lp_simplex_EXIT_SUCCESS` and
`result.status == lp_simplex_Success`.

The `_1.c` program only constructs the model manually. Its source includes allocation
checks, status checks, example-objective validation, and memory cleanup.

The solve entry requires finite matrix values, costs, objective offset, RHS and active bound endpoints, with objective sense `+1` or `-1`. Encode absent endpoints using FR/LO/UP types; inactive endpoints are ignored. Invalid input returns `CondUnsatisfied` with zero work before progress callbacks or writes to x, for dense, cached-sparse and sparse-only models. Callers still maintain agreement between matrix representations; the solver does not reconcile them.

## 3. Read MPS when the model is already in a file

The same problem is stored in [`production_planning.mps`](../data/simple_examples/production_planning.mps).
The separate reader program is [`production_planning_2.c`](../data/simple_examples/production_planning_2.c).
It replaces manual construction with one call before using the same solve API:

```c
model = lp_read_mps(argv[1]);
```

`argv[1]` is the MPS path supplied at runtime. `NULL` indicates a read failure.
On success, the objective, constraints, bounds, and sparse storage are ready;
no additional `lp_model_build_sparse` call is required. The `_2.c` program contains
only the file-input path, without a manual-modeling branch.

Here is the complete MPS file. Preserve the fixed-column field spacing when copying or editing:

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

| Model component | MPS representation |
| --- | --- |
| Objective `-3 x1 - 2 x2` | `N OBJECTIV` and objective entries in COLUMNS |
| Four `<=` constraints | `L C1` through `L C4` |
| Nonzero coefficients by variable | COLUMNS |
| Right-hand sides `9, 18, 7, 6` | RHS |
| Nonnegative variables | Default lower bound zero without BOUNDS |

Both model constructors transfer ownership to the caller: release either model
with `lp_model_free`. Allocate and free the solution vector separately. Dense
`constraints[i].coef` row pointers belong to the model and must not be freed separately.

## 4. Compile and run the two programs

Requirements are a C compiler, CMake 3.12+, BLAS, and LAPACK; SuiteSparse KLU/SPQR
are optional sparse backends. Run from the repository root:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DLP_SIMPLEX_BUILD_EXAMPLES=ON
cmake --build build -j 4

# Manual modeling: no input file.
./build/data/simple_examples/example_production_planning_manual

# MPS input: pass the matching file explicitly.
./build/data/simple_examples/example_production_planning_mps data/simple_examples/production_planning.mps
```

Both programs should produce:

```text
objective=-22.5 primal=0 iterations=2
X1=4.5
X2=4.5
```

The last residual digits and iteration count may vary with the build and algorithm.
Each program validates solve status, original-model feasibility, and its example's
reference objective before returning exit code zero.

To build only this pair:

```sh
cmake --build build --target example_production_planning_manual example_production_planning_mps
```

Alternatively, compile each source directly against the built shared library:

```sh
cc -std=c90 -I include data/simple_examples/production_planning_1.c \
  -L build -Wl,-rpath,"$PWD/build" -llp_simplex -lm -o production_manual
cc -std=c90 -I include data/simple_examples/production_planning_2.c \
  -L build -Wl,-rpath,"$PWD/build" -llp_simplex -lm -o production_mps
./production_manual
./production_mps data/simple_examples/production_planning.mps
```

## 5. Browse the complete examples

Every model has three matching files: `xxx.mps` for data, `xxx_1.c` for manual
modeling, and `xxx_2.c` for MPS input. Each C file is a standalone program using
the public API. These are usage examples and functional tests, excluded from
Netlib, PDS, and Railway performance benchmarks.

| Example | Topic | Manual `_1.c` | MPS input `_2.c` | Data |
| --- | --- | --- | --- | --- |
| `production_planning` | Production planning and negative-profit minimization | [Source](../data/simple_examples/production_planning_1.c) | [Source](../data/simple_examples/production_planning_2.c) | [MPS](../data/simple_examples/production_planning.mps) |
| `mixed_inequalities` | Mixed inequality directions | [Source](../data/simple_examples/mixed_inequalities_1.c) | [Source](../data/simple_examples/mixed_inequalities_2.c) | [MPS](../data/simple_examples/mixed_inequalities.mps) |
| `equality_system` | Equality constraints | [Source](../data/simple_examples/equality_system_1.c) | [Source](../data/simple_examples/equality_system_2.c) | [MPS](../data/simple_examples/equality_system.mps) |
| `free_variables` | Free variables | [Source](../data/simple_examples/free_variables_1.c) | [Source](../data/simple_examples/free_variables_2.c) | [MPS](../data/simple_examples/free_variables.mps) |
| `free_and_lower_bounds` | Free variables and negative lower bounds | [Source](../data/simple_examples/free_and_lower_bounds_1.c) | [Source](../data/simple_examples/free_and_lower_bounds_2.c) | [MPS](../data/simple_examples/free_and_lower_bounds.mps) |
| `variable_bounds` | One-sided and two-sided bounds | [Source](../data/simple_examples/variable_bounds_1.c) | [Source](../data/simple_examples/variable_bounds_2.c) | [MPS](../data/simple_examples/variable_bounds.mps) |
| `ranged_rows` | Ranged rows as pairs of inequalities | [Source](../data/simple_examples/ranged_rows_1.c) | [Source](../data/simple_examples/ranged_rows_2.c) | [MPS](../data/simple_examples/ranged_rows.mps) |
| `degenerate_pivot` | Degenerate pivots | [Source](../data/simple_examples/degenerate_pivot_1.c) | [Source](../data/simple_examples/degenerate_pivot_2.c) | [MPS](../data/simple_examples/degenerate_pivot.mps) |
| `energy_system` | A larger energy-system model | [Source](../data/simple_examples/energy_system_1.c) | [Source](../data/simple_examples/energy_system_2.c) | [MPS](../data/simple_examples/energy_system.mps) |

Other targets follow the same `example_<name>_manual` and `example_<name>_mps` convention:

```sh
./build/data/simple_examples/example_free_variables_manual
./build/data/simple_examples/example_free_variables_mps data/simple_examples/free_variables.mps
```

`ranged_rows_1.c` expands each ranged row into ordinary lower and upper constraints;
`ranged_rows_2.c` reads the RANGES section directly. Both represent the same feasible
region. `energy_system_1.c` uses a list of nonzero entries to keep its larger matrix readable.
With BUILD_TESTING enabled, CMake registers 18 tests for the nine manual/MPS pairs:

```sh
ctest --test-dir build -R '^example_' --output-on-failure
```

For model fields, algorithms, tolerances, statuses, diagnostics, installation, and
troubleshooting, see the [parameters and results reference](./parameters_en.md).
