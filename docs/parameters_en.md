# Parameters and Results Reference

[简体中文](./parameters.md) | English | [Usage guide](./usage_en.md) | [Documentation](./README_en.md)

This reference collects modeling parameters, solver options, statuses, diagnostics,
build options, and installation troubleshooting. Start with
the [usage guide](./usage_en.md) for an end-to-end example. Public declarations are in
[model.h](../include/lp_simplex/model.h), [solve.h](../include/lp_simplex/solve.h), and
[status.h](../include/lp_simplex/status.h).

## 1. Model parameters

`lp_model_create(m, n)` allocates `m` constraints and `n` variables; both counts must
be positive. Variables default to continuous and nonnegative, with a minimization
objective. Call `lp_model_build_sparse` after filling or changing the dense matrix.

| Field | Meaning |
| --- | --- |
| `objective[j]` | Objective coefficient in the internal minimization form |
| `objective_offset` | Constant in the same objective sense; default 0 |
| `objective_sense` | `+1` for minimization; for maximization store the negated objective and set `-1` |
| `constraints[i].coef[j]` | Coefficient of variable `j` in constraint `i`; indices start at zero |
| `constraints[i].rhs` | Right-hand side |
| `constraints[i].type` | `optm_CONS_T_LE` / `optm_CONS_T_EQ` / `optm_CONS_T_GE` mean `<=` / `=` / `>=` |
| `bounds[j].b_type` | Which bounds are active; see below |
| `bounds[j].lb` / `.ub` | Lower / upper value, used only when enabled by the bound type |
| `bounds[j].v_type` | Defaults to `optm_VAR_T_REAL`; integer and binary variables are rejected by the solver |

| Bound type | Active values | Example |
| --- | --- | --- |
| `optm_BOUND_T_FR` | None | Free variable |
| `optm_BOUND_T_LO` | `lb` | Default `x >= 0`, or a negative lower bound |
| `optm_BOUND_T_UP` | `ub` | Upper bound only; unbounded toward negative infinity |
| `optm_BOUND_T_BS` | `lb` and `ub` | Two-sided bounds; equal values fix a variable |

Use these snippets after model creation and before solving. Values must be finite
and bounds must be consistent.

```c
model->bounds[0].b_type = optm_BOUND_T_FR; /* X1 is free; its returned value is x[0]. */
model->bounds[1].b_type = optm_BOUND_T_BS;
model->bounds[1].lb = 1.;
model->bounds[1].ub = 5.;                 /* 1 <= X2 <= 5; returned as x[1]. */
model->constraints[0].type = optm_CONS_T_EQ;
model->constraints[1].type = optm_CONS_T_GE;
```

For a fixed variable use `optm_BOUND_T_BS` with `lb == ub`. Use `optm_BOUND_T_UP`
for an upper bound only, or `optm_BOUND_T_LO` for a lower bound only. The type
controls which bound values are active; changing numbers alone is insufficient.

For a programmatically constructed maximization model, negate the objective
coefficients and constant before storing them, and set `objective_sense = -1`.
The result restores the original direction. The [production-planning example](../data/simple_examples/production_planning_1.c) already stores negative-profit coefficients, so reporting
maximum profit only requires this change before solving:

```c
/* Coefficients remain -3 and -2; do not negate them again. */
model->objective_sense = -1.;
```

It reports `+22.5` with the same `X1 = X2 = 4.5`. Keep the objective sense consistent
when comparing it with the default minimization example. The MPS reader handles
`OBJSENSE MAXIMIZE` automatically; do not negate its output again.

## 2. Solver options: lp_simplex_Options

Initialize defaults, then override the required fields:

```c
lp_simplex_default_options(&options, lp_simplex_ALGORITHM_AUTO);
options.iteration_limit = 300000;
options.primal_tolerance = 1e-7;
options.dual_tolerance = 1e-7;
options.pivot_tolerance = 1e-9;
options.presolve = 1;
```

| Field | Default | Meaning and constraints |
| --- | --- | --- |
| `algorithm` | Selected by the default-options function's second argument | Supported algorithms are listed below |
| `pricing` | Depends on algorithm | Only the matching enum values below are accepted |
| `iteration_limit` | 300000 for AUTO; 100000 otherwise | Positive cumulative work budget, not seconds; shared across AUTO stages and retries |
| `primal_tolerance` | `1e-7` | Final maximum absolute original row/bound violation must not exceed this value |
| `dual_tolerance` | `1e-7` | Kernel search/optimality tolerance, not a universal absolute threshold for public dual diagnostics |
| `pivot_tolerance` | `1e-9` | Numerical pivot-related tolerance |
| `presolve` | 1 | Must be 0 or 1: 0 disables presolve, 1 enables it |

All three tolerances must be finite and strictly positive. Zero-initializing the
structure does not provide valid defaults. The public API has no wall-clock timeout;
a factorization can take substantial time. The benchmark script enforces a separate
external process deadline.

| Algorithm | `algorithm` enum | Allowed `pricing`; default in bold |
| --- | --- | --- |
| Automatic scheduling | `lp_simplex_ALGORITHM_AUTO` | **`lp_simplex_PRICING_DUAL_STEEPEST_EDGE`** |
| Dual revised simplex | `lp_simplex_ALGORITHM_DUAL_REVISED` | **`lp_simplex_PRICING_DUAL_STEEPEST_EDGE`** |
| Pan/BDA | `lp_simplex_ALGORITHM_PAN_BDA` | **`lp_simplex_PRICING_DANTZIG`**, `lp_simplex_PRICING_PAN_NORMALIZED` |
| Two-phase tableau | `lp_simplex_ALGORITHM_TABLEAU` | **`lp_simplex_PRICING_BLAND`**, `lp_simplex_PRICING_DANTZIG` |

When switching algorithms, reinitialize defaults before applying custom settings to
avoid retaining incompatible pricing. Internal experimental environment switches
are not part of this public options interface.

### AUTO stage transitions

AUTO starts with dual simplex and applies the following rules:

- Continued progress renews the same live basis and factorization, without the former 16384-iteration cutoff. One eighth of the total work budget is reserved for recovery routes.
- Feasibility restoration and original-objective cleanup track progress separately. Within an observation interval, residual reduction, a decrease in the number of violated constraints, or original-objective improvement must exceed its own noise budget and 1% threshold. Working-bound changes do not erase stagnation history.
- Material Phase-I cost repairs receive a finite recovery allowance; symbolic tie perturbations alone do not qualify.
- A negative accumulated feasibility merit triggers a recomputation of row violations. Roundoff alone cannot trigger a route change.
- An alternate initializer is considered only after dual yields without a feasible point. If finite bounds introduced by presolve substantially expand the equality form, NNLS may search the original model. A candidate must pass feasibility checks on both the original model and its mapped reduced model before entering the hybrid path.
- Original-objective cleanup trials Devex on the current basis when at least half the pivots in an observation interval have zero step and dual error decreases by less than 1%. Zero steps alone do not justify an exit.
- A near-optimal feasible point may receive a bounded multiplier correction after failing the dual audit. Failure restores the original multiplier; acceptance still requires the existing componentwise KKT and primal–dual gap checks.

Initialization, multiplier correction and recovery share `iteration_limit`. Local stagnation, numerical failure and resource limits alone prove neither infeasibility, unboundedness nor optimality. Final acceptance requirements are independent of the route.

## 3. API return values and status codes

```c
state = lp_simplex_solve(model, &options, x, &result);
if (state == lp_simplex_EXIT_SUCCESS && result.status == lp_simplex_Success) {
    /* Use the solved point and objective. */
}
```

API return values are 0 (`lp_simplex_EXIT_SUCCESS`) and -1 (`lp_simplex_EXIT_FAILURE`).
Check both the return value and `result.status`; a finite candidate on failure is
not an optimal solution.

| `result.status` | Value | Meaning and handling |
| --- | ---: | --- |
| `lp_simplex_Success` | 0 | Success under the unified acceptance contract |
| `lp_simplex_MemoryAllocError` | 1 | Allocation failure; inspect model size and dense storage |
| `lp_simplex_CondUnsatisfied` | 2 | Invalid input/options or unmet preconditions |
| `lp_simplex_ExceedIterLimit` | 3 | Work budget exhausted; does not establish infeasibility or unboundedness |
| `lp_simplex_Singularity` | 4 | Singular basis or factorization failure; inspect scaling and alternatives |
| `lp_simplex_OverDetermination` | 5 | Backend constraint-system handling failed; not by itself a proof of original-model infeasibility |
| `lp_simplex_Unboundedness` | 6 | Solver classified the objective as unbounded |
| `lp_simplex_Infeasibility` | 7 | Solver classified the model as infeasible |
| `lp_simplex_Degeneracy` | 8 | Degeneracy-related failure, not an optimal status |
| `lp_simplex_PrecisionError` | 9 | Numerical-accuracy failure; preserve input and settings for reproduction |

## 4. Result fields and the unified acceptance contract

All `lp_simplex_solve` algorithm entry points, including AUTO, dual, Pan, and tableau, share `simplex_finalize_result`. Final accuracy requirements on the original model are not relaxed for any algorithm.

| Field | Meaning |
|---|---|
| status / return code | Success must correspond to EXIT_SUCCESS; a feasible candidate does not substitute for an optimal status |
| objective | Recomputed at the returned x in the original model, including the objective constant and maximization sense; NaN when no finite point is available |
| primal_infeasibility | Maximum absolute violation of original-model row constraints and variable bounds; infinity when no finite point is available |
| dual_infeasibility | A KKT sign diagnostic in an explicitly identified model space; interpret together with dual_residual_space; NaN when unavailable |
| dual_residual_space | ORIGINAL / PRESOLVED / STANDARD / SCALED / AUXILIARY / UNAVAILABLE |
| iterations | Cumulative work for the call, including backend and original-model retries; retries do not reset the budget |

Final success requires the kernel to establish optimality for the original objective, a finite original-model objective, and an original-model primal residual no greater than `options.primal_tolerance`. Any supplied dual diagnostic must also be finite. The former final relaxation to `10 * primal_tolerance` has been removed. All tolerance parameters must be finite and strictly positive.

If a failed solve returns a finite candidate, the objective and primal diagnostic describe that candidate; the status still indicates failure. Invalid input initializes the result as unavailable rather than treating the caller's existing x as a computed solution.

After reversing row-multiplier scaling, dual simplex recomputes the dual diagnostic in the input problem's original units. Native BDA uses the same checking function; standard-form Pan recomputes its diagnostic in standard form at the final point. Sign residuals use costs normalized to minimization: reduced costs must be nonnegative at active lower bounds, nonpositive at active upper bounds, and zero for free or interior variables. Fixed variables impose no reduced-cost restriction. Bound activity is determined using the primal tolerance. For inequality rows, multiplier signs and zero multipliers on inactive rows are also checked. Standard form has nonnegative variables and equality rows.

The dual kernel uses `c-A^T*pi`, whereas native BDA uses `c+A^T*y`. The checking function receives the sign convention explicitly; the numerical meaning is the same. Kernels use their own scaling and cost-scale tolerances for search and optimality checks. The public absolute sign residual is not a relative duality gap: `dual_infeasibility <= dual_tolerance` alone must not replace kernel certificate checks. AUTO does not impose this dimensional absolute diagnostic as an additional kernel acceptance test.

Diagnostics after presolve without reconstructed original-model row multipliers are marked PRESOLVED, not ORIGINAL. If standard-form variables and bound-width rows have not been recovered, the space is STANDARD. When presolve directly solves the problem or tableau does not export multipliers, the space is UNAVAILABLE. This means the public diagnostic is missing, not zero, and does not automatically invalidate the kernel's optimality conclusion. Complete original-model dual and Farkas recovery is not yet implemented uniformly.

`primal_infeasibility` is always computed on the original model. `dual_residual_space` describes only the model space of `dual_infeasibility`.

The provenance field added to `lp_simplex_Result` requires C callers to recompile against the updated headers; binary compatibility with the previous structure layout is not guaranteed. Tests cover all four entry points, native and standard forms, scaling enabled and disabled, presolve, maximization and objective constants, iteration-limit candidates, nonfinite inputs, and strict primal accuracy.

### Dual diagnostic spaces

The complete enum names have the prefix `lp_simplex_RESIDUAL_`. Never infer the
space from solve status alone.

| Suffix | Value | Dual diagnostic space |
| --- | ---: | --- |
| `UNAVAILABLE` | 0 | Unavailable, not zero |
| `ORIGINAL` | 1 | Original model |
| `PRESOLVED` | 2 | Presolved model |
| `STANDARD` | 3 | Standard form |
| `SCALED` | 4 | Scaled model |
| `AUXILIARY` | 5 | Auxiliary problem |

[MPS diagnostic and benchmark command-line options](../data/README_en.md) are documented with the data tools.

## 5. Build, installation, and input troubleshooting

| CMake option | Default | Meaning |
| --- | --- | --- |
| `CMAKE_BUILD_TYPE` | Build configuration dependent | Use Release for timing and Debug for debugging |
| `LP_SIMPLEX_BUILD_EXAMPLES` | ON | Build 18 standalone usage programs |
| `LP_SIMPLEX_BUILD_BENCHMARKS` | ON | Build the single-model diagnostic/benchmark driver |
| `BUILD_TESTING` | ON | Register regression and example tests |
| `LP_SIMPLEX_USE_KLU` | ON | Detect and use SuiteSparse KLU when available |
| `LP_SIMPLEX_USE_SPQR` | ON | Detect and use CHOLMOD/SPQR when available |
| `CMAKE_INSTALL_PREFIX` | CMake's default prefix | Set library/header destination at configuration time |

The reader supports ROWS, COLUMNS, RHS, RANGES, BOUNDS, and OBJSENSE. Names use
a fixed-column layout; numeric tokens retain long decimal values. Arbitrary
free-format MPS is not supported. Integer declarations such as INTORG/INTEND,
BV, LI, and UI cause the solve entry point to reject the model; explicitly
construct a continuous relaxation if that is intended. SOS, quadratic objectives,
and multiple RHS/RANGES/BOUNDS vectors are unsupported.

```sh
cmake -S . -B build -DCMAKE_INSTALL_PREFIX="$PWD/install"
cmake --build build -j 4
cmake --install build
cc -std=c90 -I install/include data/simple_examples/production_planning_1.c \
  -L install/lib -Wl,-rpath,"$PWD/install/lib" -llp_simplex -lm -o solve_installed
./solve_installed
```

Set the install prefix at configuration time so the macOS library install name matches
its destination, rather than changing only the final installation prefix.
Installation provides the shared library and `include/lp_simplex/` public headers.
Run the examples from the build directory. When linking an application against an
installation, set its include/library paths and the runtime rpath or platform
library search path.

- Missing example executables: enable `LP_SIMPLEX_BUILD_EXAMPLES` and use `build/data/simple_examples/`.
- Missing shared library: check `-L` and runtime rpath, especially after changing build directories.
- MPS read failure: check fixed columns, unsupported declarations, and whether conversion actually completed.
- More small models: use `data/simple_examples/*.mps` directly. For large datasets and benchmarks, see [data instructions](../data/README_en.md).
- Netlib regression tests: run `python3 data/netlib/prepare.py --all` before configuring CMake, then configure, build, and run `ctest --test-dir build --output-on-failure` so data-dependent tests are registered.
