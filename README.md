# lp_simplex
A simple implementation of the Simplex algorithm for linear programming in C

The library solves continuous minimization problems with equality and
inequality constraints, free variables, lower bounds, upper bounds, and
two-sided bounds. It provides two independent simplex implementations:

- a two-phase tableau solver with Bland or Dantzig pricing;
- a CSC dual revised simplex solver with logical row variables, dual
  steepest-edge pricing, Harris bound flipping, and product-form basis updates.

It can read fixed-column MPS files (the `RANGES` section is not supported).

## Build and test

BLAS and LAPACK development libraries are required.
SuiteSparse KLU is optional and detected automatically; when available it is
used for sparse basis reinversion. Otherwise the repository's sparse LU backend
is used with the same FTRAN/BTRAN interface.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build
ctest --test-dir build --output-on-failure
```

Use `-DBUILD_TESTING=OFF` for a library-only build. Integer and binary variable
types may be represented by the model structure but are rejected by the LP
solver rather than silently relaxed.

## Solver API

The public solve interface uses explicit options and a structured result:

```c
struct lp_simplex_Options options;
struct lp_simplex_Result result;
double *x;

lp_simplex_default_options(&options, lp_simplex_ALGORITHM_DUAL_REVISED);
options.iteration_limit = 100000;

if (lp_simplex_solve(model, &options, x, &result) ==
    lp_simplex_EXIT_SUCCESS) {
	printf("objective = %.15g, iterations = %d\n",
	       result.objective, result.iterations);
}
```

Select `lp_simplex_ALGORITHM_TABLEAU` to use the original solver. Its default
pricing rule is Bland; set `options.pricing` to
`lp_simplex_PRICING_DANTZIG` when desired.

## Source layout

The public headers under `include/lp_simplex/` separate the model, MPS reader,
termination status, and simplex solve API. The implementation is divided by
its current responsibilities:

- `model.c` and `mps.c`: model ownership and fixed-column MPS input;
- `simplex.c`: public validation, default options, and algorithm dispatch;
- `simplex_tableau_solver.c` and `simplex_transform.c`: tableau orchestration,
  bound conversion, and recovery of original variables;
- `simplex_tableau.c`, `simplex_phase.c`, and `simplex_pivot.c`: tableau
  construction, the two-phase procedure, and pivot iterations;
- `simplex_csc.c`: immutable CSC construction and augmented-column access;
- `simplex_basis.c` and `simplex_sparse_lu.c`: the basis-factorization boundary,
  sparse LU, FTRAN/BTRAN, product-form updates, and periodic reinversion;
- `simplex_singleton_dual.c`: structure-detected dualization of equality models
  with paired positive/negative singleton residual columns;
- `simplex_dual.c`: dual crash, dual steepest-edge leaving selection, Harris
  ratio testing, bound flipping, and revised-simplex iterations;
- `linalg.c` and `utils.c`: private numerical and support routines.

Tableau and linear-algebra internals are deliberately absent from the installed
public API.

## Examples

The small example models are stored as fixed-column MPS files in
[`data/simple_examples`](./data/simple_examples). A single table-driven test,
[`test_simple_examples.c`](./test/test_simple_examples.c), loads every model
through the public MPS API and verifies its documented predicted result.

To diagnose one Netlib instance, build the optional tools (enabled by default)
and pass an MPS path to `test_netlib`:

```sh
./build/tools/test_netlib data/netlib/feasible/afiro.mps
./build/tools/test_netlib --algorithm dual-revised \
    data/netlib/feasible/afiro.mps
./build/tools/test_netlib --criteria dantzig --iterations 200000 \
    data/netlib/feasible/25fv47.mps
```

For feasible Netlib models the tool automatically reads the reference objective
and Gurobi solve time from the bundled CSV, then reports a timing ratio. Paths
under `infeasible/` are automatically expected to be infeasible. Bland remains
the conservative default; Dantzig pricing is faster on some models but may be
less stable on degenerate instances.
