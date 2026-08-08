# lp_simplex
A simple implementation of the Simplex algorithm for linear programming in C

The library solves continuous minimization problems with equality and
inequality constraints, free variables, lower bounds, upper bounds, and
two-sided bounds. It provides two independent simplex implementations:

- a two-phase tableau solver with Bland or Dantzig pricing;
- a CSC dual revised simplex solver with logical row variables, dual
  steepest-edge pricing, Harris bound flipping, and product-form basis updates.

It can read fixed-column MPS files, including ranged rows from the `RANGES`
section.

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
options.presolve = 1;

if (lp_simplex_solve(model, &options, x, &result) ==
    lp_simplex_EXIT_SUCCESS) {
	printf("objective = %.15g, iterations = %d\n",
	       result.objective, result.iterations);
}
```

Select `lp_simplex_ALGORITHM_TABLEAU` to use the original solver. Its default
pricing rule is Bland; set `options.pricing` to
`lp_simplex_PRICING_DANTZIG` when desired.
Presolve is enabled by default; set `options.presolve = 0` for differential
diagnostics or to solve the original model directly.

## Source layout

The public headers under `include/lp_simplex/` separate the model, MPS reader,
termination status, and simplex solve API. The implementation is divided by
its current responsibilities:

- `src/core/`: model ownership, MPS input, public validation and algorithm
  dispatch, plus the solver's immutable problem representation;
- `src/presolve/`: sparse reductions, activity and bound propagation,
  substitution, postsolve journaling, and solution reconstruction;
- `src/dual/`: revised-simplex orchestration, feasibility, pricing,
  degeneracy control, and structure-detected dualization;
- `src/basis/`: the opaque basis-factorization boundary, sparse LU/KLU
  backends, FTRAN/BTRAN, product-form updates, and periodic reinversion;
- `src/matrix/`: immutable CSC and sparse-vector primitives;
- `src/tableau/`: bound transformation, two-phase tableau construction,
  pricing, pivoting, and recovery of original variables;
- `src/common/`: private numerical and support routines.

Tableau and linear-algebra internals are deliberately absent from the installed
public API.

## Algorithm documentation

See the [overall algorithm guide](./docs/algorithm.md) for the end-to-end
flowcharts and the current presolve, dual revised simplex, basis-update,
degeneracy-control, tableau, and postsolve behavior.

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
