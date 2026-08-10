# lp_simplex
A simple implementation of the Simplex algorithm for linear programming in C

The library solves continuous minimization problems with equality and
inequality constraints, free variables, lower bounds, upper bounds, and
two-sided bounds. It provides three independent simplex implementations:

- a two-phase tableau solver with Bland or Dantzig pricing;
- a CSC dual revised simplex solver with logical row variables, dual
  steepest-edge pricing, Harris bound flipping, and product-form basis updates.
- Pan's generalized simplex (Pan/BDA), with reliable NNLS Phase I, a
  dynamically deficient basis, minimum-norm multipliers, and orthogonal
  recovery/certification.

It can read fixed-column MPS files, including ranged rows from the `RANGES`
section.

## Build and test

BLAS and LAPACK development libraries are required.
SuiteSparse KLU is optional and detected automatically; when available it is
used for sparse basis reinversion. Otherwise the repository's sparse LU backend
is used with the same FTRAN/BTRAN interface.
SuiteSparse CHOLMOD and SuiteSparseQR are also detected for Pan/BDA. When both
are available, Pan uses sparse normalized semi-normal equations with CHOLMOD
and reserves SPQR for orthogonal recovery/certification; otherwise it builds
the portable dense fallback. Use `-DLP_SIMPLEX_USE_SPQR=OFF` to test that
fallback explicitly.

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
Select `lp_simplex_ALGORITHM_PAN_BDA` to run the full Pan/BDA path. It is an
independent solver, not the Pan-inspired anti-stalling policy inside the dual
revised implementation. Pan defaults to Dantzig violation pricing; set
`options.pricing = lp_simplex_PRICING_PAN_NORMALIZED` for normalized violation
pricing.
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
  Pan-inspired anti-stalling control, and structure-detected dualization;
- `src/degeneracy/pan/`: full Pan/BDA standard-form conversion, NNLS Phase I,
  dynamic deficient basis, minimum-norm pricing, and numerical certification;
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
degeneracy-control, tableau, and postsolve behavior. The full implementation is
specified separately in the [Pan/BDA algorithm note](./docs/pan-bda.md).

## Examples

The small example models are stored as fixed-column MPS files in
[`data/simple_examples`](./data/simple_examples). A single table-driven test,
[`test_simple_examples.c`](./test/test_simple_examples.c), loads every model
through the public MPS API and verifies its documented predicted result.

To diagnose one Netlib instance, build the optional tools (enabled by default)
and pass an MPS path to `test_netlib`:

```sh
./build/tools/test_netlib data/netlib/feasible/afiro.mps
./build/tools/test_netlib --algorithm tableau --criteria bland \
    data/netlib/feasible/afiro.mps
./build/tools/test_netlib --algorithm pan-bda \
    data/netlib/feasible/degen2.mps
./build/tools/test_netlib --iterations 500000 \
    data/netlib/feasible/25fv47.mps
```

For feasible Netlib models the tool automatically reads the reference objective
and Gurobi solve time from the bundled CSV, then reports a timing ratio. Paths
under `infeasible/` are automatically expected to be infeasible. The tool
defaults to the dual revised solver with presolve, dual steepest-edge pricing,
Pan-inspired anti-stalling control, and a 300000-pivot limit. The tableau solver remains
available for differential diagnostics and defaults to Dantzig pricing.
