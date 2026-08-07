# lp_simplex
A simple implementation of the Simplex algorithm for linear programming in C

The library solves continuous minimization problems with equality and
inequality constraints, free variables, lower bounds, upper bounds, and
two-sided bounds. It uses a two-phase tableau simplex implementation and can
read fixed-column MPS files (the `RANGES` section is not supported).

## Build and test

BLAS and LAPACK development libraries are required.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build
ctest --test-dir build --output-on-failure
```

Use `-DBUILD_TESTING=OFF` for a library-only build. The default pivot rule is
Bland's rule; `dantzig` can be selected explicitly. Integer and binary variable
types may be represented by the model structure but are rejected by the LP
solver rather than silently relaxed.

## Source layout

The public headers under `include/lp_simplex/` separate the model, MPS reader,
termination status, and simplex solve API. The implementation is divided by
its current responsibilities:

- `model.c` and `mps.c`: model ownership and fixed-column MPS input;
- `simplex.c` and `simplex_transform.c`: solve orchestration, bound conversion,
  and recovery of the original variables;
- `simplex_tableau.c`, `simplex_phase.c`, and `simplex_pivot.c`: tableau
  construction, the two-phase procedure, and pivot iterations;
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
./build/tools/test_netlib --criteria dantzig --iterations 200000 \
    data/netlib/feasible/25fv47.mps
```

For feasible Netlib models the tool automatically reads the reference objective
and Gurobi solve time from the bundled CSV, then reports a timing ratio. Paths
under `infeasible/` are automatically expected to be infeasible. Bland remains
the conservative default; Dantzig pricing is faster on some models but may be
less stable on degenerate instances.
