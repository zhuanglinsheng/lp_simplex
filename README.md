# lp_simplex

`lp_simplex` is a C library for continuous linear programming. It supports
equality and inequality constraints, variable bounds, and fixed-column MPS
input.

## Build

Requirements are CMake 3.12 or later, BLAS, and LAPACK. SuiteSparse is
optional.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

Run the tests:

```sh
ctest --test-dir build --output-on-failure
```

Install to a chosen prefix:

```sh
cmake --install build --prefix /path/to/install
```

## Usage

The following program reads an MPS file, solves it with the default dual
revised simplex algorithm, and prints the objective value:

```c
#include <lp_simplex/lp_simplex.h>

#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv)
{
	struct lp_simplex_Options options;
	struct lp_simplex_Result result;
	struct lp_Model *model;
	double *x;
	int exit_code;

	if (argc != 2) {
		fprintf(stderr, "usage: %s MODEL.mps\n", argv[0]);
		return 2;
	}

	model = lp_read_mps(argv[1]);
	if (model == NULL)
		return 1;
	x = (double *)malloc((size_t)model->n * sizeof(double));
	if (x == NULL) {
		lp_model_free(model);
		return 1;
	}

	lp_simplex_default_options(&options,
		lp_simplex_ALGORITHM_DUAL_REVISED);
	exit_code = lp_simplex_solve(model, &options, x, &result);
	if (exit_code == lp_simplex_EXIT_SUCCESS)
		printf("objective = %.15g\n", result.objective);
	else
		fprintf(stderr, "solve failed, status = %d\n", result.status);

	free(x);
	lp_model_free(model);
	return exit_code == lp_simplex_EXIT_SUCCESS ? 0 : 1;
}
```

If the library is installed under `/path/to/install`, compile the example with:

```sh
cc solve_mps.c -I/path/to/install/include -L/path/to/install/lib \
    -llp_simplex -o solve_mps
```

To select another algorithm, change the argument passed to
`lp_simplex_default_options`:

```c
lp_simplex_default_options(&options, lp_simplex_ALGORITHM_PAN_BDA);
/* Or lp_simplex_ALGORITHM_TABLEAU. */
```

`lp_simplex_default_options` initializes every option. The iteration limit and
presolve setting can then be overridden if needed:

```c
options.iteration_limit = 100000;
options.presolve = 0;
```

The solver accepts continuous variables only. `x` must contain at least
`model->n` elements. Release the model with `lp_model_free` after use.

## Command-line example

When the command-line tools are enabled, bundled MPS instances can be solved
directly:

```sh
./build/tools/test_netlib data/netlib/feasible/afiro.mps
./build/tools/test_netlib --algorithm pan-bda \
    data/netlib/feasible/afiro.mps
./build/tools/test_netlib --algorithm tableau --criteria bland \
    data/netlib/feasible/afiro.mps
```

Show all available options:

```sh
./build/tools/test_netlib --help
```

See the [documentation index](./docs/README.md) for algorithmic and numerical
details.
