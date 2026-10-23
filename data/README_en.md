# Data and Benchmarks

[简体中文](./README.md) | English | [Project home](../README_en.md) | [Documentation](../docs/README_en.md)

This directory maintains test data, scripts for generating MPS files on demand, and performance comparisons of lp_simplex, SCIP/SoPlex, HiGHS, Clp, and GLPK. Run all commands below from the repository root.

## Datasets

| Dataset | Contents | Performance report |
|---|---|---|
| [Netlib](netlib/README_en.md) | 114 feasible LPs, plus infeasible models and presolved variants; benchmarks use only the original feasible models | [Per-instance results](../docs/benchmarks/netlib_en.md) |
| [PDS](pds/README_en.md) | 9 large LPs with approximately 106,000–505,000 variables | [Per-instance results](../docs/benchmarks/pds_en.md) |
| [Railway](railway/README_en.md) | 7 railway set-covering problems, tested as continuous LP relaxations, with approximately 47,000–1,093,000 variables | [Per-instance results](../docs/benchmarks/railway_en.md) |
| [Small examples](simple_examples/README_en.md) | 9 functional and regression test models | Excluded from benchmarks |

Netlib, PDS, and Railway store only raw inputs and files needed for conversion. Each dataset's `prepare.py` generates MPS files on demand. `simple_examples` stores small MPS files and matching C examples directly; see the [usage guide](../docs/usage_en.md) for modeling and solver API examples.

## Prepare the environment

Build the single-model driver with CMake, a C compiler, and BLAS/LAPACK:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DLP_SIMPLEX_BUILD_BENCHMARKS=ON
cmake --build build --target test_netlib -j 4
```

The executable is `build/data/netlib/test_netlib`. Running it separately requires a prepared MPS file. Generating models from raw inputs also requires Python 3.9+; Python temporarily builds the PDS converter with a C compiler.

Benchmarks additionally require SCIP, Clp, and GLPK on PATH (or paths supplied with `--scip`, `--clp`, and `--glpk`), plus NumPy, SciPy, and highspy or a compatible SciPy HiGHS binding in the Python environment. On macOS, install Clp and GLPK with `brew install clp glpk`. For example:

```sh
python3 -m pip install numpy scipy highspy
```

## Run full benchmarks

Run every benchmark instance in each dataset:

```sh
python3 data/benchmark.py --dataset netlib
python3 data/benchmark.py --dataset pds
python3 data/benchmark.py --dataset railway
```

Defaults are 3 repetitions per instance and solver. The per-run limit is 60 seconds for Netlib and 120 seconds for PDS and Railway. Netlib uses only its 114 original feasible models; `simple_examples` is excluded from performance benchmarks. For a custom build directory, use `--binary PATH` to locate `test_netlib`.

For each instance, the script:

1. Generates a temporary MPS file from raw inputs.
2. Runs all five solvers and repetitions serially, recording timings and statuses and auditing solutions against the original model.
3. Deletes the instance's temporary MPS file.
4. Saves run records and updates the dataset reports and both language summaries before proceeding to the next instance.

No dataset-wide MPS preparation is needed. Temporary MPS files live inside the run directory; manually generated MPS files are outside this cleanup process. Dataset reports maintain the latest per-instance results while retaining test information and statistics separately for each run.

For debugging, select a single instance with fewer repetitions and a shorter limit:

```sh
python3 data/benchmark.py --dataset netlib --cases afiro --repeats 1 --seconds 10
```

`--cases` takes instance names without directories or file extensions, such as `afiro` for Netlib, `pds-30` for PDS, or `rail507` for Railway. It limits only the current run's selection.

After an interruption, retain the old run directory and start a new run with `--cases` selecting unfinished instances. Keep the original algorithms, solvers, work budgets, and time limits, and use a new directory for `--run-dir`. Unselected completed instance rows are preserved, while statistics remain separate for each run. The script cannot resume in the original run directory or automatically skip completed solvers or repetitions within an instance; an interrupted instance starts again. `--report-only` rebuilds reports without continuing solver runs.

## Check a single model

[`netlib/test_netlib.c`](netlib/test_netlib.c) solves one model and checks the expected outcome. For example, generate Afiro, then solve it and print its variables:

```sh
python3 data/netlib/prepare.py feasible/afiro
./build/data/netlib/test_netlib --print-solution data/netlib/feasible/afiro.mps
```

Here `prepare.py` requires the `feasible/` category prefix to distinguish Netlib models with the same name in different categories. The generated MPS remains in the dataset directory. See the other dataset READMEs for their preparation commands.

The driver automatically looks up available Netlib reference objectives and expects infeasibility for paths containing `/infeasible/`. Use `--expected VALUE` or `--infeasible` to specify an expectation explicitly. Without a reference objective, only a successful solve is required.

Exit code 0 means the expected outcome matched, 1 means reading, solving, or validation failed, and 2 means incorrect usage. An expected infeasible result therefore also exits 0 when validated.

## View results and rebuild reports

The [benchmark summary](../docs/benchmarks/README_en.md) provides test environments, comparison statistics, statistical definitions, reproduction commands, and links to all three dataset reports. PASS requires an optimal solver status and successful original-model feasibility, objective recomputation, and available reference-objective checks. It is not a complete independent optimality certificate.

Each run stores configurations, environments, results, logs, and solution vectors under `tmp/benchmarks/<run-directory>/`, which Git ignores. Back up runs you need to retain. Use `--run-dir` to select a new run directory and `--docs-dir` to select the report directory.

Rebuild reports from existing run records without solving again. Replace the path below with the actual run directory:

```sh
python3 data/benchmark.py --report-only tmp/benchmarks/<run-directory>
```

The script updates the corresponding dataset reports and summaries in both languages. Incomplete runs retain visible pending results.

## Command-line reference

### Three-solver benchmarks: benchmark.py

| Option | Default | Meaning |
| --- | --- | --- |
| `--dataset` | `netlib` | `netlib`, `pds`, or `railway`; excludes simple_examples |
| `--cases` | Full selected dataset | One or more raw instance names |
| `--binary` | `build/data/netlib/test_netlib` | lp_simplex single-model driver |
| `--scip` | `scip` | SCIP executable name or path |
| `--clp` | `clp` | Clp executable name or path |
| `--glpk` | `glpsol` | GLPK executable name or path |
| `--solvers` | All five solvers | Solver selection; must include `lp_simplex` |
| `--repeats` | 3 | Repetitions per instance and solver |
| `--seconds` | Netlib: 60; PDS/Railway: 120 | Per-run wall-clock deadline |
| `--iterations` | 300000 | lp_simplex work budget |
| `--lp-algorithm` | `auto` | `auto`, `dual-revised`, or `pan-bda` |
| `--highs-algorithm` | `choose` | `choose`, `simplex`, or `ipm` |
| `--scip-algorithm` | `default` | `default` or `dual` |
| `--run-dir` | `tmp/benchmarks/<UTC timestamp>` | New directory for raw run records |
| `--docs-dir` | `docs/benchmarks` | Report destination |
| `--report-only DIR` | Off | Regenerate reports without solving |

### Single-model checks: test_netlib

```sh
./build/data/netlib/test_netlib [options] MODEL.mps
```

| Option | Default | Meaning |
| --- | --- | --- |
| `--algorithm` | `auto` | `auto`, `dual-revised`, `pan-bda`, or `tableau` |
| `--criteria` | `dantzig` | Use `bland/dantzig` for manual tableau, `dantzig/normalized` for Pan; AUTO/dual select pricing automatically |
| `--iterations` | 300000 | Positive work budget overriding the selected algorithm's default |
| `--no-presolve` | Off | Disable presolve |
| `--tolerance` | `1e-7` | Relative reference-objective comparison only; does not change primal/dual tolerances |
| `--expected VALUE` | Automatic reference lookup, when available | Override reference objective |
| `--infeasible` | Off | Expect infeasibility; paths containing `/infeasible/` also select this expectation |
| `--print-solution` | Off | Print original variables |
| `--save-point PATH` | No export | Export a finite point passing original primal checks, potentially at an iteration limit; not an optimality guarantee |
| `--help` | — | Show current options |

These options belong to the test tools. See the [parameters and results reference](../docs/parameters_en.md) for C API solver options and result fields. Current tool options are also available through `python3 data/benchmark.py --help` and `./build/data/netlib/test_netlib --help`.
