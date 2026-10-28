# Data and Benchmarks

[简体中文](./README.md) | English | [Project home](../README_en.md) | [Documentation](../docs/README_en.md)

This directory maintains test data, scripts for generating MPS files on demand, and performance comparisons of lp_simplex, SCIP/SoPlex, HiGHS, Clp, and GLPK. Run all commands below from the repository root.

## Datasets

| Dataset | Contents | Performance report |
|---|---|---|
| [Netlib](netlib/README_en.md) | 114 feasible LPs, plus infeasible models and presolved variants; this suite tests feasible original models only | [Per-instance results](../docs/benchmarks/netlib_en.md) |
| Netlib infeasible | 29 infeasible LPs; infeasibility decisions summarized separately | [Per-instance results](../docs/benchmarks/netlib-infeasible_en.md) |
| [PDS](pds/README_en.md) | 9 large LPs with approximately 106,000–505,000 variables | [Per-instance results](../docs/benchmarks/pds_en.md) |
| [Railway](railway/README_en.md) | 7 railway set-covering problems, tested as continuous LP relaxations, with approximately 47,000–1,093,000 variables | [Per-instance results](../docs/benchmarks/railway_en.md) |
| [Small examples](simple_examples/README_en.md) | 9 functional and regression test models | Excluded from benchmarks |

Netlib, PDS, and Railway store only raw inputs and files needed for conversion. Each dataset's `prepare.py` generates MPS files on demand. `simple_examples` stores small MPS files and matching C examples directly; see the [usage guide](../docs/usage_en.md) for modeling and solver API examples.

The command-line driver uses production numerical defaults and ignores numerical `LP_SIMPLEX_*` experiment switches. Add `--experimental-environment` explicitly to reproduce an environment-based experiment with the single-model driver; default benchmarks use production configuration. See [runtime configuration](../docs/runtime-configuration_en.md).

## Timing scopes

`benchmark.py` retains `wall_seconds`: CLI startup, reading, solving and point export, plus adapter log/point parsing, excluding adapter initialization and independent original-model auditing. CLI runs also record `cli_wall_seconds` (launch through exit, excluding later parsing), `cli_user_seconds`, `cli_system_seconds` and their sum `cli_cpu_seconds`. CPU is the isolated worker's child `getrusage` delta for the CLI and its descendants. Unsupported or invalid counters are `null`, never zero. The HiGHS Python binding has no CLI subprocess and does not provide these fields.

`compare_performance.py` keeps its legacy wall summary and adds separate summaries under `cli_timing`. Ratios require complete passing repetitions and valid measurements on both sides. Numerical PASS and timing availability are separate. Driver `time:` measures solve plus one original-model audit in process CPU time; it cannot replace full CLI timing. Hard deadlines and success criteria are unchanged. A killed child without a complete measurement has unknown CPU. Per-instance measurements, audit summaries, and reuse fingerprints are stored in the Markdown report; temporary logs and solution files may be cleaned after auditing.

## Prepare the environment

Build the single-model driver with CMake, a C compiler, and BLAS/LAPACK:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DLP_SIMPLEX_BUILD_BENCHMARKS=ON

cmake --build build --target lp_solve_mps -j 4
```

The executable is `build/data/lp_solve_mps`. Running it separately requires a prepared MPS file. Generating models from raw inputs also requires Python 3.9+; Python temporarily builds the PDS converter with a C compiler.

SCIP, Clp, and GLPK executables are required only when selected; provide them on PATH or through `--scip`, `--clp`, and `--glpk`. The Python environment requires NumPy, SciPy, and highspy or a compatible SciPy HiGHS binding. On macOS, install Clp and GLPK with `brew install clp glpk`. For example:

```sh
python3 -m pip install numpy scipy highspy
```

## Run full benchmarks

PDS and Railway run in ascending instance-number order, progressing from smaller to larger models. For PDS and Railway only, after the first GLPK timeout, all remaining GLPK runs are marked `SKIPPED`; other solvers continue. Netlib (including infeasible models) never skips subsequent instances. Skips count as neither timeouts nor failures and are excluded from timing statistics.

Run every benchmark instance in each dataset:

```sh
python3 data/benchmark.py --dataset netlib
python3 data/benchmark.py --dataset pds
python3 data/benchmark.py --dataset railway
```

By default, only lp_simplex is selected, with one measurement per instance. Matching measurements on the same machine are reused; only missing or changed entries run. The per-run limit is 60 seconds for Netlib and 120 seconds for PDS and Railway. Netlib uses only its 114 original feasible models; `simple_examples` is excluded from performance benchmarks. For a custom build directory, use `--binary PATH` to locate `lp_solve_mps`.

Instances with complete matching measurements need neither MPS generation nor solver launches. Otherwise, the script:

1. Generates a temporary MPS file from raw inputs.
2. Reuses matching measurements and runs only missing or changed entries serially, recording timings and statuses and auditing solutions against the original model.
3. Deletes the instance's temporary MPS file.
4. Saves run records and updates the dataset reports and both language summaries before proceeding to the next instance.

No dataset-wide MPS preparation is needed. Temporary MPS files live inside the run directory; manually generated MPS files are outside this cleanup process. Dataset reports retain current per-instance measurements and recompute comparisons from compatible records.

For debugging, select a single instance with fewer repetitions and a shorter limit:

```sh
python3 data/benchmark.py --dataset netlib --cases afiro --repeats 1 --seconds 10
```

`--cases` takes instance names without directories or file extensions, such as `afiro` for Netlib, `pds-30` for PDS, or `rail507` for Railway. It limits only the current run's selection.

After an interruption, use a new `--run-dir` with the original algorithms, work budgets, and time limits. Completed solver/instance/repetition entries with matching fingerprints are reused from the report, so temporary run directories are unnecessary for reuse. Changes to the machine, solver binaries or shared solver library, inputs, settings, or timing adapter invalidate the corresponding measurements; changing lp_simplex does not invalidate unchanged competitors. Use `--rerun lp_simplex` to explicitly retest a selected solver. For reports without fingerprints, a one-time `--adopt-existing` imports competitor results after checking the machine, versions, dimensions, and settings; it never imports an old lp_simplex build. `--report-only` rebuilds reports without continuing solver runs.

## Common operations

Use the current PDS report settings and run only missing or changed lp_simplex entries:

```sh
python3 -B data/benchmark.py --dataset pds --binary build/data/lp_solve_mps \
  --lp-algorithm auto --seconds 300 --iterations 300000
```

Append `--rerun lp_simplex` for an intentional retest. To update formatting, translations, or summaries only:

```sh
python3 -B data/benchmark.py --refresh-reports
```

Measurements and fingerprints live in the dataset Markdown, so deleting temporary run directories does not prevent reuse. `--report-only DIR` reads a retained raw run; use `--refresh-reports` when that directory is unavailable.

## Check a single model

[`lp_solve_mps.c`](lp_solve_mps.c) solves one model and checks the expected outcome. For example, generate Afiro, then solve it and print its variables:

```sh
python3 data/netlib/prepare.py feasible/afiro

./build/data/lp_solve_mps --expected -464.7531429 --print-solution data/netlib/feasible/afiro.mps
```

Here `prepare.py` requires the `feasible/` category prefix to distinguish Netlib models with the same name in different categories. The generated MPS remains in the dataset directory. See the other dataset READMEs for their preparation commands.

The driver supports Netlib, PDS, Railway, and other MPS models without inferring expectations from model names or directories. Use `--expected VALUE` to supply a reference objective or `--infeasible` to expect infeasibility; otherwise, only a successful solve is required. Dataset preparation, reference selection, and independent audits are managed by `benchmark.py`.

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

### Solver benchmarks: benchmark.py

| Option | Default | Meaning |
| --- | --- | --- |
| `--dataset` | `netlib` | `netlib`, `netlib-infeasible`, `pds`, or `railway`; excludes simple_examples |
| `--cases` | Full selected dataset | One or more raw instance names |
| `--binary` | `build/data/lp_solve_mps` | lp_simplex single-model driver |
| `--scip` | `scip` | SCIP executable name or path |
| `--clp` | `clp` | Clp executable name or path |
| `--glpk` | `glpsol` | GLPK executable name or path |
| `--glpk-gate-case` | Off | PDS/Railway only: test this instance first; if GLPK does not pass, skip its remaining instances. Skips are recorded separately from timeouts and failures |
| `--solvers` | `lp_simplex` | Solvers allowed to run; only selected executables are required. Compatible saved comparisons remain included |
| `--repeats` | 1 | Repetitions per instance and solver |
| `--seconds` | Netlib feasible: 60; infeasible: 90; PDS/Railway: 120 | Per-run wall-clock deadline; keep it consistent with the report for reuse |
| `--iterations` | 300000 | lp_simplex work budget |
| `--lp-algorithm` | `auto` | `auto`, `dual-revised`, `pan-bda`, or `tableau` |
| `--highs-algorithm` | `choose` | `choose`, `simplex`, or `ipm` |
| `--scip-algorithm` | `default` | `default` or `dual` |
| `--run-dir` | `tmp/benchmarks/<UTC timestamp>` | New directory for raw run records |
| `--docs-dir` | `docs/benchmarks` | Report destination |
| `--rerun SOLVER ...` | Off | Explicitly retest selected solvers instead of reusing their measurements |
| `--adopt-existing` | Off | One-time import of compatible competitor reports without fingerprints; never imports an old lp_simplex build |
| `--refresh-reports` | Off | Refresh bilingual Markdown reports and summaries without raw logs or solver launches |
| `--report-only DIR` | Off | Rebuild a report from a retained raw run directory without solving |

Routine optimization permits only lp_simplex to run, once per case. Unchanged competitor binaries and settings reuse their measurements; updating lp_simplex does not rerun competitors. `--solvers` permits execution without forcing it; add `--rerun` for an explicit retest. Compatible complete measurements contribute to comparisons on the common solved cases. The HiGHS binding, NumPy, and SciPy remain required for model loading and independent audits even when HiGHS is not selected.

### Single-model checks: lp_solve_mps

```sh
./build/data/lp_solve_mps [options] MODEL.mps
```

| Option | Default | Meaning |
| --- | --- | --- |
| `--algorithm` | `auto` | `auto`, `dual-revised`, `pan-bda`, or `tableau` |
| `--criteria` | `dantzig` | Use `bland/dantzig` for manual tableau, `dantzig/normalized` for Pan; AUTO/dual select pricing automatically |
| `--iterations` | 300000 | Positive work budget overriding the selected algorithm's default |
| `--no-presolve` | Off | Disable presolve |
| `--tolerance` | `1e-7` | Relative reference-objective comparison only; does not change primal/dual tolerances |
| `--expected VALUE` | No reference-objective check | Specify a reference objective |
| `--infeasible` | Off | Explicitly expect infeasibility |
| `--print-solution` | Off | Print original variables |
| `--progress` | Off | Flush original-model bounds and gap to stderr |
| `--progress-interval S` | 5 CPU seconds | Positive report interval; enables progress |
| `--save-point PATH` | No export | Atomically checkpoint improved feasible points during the solve; not an optimality guarantee |
| `--help` | — | Show current options |

These options belong to the test tools. See the [parameters and results reference](../docs/parameters_en.md) for C API solver options and result fields. Current tool options are also available through `python3 data/benchmark.py --help` and `./build/data/lp_solve_mps --help`.

Run the infeasibility benchmark:

```sh
python3 data/benchmark.py --dataset netlib-infeasible --binary build/data/lp_solve_mps --solvers lp_simplex scip highs clp glpk --repeats 1 --seconds 90
```

## Script structure and regression checks

| File | Responsibility |
| --- | --- |
| [benchmark.py](benchmark.py) | CLI, model preparation, solver execution, original-model audits, and measurement reuse |
| [benchmark_reports.py](benchmark_reports.py) | Statistics, Markdown persistence, bilingual reports, and reproduction commands |
| [compare_performance.py](compare_performance.py) | Paired diagnostics alternating between two builds |

Timing, reuse, and solver/report lifecycle checks live in `test/test_benchmark_timing.py`, `test/test_benchmark_cache.py`, and `test/test_benchmark_solvers.py`, respectively. Run the complete tool regression suite:

```sh
python3 -B -m unittest discover -s test -p 'test_benchmark_*.py'
```
