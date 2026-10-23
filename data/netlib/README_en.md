# Netlib Dataset

[简体中文](./README.md) | English | [Project home](../../README_en.md) | [Data overview](../README_en.md)

The collection contains 114 feasible LPs, 29 infeasible LPs, and Gurobi 10 presolved variants. Benchmarks use only the original feasible models. Reference objectives are stored in [feasible_gurobi_1e-8.csv](feasible_gurobi_1e-8.csv).

## Prepare models

Python 3.9+ is required. Run from the repository root:

```sh
python3 data/netlib/prepare.py feasible/afiro feasible/degen2
python3 data/netlib/prepare.py --all
```

Names include a category because different categories may contain the same name. Files are restored to their original paths, such as `data/netlib/feasible/afiro.mps`. Use `--output-dir DIR` to choose another output directory.

[`test_netlib.c`](test_netlib.c) in this directory solves a single model and checks its reference result. After building, run:

```sh
./build/data/netlib/test_netlib data/netlib/feasible/afiro.mps
```

[`../benchmark.py`](../benchmark.py) coordinates batch comparisons of the five solvers.

## Stored data

- `raw/{feasible,infeasible,netlib_grbp}/*.mps.gz`: lossless gzip archives of the repository's existing MPS files.
- `SHA256SUMS`: compressed-input checksums.
- Generated MPS files are ignored by Git. Existing outputs are replaced only after successful generation.

See the [original Netlib collection](https://www.netlib.org/lp/) and the [Netlib benchmark report](../../docs/benchmarks/netlib_en.md).
