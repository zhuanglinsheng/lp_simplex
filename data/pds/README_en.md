# PDS Dataset

[简体中文](./README.md) | English | [Project home](../../README_en.md) | [Data overview](../README_en.md)

The collection contains 9 LP instances from Mittelmann's PDS directory. PDS stands for Patient Distribution System, a casualty distribution and evacuation problem. The `pds-20` model matches its Netlib counterpart; the smaller `pds-02/06/10` models are in Netlib.

PDS tests and report rows use ascending numeric order, from `pds-20` through `pds-100`, with `pds-100` executed and displayed last.

## Instance sizes

| Instance | Constraints | Variables | Constraint nonzeros |
|---|---:|---:|---:|
| pds-20 | 33,874 | 105,728 | 230,200 |
| pds-30 | 49,944 | 154,998 | 337,144 |
| pds-40 | 66,844 | 212,859 | 462,128 |
| pds-50 | 83,060 | 270,095 | 585,114 |
| pds-60 | 99,431 | 329,643 | 712,779 |
| pds-70 | 114,944 | 382,311 | 825,771 |
| pds-80 | 129,181 | 426,278 | 919,524 |
| pds-90 | 142,823 | 466,671 | 1,005,359 |
| pds-100 | 156,243 | 505,360 | 1,086,785 |

Nonzeros count only constraint-matrix entries, excluding objective coefficients and explicit zeros. See `manifest.csv` and `SHA256SUMS` for sources, sizes, and checksums.

## Prepare models

Python 3.9+ and a C compiler are required. Run from the repository root:

```sh
python3 data/pds/prepare.py pds-30
python3 data/pds/prepare.py --all
```

Inputs are bzip2-compressed Netlib MPC files. The script works offline: it decompresses inputs, temporarily compiles `raw/emps.c` using `CC` (default `cc`), and converts them to MPS. Output defaults to `data/pds/mps/`; override it with `--output-dir DIR`. Intermediate files are cleaned automatically, and failures leave existing outputs unchanged.

## Sources and results

Public reference objectives for all 9 instances are stored in [reference_objectives.csv](reference_objectives.csv), with the instance name, source-reported status, objective value (`pobj`), and source URL. Values come from the entries where both solvers report optimality in the [COIN-OR CPLEX11-Barrier versus CLP comparison](https://www.coin-or.org/GAMSlinks/benchmarks/LP/allSolver_080107/2-CPLEX-1.trc_4-CBC-0.trc_sqr.htm); both solvers report the same objective values. These are numerical reference results, not exact mathematical certificates. The benchmark script uses this file to check reference objectives for all 9 PDS instances and records its SHA-256. New configurations are loaded only when a new run starts.

- [Mittelmann PDS directory](https://plato.asu.edu/ftp/lptestset/pds/) and [official format notes](https://plato.asu.edu/ftp/lptestset/pds/00README).
- [Netlib MPC-to-MPS converter](https://www.netlib.org/lp/data/emps.c), stored locally as `raw/emps.c`.
- [PDS benchmark report](../../docs/benchmarks/pds_en.md) and [shared benchmark entry point](../README_en.md).
