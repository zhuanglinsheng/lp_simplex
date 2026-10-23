# Railway Dataset

[简体中文](./README.md) | English | [Project home](../../README_en.md) | [Data overview](../README_en.md)

The collection contains all 7 OR-Library railway set-covering instances, originating from an Italian railway application. This project solves their continuous LP relaxations, with variable bounds `[0, 1]`.

Benchmarks and report rows use ascending numeric order: `rail507`, `rail516`, `rail582`, `rail2536`, `rail2586`, `rail4284`, `rail4872`.

## Instance sizes

| Instance | Constraints | Variables | Matrix nonzeros |
|---|---:|---:|---:|
| rail507 | 507 | 63,009 | 409,349 |
| rail516 | 516 | 47,311 | 314,896 |
| rail582 | 582 | 55,515 | 401,708 |
| rail2536 | 2,536 | 1,081,841 | 10,993,311 |
| rail2586 | 2,586 | 920,683 | 8,008,776 |
| rail4284 | 4,284 | 1,092,610 | 11,279,748 |
| rail4872 | 4,872 | 968,672 | 9,244,093 |

Million-scale refers to the number of variables. The number in each instance name is its constraint count. Nonzeros exclude objective coefficients. See `manifest.csv` and `SHA256SUMS` for sources, file sizes, and checksums.

## Prepare models

Python 3.9+ is required; no third-party packages are needed. Run from the repository root:

```sh
python3 data/railway/prepare.py rail507
python3 data/railway/prepare.py --all
```

The official compressed inputs, `raw/rail*.gz`, contain a Railway-specific column format, not MPS. The script reads column costs and covered rows, then generates fixed-column MPS continuous models. Output defaults to `data/railway/mps/`; override it with `--output-dir DIR`.

Large instances produce large MPS files, so prepare them as needed. Compressed inputs retain their original bytes, and generated files are ignored by Git. Failed conversions leave existing outputs unchanged.

## Sources and results

- [OR-Library dataset notes](https://people.brunel.ac.uk/~mastjjb/jeb/orlib/scpinfo.html) and [file directory](https://people.brunel.ac.uk/~mastjjb/jeb/orlib/files/).
- `raw/scpinfo.html` preserves the official format description. `manifest.csv` records each instance's source and parsing results.
- [Railway benchmark report](../../docs/benchmarks/railway_en.md) and [shared benchmark entry point](../README_en.md). All five solvers use the same continuous LP relaxation.
