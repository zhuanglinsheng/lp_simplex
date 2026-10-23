# NETLIB

114 feasible and 29 infeasible LP instances, plus Gurobi 10 presolved variants
and the reference objective values in `feasible_gurobi_1e-8.csv`.

Only compressed data is stored. `raw/{feasible,infeasible,netlib_grbp}/*.mps.gz`
are lossless gzip archives of this repository's existing MPS files, not newly
downloaded official Netlib MPC files. All model variants are preserved.

Generate selected models from the repository root (Python 3.9+):

```sh
python3 data/netlib/prepare.py feasible/afiro feasible/degen2
python3 data/netlib/prepare.py --all
```

Names include the category because feasible and infeasible models can share a
name. Outputs restore the original paths, e.g. `feasible/afiro.mps`, and are
ignored by Git. `--output-dir DIR` changes the output root. Existing selected
outputs are regenerated atomically; invalid input leaves them untouched.
`SHA256SUMS` records the compressed input checksums.

Original collection: <https://www.netlib.org/lp/>.
