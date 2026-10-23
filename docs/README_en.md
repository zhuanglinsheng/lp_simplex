# Project Documentation

[简体中文](./README.md) | English | [Project home](../README_en.md) | [Algorithm overview](./algorithm/README_en.md)

This directory provides entry points to algorithm theory, solver result interpretation, and performance tests. See the [project home](../README_en.md) for building, installation, and basic usage.

## Using the solver

The [usage guide](./usage_en.md) starts from `data/simple_examples/production_planning.mps` and follows the mathematical model through manual C modeling, MPS input, compilation, and standalone examples, with links to other complete programs.

## Algorithm theory

The  [Algorithm collection](./algorithm/README_en.md) presents mathematical models, assumptions, derivations, and optimality conditions. It covers presolve, two-phase primal simplex, dual revised simplex, anti-degeneracy, NNLS, Pan/BDA, tangent-cone recovery, PDHG, and basis linear algebra. Every article has Chinese and English versions, with language links below the title.

Start with the notation and primal–dual conditions in the overview, then read the chapters relevant to your interests.

## Parameters and results

The [parameters and results reference](./parameters_en.md) collects model fields, bounds, algorithms, pricing, budgets, tolerances, all status codes, and the unified result contract, plus CMake, installation, and troubleshooting reference.

## Performance tests

The [benchmark summary](./benchmarks/README_en.md) provides test information, comparison statistics, and reproduction instructions for lp_simplex, SCIP/SoPlex, and HiGHS. Per-instance results are available for:

- [Netlib feasible](./benchmarks/netlib_en.md)
- [PDS](./benchmarks/pds_en.md)
- [Railway continuous LP relaxations](./benchmarks/railway_en.md)

See [data and benchmark instructions](../data/README_en.md) for data preparation, parameters, and audit criteria. The benchmark summary is available in Chinese and English. Data instructions are also available in both languages. Per-run configurations, raw results, logs, and solution vectors are stored locally under `tmp/benchmarks/`, which Git ignores.
