# Production configuration and experimental environment controls

[简体中文](./runtime-configuration.md) | English | [Documentation](./README_en.md) | [Parameters/results](./parameters_en.md)

`lp_simplex_solve` is the production entry. All four public algorithms use public Options and deterministic numerical defaults, ignoring numerical `LP_SIMPLEX_*` environment controls. Default AUTO needs no environment setup and cannot inherit presolve, scaling, perturbation, pricing or backend preferences from an earlier experiment. `LP_SIMPLEX_PROFILE` remains diagnostic. Thread controls belonging to third-party BLAS/runtime libraries are outside this boundary.

## Calls and migration

Initialize Options with `lp_simplex_default_options`, set supported fields, then call `lp_simplex_solve`. To reproduce environment-based experiments, explicitly call `lp_simplex_solve_experimental` with the same arguments. It resolves numerical environment controls at the start and uses a solve-local snapshot thereafter. The solver never changes process environment; callback changes cannot alter the current numerical configuration. Callers must synchronize concurrent process-environment mutation; numerical snapshots do not make environment writes thread-safe. Experimental switches are not stable public parameters or new algorithms.

The command-line driver also defaults to production. Old environment experiments require an explicit switch:

```sh
LP_SIMPLEX_DISABLE_SCALING=1 build/data/lp_solve_mps \
  --experimental-environment --algorithm dual-revised --no-presolve MODEL.mps
```

Driver output and PROFILE root metadata identify the mode. No Options, Progress or Result fields were added; existing structure layouts are unchanged. Programs previously relying on numerical environment controls through the normal entry must migrate to the experimental entry. The ordinary Dual preset uses bounded initial cost separation and at most one later face perturbation; the explicit compatibility entry retains the historical empty-environment experiment preset.

The new exact equality-dependency reduction currently requires `LP_SIMPLEX_PRESOLVE_DEPENDENCIES=1` and the explicit experimental entry. Production defaults and the plain private `simplex_presolve_run` leave it disabled; configured presolve receives only the resolved boolean. Limited-pivot savings have not established complete-solve gains, and a same-limit regression changed qap12 from passing to timing out. This admission decision applies to every model without instance-name exceptions. See the [qap15 analysis](./development/qap15-presolve-analysis-2026-10-28.md).

## Ownership

The public root creates `simplex_SolveContext`, resolves `simplex_RunConfig` and fixes the configuration. Standalone private AUTO calls install production defaults when no configuration exists. Presolve, scaling and child models inherit the root snapshot. Dual copies normalized `simplex_DualControls` into its iteration state; hot loops read fields. Configured LU creation receives a concrete KLU preference, never an AUTO/experimental label.

Manual BDA entries use the root snapshot; AUTO's complete BDA routes retain deterministic adapter-selected configurations. Original-form NNLS initialization receives resolved numerical controls and explicitly selects standard form, avoiding an indirect environment read through legacy basis construction. Failed-start comparison uses controls actually consumed by Dual; PROFILE and unrelated PAN settings cannot change its eligibility.

Experimental AUTO may enable presolve/scaling and Dual controls while keeping its complete BDA route factories isolated from PAN environment choices. Original-form NNLS uses captured numerical settings while its owner selects the representation.

Unconfigured lower-level private compatibility entries retain legacy experiment behavior for internal research/tests. They are not the production public API. New composed workflows must bind a shared configured context or call explicit configured numerical interfaces. A valid point, basis or model identity does not authorize reloading environment.

## Legacy experiment semantics

Presence flags retain their historical semantics: the string `0` enables them, in the explicit experimental entry only. Pan integer work/epoch/block settings retain bounded parsing; invalid/negative values disable, overflow saturates and block sizes are capped. Weights must be finite and valid. Production does not parse these numerical environment controls.

Both entries retain the same validation, cumulative budgets, original-model audits and result contracts; experiments cannot lower certification standards. A snapshot is a configuration source, not runtime state, a candidate, a journal or an algorithm route. See `simplex_run_config.c/.h`, `test_environment` and the [development record](./development/README.md) for implementation, evidence and limits.


## Persistent Dual working objective

Ordinary solves now enable the existing bounded initial cost separation and disable repeated face-perturbation epochs. The separation applies to structural nonbasic variables at finite endpoints; its magnitude is `5 * dual_tolerance * (1 + abs(original_cost)) * (1 + hash_fraction)`, followed by the existing outward rounding. The factor 5 is the existing experimental policy, not fitted to a dataset. It is a working-objective device, not a relaxation of final feasibility or optimality tolerances. The original objective is restored before final certification. Symbolic repetition detection, refinement, bound restoration, and cumulative work grants remain in force.

The explicitly requested legacy environment entry keeps its historical empty-environment values: initial cost separation off and multiple epochs on. Set both `LP_SIMPLEX_INITIAL_COST_PERTURB` and `LP_SIMPLEX_SINGLE_PERTURB_EPOCH` through that entry to reproduce the ordinary Dual preset. Numerical kernels receive the resolved values and never branch on entry mode or dataset names. The two-column restoration regression deliberately reverses the working preference, then requires the original optimum and original reduced costs for both cost signs and basis backends.
