#!/usr/bin/env python3
"""Solve MPS files with reference HiGHS algorithms and audit original-model residuals.

Uses highspy when installed, otherwise SciPy's bundled private HiGHS binding
(tested with SciPy's HiGHS 1.12.0). The latter is a diagnostic compatibility
fallback, not a new dependency or solver backend of the C library.
"""
import argparse
import csv
import json
import math
from pathlib import Path
import time


def load_highs():
    try:
        import highspy as core
        return core, core.Highs, 'highspy'
    except ImportError:
        try:
            from scipy.optimize._highspy import _core as core
            return core, core._Highs, 'scipy-private-binding'
        except ImportError as exc:
            raise SystemExit('Install highspy, or use SciPy with its bundled HiGHS binding.') from exc


def json_safe(value):
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, dict):
        return {key: json_safe(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [json_safe(item) for item in value]
    return value


def solve(model, args, core, factory, binding):
    import numpy as np
    from scipy.sparse import csc_matrix
    highs = factory()
    options = dict(solver='ipm' if args.algorithm == 'ipm' else 'simplex',
                   threads=1, time_limit=args.seconds, presolve=args.presolve,
                   primal_feasibility_tolerance=args.tolerance,
                   dual_feasibility_tolerance=args.tolerance, random_seed=0,
                   run_crossover='on')
    if args.algorithm != 'ipm':
        options['simplex_strategy'] = 1 if args.algorithm == 'dual' else 4
    if args.disable_perturbation:
        options.update(dual_simplex_cost_perturbation_multiplier=0.,
                       primal_simplex_bound_perturbation_multiplier=0.)
    for name, value in options.items():
        if highs.setOptionValue(name, value) != core.HighsStatus.kOk:
            raise RuntimeError('Unsupported HiGHS option: ' + name)
    log_path = args.output / (model.stem + '.log')
    if highs.setOptionValue('log_file', str(log_path.resolve())) != core.HighsStatus.kOk:
        raise RuntimeError('Could not configure HiGHS log file')
    if highs.readModel(str(model)) != core.HighsStatus.kOk:
        raise RuntimeError('Could not read model: ' + str(model))
    start = time.monotonic()
    run_status = highs.run()
    elapsed = time.monotonic() - start
    info, lp, solution = highs.getInfo(), highs.getLp(), highs.getSolution()
    result = dict(model=str(model.resolve()), algorithm=args.algorithm,
                  binding=binding, version=highs.version(), options=options,
                  run_status=str(run_status), status=str(highs.getModelStatus()),
                  wall_seconds=elapsed, objective=info.objective_function_value,
                  simplex_iterations=info.simplex_iteration_count,
                  ipm_iterations=info.ipm_iteration_count,
                  crossover_iterations=info.crossover_iteration_count,
                  max_primal=info.max_primal_infeasibility,
                  max_dual=info.max_dual_infeasibility,
                  independent_primal=None, stationarity=None,
                  relative_stationarity=None, validated=False)
    if solution.value_valid:
        x = np.asarray(solution.col_value)
        # readModel exposes the constraint matrix in column form.
        matrix = csc_matrix((np.asarray(lp.a_matrix_.value_),
                             np.asarray(lp.a_matrix_.index_),
                             np.asarray(lp.a_matrix_.start_)),
                            shape=(lp.num_row_, lp.num_col_))
        activity = matrix @ x
        vectors = [x, activity]
        violations = [np.asarray(lp.col_lower_) - x,
                      x - np.asarray(lp.col_upper_),
                      np.asarray(lp.row_lower_) - activity,
                      activity - np.asarray(lp.row_upper_)]
        primal = max([0.] + [float(np.max(v, initial=0.)) for v in violations])
        if not all(np.all(np.isfinite(v)) for v in vectors):
            primal = float('inf')
        result['independent_primal'] = primal
        if solution.dual_valid:
            y, rc = np.asarray(solution.row_dual), np.asarray(solution.col_dual)
            priced = matrix.T @ y
            costs = np.asarray(lp.col_cost_)
            residual = costs - priced - rc
            result['stationarity'] = float(np.max(np.abs(residual), initial=0.))
            scale = 1. + np.abs(costs) + np.abs(priced) + np.abs(rc)
            result['relative_stationarity'] = float(np.max(np.abs(residual) / scale, initial=0.))
        computed_objective = float(np.sum(np.asarray(lp.col_cost_, dtype=np.longdouble) * x) + lp.offset_)
        result['independent_objective'] = computed_objective
        reference = Path(__file__).resolve().parents[1] / 'data/netlib/feasible_gurobi_1e-8.csv'
        if reference.exists():
            with reference.open() as stream:
                for row in csv.DictReader(stream):
                    if row['name'] == model.stem and row['status'] == 'OPTIMAL':
                        expected = float(row['pobj'])
                        result['reference_objective'] = expected
                        result['reference_matches'] = abs(info.objective_function_value - expected) <= args.tolerance * max(1., abs(expected))
                        break
        result['validated'] = bool(
            run_status == core.HighsStatus.kOk and
            highs.getModelStatus() == core.HighsModelStatus.kOptimal and
            math.isfinite(primal) and primal <= args.tolerance and
            math.isfinite(info.max_dual_infeasibility) and
            info.max_dual_infeasibility <= args.tolerance and
            result['stationarity'] is not None and
            math.isfinite(result['stationarity']) and
            result['relative_stationarity'] <= args.tolerance and
            math.isfinite(computed_objective) and
            abs(computed_objective - info.objective_function_value) <= args.tolerance * (1. + abs(computed_objective)) and
            result.get('reference_matches', True))
        if result['validated']:
            np.savez_compressed(args.output / (model.stem + '-solution.npz'),
                                primal=x, row_dual=solution.row_dual,
                                column_dual=solution.col_dual)
    result = json_safe(result)
    (args.output / (model.stem + '.json')).write_text(json.dumps(result, indent=2, allow_nan=False))
    print(json.dumps(result, allow_nan=False), flush=True)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('models', type=Path, nargs='+')
    parser.add_argument('--algorithm', choices=['ipm', 'dual', 'primal'], default='ipm')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--seconds', type=float, default=90.)
    parser.add_argument('--tolerance', type=float, default=1e-7)
    parser.add_argument('--presolve', choices=['on', 'off'], default='on')
    parser.add_argument('--disable-perturbation', action='store_true')
    args = parser.parse_args()
    if not math.isfinite(args.seconds) or args.seconds <= 0 or not 1e-10 <= args.tolerance <= 1e-3:
        parser.error('seconds must be finite and positive; tolerance must be in [1e-10, 1e-3]')
    args.output.mkdir(parents=True, exist_ok=True)
    core, factory, binding = load_highs()
    results = [solve(model, args, core, factory, binding) for model in args.models]
    (args.output / 'runs.json').write_text(json.dumps(results, indent=2, allow_nan=False))
    return 0 if all(item['validated'] for item in results) else 1


if __name__ == '__main__':
    raise SystemExit(main())
