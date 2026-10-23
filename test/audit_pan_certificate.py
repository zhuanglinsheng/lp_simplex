#!/usr/bin/env python3
"""Independently audit a minimization LP primal/row-dual certificate.

Uses the external library solely to read the original MPS, never to optimize.
Multipliers follow c + A^T y; nonnegative bound multipliers are reconstructed.
"""
import argparse
import hashlib
import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'data'))
from benchmark import load_highs


def audit(model, primal_file, dual_file, primal_tolerance=1e-7, dual_tolerance=1e-7):
    core, factory, binding = load_highs()
    reader = factory()
    reader.setOptionValue('output_flag', False)
    if reader.readModel(str(model)) != core.HighsStatus.kOk:
        raise ValueError(f'Cannot read {model}')
    lp = reader.getLp()
    if 'max' in str(lp.sense_).lower():
        raise ValueError('This certificate convention requires a minimization model')
    starts = list(lp.a_matrix_.start_)
    indices = list(lp.a_matrix_.index_)
    values = list(lp.a_matrix_.value_)
    costs = list(lp.col_cost_)
    lower, upper = list(lp.col_lower_), list(lp.col_upper_)
    row_lower, row_upper = list(lp.row_lower_), list(lp.row_upper_)
    x = [float(t) for t in primal_file.read_text().split()]
    y = [float(t) for t in dual_file.read_text().split()]
    if len(x) != lp.num_col_ or len(y) != lp.num_row_ or not all(map(math.isfinite, x + y)):
        raise ValueError('Invalid certificate dimension or nonfinite coordinate')
    rows = [[] for _ in y]
    reduced = []
    for j, value in enumerate(x):
        reduced.append(math.fsum([costs[j]] + [values[k] * y[indices[k]]
                                             for k in range(starts[j], starts[j + 1])]))
        for k in range(starts[j], starts[j + 1]):
            rows[indices[k]].append(values[k] * value)
    activity = [math.fsum(row) for row in rows]
    dual_terms, complementarity, stationarity, relative = [], [], [], []
    for j, g in enumerate(reduced):
        z_lower = max(g, 0.) if math.isfinite(lower[j]) else 0.
        z_upper = max(-g, 0.) if math.isfinite(upper[j]) else 0.
        if math.isfinite(lower[j]):
            dual_terms.append(lower[j] * z_lower)
            complementarity.append(abs((x[j] - lower[j]) * z_lower))
        if math.isfinite(upper[j]):
            dual_terms.append(-upper[j] * z_upper)
            complementarity.append(abs((upper[j] - x[j]) * z_upper))
        error = abs(g - z_lower + z_upper)
        stationarity.append(error)
        relative.append(error / (1. + abs(costs[j])))
    for i, value in enumerate(y):
        z_upper = max(value, 0.) if math.isfinite(row_upper[i]) else 0.
        z_lower = max(-value, 0.) if math.isfinite(row_lower[i]) else 0.
        if math.isfinite(row_upper[i]):
            dual_terms.append(-row_upper[i] * z_upper)
            if row_upper[i] != row_lower[i]:
                complementarity.append(abs((row_upper[i] - activity[i]) * z_upper))
        if math.isfinite(row_lower[i]):
            dual_terms.append(row_lower[i] * z_lower)
            if row_upper[i] != row_lower[i]:
                complementarity.append(abs((activity[i] - row_lower[i]) * z_lower))
        error = abs(value - z_upper + z_lower)
        stationarity.append(error)
        relative.append(error)
    primal = max([0.] + [max(lo - v, v - hi) for lo, v, hi in zip(lower, x, upper)] +
                 [max(lo - v, v - hi) for lo, v, hi in zip(row_lower, activity, row_upper)])
    objective = math.fsum(c * v for c, v in zip(costs, x)) + lp.offset_
    dual_objective = math.fsum(dual_terms) + lp.offset_
    result = dict(model=str(model), primal_file=str(primal_file), dual_file=str(dual_file),
                  reader=binding, optimizer_called=False, presolve_in_audit=False,
                  objective=objective, dual_objective=dual_objective, primal=primal,
                  stationarity=max([0.] + stationarity), relative_stationarity=max([0.] + relative),
                  complementarity_l1=math.fsum(complementarity), gap=abs(objective - dual_objective),
                  gap_limit=dual_tolerance * (1. + abs(objective)),
                  primal_tolerance=primal_tolerance, dual_tolerance=dual_tolerance,
                  sha256={str(p): hashlib.sha256(p.read_bytes()).hexdigest()
                          for p in (model, primal_file, dual_file)})
    result['verified'] = (primal <= primal_tolerance and
                          result['relative_stationarity'] <= dual_tolerance and
                          result['gap'] <= result['gap_limit'] and
                          result['complementarity_l1'] <= result['gap_limit'])
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('model', type=Path)
    parser.add_argument('primal', type=Path)
    parser.add_argument('dual', type=Path)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = audit(args.model, args.primal, args.dual)
    text = json.dumps(result, indent=2) + '\n'
    if args.output:
        args.output.write_text(text)
    print(text, end='')
    raise SystemExit(0 if result['verified'] else 1)
