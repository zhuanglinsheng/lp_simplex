#!/usr/bin/env python3
"""Bounded, reproducible diagnostics for primal/dual degeneracy (stdlib only).

Timeouts are recorded as incomplete, never as an infeasibility/optimality result.
Use separate build directories for comparisons; do not rebuild a running binary.
"""
import argparse
import concurrent.futures
import json
import os
from pathlib import Path
import re
import subprocess
import time


def profile_fields(text):
    result = {}
    for label, pattern in {
        'iterations': r'^iterations:\s+(\d+)',
        'cpu_seconds': r'^time:\s+([\d.]+)',
        'objective': r'^objective:\s+([^\n]+)',
        'status': r'^result:\s+([^\n]+)',
    }.items():
        match = re.search(pattern, text, re.M)
        if match:
            result[label] = match.group(1)
    for name in ('main progress=', 'main summary ', 'cleanup summary ',
                 'primal cleanup progress=', 'singleton crash rank=',
                 'recovery stage=', 'original bounds restored merit='):
        lines = [s for s in text.splitlines() if name in s]
        if lines:
            result[name.strip(' =')] = lines[-1]
    result['validated'] = '[PASS] prediction reached' in text
    return result


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--cases', nargs='+', default=[
        'ken-18', 'pds-20', 'dfl001', 'qap12', 'qap15'])
    parser.add_argument('--iterations', type=int, default=20000)
    parser.add_argument('--seconds', type=float, default=120)
    parser.add_argument('--jobs', type=int, default=1)
    parser.add_argument('--algorithm', choices=['dual-revised', 'pan-bda'],
                        default='dual-revised')
    parser.add_argument('--devex', action='store_true')
    parser.add_argument('--pan-singleton-crash', action='store_true')
    parser.add_argument('--initial-cost-perturb', action='store_true')
    parser.add_argument('--stable-klu', action='store_true')
    bounds = parser.add_mutually_exclusive_group()
    bounds.add_argument('--no-bound-perturb', action='store_true')
    bounds.add_argument('--bound-perturb', action='store_true')
    bounds.add_argument('--force-bound-perturb', action='store_true')
    args = parser.parse_args()
    if args.iterations <= 0 or args.seconds <= 0 or args.jobs <= 0:
        parser.error('iterations, seconds and jobs must be positive')
    args.output.mkdir(parents=True, exist_ok=True)
    binary = args.binary.resolve()
    env = dict(os.environ)
    # Prevent inherited experimental switches from contaminating a baseline.
    for key in ('LP_SIMPLEX_CLEANUP_DEVEX', 'LP_SIMPLEX_PAN_SINGLETON_CRASH',
                'LP_SIMPLEX_DISABLE_PAN', 'LP_SIMPLEX_SINGLE_PERTURB_EPOCH',
                'LP_SIMPLEX_INITIAL_COST_PERTURB', 'LP_SIMPLEX_STABLE_KLU',
                'LP_SIMPLEX_DISABLE_BOUND_PERTURB', 'LP_SIMPLEX_FORCE_BOUND_PERTURB',
                'LP_SIMPLEX_PRIMAL_BOUND_PERTURB'):
        env.pop(key, None)
    env.update(LP_SIMPLEX_PROFILE='1', OMP_NUM_THREADS='1',
               OPENBLAS_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='1')
    if args.devex:
        env['LP_SIMPLEX_CLEANUP_DEVEX'] = '1'
    if args.pan_singleton_crash:
        env['LP_SIMPLEX_PAN_SINGLETON_CRASH'] = '1'
    if args.initial_cost_perturb:
        env['LP_SIMPLEX_INITIAL_COST_PERTURB'] = '1'
    if args.stable_klu:
        env['LP_SIMPLEX_STABLE_KLU'] = '1'
    if args.no_bound_perturb:
        env['LP_SIMPLEX_DISABLE_BOUND_PERTURB'] = '1'
    if args.bound_perturb:
        env['LP_SIMPLEX_PRIMAL_BOUND_PERTURB'] = '1'
    if args.force_bound_perturb:
        env['LP_SIMPLEX_FORCE_BOUND_PERTURB'] = '1'

    def run(case):
        model = root / 'data/netlib/feasible' / (case + '.mps')
        command = [str(binary), '--iterations', str(args.iterations),
                   '--algorithm', args.algorithm, str(model)]
        log = args.output / (case + '.log')
        start = time.monotonic()
        timeout = False
        with log.open('w') as stream:
            process = subprocess.Popen(command, stdout=stream, stderr=stream, env=env)
            try:
                returncode = process.wait(timeout=args.seconds)
            except subprocess.TimeoutExpired:
                timeout = True
                process.terminate()
                try:
                    returncode = process.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    process.kill()
                    returncode = process.wait()
            except BaseException:
                process.kill()
                process.wait()
                raise
        result = dict(case=case, command=command, timeout=timeout,
                      returncode=returncode, wall_seconds=time.monotonic() - start,
                      switches={k: v for k, v in env.items()
                                if k.startswith('LP_SIMPLEX_')},
                      **profile_fields(log.read_text()))
        if timeout:
            result['validated'] = False
        (args.output / (case + '.json')).write_text(json.dumps(result, indent=2))
        print(json.dumps(result, ensure_ascii=False), flush=True)
        return result

    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as executor:
        results = list(executor.map(run, args.cases))
    (args.output / 'runs.json').write_text(json.dumps(results, indent=2))


if __name__ == '__main__':
    main()
