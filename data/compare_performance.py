#!/usr/bin/env python3
"""Alternate two builds on identical inputs with the existing original-model audit."""
import argparse
import datetime as dt
import math
import os
from pathlib import Path
import platform
import statistics
import sys

import benchmark as bench


def build_fingerprint(binary):
    root = bench.build_root(binary)
    files = {binary, *root.rglob('liblp_simplex.*')}
    return {str(p): bench.sha(p) for p in sorted(files) if p.is_file()}


def summarize(results, models, repeats, metric="wall_seconds"):
    rows, ratios = [], []
    for model in models:
        row = dict(instance=model['name'])
        for variant in ('baseline', 'candidate'):
            runs = [r for r in results if r['instance'] == model['name']
                    and r['variant'] == variant]
            row[variant + '_passed'] = sum(r['status'] == 'PASS' for r in runs)
            row[variant + '_seconds'] = (
                statistics.median(r[metric] for r in runs)
                if len(runs) == repeats and all(r['status'] == 'PASS'
                    and isinstance(r.get(metric), (int, float))
                    and not isinstance(r[metric], bool)
                    and math.isfinite(r[metric]) and r[metric] >= 0 for r in runs)
                else None)
        if row['baseline_seconds'] and row['candidate_seconds']:
            row['speedup'] = row['baseline_seconds'] / row['candidate_seconds']
            ratios.append(row['speedup'])
        rows.append(row)
    return dict(metric=metric, cases=rows, common_passed=len(ratios), total=len(models),
                geometric_mean_speedup=math.exp(statistics.mean(map(math.log, ratios)))
                if ratios else None,
                all_passed=all(r['status'] == 'PASS' for r in results)
                and len(results) == len(models) * repeats * 2)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline', type=Path, required=True)
    parser.add_argument('--candidate', type=Path, required=True)
    parser.add_argument('--dataset', choices=bench.DATASETS, default='netlib')
    parser.add_argument('--cases', nargs='+')
    parser.add_argument('--repeats', type=int, default=3)
    parser.add_argument('--seconds', type=float, default=90.)
    parser.add_argument('--iterations', type=int, default=300000)
    parser.add_argument('--lp-algorithm', choices=bench.LP_ALGORITHMS, default='auto')
    parser.add_argument('--run-dir', type=Path, required=True)
    args = parser.parse_args()
    if args.repeats < 1 or args.iterations < 1 or not math.isfinite(args.seconds) or args.seconds <= 0:
        parser.error('repeats, iterations and seconds must be positive')
    binaries = dict(baseline=args.baseline.resolve(), candidate=args.candidate.resolve())
    if any(not p.is_file() for p in binaries.values()):
        parser.error('both binaries must exist')
    models = bench.select_models(parser, args.cases, args.dataset)
    run_dir = args.run_dir.resolve()
    if run_dir.exists():
        parser.error('run directory already exists')
    (run_dir / 'logs').mkdir(parents=True)
    env = {k: v for k, v in os.environ.items() if not k.startswith('LP_SIMPLEX_')}
    env.update(bench.THREADS, LC_ALL='C', PYTHONDONTWRITEBYTECODE='1')
    os.environ.update(bench.THREADS)
    bench.load_highs()  # Verify the independent MPS reader before running solvers.
    fingerprints = {v: build_fingerprint(p) for v, p in binaries.items()}
    metadata = dict(created_utc=dt.datetime.now(dt.timezone.utc).isoformat(),
                    invocation=sys.argv, platform=platform.platform(),
                    python=sys.version, thread_environment=bench.THREADS,
                    builds=fingerprints, source_sha256=bench.sources(),
                    models=models, completed=False,
                    timing='legacy wall_seconds includes CLI lifetime plus log/point parsing; excludes adapter initialization and independent audit',
                    cli_timing='cli_wall_seconds: process launch through exit; child user/system CPU from getrusage when available; includes read, solve and point export; excludes adapter parsing/audit; unavailable remains null',
                    audit='original-model primal/objective; reference objective when available; no independent dual proof')
    bench.dump(run_dir / 'metadata.json', metadata)
    results = []
    try:
        for index, model in enumerate(models):
            with bench.prepared_model(model, run_dir):
                for repeat in range(1, args.repeats + 1):
                    order = ('baseline', 'candidate')
                    if (index + repeat) % 2 == 0:
                        order = order[::-1]
                    for variant in order:
                        variant_dir = run_dir / variant
                        (variant_dir / 'logs').mkdir(parents=True, exist_ok=True)
                        config = dict(binary=str(binaries[variant]), models=[model],
                                      seconds=args.seconds, tolerance=1e-7,
                                      iterations=args.iterations, lp_algorithm=args.lp_algorithm)
                        config_path = variant_dir / 'config.json'
                        bench.dump(config_path, config)
                        print(variant, end=' ', flush=True)
                        result = bench.run_one(config_path, config, 'lp_simplex', model, repeat, env)
                        result['variant'] = variant
                        results.append(result)
                        bench.dump(run_dir / 'results.json', results)
                if bench.sha(Path(model['path'])) != model['sha256']:
                    raise RuntimeError('MPS changed during comparison')
        metadata['completed'] = True
    finally:
        metadata['builds_unchanged'] = all(build_fingerprint(p) == fingerprints[v]
                                           for v, p in binaries.items())
        metadata['source_unchanged'] = bench.sources() == metadata['source_sha256']
        bench.dump(run_dir / 'metadata.json', metadata)
        summary = summarize(results, models, args.repeats)
        summary['cli_timing'] = {metric: summarize(results, models, args.repeats, metric)
                                 for metric in ('cli_wall_seconds', 'cli_cpu_seconds',
                                                'cli_user_seconds', 'cli_system_seconds')}
        summary['valid'] = metadata['completed'] and metadata['builds_unchanged'] and metadata['source_unchanged']
        bench.dump(run_dir / 'summary.json', summary)
    for row in summary['cases']:
        print(row)
    print('Geometric mean speedup:', summary['geometric_mean_speedup'])
    return 0 if summary['valid'] and summary['all_passed'] else 1


if __name__ == '__main__':
    sys.exit(main())
