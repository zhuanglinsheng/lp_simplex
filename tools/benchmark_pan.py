#!/usr/bin/env python3
"""Bounded, reproducible Pan cold/own-Phase-I/reference-point experiments.

Reference points isolate Phase II only and never count as independent solves.
"""
import argparse
import concurrent.futures
import csv
import datetime
import hashlib
import json
import os
import platform
import re
from pathlib import Path
import subprocess
import time


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--cases', nargs='+', default=['dfl001', 'qap12', 'qap15'])
    parser.add_argument('--mode', choices=['nnls', 'dual-feasible', 'reference'], default='nnls')
    parser.add_argument('--points', type=Path)
    parser.add_argument('--affine', action='store_true')
    parser.add_argument('--complete', action='store_true')
    parser.add_argument('--stabilize', action='store_true')
    parser.add_argument('--compatible', action='store_true')
    parser.add_argument('--cone', action='store_true')
    parser.add_argument('--native-bounds', action='store_true')
    parser.add_argument('--native-adaptive', action='store_true')
    parser.add_argument('--native-predict', action='store_true')
    parser.add_argument('--native-pdhg', type=int, default=0)
    parser.add_argument('--pdhg-weight', type=float)
    parser.add_argument('--pdhg-restore', action='store_true')
    parser.add_argument('--cone-backend', choices=['explicit', 'implicit'], default='explicit')
    parser.add_argument('--no-presolve', action='store_true')
    parser.add_argument('--block-crash', action='store_true')
    parser.add_argument('--full-refinement', action='store_true')
    parser.add_argument('--dual-center', action='store_true')
    parser.add_argument('--cone-batch', action='store_true')
    parser.add_argument('--nnls-block', type=int, default=1, choices=range(1, 33))
    parser.add_argument('--face-guard', action='store_true')
    parser.add_argument('--dual-polish', type=int, default=0)
    parser.add_argument('--first-order', type=int, default=0)
    parser.add_argument('--newton-qr', action='store_true')
    parser.add_argument('--center-epoch', type=int, default=0)
    parser.add_argument('--implicit-bulk', action='store_true')
    parser.add_argument('--implicit-updates', action='store_true')
    parser.add_argument('--save-points', action='store_true')
    parser.add_argument('--pricing', choices=['dantzig', 'normalized'], default='dantzig')
    parser.add_argument('--iterations', type=int, default=20000)
    parser.add_argument('--seconds', type=float, default=30.)
    parser.add_argument('--jobs', type=int, default=1)
    args = parser.parse_args()
    if args.iterations <= 0 or args.seconds <= 0 or args.jobs <= 0:
        parser.error('Budgets must be positive')
    if args.mode == 'reference' and args.points is None:
        parser.error('--mode reference requires --points')
    if args.mode == 'reference' and args.pricing != 'dantzig':
        parser.error('The reference-point probe currently uses Dantzig pricing')
    if args.cone and (args.affine or args.complete or args.stabilize or args.compatible or args.pricing != 'dantzig'):
        parser.error('Cone direction NNLS replaces classical Phase II; classical pricing/recovery flags do not apply')
    if not args.cone and args.cone_backend != 'explicit':
        parser.error('--cone-backend requires --cone')
    if args.dual_center and not args.cone:
        parser.error('--dual-center requires --cone')
    if args.native_bounds and not args.cone:
        parser.error('--native-bounds requires --cone')
    if args.native_adaptive and not args.native_bounds:
        parser.error('--native-adaptive requires --native-bounds')
    if args.native_predict and not args.native_adaptive:
        parser.error('--native-predict requires --native-adaptive')
    if args.native_pdhg < 0 or (args.native_pdhg and not args.native_bounds):
        parser.error('--native-pdhg requires --native-bounds and a positive budget')
    if args.pdhg_weight is not None and (not args.native_pdhg or not 1e-12 <= args.pdhg_weight <= 1e12):
        parser.error('--pdhg-weight requires --native-pdhg and weight in [1e-12, 1e12]')
    if args.pdhg_restore and not args.native_pdhg:
        parser.error('--pdhg-restore requires --native-pdhg')
    if args.native_bounds and (args.first_order or args.newton_qr or args.center_epoch or
                               args.face_guard or args.cone_batch or args.dual_center):
        parser.error('Native bounds currently uses its own face/centering flow; standard-cone prediction, epoch and face flags do not apply')
    if args.cone_batch and not args.cone:
        parser.error('--cone-batch requires --cone')
    if args.nnls_block != 1 and (not args.cone or args.cone_backend != 'explicit'):
        parser.error('--nnls-block requires --cone --cone-backend explicit')
    if args.face_guard and not args.cone:
        parser.error('--face-guard requires --cone')
    if args.dual_polish < 0 or (args.dual_polish and not args.cone):
        parser.error('--dual-polish requires --cone and a positive iteration budget')
    if args.first_order < 0 or (args.first_order and not args.cone):
        parser.error('--first-order requires --cone and a positive iteration budget')
    if args.newton_qr and not args.first_order:
        parser.error('--newton-qr requires --first-order')
    if args.center_epoch < 0 or (args.center_epoch and not (args.dual_center and args.dual_polish)):
        parser.error('--center-epoch requires --dual-center and --dual-polish')
    if args.implicit_bulk and (not args.cone or args.cone_backend != 'implicit'):
        parser.error('--implicit-bulk requires the implicit cone backend')
    if args.implicit_updates and (not args.cone or args.cone_backend != 'implicit'):
        parser.error('--implicit-updates requires the implicit cone backend')
    if args.mode == 'reference' and args.no_presolve:
        parser.error('The reference-point probe always starts from the original model')
    root = Path(__file__).resolve().parents[1]
    binary = args.binary.resolve()
    args.output.mkdir(parents=True, exist_ok=True)
    with (root/'data/netlib/feasible_gurobi_1e-8.csv').open() as stream:
        references = {r['name']: r['pobj'] for r in csv.DictReader(stream)}
    env = {k: v for k, v in os.environ.items() if not k.startswith('LP_SIMPLEX_')}
    env.update(LP_SIMPLEX_PROFILE='1', OMP_NUM_THREADS='1',
               VECLIB_MAXIMUM_THREADS='1', OPENBLAS_NUM_THREADS='1')
    if args.mode == 'dual-feasible':
        env['LP_SIMPLEX_PAN_DUAL_START'] = '1'
    if args.affine:
        env['LP_SIMPLEX_PAN_AFFINE_MULTIPLIER'] = '1'
    if args.complete:
        env['LP_SIMPLEX_PAN_COMPLETE_START'] = '1'
    if args.stabilize:
        env['LP_SIMPLEX_PAN_STABILIZE'] = '1'
    if args.compatible:
        env['LP_SIMPLEX_PAN_COMPATIBLE_PRICING'] = '1'
    if args.cone:
        env['LP_SIMPLEX_PAN_CONE'] = 'implicit' if args.cone_backend == 'implicit' else '1'
    if args.native_bounds:
        env['LP_SIMPLEX_PAN_NATIVE_BOUNDS'] = '1'
    if args.native_adaptive:
        env['LP_SIMPLEX_PAN_NATIVE_ADAPTIVE'] = '1'
    if args.native_predict:
        env['LP_SIMPLEX_PAN_NATIVE_PREDICT'] = '1'
    if args.native_pdhg:
        env['LP_SIMPLEX_PAN_NATIVE_PDHG'] = str(args.native_pdhg)
    if args.pdhg_weight is not None:
        env['LP_SIMPLEX_PAN_PDHG_WEIGHT'] = str(args.pdhg_weight)
    if args.pdhg_restore:
        env['LP_SIMPLEX_PAN_PDHG_RESTORE'] = '1'
    if args.block_crash:
        env['LP_SIMPLEX_PAN_BLOCK_CRASH'] = '1'
    if args.full_refinement:
        env['LP_SIMPLEX_PAN_FULL_REFINEMENT'] = '1'
    if args.dual_center:
        env['LP_SIMPLEX_PAN_DUAL_CENTER'] = '1'
    if args.cone_batch:
        env['LP_SIMPLEX_PAN_CONE_BATCH'] = '1'
    if args.nnls_block != 1:
        env['LP_SIMPLEX_PAN_NNLS_BLOCK'] = str(args.nnls_block)
    if args.face_guard:
        env['LP_SIMPLEX_PAN_FACE_GUARD'] = '1'
    if args.dual_polish:
        env['LP_SIMPLEX_PAN_DUAL_POLISH'] = str(args.dual_polish)
    if args.first_order:
        env['LP_SIMPLEX_PAN_FIRST_ORDER'] = str(args.first_order)
    if args.newton_qr:
        env['LP_SIMPLEX_PAN_NEWTON_QR'] = '1'
    if args.center_epoch:
        env['LP_SIMPLEX_PAN_CENTER_EPOCH'] = str(args.center_epoch)
    if args.implicit_bulk:
        env['LP_SIMPLEX_PAN_IMPLICIT_BULK'] = '1'
    if args.implicit_updates:
        env['LP_SIMPLEX_PAN_IMPLICIT_UPDATES'] = '1'
    tested = [binary, binary.parent/'pan_probe', binary.parent.parent/'liblp_simplex.dylib']
    hashes = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in tested if p.exists()}
    source_files = [root/'CMakeLists.txt']
    for directory in ['src', 'include', 'test', 'tools']:
        source_files.extend(p for p in (root/directory).rglob('*')
                            if p.is_file() and (p.suffix in {'.c', '.h', '.py'} or p.name == 'CMakeLists.txt'))
    manifest = dict(started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                    options=vars(args), binary_sha256=hashes,
                    platform=platform.platform(), python=platform.python_version(),
                    solver_environment={k: v for k, v in env.items() if k.startswith('LP_SIMPLEX_') or
                                        k in {'OMP_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'OPENBLAS_NUM_THREADS'}},
                    source_sha256={str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
                                   for p in sorted(source_files)})
    (args.output/'manifest.json').write_text(json.dumps(manifest, indent=2, default=str)+'\n')

    def run(case):
        model = root/'data/netlib/feasible'/f'{case}.mps'
        if args.mode == 'reference':
            command = [str(binary.parent/'pan_probe'), str(model),
                       str((args.points/f'{case}.txt').resolve()), str(args.iterations)]
            if args.save_points:
                command.extend([references.get(case, '-'), str((args.output/f'{case}-point.txt').resolve())])
            elif case in references:
                command.append(references[case])
        else:
            command = [str(binary), '--algorithm', 'pan-bda', '--iterations',
                       str(args.iterations), '--criteria', args.pricing, str(model)]
            if args.no_presolve:
                command.insert(1, '--no-presolve')
            if args.save_points:
                command[1:1] = ['--save-point', str((args.output/f'{case}-point.txt').resolve())]
        result = dict(case=case, mode=args.mode, affine=args.affine, complete=args.complete,
                      stabilize=args.stabilize, pricing=args.pricing, compatible=args.compatible, cone=args.cone,
                      cone_backend=args.cone_backend if args.cone else None, presolve=args.mode != 'reference' and not args.no_presolve,
                      native_bounds=args.native_bounds,
                      native_adaptive=args.native_adaptive,
                      native_predict=args.native_predict,
                      native_pdhg=args.native_pdhg,
                      pdhg_weight=args.pdhg_weight,
                      pdhg_restore=args.pdhg_restore,
                      block_crash=args.block_crash,
                      full_refinement=args.full_refinement,
                      dual_center=args.dual_center,
                      cone_batch=args.cone_batch,
                      nnls_block=args.nnls_block,
                      face_guard=args.face_guard,
                      dual_polish=args.dual_polish,
                      first_order=args.first_order,
                      newton_qr=args.newton_qr,
                      center_epoch=args.center_epoch,
                      implicit_bulk=args.implicit_bulk,
                      implicit_updates=args.implicit_updates,
                      model_sha256=hashlib.sha256(model.read_bytes()).hexdigest(),
                      independent_solve=args.mode != 'reference', command=command,
                      sha256=hashes, validated=False, timeout=False)
        begun = time.monotonic()
        logfile = args.output/f'{case}.log'
        with logfile.open('w') as stream:
            try:
                process = subprocess.run(command, env=env, stdout=stream,
                                         stderr=subprocess.STDOUT, timeout=args.seconds)
                result['returncode'] = process.returncode
            except subprocess.TimeoutExpired:
                result['timeout'] = True
        result['wall_seconds'] = time.monotonic()-begun
        result['binary_unchanged'] = all(hashlib.sha256(Path(p).read_bytes()).hexdigest() == digest
                                         for p, digest in hashes.items())
        text = logfile.read_text(errors='replace')
        first = re.search(r'^pan-bda first feasible: .*cpu=([0-9.eE+-]+).*$', text, re.M)
        if first:
            result['first_feasible_cpu_seconds'] = float(first.group(1))
            result['first_feasible_line'] = first.group(0)
        point = args.output/f'{case}-point.txt'
        if args.save_points and point.exists() and 'saved point:' in text:
            result['saved_point_sha256'] = hashlib.sha256(point.read_bytes()).hexdigest()
        result['validated'] = (result['binary_unchanged'] and not result['timeout'] and result.get('returncode') == 0 and
                               ('[PASS] prediction reached' in text or 'verified=1' in text))
        for label in ['pan-bda cone:', 'pan-bda cone start:', 'pan-bda cone summary:', 'pan-bda cone qr:', 'pan-bda cone timing:', 'pan-bda dual polish:', 'pan-bda first order:', 'pan-bda first-order seed:', 'pan-bda nnls blocks:', 'pan-bda support:', 'pan-bda profile:', 'pan-bda phase2:', 'pan-bda precision failure:',
                      'pan probe:', 'saved point:', 'result:', 'iterations:', 'objective:', 'residuals:']:
            lines = [line for line in text.splitlines() if line.startswith(label)]
            if lines:
                result[label.rstrip(':')] = lines[-1]
        (args.output/f'{case}.json').write_text(json.dumps(result, indent=2)+'\n')
        print(json.dumps(result), flush=True)
        return result

    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        results = list(pool.map(run, args.cases))
    (args.output/'runs.json').write_text(json.dumps(results, indent=2)+'\n')


if __name__ == '__main__':
    main()
