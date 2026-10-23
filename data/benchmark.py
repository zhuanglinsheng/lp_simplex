#!/usr/bin/env python3
"""Serial LP comparison and generated reports. Run --help for reproducible options."""
import argparse
import csv
from contextlib import contextmanager
import datetime as dt
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import re
import shutil
import shlex
import signal
import statistics
import struct
import subprocess
import sys
import tempfile
import time

ROOT = Path(__file__).resolve().parents[1]
DATASETS = ('netlib', 'pds', 'railway')
LEGACY_SOLVERS = ('lp_simplex', 'scip', 'highs')
SOLVERS = (*LEGACY_SOLVERS, 'clp', 'glpk')
LABELS = {'lp_simplex': 'lp_simplex', 'scip': 'SCIP / SoPlex', 'highs': 'HiGHS', 'clp': 'Clp', 'glpk': 'GLPK'}
THREADS = dict(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
               VECLIB_MAXIMUM_THREADS='1', MKL_NUM_THREADS='1',
               BLIS_NUM_THREADS='1', NUMEXPR_NUM_THREADS='1')


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


def dump(path, value):
    atomic_text(path, json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + '\n')


def atomic_text(path, text):
    temporary = path.with_name('.' + path.name + '.tmp')
    try:
        temporary.write_text(text)
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def artifact(prefix, suffix):
    return Path(str(prefix) + suffix)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def capture(command):
    try:
        p = subprocess.run(command, cwd=ROOT, capture_output=True, text=True, timeout=15)
        return p.stdout.strip() if p.returncode == 0 else 'unknown'
    except (OSError, subprocess.TimeoutExpired):
        return 'unknown'


def sources():
    paths = [ROOT / 'CMakeLists.txt', *ROOT.glob('src/**/*.[ch]'),
             *ROOT.glob('include/**/*.h'), ROOT / 'data/CMakeLists.txt',
             ROOT / 'data/netlib/CMakeLists.txt',
             *ROOT.glob('data/*.[ch]'), *ROOT.glob('data/*.py'),
             *ROOT.glob('data/*/prepare.py'), *ROOT.glob('data/*/raw/*.c')]
    return {str(p.relative_to(ROOT)): sha(p) for p in sorted(paths)}


def build_root(binary):
    """Find the CMake build tree even when the driver is in a nested directory."""
    binary = Path(binary).resolve()
    return next((parent for parent in binary.parents
                 if (parent / 'CMakeCache.txt').is_file()), binary.parent.parent)


def environment(binary, scip, clp, glpk):
    import numpy
    import scipy
    core, factory, binding = load_highs()
    h = factory()
    if sys.platform == 'darwin':
        cpu = capture(['sysctl', '-n', 'machdep.cpu.brand_string'])
        memory = capture(['sysctl', '-n', 'hw.memsize'])
        linkage = capture(['otool', '-L', str(binary)])
        # sysctl can be denied in a sandbox; retain only non-identifying hardware fields.
        if cpu == 'unknown' or memory == 'unknown':
            try:
                hardware = json.loads(capture(['system_profiler', 'SPHardwareDataType', '-json']))['SPHardwareDataType'][0]
                cpu = hardware.get('chip_type', hardware.get('cpu_type', cpu))
                amount = re.fullmatch(r'([\d.]+) (GB|TB)', hardware.get('physical_memory', ''))
                if amount:
                    memory = str(int(float(amount[1]) * 1024 ** (3 if amount[2] == 'GB' else 4)))
            except (ValueError, KeyError, IndexError):
                pass
    else:
        cpu = platform.processor() or 'unknown'
        cpuinfo = Path('/proc/cpuinfo')
        if cpuinfo.exists():
            match = re.search(r'^model name\s*:\s*(.+)$', cpuinfo.read_text(), re.M)
            if match:
                cpu = match.group(1)
        try:
            memory = str(os.sysconf('SC_PHYS_PAGES') * os.sysconf('SC_PAGE_SIZE'))
        except (ValueError, OSError):
            memory = 'unknown'
        linkage = capture(['ldd', str(binary)])
    # lp_solve_mps links the actual shared solver. Record it, not just the launcher.
    libraries = {}
    for name in re.findall(r'(/[^\s()]*liblp_simplex[^\s()]*)', linkage):
        p = Path(name)
        if p.is_file():
            libraries[str(p)] = sha(p)
    build_dir = build_root(binary)
    for p in build_dir.glob('liblp_simplex.*'):
        if p.is_file():
            libraries[str(p.resolve())] = sha(p)
    cache = build_dir / 'CMakeCache.txt'
    config = {}
    if cache.exists():
        for line in cache.read_text().splitlines():
            if re.match(r'(CMAKE_(BUILD_TYPE|C_COMPILER[^:]*|C_FLAGS[^:]*)|LP_SIMPLEX_[^:]+|.*(?:BLAS|LAPACK|SPQR|KLU|CHOLMOD).*):', line):
                if '=' in line:
                    k, v = line.split('=', 1)
                    config[k] = v
    compiler = config.get('CMAKE_C_COMPILER:FILEPATH', 'cc')
    versions = {
        'lp_simplex': {'git_head': capture(['git', 'rev-parse', 'HEAD']),
                       'git_status': capture(['git', 'status', '--short']),
                       'binary': str(binary), 'binary_sha256': sha(binary),
                       'solver_libraries_sha256': libraries, 'linkage': linkage,
                       'cmake_cache': config,
                       'compiler': capture([compiler, '--version']).splitlines()[0]},
        'scip': {'binary': scip, 'binary_sha256': sha(scip),
                 'version_output': capture([scip, '--version'])},
        'highs': {'version': h.version(), 'binding': binding,
                  'binding_file': core.__file__, 'binding_sha256': sha(core.__file__)}}
    for name, executable in (('clp', clp), ('glpk', glpk)):
        if not executable:
            continue
        versions[name] = dict(binary=executable, binary_sha256=sha(executable),
                              version_output=capture([executable, '-quit' if name == 'clp' else '--version']))
    return dict(os=platform.platform(), architecture=platform.machine(), cpu=cpu,
                logical_cpus=os.cpu_count(), memory_bytes=memory,
                python=sys.version, python_executable=sys.executable,
                numpy=numpy.__version__, scipy=scipy.__version__,
                versions=versions, thread_environment=THREADS)


def audit(lp, x, reported, reference, tolerance):
    """Original-model primal and objective audit; NOT an independent dual proof."""
    import numpy as np
    from scipy.sparse import csc_matrix
    x = np.asarray(x, dtype=np.longdouble)
    if len(x) != lp.num_col_ or not np.all(np.isfinite(x)):
        return dict(audit_pass=False, primal_residual=None, objective=None,
                    reason='invalid primal vector')
    matrix = csc_matrix((np.asarray(lp.a_matrix_.value_, dtype=np.longdouble),
                         lp.a_matrix_.index_, lp.a_matrix_.start_),
                        shape=(lp.num_row_, lp.num_col_))
    activity = matrix @ x
    if not np.all(np.isfinite(activity)):
        return dict(audit_pass=False, primal_residual=None, objective=None,
                    reason='nonfinite row activity')
    residual = max(0., *(float(np.max(v, initial=0.)) for v in (
        np.asarray(lp.col_lower_) - x, x - np.asarray(lp.col_upper_),
        np.asarray(lp.row_lower_) - activity, activity - np.asarray(lp.row_upper_))))
    objective = float(np.sum(np.asarray(lp.col_cost_, dtype=np.longdouble) * x) + lp.offset_)
    valid = (math.isfinite(objective) and math.isfinite(reported) and
             residual <= tolerance and
             abs(objective - reported) <= tolerance * (1. + abs(objective)) and
             (reference is None or abs(objective - reference) <= tolerance * (1. + abs(reference))))
    return dict(audit_pass=valid, primal_residual=residual,
                objective=objective if math.isfinite(objective) else None,
                reference_checked=reference is not None,
                reference_error=abs(objective - reference) if reference is not None and math.isfinite(objective) else None)


def parse_scip_solution(path, names):
    text = path.read_text()
    match = re.search(r'^objective value:\s*(\S+)', text, re.M)
    if not match:
        raise ValueError('SCIP solution has no objective')
    positions = {name: j for j, name in enumerate(names)}
    x = [0.] * len(names)  # SCIP omits zero coordinates.
    seen = set()
    for line in text.splitlines()[2:]:
        parts = line.split()
        if not parts:
            continue
        if len(parts) < 2 or parts[0] not in positions or parts[0] in seen:
            raise ValueError('Unknown/duplicate SCIP solution coordinate: ' + line)
        seen.add(parts[0])
        x[positions[parts[0]]] = float(parts[1])
    return x, float(match.group(1))


def parse_clp_solution(path, rows, columns):
    """Read Clp saveSolution's native binary format without decimal rounding."""
    data = path.read_bytes()
    if len(data) < 16:
        raise ValueError('Truncated Clp solution')
    nr, nc, objective = struct.unpack_from('=iid', data)
    if (nr, nc) != (rows, columns) or len(data) != 16 + 16 * (nr + nc):
        raise ValueError('Clp solution dimensions/length mismatch')
    return list(struct.unpack_from(f'={nc}d', data, 16 + 16 * nr)), objective


def parse_glpk_solution(path, columns):
    lines = [line.split() for line in path.read_text().splitlines()]
    header = next((line for line in lines if line[:2] == ['s', 'bas']), None)
    if header is None or int(header[3]) != columns or header[4:6] != ['f', 'f']:
        raise ValueError('GLPK solution is not optimal or has wrong dimensions')
    x, seen = [0.] * columns, set()
    for line in lines:
        if line[:1] != ['j']:
            continue
        j = int(line[1]) - 1
        if not 0 <= j < columns or j in seen:
            raise ValueError('Unknown/duplicate GLPK coordinate')
        seen.add(j)
        x[j] = float(line[3])
    if len(seen) != columns:
        raise ValueError('Incomplete GLPK solution')
    return x, float(header[6])


def glpk_input(model):
    """GLPK's free MPS reader rejects OBJSENSE; CLI --min/--max carries it."""
    needs_copy = False
    with model.open() as stream:
        for line in stream:
            if line.strip().split()[:1] == ['OBJSENSE']:
                needs_copy = True
                break
            if line.strip() == 'ROWS':
                break
    if not needs_copy:
        return model
    target = model.with_name(model.stem + '-glpk.mps')
    skip_sense = False
    with model.open() as source, target.open('w') as output:
        for line in source:
            fields = line.split()
            if fields[:1] == ['OBJSENSE']:
                skip_sense = len(fields) == 1
                continue
            if skip_sense and fields and not line.startswith('*'):
                if fields[0].upper() not in ('MIN', 'MINIMIZE', 'MAX', 'MAXIMIZE'):
                    raise ValueError('Malformed OBJSENSE')
                skip_sense = False
                continue
            output.write(line)
    return target


def worker(config, solver, model_info, prefix):
    core, factory, _ = load_highs()
    reader = factory()
    reader.setOptionValue('output_flag', False)
    model = ROOT / model_info['path']
    if reader.readModel(str(model)) != core.HighsStatus.kOk:
        raise ValueError('Cannot read original MPS')
    lp = reader.getLp()
    if lp.num_col_ != len(set(lp.col_names_)):
        raise ValueError('Original MPS column names must be unique')
    tol, seconds = config['tolerance'], config['seconds']
    point, log = artifact(prefix, '.point'), artifact(prefix, '.log')
    record = dict(rows=lp.num_row_, columns=lp.num_col_, nonzeros=len(lp.a_matrix_.value_))
    dump(artifact(prefix, '.model.json'), record)
    reported, x, optimal = math.nan, None, False
    if solver == 'highs':
        h = factory()
        options = dict(output_flag=True, log_to_console=False, log_file=str(log), threads=1, presolve='on',
                       solver=config['highs_algorithm'], time_limit=seconds,
                       primal_feasibility_tolerance=tol, dual_feasibility_tolerance=tol,
                       ipm_optimality_tolerance=tol, random_seed=0, run_crossover='on')
        if config['highs_algorithm'] == 'simplex':
            options['simplex_strategy'] = 1
        for k, v in options.items():
            if h.setOptionValue(k, v) != core.HighsStatus.kOk:
                raise ValueError('Unsupported HiGHS option: ' + k)
        record['options'] = options
    else:
        if solver == 'lp_simplex':
            command = [config['binary'], '--algorithm', config['lp_algorithm'],
                       '--iterations', str(config['iterations']),
                       '--save-point', str(point), str(model)]
            if model_info['reference'] is not None:
                command[1:1] = ['--expected', str(model_info['reference'])]
        elif solver == 'clp':
            command = [config['clp'], str(model),
                       '-maximize' if lp.sense_ == core.ObjSense.kMaximize else '-minimize',
                       '-primalTolerance', str(tol),
                       '-dualTolerance', str(tol), '-seconds', str(seconds),
                       '-presolve', 'on', '-dualsimplex', '-saveSolution', str(point), '-quit']
        elif solver == 'glpk':
            compatible = glpk_input(model)
            command = [config['glpk'], '--freemps', str(compatible),
                       '--max' if lp.sense_ == core.ObjSense.kMaximize else '--min',
                       '--nomip', '--simplex',
                       '--tmlim', str(math.ceil(seconds)), '--write', str(point)]
        else:
            command = [config['scip']]
            commands = ['set parallel maxnthreads 1', 'set lp advanced threads 1',
                        'set presolving milp threads 1', 'set randomization randomseedshift 0',
                        f'set numerics feastol {tol}', f'set numerics dualfeastol {tol}',
                        f'set limits time {seconds}', f'read "{model}"', 'optimize',
                        f'write solution "{point}"', 'display statistics', 'quit']
            if config['scip_algorithm'] == 'dual':
                commands[0:0] = ['set lp initalgorithm d', 'set lp resolvealgorithm d']
            for item in commands:
                command.extend(['-c', item])
        record['command'] = command
    # Keep commands/options even when the worker is killed at its deadline.
    dump(artifact(prefix, '.request.json'), record)
    # Parent applies the same hard deadline to reading + solving + point export.
    dump(artifact(prefix, '.clock.json'), {'start': time.monotonic(), 'done': False})
    start = time.monotonic()
    if solver == 'highs':
        read_status = h.readModel(str(model))
        if read_status != core.HighsStatus.kOk:
            raise ValueError('HiGHS read failed')
        run_status = h.run()
        solution, info, status = h.getSolution(), h.getInfo(), h.getModelStatus()
        optimal = run_status == core.HighsStatus.kOk and status == core.HighsModelStatus.kOptimal
        record['solver_status'] = str(status)
        reported = info.objective_function_value
        if solution.value_valid:
            x = list(solution.col_value)
            point.write_text(''.join(f'{v:.17g}\n' for v in x))
        if status == core.HighsModelStatus.kTimeLimit:
            record['status'] = 'TIMEOUT'
    else:
        with log.open('w') as stream:
            process = subprocess.run(command, cwd=ROOT, stdout=stream,
                                     stderr=subprocess.STDOUT, check=False)
        record['returncode'] = process.returncode
        text = log.read_text(errors='replace')
        if solver == 'lp_simplex':
            optimal = process.returncode == 0 and bool(re.search(r'^result:.*code=0 \(success\)', text, re.M))
            found = re.search(r'^objective:\s*(\S+)', text, re.M)
            if found:
                reported = float(found.group(1))
            record['solver_status'] = next((s for s in text.splitlines() if s.startswith('result:')), 'ERROR')
            if point.exists():
                x = [float(v) for v in point.read_text().split()]
        elif solver in ('clp', 'glpk'):
            optimal = process.returncode == 0 and point.exists() and (
                bool(re.search(r'^Optimal(?: objective| - objective)', text, re.M)) if solver == 'clp' else
                bool(re.search(r'^c Status:\s+OPTIMAL$', point.read_text(), re.M)))
            record['solver_status'] = 'optimal' if optimal else 'not optimal'
            if re.search(r'(?i)(time limit|stopped on time|maximum time)', text):
                record['status'] = 'TIMEOUT'
            if optimal:
                x, reported = (parse_clp_solution(point, lp.num_row_, lp.num_col_) if solver == 'clp'
                               else parse_glpk_solution(point, lp.num_col_))
                if solver == 'glpk':
                    record['native_reported_objective'] = reported
                    record['objective_offset_adjustment'] = 2 * lp.offset_
                    reported += 2 * lp.offset_
        else:
            if re.search(r'(?i)(error|invalid parameter|unknown parameter)', text):
                raise ValueError('SCIP configuration/runtime error; see log')
            found = re.search(r'^SCIP Status\s*:\s*(.+)$', text, re.M)
            record['solver_status'] = found.group(1) if found else 'ERROR'
            optimal = process.returncode == 0 and 'optimal solution found' in record['solver_status']
            if 'time limit' in record['solver_status']:
                record['status'] = 'TIMEOUT'
            if point.exists():
                x, reported = parse_scip_solution(point, list(lp.col_names_))
    record['wall_seconds'] = time.monotonic() - start
    dump(artifact(prefix, '.clock.json'), {'done': True})
    if record['wall_seconds'] > seconds:
        record['status'] = 'TIMEOUT'
    record.setdefault('status', 'SOLVER_FAILED')
    if x is not None:
        record.update(audit(lp, x, reported, model_info['reference'], tol))
    if record['status'] != 'TIMEOUT' and optimal:
        record['status'] = 'PASS' if record.get('audit_pass') else 'AUDIT_FAILED'
    if point.exists():
        with point.open('rb') as src, gzip.open(str(point) + '.gz', 'wb') as dest:
            shutil.copyfileobj(src, dest)
        point.unlink()
    dump(artifact(prefix, '.json'), record)


def kill(process):
    try:
        os.killpg(process.pid, signal.SIGKILL)
    except ProcessLookupError:
        pass
    process.wait()


def run_one(config_path, config, solver, model, repeat, env):
    prefix = config_path.parent / 'logs' / f"{model['name']}-{solver}-{repeat}"
    command = [sys.executable, str(Path(__file__).resolve()), '--worker',
               str(config_path), solver, model['name'], str(prefix)]
    start = time.monotonic()
    timed_out = False
    with artifact(prefix, '.adapter.log').open('w') as stream:
        p = subprocess.Popen(command, cwd=ROOT, env=env, stdout=stream,
                             stderr=subprocess.STDOUT, start_new_session=True)
        try:
            while p.poll() is None:
                clock = artifact(prefix, '.clock.json')
                if clock.exists():
                    try:
                        state = json.loads(clock.read_text())
                    except json.JSONDecodeError:
                        state = {}
                    if not state.get('done', True) and time.monotonic() - state['start'] > config['seconds']:
                        timed_out = True
                        kill(p)
                        break
                if time.monotonic() - start > config['seconds'] + 120:
                    kill(p)
                    break
                time.sleep(.02)
        except BaseException:
            kill(p)
            raise
    result_path = artifact(prefix, '.json')
    if result_path.exists() and not timed_out and p.returncode == 0:
        result = json.loads(result_path.read_text())
    else:
        result = json.loads(artifact(prefix, '.request.json').read_text()) if artifact(prefix, '.request.json').exists() else {}
        result.update(status='TIMEOUT' if timed_out else 'ERROR',
                      wall_seconds=config['seconds'] if timed_out else None,
                      returncode=p.returncode)
    if artifact(prefix, '.model.json').exists():
        result.update(json.loads(artifact(prefix, '.model.json').read_text()))
    result.update(instance=model['name'], solver=solver, repeat=repeat,
                  log=str(artifact(prefix, '.log').relative_to(config_path.parent)),
                  adapter_log=str(artifact(prefix, '.adapter.log').relative_to(config_path.parent)))
    print(f"{model['name']} {solver} #{repeat}: {result['status']} "
          f"{result.get('wall_seconds')}s", flush=True)
    return result


def aggregate(results, names, repeats, seconds, solvers=SOLVERS):
    rows = []
    for name in names:
        row = dict(instance=name)
        for solver in solvers:
            runs = [r for r in results if r['instance'] == name and r['solver'] == solver]
            passed = len(runs) == repeats and all(r['status'] == 'PASS' for r in runs)
            row[solver] = dict(passed=passed,
                               status='PASS' if passed else ('/'.join(sorted({r['status'] for r in runs}))
                                       if len(runs) == repeats else 'PENDING' if not runs else 'INCOMPLETE'),
                               median=statistics.median(r['wall_seconds'] for r in runs) if passed else None,
                               penalized=statistics.median(r['wall_seconds'] if r['status'] == 'PASS'
                                                          else 10 * seconds for r in runs) if len(runs) == repeats else None)
        rows.append(row)
    return rows


def report(run_dir, docs_dir, allow_partial=False):
    m = json.loads((run_dir / 'metadata.json').read_text())
    results = json.loads((run_dir / 'results.json').read_text())
    c, e = m['config'], m['environment']
    solvers = tuple(c.get('solvers', LEGACY_SOLVERS))
    if (not m.get('completed') and not allow_partial) or any(m.get(k, True) for k in (
            'source_changed_during_run', 'binary_changed_during_run',
            'library_changed_during_run', 'model_changed_during_run',
            'competitor_changed_during_run')):
        raise ValueError('Incomplete run or changed implementation/models; cannot publish')
    expected = {(model['name'], solver, repeat) for model in c['models']
                for solver in solvers for repeat in range(1, c['repeats'] + 1)}
    actual = [(r['instance'], r['solver'], r['repeat']) for r in results]
    if (len(actual) != len(set(actual)) or not set(actual) <= expected or
            (m.get('completed') and set(actual) != expected)):
        raise ValueError('Missing or duplicate runs; cannot publish')
    dataset = c.get('dataset', 'netlib')
    if dataset not in DATASETS:
        raise ValueError('Unknown benchmark dataset: ' + dataset)
    ensure_reports(docs_dir)
    target = docs_dir / (dataset + '.md')
    stored = read_table(target)
    rows = aggregate(results, [v['name'] for v in c['models']], c['repeats'], c['seconds'], solvers)
    by_name = {model['name']: model for model in c['models']}
    for row in rows:
        name = row['instance']
        runs = [r for r in results if r['instance'] == name]
        dims = next((r for r in runs if 'rows' in r), by_name[name])
        shape = [str(dims.get(k, '—')) for k in ('rows', 'columns', 'nonzeros')]
        def cell(solver):
            value = row[solver]
            elapsed = f"{value['median']:.3f}" if value['passed'] else '—'
            objectives = [r['objective'] for r in runs if r['solver'] == solver
                          and r.get('objective') is not None and math.isfinite(r['objective'])]
            objective = f'{statistics.median(objectives):.12g}' if objectives else '—'
            return f"{value['status']} / {elapsed} / {objective}"
        model = by_name[name]
        reference = f"{model['reference']:.12g}" if model.get('reference') is not None else '—'
        prior = stored.get(name, '')
        previous = [v.strip() for v in prior.strip('|').split('|')] if prior else []
        values = [name, ' / '.join(shape), reference,
                  *(cell(solver) if solver in solvers else
                    previous[3 + SOLVERS.index(solver)] if previous else 'NOT_RUN / — / —'
                    for solver in SOLVERS)]
        stored[name] = '| ' + ' | '.join(values) + ' |'
    blocks = {run_dir.name: run_summary(run_dir, m, results)}
    atomic_text(target, render_report(dataset, stored, blocks))
    overview(docs_dir)


TABLE_HEADER = ('| 实例 | 规模（行 / 列 / 非零元） | 参考目标 | ' +
                ' | '.join(f'{LABELS[k]}（状态 / 秒 / 目标）' for k in SOLVERS) + ' |\n' +
                '| --- | --- | ---: | ' + ' | '.join(['---'] * len(SOLVERS)) + ' |\n')


SUMMARY_HEADER = ('| 求解器 | 成功案例 | 超时次数 | 其他未通过次数 | 与 lp_simplex 共同成功数 | '
                  '共同案例几何平均/秒（lp_simplex / 本求解器） | 相对耗时（lp_simplex / 各求解器） |')
SUMMARY_SEPARATOR = '| --- | ---: | ---: | ---: | ---: | ---: | ---: |'


def read_table(path):
    rows = {}
    if path.exists():
        inside, legacy = False, False
        for line in path.read_text().splitlines():
            if line.startswith('| 实例 |'):
                inside, legacy = True, line.startswith('| 实例 | 行 |')
                continue
            if inside and not line.startswith('|'):
                break
            if inside and line.startswith('| ') and not line.startswith('| ---'):
                cells = [v.strip() for v in line.strip('|').split('|')]
                if legacy:
                    cells = [cells[0], ' / '.join(cells[1:4]), cells[4], *cells[5:8],
                             'NOT_RUN / — / —', 'NOT_RUN / — / —']
                rows[cells[0]] = '| ' + ' | '.join(cells) + ' |'
    return rows


def read_run_blocks(path):
    if not path.exists():
        return {}
    return {name: re.sub(r'^### (?:运行|Run) `[^`\n]+`\n\n', '', block, flags=re.M)
            for name, block in re.findall(
        r'<!-- benchmark-run:([^\n]+) -->\n(.*?)<!-- /benchmark-run -->',
        path.read_text(), re.S)}


def geometric_mean(values):
    if not values:
        return None
    if any(not math.isfinite(v) or v < 0 for v in values):
        raise ValueError('Invalid elapsed time')
    return 0. if 0. in values else math.exp(statistics.mean(math.log(v) for v in values))


def success_statistics(results, names, repeats, seconds, solvers=SOLVERS):
    rows = aggregate(results, names, repeats, seconds, solvers)
    solved = {solver: {row['instance']: row[solver]['median'] for row in rows
                       if row[solver]['passed'] and all(
                           r.get('wall_seconds') is not None and math.isfinite(r['wall_seconds'])
                           and 0 <= r['wall_seconds'] <= seconds
                           for r in results if r['instance'] == row['instance'] and r['solver'] == solver)}
              for solver in solvers}
    output = {}
    for solver in solvers:
        common = set(solved['lp_simplex']) & set(solved[solver])
        runs = [r for r in results if r['solver'] == solver]
        output[solver] = dict(successes=len(solved[solver]),
            timeouts=sum(r['status'] == 'TIMEOUT' for r in runs),
            failures=sum(r['status'] not in ('PASS', 'TIMEOUT', 'SKIPPED') for r in runs),
            skipped=sum(r['status'] == 'SKIPPED' for r in runs),
            pending=len(names) * repeats - len(runs),
            lp_common_mean=geometric_mean([solved['lp_simplex'][name] for name in sorted(common)]),
            common_mean=geometric_mean([solved[solver][name] for name in sorted(common)]),
            common_count=len(common))
    return output


def run_summary(run_dir, metadata, results):
    c, e = metadata['config'], metadata['environment']
    solvers = tuple(c.get('solvers', LEGACY_SOLVERS))
    names = [model['name'] for model in c['models']]
    stats = success_statistics(results, names, c['repeats'], c['seconds'], solvers)
    state = '已完成' if metadata.get('completed') else metadata.get('state', '进行中')
    done = len(metadata.get('finished_cases', names if metadata.get('completed') else []))
    versions = e.get('versions', {})
    lp, scip, highs = (versions.get(key, {}) for key in LEGACY_SOLVERS)
    def clean(value):
        return str(value).replace('|', '/').replace('\n', '；')
    info = [
        ('状态 / 进度', f'{state}；{done}/{len(names)} 个案例'),
        ('时间（UTC）', metadata.get('created_utc', '未知')),
        ('参与求解器', ', '.join(LABELS[key] for key in solvers)),
        ('案例范围', c.get('suite_label', c.get('dataset', 'netlib')) + '；' + ', '.join(ordered_cases(c.get('dataset', 'netlib'), names))),
        ('系统 / CPU', f"{e.get('os', '未知')} / {e.get('cpu', '未知')}"),
        ('内存 / 逻辑 CPU', f"{e.get('memory_bytes', '未知')} bytes / {e.get('logical_cpus', '未知')}"),
        ('Python / NumPy / SciPy', f"{e.get('python', '未知').split()[0]} / {e.get('numpy', '未知')} / {e.get('scipy', '未知')}"),
        ('lp_simplex 版本', lp.get('git_head', '未知') + ('；工作区有改动' if lp.get('git_status') else '')),
        ('构建 / 编译器', f"{lp.get('cmake_cache', {}).get('CMAKE_BUILD_TYPE:STRING', '未知')} / {lp.get('compiler', '未知')}"),
        ('SCIP / LP 后端', (scip.get('version_output') or '未知').splitlines()[0]),
        ('HiGHS', f"{highs.get('version', '未知')} / {highs.get('binding', '未知')}"),
        ('算法（lp / scip / highs）', ' / '.join(c[key] for key in ('lp_algorithm', 'scip_algorithm', 'highs_algorithm'))),
        ('重复 / 单次限时 / 容差 / 工作预算', f"{c['repeats']} / {c['seconds']:g} 秒 / {c['tolerance']:g} / {c['iterations']}"),
        ('资源', '串行；每个求解器单线程；未限制内存，不保证机器独占'),
        ('原始记录', os.path.relpath(run_dir, ROOT) + '（配置、环境、SHA-256、命令、日志、解和审核残差）')]
    for solver in ('clp', 'glpk'):
        if solver in solvers:
            info.append((LABELS[solver], versions.get(solver, {}).get('version_output', '未知').splitlines()[0]))
    if 'clp' in solvers:
        info.append(('Clp / GLPK 算法', 'dual simplex / primal simplex；presolve on'))
    if c.get('glpk_stop_on_timeout'):
        info.append(('GLPK 跳过策略', f"任一案例超时后跳过后续测试；已跳过 {stats['glpk']['skipped']} 次"))
    elif c.get('glpk_gate_case'):
        info.append(('GLPK 跳过策略', f"先测 {c['glpk_gate_case']}；未通过则跳过其余案例；已跳过 {stats['glpk']['skipped']} 次"))
    lines = ['| 测试信息 | 值 |', '| --- | --- |']
    lines.extend(f'| {key} | {clean(value)} |' for key, value in info)
    lines += ['', SUMMARY_HEADER, SUMMARY_SEPARATOR]
    def number(value):
        return '—' if value is None else f'{value:.6g}'
    for solver in solvers:
        value = stats[solver]
        if solver == 'lp_simplex':
            count, means, relative = '—', '—', '—'
        else:
            baseline = value['lp_common_mean']
            ratio = baseline / value['common_mean'] if baseline is not None and value['common_mean'] is not None and value['common_mean'] > 0. else None
            count = str(value['common_count'])
            means = f"{number(baseline)} / {number(value['common_mean'])}"
            relative = number(ratio)
        lines.append(f"| {LABELS[solver]} | {value['successes']}/{len(names)} | {value['timeouts']} | "
                     f"{value['failures']} | {count} | {means} | {relative} |")
    lines += ['', '复现本轮案例与参数（新建运行目录）：', '', '```sh', reproduction_command(c), '```']
    return '\n'.join(lines) + '\n\n'


def reproduction_command(config):
    command = ['python3', 'data/benchmark.py', '--dataset', config.get('dataset', 'netlib'),
               '--cases', *(model['name'] for model in config['models'])]
    for option, key in (('--binary', 'binary'), ('--scip', 'scip'), ('--clp', 'clp'), ('--glpk', 'glpk'), ('--repeats', 'repeats'),
                        ('--seconds', 'seconds'), ('--iterations', 'iterations'),
                        ('--lp-algorithm', 'lp_algorithm'), ('--scip-algorithm', 'scip_algorithm'),
                        ('--highs-algorithm', 'highs_algorithm')):
        if key in config:
            command.extend([option, str(config[key])])
    command.extend(['--solvers', *config.get('solvers', LEGACY_SOLVERS)])
    if config.get('glpk_gate_case'):
        command.extend(['--glpk-gate-case', config['glpk_gate_case']])
    return shlex.join(command)


def english_result_row(row):
    return row


def ordered_cases(dataset, names):
    pattern = {'pds': r'pds-(\d+)', 'railway': r'rail(\d+)'}.get(dataset)
    def key(name):
        match = re.fullmatch(pattern, name) if pattern else None
        return (0, int(match[1]), name) if match else (1, 0, name)
    return sorted(names, key=key)


def render_report(dataset, stored, blocks, english=False):
    def tr(chinese, translated):
        return translated if english else chinese
    title = {'netlib': 'Netlib feasible', 'pds': 'PDS',
             'railway': tr('Railway 连续 LP 松弛', 'Railway continuous LP relaxations')}[dataset]
    navigation = tr(
        f'简体中文 | [English](./{dataset}_en.md) | [项目主页](../../README.md) | [基准测试汇总](./README.md)',
        f'[简体中文](./{dataset}.md) | English | [Project home](../../README_en.md) | [Benchmark summary](./README_en.md)')
    lines = [tr(f'# {title} 基准测试', f'# {title} Benchmark'), '', navigation, '',
        tr('由 `data/benchmark.py` 逐案例生成 MPS、执行求解器对比测试、清理 MPS 后更新。',
           'Updated by `data/benchmark.py` after each instance: generate an MPS, run the configured solvers, and clean up the MPS.'), '',
        tr('## 测试与统计口径', '## Test and statistical criteria'), '',
        tr('- 计时为读取、求解和解导出的墙钟时间；不含 MPS 生成/删除、格式兼容处理、适配器公共初始化和独立审核。CLI 含进程启动，HiGHS 使用 Python binding，小案例开销有差异。',
           '- Timing is wall-clock time for model loading, solving, and solution export. It excludes MPS generation and deletion, format compatibility processing, shared adapter initialization, and independent audits. CLI timing includes process startup; HiGHS uses a Python binding, so overhead differs on small instances.'),
        tr('- PASS 要求求解器报告最优，原模型行约束、变量界和目标复算通过；有参考目标时再核对参考值（误差限为容差乘以 `1 + abs(reference)`）。无参考值显示 `—`，不代表独立最优性认证。',
           '- PASS requires an optimal solver status and successful audits of original-model rows, bounds, and objective values. Reference objectives are also checked when available, within tolerance times `1 + abs(reference)`. Missing references are shown as `—`; their absence does not establish independent optimality certification.'),
        tr('- 成功案例须全部重复为 PASS 且每次均在限时内；先取各案例重复耗时中位数，再计算普通几何平均（不移位、不对失败加惩罚）。空集合显示 `—`。',
           '- An instance is solved only when every repetition is PASS within the time limit. Take each instance\'s median time, then an ordinary geometric mean, with no shift or failure penalty. Empty sets are shown as `—`.'),
        tr('- 每个求解器分别与 lp_simplex 比较，仅使用双方共同成功的案例。几何平均耗时按“lp_simplex / 本求解器”列出；各行集合可能不同。相对耗时为 lp_simplex 耗时除以本求解器耗时，大于 1 表示本求解器更快；lp_simplex 行不列自比较统计。',
           '- Each solver is compared separately with lp_simplex on their pairwise intersection of solved instances. Geometric mean times are listed as lp_simplex / solver; sets may differ between rows. Relative time is lp_simplex time divided by the solver time; above 1 means this solver is faster. The lp_simplex row omits self-comparisons.'),
        tr('- 成功数按案例计，超时与其他未通过数按重复运行计；统计仅属于对应运行的明确案例范围，进行中结果不构成全集结论。',
           '- Successes count instances; timeouts and other unsuccessful results count individual runs. Statistics apply only to the specified run and case selection. In-progress results do not establish full-suite conclusions.'),
        tr('- NOT_RUN / PENDING / INCOMPLETE 表示从未测试 / 本轮尚未测试 / 未完成全部重复。逐例单元格为“状态 / 秒数中位数 / 复算目标中位数”；非 PASS 的目标值仍须结合状态阅读。',
           '- NOT_RUN / PENDING / INCOMPLETE mean never tested / not yet tested in this run / repetitions incomplete. Instance cells show status / median seconds / median audited objective. Interpret non-PASS objectives together with their status.'),
        tr('- SKIPPED 表示按预设规则跳过，不计为超时或其他未通过，也不参与耗时统计。',
           '- SKIPPED means omitted under a predefined rule; it counts as neither a timeout nor another unsuccessful run and is excluded from timing statistics.'),
        tr('- 每行保留该案例最新记录，未选中案例继续保留。不同运行的环境与配置可能不同，不合并计算跨轮速度排名。',
           '- Each row retains the latest record for that instance. Unselected instances are preserved. Environments and settings may differ between runs, so speed rankings are not pooled across runs.'), '',
        tr('## 测试信息与统计', '## Test information and statistics'), '',
        tr('仅展示最近一次报告更新对应的运行；历史记录保存在各轮原始运行目录中。',
           'Only the run from the most recent report update is shown; historical records remain in their original run directories.'), '']
    if not blocks:
        empty = ['尚未运行；环境、版本、算法、限时和重复次数尚未记录。', '',
                 SUMMARY_HEADER, SUMMARY_SEPARATOR]
        empty.extend(f"| {LABELS[solver]} | 0/{len(stored)} | 0 | 0 | "
                     + ('— | — | — |' if solver == 'lp_simplex' else '0 | — / — | — |')
                     for solver in SOLVERS)
        lines += [english_summary('\n'.join(empty)) if english else '\n'.join(empty), '']
    # Older reports contain multiple run blocks. Keep only the last saved run.
    for name, block in list(blocks.items())[-1:]:
        block = re.sub(r'^(\| 案例范围 \| [^；\n]*；)([^|\n]+)\|$',
                       lambda match: match[1] + ', '.join(ordered_cases(
                           dataset, [case.strip() for case in match[2].split(',')])) + ' |',
                       block, flags=re.M)
        lines += [f'<!-- benchmark-run:{name} -->',
                  english_summary(block.rstrip()) if english else block.rstrip(),
                  '<!-- /benchmark-run -->', '']
    header = ('| Instance | Size (rows / columns / nonzeros) | Reference objective | ' +
              ' | '.join(f'{LABELS[k]} (status / s / objective)' for k in SOLVERS) + ' |\n' +
              TABLE_HEADER.splitlines()[1]) if english else TABLE_HEADER.rstrip()
    lines += [tr('## 逐例最新结果', '## Latest per-instance results'), '', header]
    lines.extend(english_result_row(stored[name]) if english else stored[name]
                 for name in ordered_cases(dataset, stored))
    return '\n'.join(lines) + '\n'


def raw_inputs(dataset):
    directory = ROOT / 'data' / dataset / 'raw'
    if dataset == 'netlib':
        return {p.name[:-7]: p for p in sorted((directory / 'feasible').glob('*.mps.gz'))}
    pattern = '*.bz2' if dataset == 'pds' else 'rail*.gz'
    return {p.stem: p for p in sorted(directory.glob(pattern))}


def ensure_reports(docs_dir):
    """Create all three tables without resetting another dataset's experiment rows."""
    docs_dir.mkdir(parents=True, exist_ok=True)
    for dataset in DATASETS:
        target = docs_dir / (dataset + '.md')
        if not target.exists():
            pending = ['| ' + ' | '.join([name, '— / — / —', '—', *(['NOT_RUN / — / —'] * len(SOLVERS))]) + ' |'
                       for name in raw_inputs(dataset)]
            stored = {line.split('|')[1].strip(): line for line in pending}
            atomic_text(target, render_report(dataset, stored, {}))


def english_summary(summary):
    """Localize saved summary prose without changing commands or numeric cells."""
    labels = {
        '测试信息': 'Test information', '值': 'Value', '参与求解器': 'Solvers in this run',
        '状态 / 进度': 'State / progress', '时间（UTC）': 'Time (UTC)',
        '案例范围': 'Case selection', '系统 / CPU': 'System / CPU',
        '内存 / 逻辑 CPU': 'Memory / logical CPUs',
        'lp_simplex 版本': 'lp_simplex version', '构建 / 编译器': 'Build / compiler',
        'SCIP / LP 后端': 'SCIP / LP backend',
        '算法（lp / scip / highs）': 'Algorithms (lp / scip / highs)',
        '重复 / 单次限时 / 容差 / 工作预算': 'Repeats / time limit / tolerance / work budget',
        '资源': 'Resources', '原始记录': 'Raw records', '求解器': 'Solver',
        '成功案例': 'Solved instances', '超时次数': 'Timeouts',
        '其他未通过次数': 'Other unsuccessful runs', '待测次数': 'Pending runs',
        '各自成功案例几何平均/秒': 'Geometric mean on own solved set / s',
        '成功案例几何平均/秒': 'Geometric mean on solved set / s',
        '三方共同成功数': 'Common solved instances',
        '与 lp_simplex 共同成功数': 'Instances solved by both this solver and lp_simplex',
        '本轮共同成功数': 'Common solved instances in this run',
        'Clp / GLPK 算法': 'Clp / GLPK algorithms', 'Clp': 'Clp', 'GLPK': 'GLPK',
        'GLPK 跳过策略': 'GLPK skip policy',
        '共同案例几何平均/秒（lp_simplex / 本求解器）': 'Geometric means on common set / s (lp_simplex / solver)',
        '共同案例几何平均/秒': 'Geometric mean on common set / s',
        '相对耗时（lp_simplex / 各求解器）': 'Relative time (lp_simplex / solver)',
    }
    def value(key, text):
        if key == '状态 / 进度':
            text = text.replace('已完成', 'Completed').replace('进行中', 'Running')
            text = text.replace('已中断', 'Interrupted').replace('个案例', 'instances')
        elif key == '案例范围':
            text = text.replace('显式子集', 'explicit subset').replace('全集', 'full suite')
        elif key == 'lp_simplex 版本':
            text = text.replace('工作区有改动', 'working tree modified')
        elif key == '重复 / 单次限时 / 容差 / 工作预算':
            text = text.replace(' 秒', ' s')
        elif key == 'Clp / GLPK 算法':
            text = text.replace('；', '; ')
        elif key == 'GLPK 跳过策略':
            text = text.replace('任一案例超时后跳过后续测试', 'Skip remaining GLPK runs after its first timeout').replace('先测 ', 'Test first: ').replace('未通过则跳过其余案例', 'skip remaining cases if it does not pass').replace('已跳过 ', 'Skipped runs: ').replace(' 次', '')
        elif key == '资源':
            return 'Serial; one thread per solver; no memory cap or exclusive machine access'
        elif key == '原始记录':
            text = text.replace('（配置、环境、SHA-256、命令、日志、解和审核残差）',
                                ' (configuration, environment, SHA-256, commands, logs, solutions, and audit residuals)')
        return text.replace('未知', 'unknown').replace('；', '; ')
    lines, in_code = [], False
    for line in summary.splitlines():
        if line.startswith('```'):
            in_code = not in_code
            lines.append(line)
        elif in_code:
            lines.append(line)
        elif line.startswith(('### 运行 ', '### Run ')):
            continue
        elif line.startswith('| '):
            cells = [cell.strip() for cell in line.strip('|').split('|')]
            key = cells[0]
            if key in labels and len(cells) == 2:
                cells[1] = value(key, cells[1])
            lines.append('| ' + ' | '.join(labels.get(cell, cell) for cell in cells) + ' |')
        elif line == '复现本轮案例与参数（新建运行目录）：':
            lines.append('Reproduce these cases and parameters in a new run directory:')
        elif line == '尚未运行；环境、版本、算法、限时和重复次数尚未记录。':
            lines.append('Not yet run; environment, versions, algorithms, time limits, and repeats are not yet recorded.')
        elif line == '仅展示最近一次报告更新对应的运行；历史记录保存在各轮原始运行目录中。':
            lines.append('Only the run from the most recent report update is shown; historical records remain in their original run directories.')
        else:
            lines.append(line)
    return '\n'.join(lines)


def overview(docs_dir):
    """Rebuild both language indexes from reports, without requiring raw run logs."""
    summaries = {}
    for dataset in DATASETS:
        path = docs_dir / (dataset + '.md')
        if path.exists():
            if '## 测试信息与统计\n\n' not in path.read_text():
                raise ValueError('Dataset report has no summary section: ' + str(path))
            stored, blocks = read_table(path), read_run_blocks(path)
            report_text = render_report(dataset, stored, blocks)
            atomic_text(path, report_text)
            summaries[dataset] = report_text.partition('## 测试信息与统计\n\n')[2].partition('## 逐例最新结果')[0].strip()
            atomic_text(docs_dir / (dataset + '_en.md'),
                        render_report(dataset, stored, blocks, english=True))
    for english in (False, True):
        def tr(chinese, translated):
            return translated if english else chinese
        navigation = tr(
            '简体中文 | [English](./README_en.md) | [项目主页](../../README.md) | [文档](../README.md)',
            '[简体中文](./README.md) | English | [Project home](../../README_en.md) | [Documentation](../README_en.md)')
        lines = [tr('# 基准测试汇总', '# Benchmark Summary'), '', navigation, '',
            tr('由 `data/benchmark.py` 自动生成；仅包含测试信息、对比统计及复现说明。',
               'Generated by `data/benchmark.py`; includes test information, comparison statistics, and reproduction instructions.'), '',
            tr('## 测试信息与统计口径', '## Test information and statistical criteria'), '',
            tr('逐案例生成临时 MPS，本轮所有求解器完成全部重复后清理，再同步更新数据集报告和本汇总。',
               'A temporary MPS is generated for each instance and removed after all configured solvers finish their repetitions. Dataset reports and this summary are then updated together.'), '',
            tr('成功案例须全部重复 PASS 且每次均在限时内。先取重复耗时中位数，再计算普通几何平均，'
               '不移位、不对失败加惩罚；空集合为 `—`。各求解器分别与 lp_simplex 使用双方共同成功案例比较；'
               '几何平均按“lp_simplex / 本求解器”列出，各行集合可能不同；lp_simplex 行不列自比较统计。'
               '相对耗时为 lp_simplex 耗时除以本求解器耗时，大于 1 表示本求解器更快。成功数按案例计，超时与其他未通过数按重复运行计。',
               'An instance is solved only if every repetition passes within the time limit. Times are medians across repetitions, followed by an ordinary geometric mean, with no shift or failure penalty. Empty sets are shown as `—`. Each solver is compared separately with lp_simplex on their pairwise intersection of solved instances. Geometric means are listed as lp_simplex / solver; sets may differ between rows. The lp_simplex row omits self-comparisons. Relative time is lp_simplex time divided by the solver time; above 1 means this solver is faster. Successes count instances; timeouts and other unsuccessful results count individual runs.'), '',
            tr('计时含模型读取、求解和解导出，不含 MPS 生成/删除、格式兼容处理、公共适配器初始化和独立审核。'
               '串行、单线程，不保证机器独占；CLI 与 Python binding 开销有差异。',
               'Timing includes model loading, solving, and solution export. It excludes MPS preparation and deletion, format compatibility processing, shared adapter initialization, and independent audits. Runs are serial and single-threaded, without guaranteed exclusive machine access. CLI and Python-binding overheads differ.'), '',
            tr('PASS 要求求解器报告最优且原模型可行性、目标复算通过；有参考目标时额外核对。'
               '缺少参考值不代表独立最优性认证。进行中或子集结果不构成全集结论，'
               '不同运行分别统计，不混合环境、配置或案例集合。',
               'PASS requires an optimal solver status and successful original-model feasibility and objective audits. Reference objectives are checked when available. A pass without a reference is not an independent optimality certificate. Partial or subset results do not establish full-suite conclusions. Runs are summarized separately, without mixing environments, settings, or case selections.'), '']
        for dataset, title in (('netlib', 'Netlib feasible'), ('pds', 'PDS'),
                               ('railway', tr('Railway 连续 LP 松弛', 'Railway continuous LP relaxations'))):
            lines += [f'## {title}', '',
                      f'[{tr("逐案例结果与完整报告", "Per-instance results and full report")}]({dataset}{"_en" if english else ""}.md)', '']
            if dataset in summaries:
                lines += [english_summary(summaries[dataset]) if english else summaries[dataset], '']
            else:
                lines += [tr('尚未生成测试信息与统计。', 'Test information and statistics have not yet been generated.'), '']
        lines += [tr('## 复现说明', '## Reproduction'), '',
            tr('在仓库根目录运行。需 CMake、C 编译器、BLAS/LAPACK，并安装 SCIP、Clp、GLPK、NumPy、SciPy，'
               '以及 highspy 或兼容的 SciPy HiGHS 绑定。PDS 转换器由脚本临时编译。',
               'Run from the repository root. Requirements are CMake, a C compiler, BLAS/LAPACK, SCIP, Clp, GLPK, NumPy, SciPy, and highspy or a compatible SciPy HiGHS binding. The script compiles the PDS converter temporarily.'), '',
            '```sh', 'cmake -S . -B build/benchmark -DCMAKE_BUILD_TYPE=Release -DLP_SIMPLEX_BUILD_BENCHMARKS=ON',
            'cmake --build build/benchmark -j 4',
            'python3 data/benchmark.py --dataset netlib --binary build/benchmark/data/lp_solve_mps',
            'python3 data/benchmark.py --dataset pds --binary build/benchmark/data/lp_solve_mps',
            'python3 data/benchmark.py --dataset railway --binary build/benchmark/data/lp_solve_mps', '```', '',
            tr('省略 `--cases` 测试所选数据集全集；短程验证可加 `--cases <实例名> --repeats 1 --seconds 10`。'
               'simple_examples 仅用于功能测试，不参与基准。上方各轮的复现命令保留其案例与具体参数，'
               '二进制路径可换为新构建结果；复核版本及 SHA-256 后再比较耗时。',
               'Omit `--cases` to test the full selected dataset. For a short check, add `--cases <instance> --repeats 1 --seconds 10`. The simple_examples models are functional tests, not benchmarks. Commands above preserve each run\'s cases and parameters. You may replace binary paths with new builds; verify versions and SHA-256 before comparing timings.'), '',
            tr('默认每轮新建 `tmp/benchmarks/<UTC 时间>/`，保存配置、环境、命令、结果、日志与审核信息。'
               '正常结束、异常、Ctrl+C 或 SIGTERM 清理临时 MPS；SIGKILL/断电残留可随运行目录删除。',
               'Each run creates `tmp/benchmarks/<UTC timestamp>/` for configuration, environment, commands, results, logs, and audit information. Temporary MPS files are cleaned on normal exit, exceptions, Ctrl+C, or SIGTERM. Files left by SIGKILL or power loss can be removed with the run directory.'), '',
            tr('仅重建报告和本汇总，不需要 MPS 或再次求解：',
               'Rebuild reports and both language summaries without MPS files or new solver runs:'), '', '```sh',
            tr('python3 data/benchmark.py --report-only tmp/benchmarks/<运行目录>',
               'python3 data/benchmark.py --report-only tmp/benchmarks/<run-directory>'), '```', '',
            tr('`--docs-dir` 可指定报告目录；重建汇总会保留其他数据集已保存的测试信息和统计。'
               '更多参数与审核说明见 [数据与基准说明](../../data/README.md)。',
               'Use `--docs-dir` to choose the report directory. Rebuilding summaries preserves saved test information and statistics for other datasets. See [data and benchmark instructions](../../data/README_en.md) for additional parameters and audit details.'), '']
        atomic_text(docs_dir / ('README_en.md' if english else 'README.md'), '\n'.join(lines))


@contextmanager
def prepared_model(model, run_dir):
    """Never reuse, overwrite or delete a user's pre-generated MPS."""
    with tempfile.TemporaryDirectory(prefix='mps-', dir=run_dir) as temporary:
        log = run_dir / 'logs' / (model['name'] + '-prepare.log')
        with log.open('w') as stream:
            dataset = model.get('dataset', 'netlib')
            name = ('feasible/' if dataset == 'netlib' else '') + model['name']
            subprocess.run([sys.executable, str(ROOT / 'data' / dataset / 'prepare.py'),
                            name, '--output-dir', temporary],
                           stdout=stream, stderr=subprocess.STDOUT, check=True)
        path = Path(temporary) / ('feasible' if dataset == 'netlib' else '') / (model['name'] + '.mps')
        digest = sha(path)
        if model.get('expected_mps_sha256') and digest != model['expected_mps_sha256']:
            raise ValueError('Generated MPS checksum differs from manifest: ' + model['name'])
        model.update(path=str(path.resolve()), sha256=digest)
        yield path


def glpk_skip_reason(config, solver, model, results):
    if solver == 'glpk' and config.get('glpk_stop_on_timeout'):
        timeout = next((row for row in results if row['solver'] == 'glpk' and row['status'] == 'TIMEOUT'), None)
        if timeout:
            return 'GLPK stopped after timeout: ' + timeout['instance']
    gate = config.get('glpk_gate_case')
    if solver != 'glpk' or not gate or model['name'] == gate:
        return None
    gate_runs = [row for row in results if row['solver'] == 'glpk' and row['instance'] == gate]
    if len(gate_runs) != config['repeats']:
        raise ValueError('GLPK gate case must finish before other cases')
    if any(row['status'] != 'PASS' for row in gate_runs):
        return 'GLPK gate case did not pass: ' + gate
    return None



def benchmark_case(config_path, config, model, case_index, env, results, metadata):
    """Keep only this case's MPS alive; persist results before cleaning it up."""
    run_dir = config_path.parent
    prepared = False
    solvers = tuple(config.get('solvers', SOLVERS))
    try:
        with prepared_model(model, run_dir) as path:
            prepared = True
            dump(config_path, config)
            try:
                for repeat in range(config['repeats']):
                    offset = (case_index + repeat) % len(solvers)
                    for solver in solvers[offset:] + solvers[:offset]:
                        try:
                            reason = glpk_skip_reason(config, solver, model, results)
                            if reason:
                                row = dict(instance=model['name'], solver=solver, repeat=repeat + 1,
                                           status='SKIPPED', wall_seconds=None, reason=reason)
                                print(f"{model['name']} {solver} #{repeat + 1}: SKIPPED ({reason})", flush=True)
                            else:
                                row = run_one(config_path, config, solver, model, repeat + 1, env)
                        except Exception as exc:
                            row = dict(instance=model['name'], solver=solver, repeat=repeat + 1,
                                       status='ERROR', wall_seconds=None, error=str(exc))
                        results.append(row)
                        dump(run_dir / 'results.json', results)
            finally:
                if not path.is_file() or sha(path) != model['sha256']:
                    metadata['model_changed_during_run'] = True
    except (OSError, ValueError, subprocess.CalledProcessError) as exc:
        if prepared:
            raise
        # A bad raw input still produces a visible row for every requested run.
        results.extend(dict(instance=model['name'], solver=solver, repeat=repeat,
                            status='PREPARE_FAILED', wall_seconds=None, error=str(exc),
                            log='logs/' + model['name'] + '-prepare.log')
                       for repeat in range(1, config['repeats'] + 1) for solver in solvers)
        dump(run_dir / 'results.json', results)


def reference_path(dataset):
    if dataset == 'netlib':
        return ROOT / 'data/netlib/feasible_gurobi_1e-8.csv'
    if dataset == 'pds':
        return ROOT / 'data/pds/reference_objectives.csv'
    return None


def select_models(parser, cases, dataset='netlib'):
    available = raw_inputs(dataset)
    if not available:
        parser.error('No raw inputs found for ' + dataset)
    names = cases if cases is not None else list(available)
    if len(names) != len(set(names)):
        parser.error('duplicate cases')
    unknown = set(names) - available.keys()
    if unknown:
        parser.error('unknown cases: ' + ', '.join(sorted(unknown)))
    if dataset in ('pds', 'railway'):
        names = ordered_cases(dataset, names)
    references = {}
    reference_file = reference_path(dataset)
    if reference_file is not None:
        try:
            with reference_file.open() as stream:
                seen = set()
                for row in csv.DictReader(stream):
                    name = row['name']
                    if name in seen:
                        parser.error('duplicate reference objective: ' + name)
                    seen.add(name)
                    if row['status'] == 'OPTIMAL':
                        references[name] = float(row['pobj'])
        except (OSError, KeyError, TypeError, ValueError) as exc:
            parser.error(f'cannot read reference objectives from {reference_file}: {exc}')
        for name in names:
            if name not in references or not math.isfinite(references[name]):
                parser.error('missing finite reference objective: ' + name)
    manifest = ROOT / 'data' / dataset / 'manifest.csv'
    info = {}
    if manifest.exists():
        with manifest.open() as stream:
            info = {row['instance']: row for row in csv.DictReader(stream)}
    models = []
    for name in names:
        row = info.get(name, {})
        model = dict(name=name, dataset=dataset, raw_path=str(available[name].relative_to(ROOT)),
                     raw_sha256=sha(available[name]), reference=references.get(name))
        for key, source in (('rows', 'constraints' if dataset == 'pds' else 'rows'),
                            ('columns', 'variables' if dataset == 'pds' else 'columns'),
                            ('nonzeros', 'constraint_nonzeros' if dataset == 'pds' else 'nonzeros')):
            if row.get(source):
                model[key] = int(row[source])
        if row.get('mps_sha256'):
            model['expected_mps_sha256'] = row['mps_sha256']
        models.append(model)
    return models


def input_metadata(dataset):
    paths = [ROOT / 'data' / dataset / 'manifest.csv']
    reference_file = reference_path(dataset)
    if reference_file is not None:
        paths.append(reference_file)
    return {str(p.relative_to(ROOT)): sha(p) for p in paths if p.exists()}


def checkpoint(run_dir, config, metadata, results, binary, scip, docs_dir):
    """Update progress only after temporary MPS cleanup and integrity checks."""
    detected = metadata['environment']
    metadata['source_changed_during_run'] = sources() != metadata['source_sha256']
    metadata['binary_changed_during_run'] = not binary.is_file() or sha(binary) != detected['versions']['lp_simplex']['binary_sha256']
    metadata['library_changed_during_run'] = any(not Path(p).is_file() or sha(p) != value
        for p, value in detected['versions']['lp_simplex']['solver_libraries_sha256'].items())
    metadata['model_changed_during_run'] |= any(not (ROOT / name).is_file() or sha(ROOT / name) != value
        for name, value in config['input_metadata_sha256'].items()) or any(not (ROOT / p['raw_path']).is_file() or
        sha(ROOT / p['raw_path']) != p['raw_sha256'] for p in config['models'])
    binding = Path(detected['versions']['highs']['binding_file'])
    metadata['competitor_changed_during_run'] = any(
        not Path(detected['versions'][key]['binary']).is_file() or
        sha(detected['versions'][key]['binary']) != detected['versions'][key]['binary_sha256']
        for key in ('clp', 'glpk') if key in detected['versions']) or (
        not Path(scip).is_file() or sha(scip) != detected['versions']['scip']['binary_sha256'] or
        not binding.is_file() or sha(binding) != detected['versions']['highs']['binding_sha256'])
    dump(run_dir / 'metadata.json', metadata)
    keys = sorted({k for row in results for k in row if k not in ('command', 'options')})
    with (run_dir / 'results.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=keys, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(results)
    report(run_dir, docs_dir, allow_partial=True)


def main():
    if len(sys.argv) > 1 and sys.argv[1] == '--worker':
        config_path, solver, name, prefix = sys.argv[2:]
        config = json.loads(Path(config_path).read_text())
        model = next(m for m in config['models'] if m['name'] == name)
        worker(config, solver, model, Path(prefix))
        return
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path, default=ROOT / 'build/data/lp_solve_mps')
    parser.add_argument('--scip', default='scip')
    parser.add_argument('--solvers', nargs='+', choices=SOLVERS, default=list(SOLVERS))
    parser.add_argument('--clp', default='clp')
    parser.add_argument('--glpk', default='glpsol')
    parser.add_argument('--glpk-gate-case', help='Test this case first; skip subsequent GLPK runs if it does not pass')
    parser.add_argument('--run-dir', type=Path, help='raw artifacts directory (default: tmp/benchmarks/<UTC timestamp>)')
    parser.add_argument('--report-only', type=Path)
    parser.add_argument('--docs-dir', type=Path, default=ROOT / 'docs/benchmarks')
    parser.add_argument('--dataset', choices=DATASETS, default='netlib', help='Dataset (default: netlib)')
    parser.add_argument('--cases', nargs='+', help='Instance names; omitted means the selected dataset in full')
    parser.add_argument('--repeats', type=int, default=3)
    parser.add_argument('--seconds', type=float, default=None,
                        help='wall-clock limit per run (default: Netlib 60 s; PDS/Railway 120 s)')
    parser.add_argument('--iterations', type=int, default=300000)
    parser.add_argument('--lp-algorithm', choices=['auto', 'dual-revised', 'pan-bda'], default='auto')
    parser.add_argument('--highs-algorithm', choices=['choose', 'simplex', 'ipm'], default='choose')
    parser.add_argument('--scip-algorithm', choices=['default', 'dual'], default='default')
    args = parser.parse_args()
    if args.report_only:
        report(args.report_only.resolve(), args.docs_dir.resolve(), allow_partial=True)
        return
    if args.seconds is None:
        args.seconds = 60. if args.dataset == 'netlib' else 120.
    if args.repeats < 1 or not math.isfinite(args.seconds) or args.seconds <= 0 or args.iterations < 1:
        parser.error('repeats, seconds, and iterations must be positive')
    models = select_models(parser, args.cases, args.dataset)
    if args.glpk_gate_case:
        if 'glpk' not in args.solvers or args.glpk_gate_case not in {m['name'] for m in models}:
            parser.error('--glpk-gate-case requires GLPK and a case in the selected suite')
        models.sort(key=lambda model: model['name'] != args.glpk_gate_case)
    binary = args.binary.resolve()
    scip = shutil.which(args.scip)
    clp, glpk = shutil.which(args.clp), shutil.which(args.glpk)
    if len(args.solvers) != len(set(args.solvers)) or 'lp_simplex' not in args.solvers:
        parser.error('solvers must be unique and include lp_simplex')
    if not binary.is_file() or not scip or ('clp' in args.solvers and not clp) or ('glpk' in args.solvers and not glpk):
        parser.error('lp_simplex binary, SCIP, Clp and GLPK (glpsol) executables are required')
    config = dict(dataset=args.dataset, binary=str(binary), scip=scip, clp=clp, glpk=glpk, solvers=args.solvers, repeats=args.repeats, seconds=args.seconds,
                  tolerance=1e-7, iterations=args.iterations, lp_algorithm=args.lp_algorithm,
                  highs_algorithm=args.highs_algorithm, scip_algorithm=args.scip_algorithm,
                  input_metadata_sha256=input_metadata(args.dataset),
                  suite_label=args.dataset + (' 全集' if not args.cases else ' 显式子集'),
                  models=models, glpk_stop_on_timeout=True)
    if args.glpk_gate_case:
        config['glpk_gate_case'] = args.glpk_gate_case
    env = {k: v for k, v in os.environ.items() if not k.startswith('LP_SIMPLEX_')}
    env.update(THREADS, LC_ALL='C', PYTHONDONTWRITEBYTECODE='1')
    os.environ.update(THREADS)
    detected = environment(binary, scip, clp, glpk)
    run_dir = (args.run_dir or ROOT / 'tmp/benchmarks' / dt.datetime.now(dt.timezone.utc).strftime('%Y%m%dT%H%M%SZ')).resolve()
    if run_dir.exists():
        parser.error('run directory exists; select a new directory')
    (run_dir / 'logs').mkdir(parents=True)
    config_path = run_dir / 'config.json'
    dump(config_path, config)
    metadata = dict(created_utc=dt.datetime.now(dt.timezone.utc).isoformat(), config=config,
                    invocation=sys.argv, environment=detected, source_sha256=sources(), completed=False,
                    state='进行中', finished_cases=[], source_changed_during_run=False,
                    binary_changed_during_run=False, library_changed_during_run=False,
                    model_changed_during_run=False, competitor_changed_during_run=False)
    results = []
    dump(run_dir / 'results.json', results)
    checkpoint(run_dir, config, metadata, results, binary, scip, args.docs_dir.resolve())
    try:
        for index, model in enumerate(models):
            benchmark_case(config_path, config, model, index, env, results, metadata)
            metadata['finished_cases'].append(model['name'])
            metadata['completed'] = len(metadata['finished_cases']) == len(models)
            metadata['state'] = '已完成' if metadata['completed'] else '进行中'
            checkpoint(run_dir, config, metadata, results, binary, scip, args.docs_dir.resolve())
    except BaseException:
        metadata['completed'] = False
        metadata['state'] = '已中断'
        try:
            checkpoint(run_dir, config, metadata, results, binary, scip, args.docs_dir.resolve())
        except Exception as exc:
            print('Progress report not updated: ' + str(exc), file=sys.stderr)
        raise
    print(f'Report: {args.docs_dir / (args.dataset + ".md")}', flush=True)


if __name__ == '__main__':
    def terminate(signum, frame):
        raise KeyboardInterrupt('Received SIGTERM')
    signal.signal(signal.SIGTERM, terminate)
    main()
