#!/usr/bin/env python3
"""Serial LP comparison and generated reports. Run --help for reproducible options."""
import argparse
import csv
from contextlib import contextmanager
import datetime as dt
import gzip
import hashlib
import inspect
import json
import math
import os
from pathlib import Path
import platform
import re
import shutil
import signal
import struct
import subprocess
import sys
import tempfile
import time

from benchmark_reports import (
    ROOT, DATASETS, LEGACY_SOLVERS, SOLVERS, LABELS,
    TABLE_HEADER, SUMMARY_HEADER, SUMMARY_SEPARATOR,
    aggregate, atomic_text, geometric_mean, success_statistics,
    read_table, read_run_blocks, read_result_cache, run_summary,
    reproduction_command, render_report, ordered_cases, overview, update_table_rows,
)

try:
    import resource
except ImportError:  # Windows has no getrusage; absence is not zero CPU.
    resource = None

LP_ALGORITHMS = ('auto', 'dual-revised', 'pan-bda', 'tableau')
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


def compiler_version(compiler):
    if Path(compiler).name.lower() != 'cl.exe':
        return capture([compiler, '--version']).splitlines()[0]
    try:
        process = subprocess.run([compiler, '/?'], cwd=ROOT, capture_output=True,
                                 text=True, errors='replace', timeout=15)
        # MSVC writes its version banner to stderr, and help text to stdout.
        lines = (process.stderr.strip() or process.stdout.strip()).splitlines()
        return lines[0] if lines else 'unknown'
    except (OSError, subprocess.TimeoutExpired):
        return 'unknown'


def sources():
    paths = [ROOT / 'CMakeLists.txt', *ROOT.glob('src/**/*.[ch]'),
             *ROOT.glob('include/**/*.h'), ROOT / 'data/CMakeLists.txt',
             ROOT / 'data/netlib/CMakeLists.txt',
             *ROOT.glob('data/*.[ch]'), *ROOT.glob('data/*.py'),
             *ROOT.glob('data/*/prepare.py'), *ROOT.glob('data/*/raw/*.c')]
    return {p.relative_to(ROOT).as_posix(): sha(p) for p in sorted(paths)}


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
        linkage = capture(['otool', '-L', str(binary)]) if binary else ''
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
            match = re.search(r'^model name\s*:\s*(.+)$', cpuinfo.read_text(encoding='utf-8'), re.M)
            if match:
                cpu = match.group(1)
        try:
            memory = str(os.sysconf('SC_PHYS_PAGES') * os.sysconf('SC_PAGE_SIZE'))
        except (AttributeError, ValueError, OSError):
            memory = 'unknown'
        linkage = capture(['ldd', str(binary)]) if binary else ''
    versions = {'highs': {'version': h.version(), 'binding': binding,
                          'binding_file': core.__file__, 'binding_sha256': sha(core.__file__)}}
    if binary:
        # Record the actual shared solver as well as its launcher.
        libraries = {}
        for name in re.findall(r'(/[^\s()]*liblp_simplex[^\s()]*)', linkage):
            p = Path(name)
            if p.is_file():
                libraries[str(p)] = sha(p)
        build_dir = build_root(binary)
        for p in build_dir.glob('liblp_simplex.*'):
            if p.is_file():
                libraries[str(p.resolve())] = sha(p)
        for p in build_dir.glob('lp_simplex.dll'):
            libraries[str(p.resolve())] = sha(p)
        cache = build_dir / 'CMakeCache.txt'
        config = {}
        if cache.exists():
            for line in cache.read_text(encoding='utf-8').splitlines():
                if re.match(r'(CMAKE_(BUILD_TYPE|C_COMPILER[^:]*|C_FLAGS[^:]*)|LP_SIMPLEX_[^:]+|.*(?:BLAS|LAPACK|SPQR|KLU|CHOLMOD).*):', line):
                    if '=' in line:
                        k, v = line.split('=', 1)
                        config[k] = v
        compiler = config.get('CMAKE_C_COMPILER:FILEPATH',
                              config.get('CMAKE_C_COMPILER:STRING', 'cc'))
        versions['lp_simplex'] = dict(
            git_head=capture(['git', 'rev-parse', 'HEAD']),
            git_status=capture(['git', 'status', '--short']),
            binary=str(binary), binary_sha256=sha(binary),
            solver_libraries_sha256=libraries, linkage=linkage, cmake_cache=config,
            compiler=compiler_version(compiler))
    for name, executable in (('scip', scip), ('clp', clp), ('glpk', glpk)):
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
    text = path.read_text(encoding='utf-8')
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
    lines = [line.split() for line in path.read_text(encoding='utf-8').splitlines()]
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
    with model.open(encoding='utf-8') as stream:
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
    with model.open(encoding='utf-8') as source, target.open('w', encoding='utf-8') as output:
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


def infeasibility_status(solver, text='', native_status='', returncode=0):
    """Require an explicit terminal decision; ambiguous statuses never pass."""
    if returncode != 0 and not (solver == 'lp_simplex' and returncode == 1):
        return None
    if solver == 'highs':
        if native_status.endswith('kInfeasible'):
            return 'infeasible'
        if native_status.endswith('kUnboundedOrInfeasible'):
            return 'ambiguous'
        if native_status.endswith('kUnbounded'):
            return 'unbounded'
    elif solver == 'lp_simplex':
        match = re.search(r'^result:.*code=\d+ \(([^)]+)\)', text, re.M)
        if match and match[1] in ('infeasible', 'unbounded', 'success'):
            if match[1] == 'infeasible' and returncode != 0:
                return None
            return 'optimal' if match[1] == 'success' else match[1]
    elif solver == 'scip':
        status = native_status.lower()
        if 'infeasible or unbounded' in status or 'unbounded or infeasible' in status:
            return 'ambiguous'
        if 'infeasible' in status:
            return 'infeasible'
        if 'unbounded' in status:
            return 'unbounded'
    elif solver == 'clp':
        if re.search(r'^Primal infeasible(?:\s|$)', text, re.M):
            return 'infeasible'
        if re.search(r'^Dual infeasible(?:\s|$)', text, re.M):
            return 'unbounded'
    elif solver == 'glpk':
        if re.search(r'PROBLEM HAS NO PRIMAL FEASIBLE SOLUTION|LP HAS NO PRIMAL FEASIBLE SOLUTION', text):
            return 'infeasible'
        if 'PROBLEM HAS UNBOUNDED SOLUTION' in text:
            return 'unbounded'
    return None


def child_cpu():
    """The worker launches one CLI synchronously; include its descendants."""
    if resource is None:
        return None
    try:
        usage = resource.getrusage(resource.RUSAGE_CHILDREN)
        values = (usage.ru_utime, usage.ru_stime)
        return values if all(math.isfinite(v) and v >= 0 for v in values) else None
    except (OSError, ValueError):
        return None


def run_cli(command, stream):
    """Measure the CLI lifetime, excluding adapter parsing/auditing afterwards."""
    cpu_before = child_cpu()
    started = time.monotonic()
    process = subprocess.run(command, cwd=ROOT, stdout=stream,
                             stderr=subprocess.STDOUT, check=False)
    wall = time.monotonic() - started
    cpu_after = child_cpu()
    user = system = cpu = None
    if cpu_before is not None and cpu_after is not None:
        delta = tuple(b - a for a, b in zip(cpu_before, cpu_after))
        if all(math.isfinite(v) and v >= 0 for v in delta):
            user, system = delta
            cpu = user + system
    return process, dict(cli_wall_seconds=wall, cli_cpu_seconds=cpu,
                         cli_user_seconds=user, cli_system_seconds=system)


def worker(config, solver, model_info, prefix):
    core, factory, _ = load_highs()
    reader = factory()
    reader.setOptionValue('output_flag', False)
    model = ROOT / model_info['path']
    reader_status = reader.readModel(str(model))
    if reader_status not in (core.HighsStatus.kOk, core.HighsStatus.kWarning):
        raise ValueError('Cannot read original MPS')
    lp = reader.getLp()
    if lp.num_col_ != len(set(lp.col_names_)):
        raise ValueError('Original MPS column names must be unique')
    tol, seconds = config['tolerance'], config['seconds']
    point, log = artifact(prefix, '.point'), artifact(prefix, '.log')
    record = dict(rows=lp.num_row_, columns=lp.num_col_, nonzeros=len(lp.a_matrix_.value_),
                  reader_status=str(reader_status))
    dump(artifact(prefix, '.model.json'), record)
    reported, x, optimal = math.nan, None, False
    expects_infeasible = model_info.get('expected_status') == 'infeasible'
    decision, text = None, ''
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
            if expects_infeasible:
                command.insert(1, '--infeasible')
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
        if read_status not in (core.HighsStatus.kOk, core.HighsStatus.kWarning):
            raise ValueError('HiGHS read failed')
        run_status = h.run()
        solution, info, status = h.getSolution(), h.getInfo(), h.getModelStatus()
        optimal = run_status == core.HighsStatus.kOk and status == core.HighsModelStatus.kOptimal
        record['solver_status'] = str(status)
        decision = infeasibility_status(solver, native_status=str(status))
        reported = info.objective_function_value
        if solution.value_valid:
            x = list(solution.col_value)
            point.write_text(''.join(f'{v:.17g}\n' for v in x), encoding='utf-8')
        if status == core.HighsModelStatus.kTimeLimit:
            record['status'] = 'TIMEOUT'
    else:
        with log.open('w', encoding='utf-8') as stream:
            process, timing = run_cli(command, stream)
        record.update(timing)
        record['returncode'] = process.returncode
        text = log.read_text(encoding='utf-8', errors='replace')
        if solver == 'lp_simplex':
            optimal = process.returncode == 0 and bool(re.search(r'^result:.*code=0 \(success\)', text, re.M))
            found = re.search(r'^objective:\s*(\S+)', text, re.M)
            if found:
                reported = float(found.group(1))
            record['solver_status'] = next((s for s in text.splitlines() if s.startswith('result:')), 'ERROR')
            if point.exists():
                x = [float(v) for v in point.read_text(encoding='utf-8').split()]
        elif solver in ('clp', 'glpk'):
            optimal = process.returncode == 0 and point.exists() and (
                bool(re.search(r'^Optimal(?: objective| - objective)', text, re.M)) if solver == 'clp' else
                bool(re.search(r'^c Status:\s+OPTIMAL$', point.read_text(encoding='utf-8'), re.M)))
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
            if point.exists() and optimal:
                x, reported = parse_scip_solution(point, list(lp.col_names_))
        decision = infeasibility_status(solver, text, record['solver_status'], process.returncode)
    record['wall_seconds'] = time.monotonic() - start
    dump(artifact(prefix, '.clock.json'), {'done': True})
    if record['wall_seconds'] > seconds:
        record['status'] = 'TIMEOUT'
    record.setdefault('status', 'SOLVER_FAILED')
    if x is not None and not expects_infeasible:
        record.update(audit(lp, x, reported, model_info['reference'], tol))
    if expects_infeasible:
        record.update(expected_status='infeasible', reported_status=decision,
                      certificate_verified=False)
        if record['status'] != 'TIMEOUT':
            record['status'] = ('PASS' if decision == 'infeasible' else
                                'AMBIGUOUS' if decision == 'ambiguous' else
                                'WRONG_STATUS' if optimal or decision in ('optimal', 'unbounded') else
                                'SOLVER_FAILED')
    elif record['status'] != 'TIMEOUT' and optimal:
        record['status'] = 'PASS' if record.get('audit_pass') else 'AUDIT_FAILED'
    if point.exists():
        with point.open('rb') as src, gzip.open(str(point) + '.gz', 'wb') as dest:
            shutil.copyfileobj(src, dest)
        point.unlink()
    dump(artifact(prefix, '.json'), record)


def kill(process):
    if os.name == 'nt':
        subprocess.run(['taskkill', '/PID', str(process.pid), '/T', '/F'],
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                       check=False)
        if process.poll() is None:
            process.kill()
        process.wait()
        return
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
    with artifact(prefix, '.adapter.log').open('w', encoding='utf-8') as stream:
        p = subprocess.Popen(command, cwd=ROOT, env=env, stdout=stream,
                             stderr=subprocess.STDOUT, start_new_session=True)
        try:
            while p.poll() is None:
                clock = artifact(prefix, '.clock.json')
                if clock.exists():
                    try:
                        state = json.loads(clock.read_text(encoding='utf-8'))
                    except (json.JSONDecodeError, PermissionError):
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
        result = json.loads(result_path.read_text(encoding='utf-8'))
    else:
        result = json.loads(artifact(prefix, '.request.json').read_text(encoding='utf-8')) if artifact(prefix, '.request.json').exists() else {}
        result.update(status='TIMEOUT' if timed_out else 'ERROR',
                      wall_seconds=config['seconds'] if timed_out else None,
                      returncode=p.returncode)
    if artifact(prefix, '.model.json').exists():
        result.update(json.loads(artifact(prefix, '.model.json').read_text(encoding='utf-8')))
    result.update(instance=model['name'], solver=solver, repeat=repeat,
                  log=str(artifact(prefix, '.log').relative_to(config_path.parent)),
                  adapter_log=str(artifact(prefix, '.adapter.log').relative_to(config_path.parent)))
    if model.get('expected_status'):
        result['expected_status'] = model['expected_status']
    print(f"{model['name']} {solver} #{repeat}: {result['status']} "
          f"{result.get('wall_seconds')}s", flush=True)
    return result


def validate_report_run(metadata, results, allow_partial):
    """Reject stale, duplicate, or incomplete measurements before publication."""
    config = metadata['config']
    solvers = tuple(config.get('solvers', LEGACY_SOLVERS))
    if (not metadata.get('completed') and not allow_partial) or any(metadata.get(k, True) for k in (
            'source_changed_during_run', 'binary_changed_during_run',
            'library_changed_during_run', 'model_changed_during_run',
            'competitor_changed_during_run')):
        raise ValueError('Incomplete run or changed implementation/models; cannot publish')
    expected = {(model['name'], solver, repeat) for model in config['models']
                for solver in solvers for repeat in range(1, config['repeats'] + 1)}
    actual = [(r['instance'], r['solver'], r['repeat']) for r in results]
    if (len(actual) != len(set(actual)) or not set(actual) <= expected or
            (metadata.get('completed') and set(actual) != expected)):
        raise ValueError('Missing or duplicate runs; cannot publish')
    dataset = config.get('dataset', 'netlib')
    if dataset not in DATASETS:
        raise ValueError('Unknown benchmark dataset: ' + dataset)
    return dataset


def comparison_results(metadata, results, cache):
    """Include compatible records without requiring another solver run."""
    config, environment = metadata['config'], metadata['environment']
    solvers = tuple(config.get('solvers', LEGACY_SOLVERS))
    dataset = config.get('dataset', 'netlib')
    compared = list(results)
    compared_solvers = list(solvers)
    for solver in SOLVERS:
        if solver in solvers:
            continue
        saved = [cached_measurement(cache, config, environment, model, solver, repeat)
                 for model in config['models'] for repeat in range(1, config['repeats'] + 1)]
        # Policy skips are report cells, not measurements. Reconstruct them
        # only from the current GLPK policy and the compatible timeout/gate.
        if solver == 'glpk':
            prior = [row for row in saved if row]
            for index, row in enumerate(saved):
                if row is None:
                    model = config['models'][index // config['repeats']]
                    policy = dict(config, glpk_stop_on_timeout=dataset in ('pds', 'railway'))
                    if glpk_skip_reason(policy, solver, model, prior):
                        saved[index] = dict(instance=model['name'], solver=solver,
                                            repeat=index % config['repeats'] + 1,
                                            status='SKIPPED', wall_seconds=None)
        if saved and all(saved):
            compared.extend(saved)
            compared_solvers.append(solver)
    compared_config = dict(
        config,
        solvers=[solver for solver in SOLVERS if solver in compared_solvers],
        glpk_stop_on_timeout=dataset in ('pds', 'railway') and 'glpk' in compared_solvers,
    )
    view = dict(metadata, config=compared_config)
    return view, compared


def report(run_dir, docs_dir, allow_partial=False):
    metadata = json.loads((run_dir / 'metadata.json').read_text(encoding='utf-8'))
    results = json.loads((run_dir / 'results.json').read_text(encoding='utf-8'))
    config, environment = metadata['config'], metadata['environment']
    dataset = validate_report_run(metadata, results, allow_partial)
    ensure_reports(docs_dir)
    target = docs_dir / (dataset + '.md')
    cache = read_result_cache(target) or metadata.get('result_cache', [])
    cache = update_result_cache(cache, config, environment, results)
    stored = update_table_rows(read_table(target), config, results)
    view, compared = comparison_results(metadata, results, cache)
    blocks = {run_dir.name: run_summary(run_dir, view, compared)}
    atomic_text(target, render_report(dataset, stored, blocks, cache=cache))
    overview(docs_dir)


def measurement_key(config, environment, model, solver):
    """Changing our solver must not invalidate another solver's measurements."""
    version = environment.get('versions', {}).get(solver, {})
    executable = version.get('binding_sha256' if solver == 'highs' else 'binary_sha256')
    if not executable or not model.get('raw_sha256') or environment.get('cpu') in (None, '', 'unknown'):
        return None
    identity = dict(executable=executable, version=version.get('version_output', version.get('version')))
    if solver == 'lp_simplex':
        identity['libraries'] = sorted(set(version.get('solver_libraries_sha256', {}).values()))
        identity['source'] = {name: value for name, value in sources().items()
                              if name.startswith(('src/', 'include/')) or name == 'CMakeLists.txt'}
    adapter = '\n'.join(inspect.getsource(function) for function in (
        audit, parse_scip_solution, parse_clp_solution, parse_glpk_solution,
        glpk_input, infeasibility_status, run_cli, worker))
    settings = {key: config.get(key) for key in ('seconds', 'tolerance')}
    if solver == 'lp_simplex':
        settings['iterations'] = config.get('iterations')
    for key in {'lp_simplex': ('lp_algorithm',), 'highs': ('highs_algorithm',),
                'scip': ('scip_algorithm',)}.get(solver, ()):
        settings[key] = config.get(key)
    machine = {key: environment.get(key) for key in (
        'os', 'architecture', 'cpu', 'memory_bytes', 'logical_cpus', 'python', 'numpy', 'scipy')}
    machine['host'] = hashlib.sha256(platform.node().encode()).hexdigest()
    converters = [ROOT / 'data' / config.get('dataset', 'netlib') / 'prepare.py']
    converters += sorted((converters[0].parent / 'raw').glob('*.c'))
    value = dict(machine=machine, threads=THREADS, identity=identity, settings=settings,
                 model={key: model.get(key) for key in (
                     'name', 'dataset', 'raw_sha256', 'expected_mps_sha256', 'reference', 'expected_status')},
                 converters=[sha(path) for path in converters if path.is_file()],
                 adapter=hashlib.sha256(adapter.encode()).hexdigest())
    return hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()


def cached_measurement(entries, config, environment, model, solver, repeat):
    key = measurement_key(config, environment, model, solver)
    if key is None:
        return None
    for entry in entries:
        row = entry['result']
        if (entry['key'] == key and row.get('instance') == model['name']
                and row.get('solver') == solver and row.get('repeat') == repeat
                and row.get('status') in ('PASS', 'TIMEOUT', 'SOLVER_FAILED', 'AUDIT_FAILED', 'WRONG_STATUS', 'AMBIGUOUS')):
            elapsed = row.get('wall_seconds')
            if elapsed is None and row['status'] not in ('PASS', 'TIMEOUT'):
                return dict(row)
            if not isinstance(elapsed, bool) and isinstance(elapsed, (int, float)) and math.isfinite(elapsed) and elapsed >= 0.:
                return dict(row)
    return None


def update_result_cache(entries, config, environment, results):
    by_slot = {(entry['result'].get('instance'), entry['result'].get('solver'),
                entry['result'].get('repeat')): entry for entry in entries}
    models = {model['name']: model for model in config['models']}
    fields = ('instance', 'solver', 'repeat', 'status', 'wall_seconds', 'objective',
              'rows', 'columns', 'nonzeros', 'audit_pass', 'primal_residual',
              'objective_error', 'reference_error', 'solver_status', 'expected_status',
              'cli_wall_seconds', 'cli_cpu_seconds', 'cli_user_seconds', 'cli_system_seconds')
    for row in results:
        key = measurement_key(config, environment, models[row['instance']], row['solver'])
        if key is not None:
            by_slot[row['instance'], row['solver'], row['repeat']] = dict(
                key=key, result={field: row[field] for field in fields if field in row})
    return list(by_slot.values())


def adopt_report_results(path, config, environment):
    """One-time import of existing competitor measurements explicitly approved by the user."""
    if config['repeats'] != 1:
        raise ValueError('Existing report medians cannot reconstruct multiple repetitions')
    blocks = read_run_blocks(path)
    if not blocks:
        return []
    info = dict(re.findall(r'^\| ([^|\n]+) \| ([^|\n]+) \|$', list(blocks.values())[-1], re.M))
    expected_machine = f"{environment.get('os')} / {environment.get('cpu')}"
    if info.get('系统 / CPU') != expected_machine:
        raise ValueError('Existing report belongs to a different machine/environment')
    settings = info.get('重复 / 单次限时 / 容差 / 工作预算', '').split(' / ')
    if len(settings) != 4 or [float(settings[0]), float(settings[1].split()[0]), float(settings[2])] != [
            1., config['seconds'], config['tolerance']]:
        raise ValueError('Existing report has different test limits/settings')
    packages = ' / '.join([environment.get('python', '未知').split()[0],
                           environment.get('numpy', '未知'), environment.get('scipy', '未知')])
    if info.get('Python / NumPy / SciPy') != packages:
        raise ValueError('Existing report has a different Python/NumPy/SciPy environment')
    algorithms = info.get('算法（lp / scip / highs）', '').split(' / ')
    if len(algorithms) != 3 or algorithms[1:] != [config['scip_algorithm'], config['highs_algorithm']]:
        raise ValueError('Existing competitor algorithm settings differ')
    versions = environment.get('versions', {})
    for solver, field in [('scip', 'SCIP / LP 后端'), ('clp', 'Clp'), ('glpk', 'GLPK')]:
        actual = versions.get(solver, {}).get('version_output', '').splitlines()
        if not actual or info.get(field) != actual[0]:
            raise ValueError('Existing report has a different ' + solver + ' version')
    highs = versions.get('highs', {})
    if info.get('HiGHS') != f"{highs.get('version')} / {highs.get('binding')}":
        raise ValueError('Existing report has a different HiGHS version')
    stored, results = read_table(path), []
    for model in config['models']:
        if model['name'] not in stored:
            continue
        cells = [cell.strip() for cell in stored[model['name']].strip('|').split('|')]
        shape = cells[1].split(' / ')
        if len(shape) != 3 or any(str(model[field]) != value for field, value in
                zip(('rows', 'columns', 'nonzeros'), shape) if field in model):
            raise ValueError('Existing report has different model dimensions: ' + model['name'])
        if model.get('reference') is not None and (cells[2] == '—' or
                abs(float(cells[2]) - model['reference']) > config['tolerance'] * (1. + abs(model['reference']))):
            raise ValueError('Existing report has a different objective reference: ' + model['name'])
        for solver in SOLVERS[1:]:
            status, seconds, objective = [value.strip() for value in cells[3 + SOLVERS.index(solver)].split(' / ')]
            if status not in ('PASS', 'TIMEOUT', 'SKIPPED', 'SOLVER_FAILED', 'AUDIT_FAILED', 'WRONG_STATUS', 'AMBIGUOUS'):
                continue
            row = dict(instance=model['name'], solver=solver, repeat=1, status=status,
                       wall_seconds=float(seconds) if status == 'PASS' else config['seconds'] if status == 'TIMEOUT' else None)
            if objective != '—':
                row['objective'] = float(objective)
            results.append(row)
    return update_result_cache([], config, environment, results)


def raw_inputs(dataset):
    directory = ROOT / 'data' / ('netlib' if dataset == 'netlib-infeasible' else dataset) / 'raw'
    if dataset in ('netlib', 'netlib-infeasible'):
        category = 'infeasible' if dataset == 'netlib-infeasible' else 'feasible'
        return {p.name[:-7]: p for p in sorted((directory / category).glob('*.mps.gz'))}
    pattern = '*.bz2' if dataset == 'pds' else 'rail*.gz'
    return {p.stem: p for p in sorted(directory.glob(pattern))}


def ensure_reports(docs_dir):
    """Create missing tables without resetting other datasets' experiment rows."""
    docs_dir.mkdir(parents=True, exist_ok=True)
    for dataset in DATASETS:
        target = docs_dir / (dataset + '.md')
        if not target.exists():
            pending = ['| ' + ' | '.join([name, '— / — / —', '—', *(['NOT_RUN / — / —'] * len(SOLVERS))]) + ' |'
                       for name in raw_inputs(dataset)]
            if dataset == 'netlib-infeasible':
                pending = [line.replace('| — / — / — | — |', '| — / — / — | INFEASIBLE |') for line in pending]
            stored = {line.split('|')[1].strip(): line for line in pending}
            atomic_text(target, render_report(dataset, stored, {}))


@contextmanager
def prepared_model(model, run_dir):
    """Never reuse, overwrite or delete a user's pre-generated MPS."""
    with tempfile.TemporaryDirectory(prefix='mps-', dir=run_dir) as temporary:
        log = run_dir / 'logs' / (model['name'] + '-prepare.log')
        with log.open('w', encoding='utf-8') as stream:
            dataset = model.get('dataset', 'netlib')
            category = ('infeasible' if dataset == 'netlib-infeasible' else 'feasible' if dataset == 'netlib' else '')
            name = (category + '/' if category else '') + model['name']
            source_dataset = 'netlib' if dataset == 'netlib-infeasible' else dataset
            subprocess.run([sys.executable, str(ROOT / 'data' / source_dataset / 'prepare.py'),
                            name, '--output-dir', temporary],
                           stdout=stream, stderr=subprocess.STDOUT, check=True)
        path = Path(temporary) / category / (model['name'] + '.mps')
        digest = sha(path)
        if model.get('expected_mps_sha256') and digest != model['expected_mps_sha256']:
            raise ValueError('Generated MPS checksum differs from manifest: ' + model['name'])
        model.update(path=str(path.resolve()), sha256=digest)
        yield path


def glpk_skip_reason(config, solver, model, results):
    if config.get('dataset') not in ('pds', 'railway'):
        return None
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
    cache = metadata.get('result_cache', [])
    saved = {(solver, repeat): cached_measurement(cache, config, metadata['environment'], model, solver, repeat)
             for solver in solvers for repeat in range(1, config['repeats'] + 1)
             if solver not in config.get('rerun', [])} if cache else {}
    for slot, row in saved.items():
        if row is None:
            reason = glpk_skip_reason(config, slot[0], model, results)
            if reason:
                saved[slot] = dict(instance=model['name'], solver=slot[0], repeat=slot[1],
                                   status='SKIPPED', wall_seconds=None, reason=reason)
    if saved and len(saved) == len(solvers) * config['repeats'] and all(saved.values()):
        results.extend(saved.values())
        dump(run_dir / 'results.json', results)
        return
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
                            if saved.get((solver, repeat + 1)):
                                row = saved[solver, repeat + 1]
                            elif reason:
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
            with reference_file.open(encoding='utf-8') as stream:
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
        with manifest.open(encoding='utf-8') as stream:
            info = {row['instance']: row for row in csv.DictReader(stream)}
    models = []
    for name in names:
        row = info.get(name, {})
        model = dict(name=name, dataset=dataset, raw_path=available[name].relative_to(ROOT).as_posix(),
                     raw_sha256=sha(available[name]), reference=references.get(name))
        if dataset == 'netlib-infeasible':
            model['expected_status'] = 'infeasible'
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
    if dataset == 'netlib-infeasible':
        paths += [ROOT / 'data/netlib/SHA256SUMS', ROOT / 'data/netlib/README.md']
    reference_file = reference_path(dataset)
    if reference_file is not None:
        paths.append(reference_file)
    return {p.relative_to(ROOT).as_posix(): sha(p) for p in paths if p.exists()}


def checkpoint(run_dir, config, metadata, results, binary, scip, docs_dir):
    """Update progress only after temporary MPS cleanup and integrity checks."""
    detected = metadata['environment']
    metadata['source_changed_during_run'] = sources() != metadata['source_sha256']
    native = detected['versions'].get('lp_simplex', {})
    metadata['binary_changed_during_run'] = bool(native) and (
        not binary.is_file() or sha(binary) != native['binary_sha256'])
    metadata['library_changed_during_run'] = any(not Path(p).is_file() or sha(p) != value
        for p, value in native.get('solver_libraries_sha256', {}).items())
    metadata['model_changed_during_run'] |= any(not (ROOT / name).is_file() or sha(ROOT / name) != value
        for name, value in config['input_metadata_sha256'].items()) or any(not (ROOT / p['raw_path']).is_file() or
        sha(ROOT / p['raw_path']) != p['raw_sha256'] for p in config['models'])
    metadata['competitor_changed_during_run'] = any(
        not Path(version['binary']).is_file() or sha(version['binary']) != version['binary_sha256']
        for key, version in detected['versions'].items() if key in ('scip', 'clp', 'glpk'))
    binding = detected['versions'].get('highs', {})
    if binding:
        path = Path(binding['binding_file'])
        metadata['competitor_changed_during_run'] |= (
            not path.is_file() or sha(path) != binding['binding_sha256'])
    dump(run_dir / 'metadata.json', metadata)
    keys = sorted({k for row in results for k in row if k not in ('command', 'options')})
    with (run_dir / 'results.csv').open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=keys, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(results)
    report(run_dir, docs_dir, allow_partial=True)


def main():
    if len(sys.argv) > 1 and sys.argv[1] == '--worker':
        config_path, solver, name, prefix = sys.argv[2:]
        config = json.loads(Path(config_path).read_text(encoding='utf-8'))
        model = next(m for m in config['models'] if m['name'] == name)
        worker(config, solver, model, Path(prefix))
        return
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path, default=ROOT / 'build/data/lp_solve_mps')
    parser.add_argument('--scip', default='scip')
    parser.add_argument('--solvers', nargs='+', choices=SOLVERS, default=['lp_simplex'],
                        help='Solvers to execute (default: lp_simplex); unselected report cells are retained')
    parser.add_argument('--clp', default='clp')
    parser.add_argument('--glpk', default='glpsol')
    parser.add_argument('--glpk-gate-case', help='PDS/Railway only: test this case first; skip subsequent GLPK runs if it does not pass')
    parser.add_argument('--run-dir', type=Path, help='raw artifacts directory (default: tmp/benchmarks/<UTC timestamp>)')
    reporting = parser.add_mutually_exclusive_group()
    reporting.add_argument('--report-only', type=Path,
                           help='Rebuild a report from a saved raw run directory')
    reporting.add_argument('--refresh-reports', action='store_true',
                           help='Refresh bilingual reports from saved Markdown without running solvers')
    parser.add_argument('--rerun', nargs='+', choices=SOLVERS, default=[],
                        help='Explicitly retest these solvers even when compatible measurements exist')
    parser.add_argument('--adopt-existing', action='store_true',
                        help='One-time import of existing competitor report cells on this machine')
    parser.add_argument('--docs-dir', type=Path, default=ROOT / 'docs/benchmarks')
    parser.add_argument('--dataset', choices=DATASETS, default='netlib', help='Dataset (default: netlib)')
    parser.add_argument('--cases', nargs='+', help='Instance names; omitted means the selected dataset in full')
    parser.add_argument('--repeats', type=int, default=1)
    parser.add_argument('--seconds', type=float, default=None,
                        help='wall-clock limit per run (default: Netlib 60 s; Netlib infeasible 90 s; PDS/Railway 120 s)')
    parser.add_argument('--iterations', type=int, default=300000)
    parser.add_argument('--lp-algorithm', choices=LP_ALGORITHMS, default='auto')
    parser.add_argument('--highs-algorithm', choices=['choose', 'simplex', 'ipm'], default='choose')
    parser.add_argument('--scip-algorithm', choices=['default', 'dual'], default='default')
    args = parser.parse_args()
    if args.refresh_reports:
        overview(args.docs_dir.resolve())
        return
    if args.report_only:
        report(args.report_only.resolve(), args.docs_dir.resolve(), allow_partial=True)
        return
    if args.seconds is None:
        args.seconds = 90. if args.dataset == 'netlib-infeasible' else 60. if args.dataset == 'netlib' else 120.
    if args.repeats < 1 or not math.isfinite(args.seconds) or args.seconds <= 0 or args.iterations < 1:
        parser.error('repeats, seconds, and iterations must be positive')
    models = select_models(parser, args.cases, args.dataset)
    cache_path = args.docs_dir.resolve() / (args.dataset + '.md')
    cache = read_result_cache(cache_path)
    if args.glpk_gate_case:
        if args.dataset not in ('pds', 'railway'):
            parser.error('--glpk-gate-case is only supported for PDS and Railway')
        if 'glpk' not in args.solvers or args.glpk_gate_case not in {m['name'] for m in models}:
            parser.error('--glpk-gate-case requires GLPK and a case in the selected suite')
        models.sort(key=lambda model: model['name'] != args.glpk_gate_case)
    if len(args.solvers) != len(set(args.solvers)):
        parser.error('solvers must be unique')
    if not set(args.rerun) <= set(args.solvers):
        parser.error('--rerun requires those solvers to be selected by --solvers')
    native_cached = any(entry['result'].get('solver') == 'lp_simplex' for entry in cache)
    binary = args.binary.resolve() if 'lp_simplex' in args.solvers or native_cached and args.binary.is_file() else None
    if binary is not None and not binary.is_file():
        parser.error('lp_simplex executable not found: ' + str(binary))
    executables = {}
    for solver, requested in (('scip', args.scip), ('clp', args.clp), ('glpk', args.glpk)):
        known = any(entry['result'].get('solver') == solver for entry in cache)
        executable = shutil.which(requested) if solver in args.solvers or known or args.adopt_existing else None
        if solver in args.solvers and not executable:
            parser.error(solver + ' executable not found: ' + requested)
        executables[solver] = executable
    scip, clp, glpk = (executables[key] for key in ('scip', 'clp', 'glpk'))
    config = dict(dataset=args.dataset, binary=str(binary) if binary else None, scip=scip, clp=clp, glpk=glpk, solvers=args.solvers, repeats=args.repeats, seconds=args.seconds,
                  tolerance=1e-7, iterations=args.iterations, lp_algorithm=args.lp_algorithm,
                  highs_algorithm=args.highs_algorithm, scip_algorithm=args.scip_algorithm,
                  input_metadata_sha256=input_metadata(args.dataset),
                  suite_label=args.dataset + (' 全集' if not args.cases else ' 显式子集'),
                  models=models, glpk_stop_on_timeout=(args.dataset in ('pds', 'railway') and 'glpk' in args.solvers))
    config['rerun'] = args.rerun
    if args.glpk_gate_case:
        config['glpk_gate_case'] = args.glpk_gate_case
    env = {k: v for k, v in os.environ.items() if not k.startswith('LP_SIMPLEX_')}
    env.update(THREADS, LC_ALL='C', PYTHONDONTWRITEBYTECODE='1')
    os.environ.update(THREADS)
    detected = environment(binary, scip, clp, glpk)
    if args.adopt_existing:
        imported = adopt_report_results(cache_path, config, detected)
        slots = {(entry['result']['instance'], entry['result']['solver'], entry['result']['repeat']) for entry in cache}
        cache += [entry for entry in imported if (entry['result']['instance'], entry['result']['solver'], entry['result']['repeat']) not in slots]
    run_dir = (args.run_dir or ROOT / 'tmp/benchmarks' / dt.datetime.now(dt.timezone.utc).strftime('%Y%m%dT%H%M%SZ')).resolve()
    if run_dir.exists():
        parser.error('run directory exists; select a new directory')
    (run_dir / 'logs').mkdir(parents=True)
    config_path = run_dir / 'config.json'
    dump(config_path, config)
    metadata = dict(created_utc=dt.datetime.now(dt.timezone.utc).isoformat(), config=config,
                    invocation=sys.argv, environment=detected, source_sha256=sources(), completed=False,
                    result_cache=cache,
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
