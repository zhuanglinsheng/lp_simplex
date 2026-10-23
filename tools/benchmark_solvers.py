#!/usr/bin/env python3
"""Serial LP comparison and generated reports. Run --help for reproducible options."""
import argparse
import csv
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
import signal
import statistics
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
SOLVERS = ('lp_simplex', 'scip', 'highs')
LABELS = {'lp_simplex': 'lp_simplex', 'scip': 'SCIP / SoPlex', 'highs': 'HiGHS'}
THREADS = dict(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
               VECLIB_MAXIMUM_THREADS='1', MKL_NUM_THREADS='1',
               BLIS_NUM_THREADS='1', NUMEXPR_NUM_THREADS='1')


def dump(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + '\n')


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
             *ROOT.glob('include/**/*.h'), *ROOT.glob('tools/*.py'),
             *ROOT.glob('tools/*.c')]
    return {str(p.relative_to(ROOT)): sha(p) for p in sorted(paths)}


def environment(binary, scip):
    from benchmark_highs import load_highs
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
    # test_netlib links the actual shared solver. Record it, not just the launcher.
    libraries = {}
    for name in re.findall(r'(/[^\s()]*liblp_simplex[^\s()]*)', linkage):
        p = Path(name)
        if p.is_file():
            libraries[str(p)] = sha(p)
    for p in binary.parent.parent.glob('liblp_simplex.*'):
        if p.is_file():
            libraries[str(p.resolve())] = sha(p)
    cache = binary.parent.parent / 'CMakeCache.txt'
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
             abs(objective - reference) <= tolerance * (1. + abs(reference)))
    return dict(audit_pass=valid, primal_residual=residual,
                objective=objective if math.isfinite(objective) else None,
                reference_error=abs(objective - reference) if math.isfinite(objective) else None)


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


def worker(config, solver, model_info, prefix):
    from benchmark_highs import load_highs
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
                       '--expected', str(model_info['reference']),
                       '--save-point', str(point), str(model)]
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


def aggregate(results, names, repeats, seconds):
    rows = []
    for name in names:
        row = dict(instance=name)
        for solver in SOLVERS:
            runs = [r for r in results if r['instance'] == name and r['solver'] == solver]
            passed = len(runs) == repeats and all(r['status'] == 'PASS' for r in runs)
            row[solver] = dict(passed=passed,
                               status='PASS' if passed else '/'.join(sorted({r['status'] for r in runs})),
                               median=statistics.median(r['wall_seconds'] for r in runs) if passed else None,
                               penalized=statistics.median(r['wall_seconds'] if r['status'] == 'PASS'
                                                          else 10 * seconds for r in runs))
        rows.append(row)
    return rows


def report(run_dir, docs_dir):
    m = json.loads((run_dir / 'metadata.json').read_text())
    results = json.loads((run_dir / 'results.json').read_text())
    c, e = m['config'], m['environment']
    if not m.get('completed') or any(m.get(k, True) for k in (
            'source_changed_during_run', 'binary_changed_during_run',
            'library_changed_during_run', 'model_changed_during_run',
            'competitor_changed_during_run')):
        raise ValueError('Incomplete run or changed implementation/models; cannot publish')
    expected = {(model['name'], solver, repeat) for model in c['models']
                for solver in SOLVERS for repeat in range(1, c['repeats'] + 1)}
    actual = [(r['instance'], r['solver'], r['repeat']) for r in results]
    if len(actual) != len(expected) or set(actual) != expected:
        raise ValueError('Missing or duplicate runs; cannot publish')
    rows = aggregate(results, [v['name'] for v in c['models']], c['repeats'], c['seconds'])
    docs_dir.mkdir(parents=True, exist_ok=True)
    link = Path(os.path.relpath(run_dir, docs_dir)).as_posix()
    def cell(value):
        return f"{value['median']:.3f}" if value['passed'] else value['status']
    scip_version = e['versions']['scip']['version_output'].splitlines()[0]
    lines = ['# 求解器性能对比', '', '> 由 `tools/benchmark_solvers.py` 自动生成，请勿手工修改结果。', '',
             f"运行：`{run_dir.name}`；范围：**{c['suite_label']}**，共 {len(rows)} 个可行 LP，每项重复 {c['repeats']} 次。", '',
             '## 环境与配置', '', '| 项目 | 自动检测值 |', '| --- | --- |',
             f"| 系统 / CPU | {e['os']} / {e['cpu']} |",
             f"| 内存 / 逻辑 CPU | {float(e['memory_bytes']) / 2**30:.1f} GiB / {e['logical_cpus']} |" if str(e['memory_bytes']).isdigit() else f"| 内存 / 逻辑 CPU | 未知 / {e['logical_cpus']} |",
             f"| Python / NumPy / SciPy | {e['python'].split()[0]} / {e['numpy']} / {e['scipy']} |",
             f"| lp_simplex | `{e['versions']['lp_simplex']['git_head'][:12]}`；工作区{'有改动' if e['versions']['lp_simplex']['git_status'] else '干净'}；{c['lp_algorithm']} |",
             f"| 编译器 | {e['versions']['lp_simplex']['compiler']} |",
             f"| 构建 | {e['versions']['lp_simplex']['cmake_cache'].get('CMAKE_BUILD_TYPE:STRING', 'unknown')}；依赖与链接库详见 metadata.json |",
             f"| SCIP / LP 后端 | {scip_version}；{c['scip_algorithm']} |",
             f"| HiGHS | {e['versions']['highs']['version']}；{e['versions']['highs']['binding']}；{c['highs_algorithm']} |",
             f"| 限制 | 1 线程；每次 {c['seconds']:g} 秒；容差 {c['tolerance']:g}；lp_simplex 工作预算 {c['iterations']} |", '',
             '## 汇总', '', '| 求解器 | 全部重复通过的实例 | 超时次数 | 其他未通过次数 | 惩罚移位几何平均 / 秒 |', '| --- | ---: | ---: | ---: | ---: |']
    for solver in SOLVERS:
        selected = [r for r in results if r['solver'] == solver]
        geom = math.expm1(statistics.mean(math.log1p(row[solver]['penalized']) for row in rows))
        lines.append(f"| {LABELS[solver]} | {sum(row[solver]['passed'] for row in rows)}/{len(rows)} | "
                     f"{sum(r['status'] == 'TIMEOUT' for r in selected)} | "
                     f"{sum(r['status'] not in ('PASS', 'TIMEOUT') for r in selected)} | {geom:.3f} |")
    lines += ['', '## 逐例结果', '', '时间为读取、求解及结果导出的墙钟秒数中位数；审计与 Python 适配器初始化不计入。', '',
              '| 实例 | 行 / 列 / 非零元 | lp_simplex | SCIP / SoPlex | HiGHS |', '| --- | ---: | ---: | ---: | ---: |']
    for row in rows:
        dims = next((r for r in results if r['instance'] == row['instance'] and 'rows' in r), {})
        shape = ' / '.join(str(dims.get(k, '?')) for k in ('rows', 'columns', 'nonzeros'))
        lines.append(f"| {row['instance']} | {shape} | " + ' | '.join(cell(row[s]) for s in SOLVERS) + ' |')
    lines += ['', 'PASS = 求解器报告最优，且独立原模型可行性、目标复算与参考目标检查通过；不等同于完整对偶证书。',
              '汇总对未通过的每次运行计入 10 倍限时惩罚，再按实例取中位数，以 1 秒移位计算几何平均；越低越好。',
              '单次或子集结果仅反映本轮观测，不代表全量排名。', '',
              f"原始记录：[环境与配置]({link}/metadata.json)、[结果 CSV]({link}/results.csv)、[逐次结果与日志索引]({link}/results.json)。", '',
              '[测试规范与复现](methodology.md)', '']
    (docs_dir / 'README.md').write_text('\n'.join(lines))
    method = f'''# 性能测试规范

由 `tools/benchmark_solvers.py` 维护。仅比较连续可行 LP；不可行与无界识别另作正确性测试。

| 项目 | 统一口径 |
| --- | --- |
| 数据集 | 默认 `data/netlib/feasible/*.mps` 全集；`--cases` 为显式子集，报告必须标明范围 |
| 算法 | 默认 lp_simplex AUTO、SCIP 默认 LP 配置、HiGHS choose；可显式指定并记录，不能混用为同一组排名 |
| 资源 | 逐题串行、求解器轮换顺序、单线程；每次读取原始模型并冷启动；不自动限制内存，不保证机器独占 |
| 计时 | 读取原模型至结果导出的墙钟时间；包含 CLI 进程启动；HiGHS 使用 Python binding，公共初始化及独立审核排除。CLI 与 binding 开销有差异，小题慎比 |
| 限时 | 父进程对上述计时区间统一强制限时，终止整个进程组；SCIP/HiGHS 另设相同内部限时 |
| 精度 | 默认原始／对偶容差 1e-7；lp_simplex 当前固定使用此默认精度。参考目标采用 `tol * (1 + abs(reference))` |
| 审核 | HiGHS 仅作为统一 MPS 读取器，用 NumPy longdouble 复算原始行／界违反与目标（有效精度取决于平台，部分平台等同 double）；报告最优且最大绝对违反不超过 tol、报告目标与参考目标均匹配才记 PASS |
| 证据边界 | 原始审核加求解器最优状态不构成完整独立最优性证明；历史参考文件只用于目标核对，不使用其计时 |
| 统计 | 默认重复 3 次；实例全部重复 PASS 才列通过耗时，否则显示状态；汇总按每次失败 10 倍限时惩罚、逐例取中位数、1 秒移位几何平均 |
| 追溯 | 每轮新建目录，保存模型／源码／二进制 SHA-256、Git 状态、构建配置、版本、命令、结果、日志及压缩解；运行后复核源码和二进制是否变化 |

## 复现

先完成 Release 构建，再运行（需要 NumPy、SciPy、SCIP，以及 highspy 或 SciPy 自带 HiGHS binding）：

```sh
python3 tools/benchmark_solvers.py \\
  --binary build/pan-controller-release/tools/test_netlib
```

短程子集验证加 `--cases degen2 qap12 dfl001 --repeats 1 --seconds 60`。
默认按 UTC 时间自动命名运行目录，并生成本目录两份文档；可用 `--docs-dir` 改变目的地。输出目录必须不存在，旧轮次不覆盖。
仅重建文档：

```sh
python3 tools/benchmark_solvers.py --report-only {os.path.relpath(run_dir, ROOT)}
```

所有机器、版本及配置以每轮 `metadata.json` 自动探测值为准；缺失值显示未知。
默认不安装依赖、不构建或修改求解器。异常中断会保留已完成结果，但不更新主页。
'''
    (docs_dir / 'methodology.md').write_text(method)


def main():
    if len(sys.argv) > 1 and sys.argv[1] == '--worker':
        config_path, solver, name, prefix = sys.argv[2:]
        config = json.loads(Path(config_path).read_text())
        model = next(m for m in config['models'] if m['name'] == name)
        worker(config, solver, model, Path(prefix))
        return
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path, default=ROOT / 'build/pan-controller-release/tools/test_netlib')
    parser.add_argument('--scip', default='scip')
    parser.add_argument('--run-dir', type=Path)
    parser.add_argument('--report-only', type=Path)
    parser.add_argument('--docs-dir', type=Path, default=ROOT / 'docs/benchmarks')
    parser.add_argument('--cases', nargs='+', help='Netlib feasible names; omitted means the full suite')
    parser.add_argument('--repeats', type=int, default=3)
    parser.add_argument('--seconds', type=float, default=120.)
    parser.add_argument('--iterations', type=int, default=300000)
    parser.add_argument('--lp-algorithm', choices=['auto', 'dual-revised', 'pan-bda'], default='auto')
    parser.add_argument('--highs-algorithm', choices=['choose', 'simplex', 'ipm'], default='choose')
    parser.add_argument('--scip-algorithm', choices=['default', 'dual'], default='default')
    args = parser.parse_args()
    if args.report_only:
        report(args.report_only.resolve(), args.docs_dir.resolve())
        return
    if args.repeats < 1 or not math.isfinite(args.seconds) or args.seconds <= 0 or args.iterations < 1:
        parser.error('repeats, seconds, and iterations must be positive')
    binary = args.binary.resolve()
    scip = shutil.which(args.scip)
    if not binary.is_file() or not scip:
        parser.error('lp_simplex binary and SCIP executable are required')
    references = {}
    with (ROOT / 'data/netlib/feasible_gurobi_1e-8.csv').open() as stream:
        for r in csv.DictReader(stream):
            if r['status'] == 'OPTIMAL':
                references[r['name']] = float(r['pobj'])
    models = sorted((ROOT / 'data/netlib/feasible').glob('*.mps'))
    if args.cases:
        if len(set(args.cases)) != len(args.cases):
            parser.error('duplicate cases')
        known = {p.stem: p for p in models}
        if set(args.cases) - known.keys():
            parser.error('unknown cases: ' + str(set(args.cases) - known.keys()))
        models = [known[name] for name in args.cases]
    for p in models:
        if p.stem not in references or not math.isfinite(references[p.stem]):
            parser.error('missing finite reference objective: ' + p.stem)
    config = dict(binary=str(binary), scip=scip, repeats=args.repeats, seconds=args.seconds,
                  tolerance=1e-7, iterations=args.iterations, lp_algorithm=args.lp_algorithm,
                  highs_algorithm=args.highs_algorithm, scip_algorithm=args.scip_algorithm,
                  reference_sha256=sha(ROOT / 'data/netlib/feasible_gurobi_1e-8.csv'),
                  suite_label='Netlib feasible 全集' if not args.cases else 'Netlib feasible 显式子集',
                  models=[dict(name=p.stem, path=str(p.relative_to(ROOT)), sha256=sha(p),
                               reference=references[p.stem]) for p in models])
    env = {k: v for k, v in os.environ.items() if not k.startswith('LP_SIMPLEX_')}
    env.update(THREADS, LC_ALL='C')
    os.environ.update(THREADS)
    detected = environment(binary, scip)
    run_dir = (args.run_dir or ROOT / 'docs/benchmarks/runs' / dt.datetime.now(dt.timezone.utc).strftime('%Y%m%dT%H%M%SZ')).resolve()
    if run_dir.exists():
        parser.error('run directory exists; select a new directory')
    (run_dir / 'logs').mkdir(parents=True)
    config_path = run_dir / 'config.json'
    dump(config_path, config)
    metadata = dict(created_utc=dt.datetime.now(dt.timezone.utc).isoformat(), config=config,
                    invocation=sys.argv, environment=detected, source_sha256=sources(), completed=False)
    dump(run_dir / 'metadata.json', metadata)
    results = []
    for case_index, model in enumerate(config['models']):
        for repeat in range(args.repeats):
            offset = (case_index + repeat) % len(SOLVERS)
            for solver in SOLVERS[offset:] + SOLVERS[:offset]:
                results.append(run_one(config_path, config, solver, model, repeat + 1, env))
                dump(run_dir / 'results.json', results)
    metadata['source_changed_during_run'] = sources() != metadata['source_sha256']
    metadata['binary_changed_during_run'] = sha(binary) != detected['versions']['lp_simplex']['binary_sha256']
    metadata['library_changed_during_run'] = any(not Path(p).exists() or sha(p) != value
        for p, value in detected['versions']['lp_simplex']['solver_libraries_sha256'].items())
    metadata['model_changed_during_run'] = any(sha(ROOT / p['path']) != p['sha256'] for p in config['models'])
    metadata['competitor_changed_during_run'] = (
        sha(scip) != detected['versions']['scip']['binary_sha256'] or
        sha(detected['versions']['highs']['binding_file']) != detected['versions']['highs']['binding_sha256'])
    metadata['completed'] = True
    dump(run_dir / 'metadata.json', metadata)
    keys = sorted({k for row in results for k in row if k not in ('command', 'options')})
    with (run_dir / 'results.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=keys, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(results)
    if any(metadata[k] for k in ('source_changed_during_run', 'binary_changed_during_run', 'library_changed_during_run',
                                  'model_changed_during_run', 'competitor_changed_during_run')):
        raise SystemExit('Implementation changed during run; raw results saved, report not published')
    report(run_dir, args.docs_dir.resolve())
    print(f'Report: {args.docs_dir / "README.md"}', flush=True)


if __name__ == '__main__':
    main()
