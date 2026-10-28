"""Benchmark statistics, Markdown persistence, and bilingual report rendering.

This module does not prepare models or launch solvers. The CLI coordinates
measurements and supplies the current compatible records for publication.
"""
import json
import math
import os
from pathlib import Path
import re
import shlex
import statistics
import time


ROOT = Path(__file__).resolve().parents[1]
DATASETS = ('netlib', 'netlib-infeasible', 'pds', 'railway')
LEGACY_SOLVERS = ('lp_simplex', 'scip', 'highs')
SOLVERS = (*LEGACY_SOLVERS, 'clp', 'glpk')
LABELS = {'lp_simplex': 'lp_simplex', 'scip': 'SCIP / SoPlex', 'highs': 'HiGHS', 'clp': 'Clp', 'glpk': 'GLPK'}
TABLE_HEADER = ('| 实例 | 规模（行 / 列 / 非零元） | 参考目标 | ' +
                ' | '.join(f'{LABELS[k]}（状态 / 秒 / 目标）' for k in SOLVERS) + ' |\n' +
                '| --- | --- | ---: | ' + ' | '.join(['---'] * len(SOLVERS)) + ' |\n')
SUMMARY_HEADER = ('| 求解器 | 成功案例 | 超时次数 | 其他未通过次数 | 与 lp_simplex 共同成功数 | '
                  '共同案例几何平均/秒（lp_simplex / 本求解器） | 相对耗时（lp_simplex / 各求解器） |')
SUMMARY_SEPARATOR = '| --- | ---: | ---: | ---: | ---: | ---: | ---: |'


def atomic_text(path, text):
    temporary = path.with_name('.' + path.name + '.tmp')
    try:
        temporary.write_text(text, encoding='utf-8')
        for attempt in range(50):
            try:
                temporary.replace(path)
                break
            except PermissionError:
                # Windows readers can briefly deny replacement of an open file.
                if os.name != 'nt' or attempt == 49:
                    raise
                time.sleep(.01)
    finally:
        temporary.unlink(missing_ok=True)


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


def update_table_rows(stored, config, results):
    """Update selected result cells while preserving all other measurements."""
    solvers = tuple(config.get('solvers', LEGACY_SOLVERS))
    rows = aggregate(results, [v['name'] for v in config['models']], config['repeats'], config['seconds'], solvers)
    by_name = {model['name']: model for model in config['models']}
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
        reference = ('INFEASIBLE' if model.get('expected_status') == 'infeasible' else
                     f"{model['reference']:.12g}" if model.get('reference') is not None else '—')
        prior = stored.get(name, '')
        previous = [v.strip() for v in prior.strip('|').split('|')] if prior else []
        values = [name, ' / '.join(shape), reference,
                  *(cell(solver) if solver in solvers else
                    previous[3 + SOLVERS.index(solver)] if previous else 'NOT_RUN / — / —'
                    for solver in SOLVERS)]
        stored[name] = '| ' + ' | '.join(values) + ' |'
    return stored


def read_table(path):
    rows = {}
    if path.exists():
        inside, legacy = False, False
        for line in path.read_text(encoding='utf-8').splitlines():
            if line.startswith('| 实例 |') and ('规模' in line or '| 行 |' in line):
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
        path.read_text(encoding='utf-8'), re.S)}


def read_result_cache(path):
    """Keep reusable measurements in the report, independent of temporary logs."""
    if not path.exists():
        return []
    found = re.search(r'<!-- benchmark-cache:v1\n(.*?)\n-->',
                      path.read_text(encoding='utf-8'), re.S)
    if not found:
        return []
    try:
        entries = json.loads(found[1])
        if not isinstance(entries, list) or any(
                not isinstance(entry, dict) or not isinstance(entry.get('result'), dict)
                or not isinstance(entry.get('key'), str) for entry in entries):
            return []
        return entries
    except (ValueError, TypeError):
        return []


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
        common = set(solved.get('lp_simplex', {})) & set(solved[solver])
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
        ('对比求解器', ', '.join(LABELS[key] for key in solvers)),
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
        ('原始记录', '逐例测量与指纹保存在本报告中；临时文件已清理' if metadata.get('artifacts_removed') else
         os.path.relpath(run_dir, ROOT) + '（配置、环境、SHA-256、命令、日志、解和审核残差）')]
    for solver in ('clp', 'glpk'):
        if solver in solvers:
            info.append((LABELS[solver], versions.get(solver, {}).get('version_output', '未知').splitlines()[0]))
    if c.get('supplemental_run'):
        info.append(('补测来源', c['supplemental_run']))
    if 'clp' in solvers:
        info.append(('Clp / GLPK 算法', 'dual simplex / primal simplex；presolve on'))
    if 'glpk' in solvers and c.get('glpk_stop_on_timeout'):
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
    if c.get('dataset') == 'netlib-infeasible':
        lines += ['', '| 求解器 | 错误判定次数 | 判定不明确次数 | 运行或数值失败次数 | 跳过次数 |',
                  '| --- | ---: | ---: | ---: | ---: |']
        for solver in solvers:
            rows = [r for r in results if r['solver'] == solver]
            counts = [sum(r['status'] in statuses for r in rows) for statuses in (
                ('WRONG_STATUS',), ('AMBIGUOUS',), ('SOLVER_FAILED', 'ERROR', 'AUDIT_FAILED', 'PREPARE_FAILED'), ('SKIPPED',))]
            lines.append('| ' + ' | '.join([LABELS[solver], *map(str, counts)]) + ' |')
        unsuccessful = [r for r in results if r['status'] not in ('PASS', 'SKIPPED')]
        if unsuccessful:
            lines += ['', '| 实例 | 求解器 | 结果 | 终止信息 |', '| --- | --- | --- | --- |']
            for row in unsuccessful:
                lines.append('| ' + ' | '.join(clean(value) for value in (
                    row['instance'], LABELS[row['solver']], row['status'],
                    row.get('solver_status', row['status']))) + ' |')
    lines += ['', '复现本轮案例与参数（新建运行目录）：', '', '```sh',
              reproduction_script(reproduction_command(c)), '```']
    return '\n'.join(lines) + '\n\n'


def format_reproduction_command(command):
    """Use repository-relative paths and PATH names for external solvers."""
    command = list(command)
    for index in range(1, len(command)):
        option, value = command[index - 1], command[index]
        if option not in ('--binary', '--scip', '--clp', '--glpk'):
            continue
        path = Path(value)
        # Saved reports may contain POSIX absolute paths even on Windows.
        if not (path.is_absolute() or value.startswith('/')):
            continue
        if option != '--binary' and not path.is_relative_to(ROOT):
            command[index] = path.name
        else:
            relative = Path(os.path.relpath(path, ROOT)).as_posix()
            command[index] = relative if '/' in relative else './' + relative
    return shlex.join(command)


def reproduction_command(config):
    command = ['python3', 'data/benchmark.py', '--dataset', config.get('dataset', 'netlib'),
               '--cases', *(model['name'] for model in config['models'])]
    for option, key in (('--binary', 'binary'), ('--scip', 'scip'), ('--clp', 'clp'), ('--glpk', 'glpk'), ('--repeats', 'repeats'),
                        ('--seconds', 'seconds'), ('--iterations', 'iterations'),
                        ('--lp-algorithm', 'lp_algorithm'), ('--scip-algorithm', 'scip_algorithm'),
                        ('--highs-algorithm', 'highs_algorithm')):
        if config.get(key) is not None:
            command.extend([option, str(config[key])])
    command.extend(['--solvers', *config.get('solvers', LEGACY_SOLVERS)])
    if config.get('glpk_gate_case'):
        command.extend(['--glpk-gate-case', config['glpk_gate_case']])
    if config.get('rerun'):
        command.extend(['--rerun', *config['rerun']])
    return format_reproduction_command(command)


def reproduction_script(command):
    """Build the driver from the repository root before running the benchmark."""
    args = shlex.split(command)
    binary = 'build/data/lp_solve_mps'
    if '--binary' in args:
        args[args.index('--binary') + 1] = binary
    else:
        args.extend(['--binary', binary])
    return '\n\n'.join([
        'cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DLP_SIMPLEX_BUILD_BENCHMARKS=ON',
        'cmake --build build --target lp_solve_mps -j 4',
        format_reproduction_command(args)])


def refresh_reproduction_script(match):
    commands = [line for line in match[1].splitlines()
                if line.startswith('python3 data/benchmark.py ') and '--dataset ' in line]
    if len(commands) != 1:
        return match[0]
    return '```sh\n' + reproduction_script(commands[0]) + '\n```'


def ordered_cases(dataset, names):
    pattern = {'pds': r'pds-(\d+)', 'railway': r'rail(\d+)'}.get(dataset)
    def key(name):
        match = re.fullmatch(pattern, name) if pattern else None
        return (0, int(match[1]), name) if match else (1, 0, name)
    return sorted(names, key=key)


def render_report(dataset, stored, blocks, english=False, cache=None):
    def tr(chinese, translated):
        return translated if english else chinese
    title = {'netlib': 'Netlib feasible', 'netlib-infeasible': 'Netlib infeasible', 'pds': 'PDS',
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
        tr('- 同一机器、求解器、输入与参数不变时直接复用逐例测量；仅补测缺失或发生变化的项。对比统计按兼容的完整记录重新计算。',
           '- Per-instance measurements are reused when the machine, solver, inputs, and settings match. Only missing or changed entries are tested; comparisons are recomputed from compatible complete records.'), '',
        tr('## 测试信息与统计', '## Test information and statistics'), '',
        tr('展示当前版本的逐例结果及对比统计。',
           'Current per-instance results and comparison statistics.'), '']
    if dataset == 'netlib-infeasible':
        lines = [line for line in lines if not line.startswith(('- PASS ', '- PASS requires'))]
        insertion = lines.index(tr('## 测试信息与统计', '## Test information and statistics'))
        lines[insertion:insertion] = [tr('- PASS 表示明确报告不可行且在限时内结束；预期状态来自 Netlib infeasible 分类。未独立验证 Farkas 证书，PASS 不代表证书审核通过。',
                     '- PASS means an explicit infeasible status within the time limit, against the Netlib infeasible classification. Farkas certificates are not independently verified; PASS is a status check, not a certificate audit.'),
                  tr('- WRONG_STATUS 表示报告最优或无界；AMBIGUOUS 表示不可行或无界；SOLVER_FAILED / ERROR 表示未得到明确判定或运行失败。TIMEOUT 优先于其他状态。',
                     '- WRONG_STATUS means optimal or unbounded; AMBIGUOUS means infeasible or unbounded; SOLVER_FAILED / ERROR means no definitive decision or an execution failure. TIMEOUT takes precedence.'),
                  tr('- 原始 MPS 按原样交给求解器。greenbea 含多个 BOUNDS 向量及重复界，pang 存在读取器兼容性差异；相关失败包含输入兼容性影响。HiGHS 的读取警告保存在 reader_status，实际求解日志另行保留。',
                     '- Original MPS files are passed to the solvers unchanged. greenbea contains multiple BOUNDS vectors and duplicate bounds; pang also exposes reader compatibility differences. Related failures include input compatibility effects. HiGHS read warnings are recorded in reader_status; solver logs are retained separately.'), '']
    if not blocks:
        empty = ['尚未运行；环境、版本、算法、限时和重复次数尚未记录。', '',
                 SUMMARY_HEADER, SUMMARY_SEPARATOR]
        empty.extend(f"| {LABELS[solver]} | 0/{len(stored)} | 0 | 0 | "
                     + ('— | — | — |' if solver == 'lp_simplex' else '0 | — / — | — |')
                     for solver in SOLVERS)
        lines += [english_summary('\n'.join(empty)) if english else '\n'.join(empty), '']
    # Older reports contain multiple run blocks. Keep only the last saved run.
    for name, block in list(blocks.items())[-1:]:
        if dataset in ('netlib', 'netlib-infeasible'):
            block = re.sub(r'^\| GLPK 跳过策略 \|[^\n]*\n', '', block, flags=re.M)
        block = re.sub(r'^(\| 案例范围 \| [^；\n]*；)([^|\n]+)\|$',
                       lambda match: match[1] + ', '.join(ordered_cases(
                           dataset, [case.strip() for case in match[2].split(',')])) + ' |',
                       block, flags=re.M)
        # Upgrade saved instructions too, without requiring original run logs.
        block = re.sub(r'```sh\n(.*?)\n```', refresh_reproduction_script,
                       block, flags=re.S)
        lines += [f'<!-- benchmark-run:{name} -->',
                  english_summary(block.rstrip()) if english else block.rstrip(),
                  '<!-- /benchmark-run -->', '']
    header = ('| Instance | Size (rows / columns / nonzeros) | Reference objective | ' +
              ' | '.join(f'{LABELS[k]} (status / s / objective)' for k in SOLVERS) + ' |\n' +
              TABLE_HEADER.splitlines()[1]) if english else TABLE_HEADER.rstrip()
    if dataset == 'netlib-infeasible':
        header = header.replace('参考目标', '预期状态').replace('Reference objective', 'Expected status')
    lines += [tr('## 逐例最新结果', '## Latest per-instance results'), '', header]
    lines.extend(stored[name] for name in ordered_cases(dataset, stored))
    if cache:
        lines += ['', '<!-- benchmark-cache:v1',
                  json.dumps(cache, ensure_ascii=False, allow_nan=False, separators=(',', ':')), '-->']
    return '\n'.join(lines) + '\n'


def english_summary(summary):
    """Localize saved summary prose without changing commands or numeric cells."""
    labels = {
        '测试信息': 'Test information', '值': 'Value', '参与求解器': 'Solvers in this run',
        '对比求解器': 'Compared solvers',
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
        '补测来源': 'Supplemental run',
        '错误判定次数': 'Wrong decisions', '判定不明确次数': 'Ambiguous decisions',
        '运行或数值失败次数': 'Execution or numerical failures', '跳过次数': 'Skipped runs',
        '实例': 'Instance', '结果': 'Result', '终止信息': 'Termination information',
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
            if text == '逐例测量与指纹保存在本报告中；临时文件已清理':
                return 'Per-instance measurements and fingerprints are stored in this report; temporary files removed'
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
        elif line == '展示当前版本的逐例结果及对比统计。':
            lines.append('Current per-instance results and comparison statistics.')
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
            if '## 测试信息与统计\n\n' not in path.read_text(encoding='utf-8'):
                raise ValueError('Dataset report has no summary section: ' + str(path))
            stored, blocks = read_table(path), read_run_blocks(path)
            cache = read_result_cache(path)
            report_text = render_report(dataset, stored, blocks, cache=cache)
            atomic_text(path, report_text)
            summaries[dataset] = report_text.partition('## 测试信息与统计\n\n')[2].partition('## 逐例最新结果')[0].strip()
            atomic_text(docs_dir / (dataset + '_en.md'),
                        render_report(dataset, stored, blocks, english=True, cache=cache))
    for english in (False, True):
        def tr(chinese, translated):
            return translated if english else chinese
        navigation = tr(
            '简体中文 | [English](./README_en.md) | [项目主页](../../README.md) | [文档](../README.md) | [历史性能](./history/README.md)',
            '[简体中文](./README.md) | English | [Project home](../../README_en.md) | [Documentation](../README_en.md) | [Performance history](./history/README_en.md)')
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
            tr('可行问题的 PASS 要求求解器报告最优且原模型可行性、目标复算通过；有参考目标时额外核对。'
               '缺少参考值不代表独立最优性认证。进行中或子集结果不构成全集结论，'
               '同一机器、求解器、输入和参数相同的测量直接复用，并按共同成功案例重算比较。',
               'For feasible suites, PASS requires an optimal solver status and successful original-model feasibility and objective audits. Reference objectives are checked when available. A pass without a reference is not an independent optimality certificate. Partial or subset results do not establish full-suite conclusions. Measurements with matching machine, solver, inputs, and settings are reused, and comparisons are recomputed on the common solved cases.'), '']
        lines += [tr('Netlib infeasible 单独统计：PASS 只核对明确的不可行状态，不代表独立不可行证书验证。',
                     'Netlib infeasible is summarized separately: PASS checks an explicit infeasible status and does not constitute independent certificate verification.'), '']
        for dataset, title in (('netlib', 'Netlib feasible'), ('netlib-infeasible', 'Netlib infeasible'), ('pds', 'PDS'),
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
            '```sh', 'cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DLP_SIMPLEX_BUILD_BENCHMARKS=ON',
            '', 'cmake --build build --target lp_solve_mps -j 4', '',
            'python3 data/benchmark.py --dataset netlib --binary build/data/lp_solve_mps', '',
            'python3 data/benchmark.py --dataset pds --binary build/data/lp_solve_mps', '',
            'python3 data/benchmark.py --dataset railway --binary build/data/lp_solve_mps', '```', '',
            tr('省略 `--cases` 测试所选数据集全集；短程验证可加 `--cases <实例名> --repeats 1 --seconds 10`。'
               'simple_examples 仅用于功能测试，不参与基准。上方各轮的复现命令保留其案例与具体参数，'
               '二进制路径可换为新构建结果；复核版本及 SHA-256 后再比较耗时。',
               'Omit `--cases` to test the full selected dataset. For a short check, add `--cases <instance> --repeats 1 --seconds 10`. The simple_examples models are functional tests, not benchmarks. Commands above preserve each run\'s cases and parameters. You may replace binary paths with new builds; verify versions and SHA-256 before comparing timings.'), '',
            tr('默认每轮新建 `tmp/benchmarks/<UTC 时间>/`，保存配置、环境、命令、结果、日志与审核信息。'
               '正常结束、异常、Ctrl+C 或 SIGTERM 清理临时 MPS；SIGKILL/断电残留可随运行目录删除。',
               'Each run creates `tmp/benchmarks/<UTC timestamp>/` for configuration, environment, commands, results, logs, and audit information. Temporary MPS files are cleaned on normal exit, exceptions, Ctrl+C, or SIGTERM. Files left by SIGKILL or power loss can be removed with the run directory.'), '',
            tr('直接从已保存的 Markdown 更新报告和汇总，不需要临时日志或再次求解：',
               'Refresh reports and summaries from saved Markdown without temporary logs or solver runs:'),
            '', '```sh', 'python3 data/benchmark.py --refresh-reports', '```', '',
            tr('如仍保留原始运行目录，也可从中重建报告：',
               'A retained raw run directory can also be used to rebuild a report:'), '', '```sh',
            tr('python3 data/benchmark.py --report-only tmp/benchmarks/<运行目录>',
               'python3 data/benchmark.py --report-only tmp/benchmarks/<run-directory>'), '```', '',
            tr('`--docs-dir` 可指定报告目录；重建汇总会保留其他数据集已保存的测试信息和统计。'
               '更多参数与审核说明见 [数据与基准说明](../../data/README.md)。',
               'Use `--docs-dir` to choose the report directory. Rebuilding summaries preserves saved test information and statistics for other datasets. See [data and benchmark instructions](../../data/README_en.md) for additional parameters and audit details.'), '']
        atomic_text(docs_dir / ('README_en.md' if english else 'README.md'), '\n'.join(lines))
