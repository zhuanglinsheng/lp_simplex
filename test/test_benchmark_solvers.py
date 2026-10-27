"""Regression checks for audit failures and benchmark accounting."""
import argparse
import bz2
import gzip
import json
import shutil
import struct
import sys
import math
import re
import shlex
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'data'))
import benchmark
import compare_performance
from benchmark import aggregate, artifact, audit, parse_scip_solution, report


class CliTimingTests(unittest.TestCase):
    def test_msvc_version_uses_stderr_banner(self):
        process = SimpleNamespace(stdout='C/C++ compiler options\n',
                                  stderr='Microsoft C/C++ Compiler Version 19.50\nCopyright\n')
        with patch.object(benchmark.subprocess, 'run', return_value=process):
            self.assertEqual(benchmark.compiler_version('cl.exe'),
                             'Microsoft C/C++ Compiler Version 19.50')

    def test_atomic_write_uses_utf8(self):
        with TemporaryDirectory() as folder:
            target = Path(folder) / 'metadata.json'
            benchmark.dump(target, {'state': '已完成'})
            self.assertEqual(json.loads(target.read_bytes().decode('utf-8')),
                             {'state': '已完成'})

    def test_atomic_write_retries_transient_windows_sharing_violation(self):
        with TemporaryDirectory() as folder:
            target = Path(folder) / 'clock.json'
            with patch.object(benchmark.os, 'name', 'nt'), \
                 patch.object(Path, 'replace', side_effect=[PermissionError(), None]) as replace, \
                 patch.object(benchmark.time, 'sleep') as sleep:
                benchmark.atomic_text(target, '{}')
            self.assertEqual(replace.call_count, 2)
            sleep.assert_called_once_with(.01)

    def test_windows_timeout_terminates_worker_tree(self):
        process = unittest.mock.Mock(pid=123)
        process.poll.return_value = 0
        with patch.object(benchmark.os, 'name', 'nt'), \
             patch.object(benchmark.subprocess, 'run') as run:
            benchmark.kill(process)
        self.assertEqual(run.call_args.args[0], ['taskkill', '/PID', '123', '/T', '/F'])
        process.wait.assert_called_once()
        process.kill.assert_not_called()

    def test_cli_boundary_and_failed_exit_keep_measured_cost(self):
        process = SimpleNamespace(returncode=7)
        with patch.object(benchmark, 'child_cpu', side_effect=[(1., 2.), (4., 2.5)]), \
             patch.object(benchmark.time, 'monotonic', side_effect=[10., 12.]), \
             patch.object(benchmark.subprocess, 'run', return_value=process):
            actual, timing = benchmark.run_cli(['solver'], None)
        self.assertIs(actual, process)
        self.assertEqual(timing, dict(cli_wall_seconds=2., cli_cpu_seconds=3.5,
                                      cli_user_seconds=3., cli_system_seconds=.5))

    def test_missing_or_invalid_cpu_is_unknown(self):
        for samples in [(None, None), ((1., 1.), None),
                        ((2., 1.), (1., 1.)), ((0., 0.), (math.nan, 1.))]:
            with self.subTest(samples=samples), \
                 patch.object(benchmark, 'child_cpu', side_effect=samples), \
                 patch.object(benchmark.subprocess, 'run', return_value=SimpleNamespace(returncode=0)):
                _, timing = benchmark.run_cli(['solver'], None)
            self.assertGreaterEqual(timing['cli_wall_seconds'], 0.)
            self.assertIsNone(timing['cli_cpu_seconds'])
            self.assertIsNone(timing['cli_user_seconds'])
            self.assertIsNone(timing['cli_system_seconds'])

    def test_unsupported_platform_is_unknown_cpu(self):
        with patch.object(benchmark, 'resource', None):
            self.assertIsNone(benchmark.child_cpu())

    def test_real_cli_cost_and_output(self):
        with TemporaryDirectory() as folder:
            path = Path(folder) / 'output'
            with path.open('w') as stream:
                process, timing = benchmark.run_cli(
                    [sys.executable, '-c', 'print(sum(i*i for i in range(20000)))'], stream)
            self.assertEqual(process.returncode, 0)
            self.assertEqual(int(path.read_text()), sum(i*i for i in range(20000)))
        self.assertGreater(timing['cli_wall_seconds'], 0.)
        if timing['cli_cpu_seconds'] is not None:
            self.assertGreaterEqual(timing['cli_user_seconds'], 0.)
            self.assertGreaterEqual(timing['cli_system_seconds'], 0.)
            self.assertEqual(timing['cli_cpu_seconds'],
                             timing['cli_user_seconds'] + timing['cli_system_seconds'])

    def test_paired_metrics_do_not_substitute_legacy_wall(self):
        rows = [dict(instance='m', variant=v, status='PASS', wall_seconds=w,
                     cli_wall_seconds=c, cli_cpu_seconds=None)
                for v, w, c in [('baseline', 10., 2.), ('candidate', 5., 4.)]]
        models = [dict(name='m')]
        legacy = compare_performance.summarize(rows, models, 1)
        cli = compare_performance.summarize(rows, models, 1, 'cli_wall_seconds')
        unknown = compare_performance.summarize(rows, models, 1, 'cli_cpu_seconds')
        self.assertEqual(legacy['geometric_mean_speedup'], 2.)
        self.assertEqual(cli['geometric_mean_speedup'], .5)
        self.assertIsNone(unknown['geometric_mean_speedup'])
        self.assertTrue(unknown['all_passed'])  # Math PASS is separate from timing availability.
        self.assertEqual(unknown['common_passed'], 0)

    def test_missing_failed_nonfinite_and_zero_metrics_are_not_ratios(self):
        models = [dict(name='m')]
        for metric in [None, math.nan, math.inf, -1., True, 0.]:
            rows = [dict(instance='m', variant=v, status='PASS', cli_cpu_seconds=metric)
                    for v in ('baseline', 'candidate')]
            summary = compare_performance.summarize(rows, models, 1, 'cli_cpu_seconds')
            self.assertIsNone(summary['geometric_mean_speedup'])
        rows = [dict(instance='m', variant=v, status=state, cli_cpu_seconds=1.)
                for v, state in [('baseline', 'PASS'), ('candidate', 'TIMEOUT')]]
        self.assertIsNone(compare_performance.summarize(rows, models, 1,
                                                       'cli_cpu_seconds')['geometric_mean_speedup'])


class BenchmarkTests(unittest.TestCase):
    def test_infeasible_status_requires_explicit_terminal_decision(self):
        cases = [
            ('lp_simplex', 'result: state=-1, code=7 (infeasible)', '', 'infeasible'),
            ('lp_simplex', 'result: state=-1, code=9 (precision error)', '', None),
            ('scip', '', 'problem is solved [infeasible]', 'infeasible'),
            ('scip', '', 'problem is solved [infeasible or unbounded]', 'ambiguous'),
            ('highs', '', 'HighsModelStatus.kInfeasible', 'infeasible'),
            ('highs', '', 'HighsModelStatus.kUnboundedOrInfeasible', 'ambiguous'),
            ('clp', 'Primal infeasible - objective value 12', '', 'infeasible'),
            ('clp', 'Presolve thinks problem is infeasible', '', None),
            ('glpk', 'PROBLEM HAS NO PRIMAL FEASIBLE SOLUTION', '', 'infeasible'),
            ('glpk', 'c Status: UNDEFINED', '', None),
        ]
        for solver, text, native, expected in cases:
            with self.subTest(solver=solver, text=text, native=native):
                self.assertEqual(benchmark.infeasibility_status(solver, text, native), expected)
                self.assertIsNone(benchmark.infeasibility_status(solver, text, native, 2))

    def test_infeasible_models_are_separate_and_prepared_by_category(self):
        import argparse
        models = benchmark.select_models(argparse.ArgumentParser(), ['greenbea'], 'netlib-infeasible')
        self.assertEqual(models[0]['expected_status'], 'infeasible')
        self.assertIn('/infeasible/', models[0]['raw_path'])
        self.assertIsNone(models[0]['reference'])
        with TemporaryDirectory() as folder:
            run = Path(folder)
            (run / 'logs').mkdir()
            with benchmark.prepared_model(models[0], run) as path:
                self.assertEqual(path.parent.name, 'infeasible')
                self.assertTrue(path.is_file())
            self.assertFalse(path.exists())

    def test_infeasible_report_does_not_claim_solution_or_certificate_audit(self):
        text = benchmark.render_report('netlib-infeasible', {}, {})
        english = benchmark.render_report('netlib-infeasible', {}, {}, english=True)
        self.assertIn('预期状态', text)
        self.assertIn('未独立验证 Farkas', text)
        self.assertNotIn('- PASS 要求求解器报告最优', text)
        self.assertIn('Expected status', english)
        self.assertIn('not independently verified', english)

    def test_report_reader_ignores_unsuccessful_status_detail_table(self):
        with TemporaryDirectory() as folder:
            path = Path(folder) / 'report.md'
            row = '| pang | 361 / 459 / 2000 | INFEASIBLE | WRONG_STATUS / — / — |'
            path.write_text('| 实例 | 求解器 | 结果 | 终止信息 |\n'
                            '| --- | --- | --- | --- |\n'
                            '| pang | lp_simplex | WRONG_STATUS | unbounded |\n\n'
                            + benchmark.TABLE_HEADER + row + '\n')
            self.assertEqual(benchmark.read_table(path), {'pang': row})

    def test_solver_selection_defaults_and_optional_executables(self):
        for selected, algorithm in (([], 'auto'), ([], 'tableau'), (['--solvers', 'highs'], 'auto')):
            with self.subTest(selected=selected), TemporaryDirectory() as folder:
                root = Path(folder)
                binary = root / 'lp_solve_mps'
                if not selected:
                    binary.touch()
                argv = ['benchmark.py', '--binary', str(binary), '--run-dir', str(root / 'run'),
                        '--docs-dir', str(root / 'docs'), '--lp-algorithm', algorithm, *selected]
                with patch.object(sys, 'argv', argv), \
                     patch.object(benchmark, 'select_models', return_value=[dict(name='afiro')]), \
                     patch.object(benchmark.shutil, 'which', side_effect=AssertionError('unselected executable lookup')), \
                     patch.object(benchmark, 'environment', return_value={}) as environment, \
                     patch.object(benchmark, 'sources', return_value={}), \
                     patch.object(benchmark, 'input_metadata', return_value={}), \
                     patch.object(benchmark, 'checkpoint'), patch.object(benchmark, 'benchmark_case'):
                    benchmark.main()
                config = json.loads((root / 'run/config.json').read_text())
                self.assertEqual(config['solvers'], ['highs'] if selected else ['lp_simplex'])
                self.assertEqual(config['lp_algorithm'], algorithm)
                self.assertEqual(config['repeats'], 1)
                self.assertFalse(config['glpk_stop_on_timeout'])
                self.assertEqual(environment.call_args.args, (None if selected else binary.resolve(), None, None, None))

    def test_lp_only_summary_does_not_require_glpk(self):
        config = dict(models=[dict(name='afiro')], solvers=['lp_simplex'], repeats=1,
                      seconds=10., tolerance=1e-7, iterations=100, lp_algorithm='auto',
                      scip_algorithm='default', highs_algorithm='choose', glpk_stop_on_timeout=True)
        summary = benchmark.run_summary(Path('run'), dict(config=config, environment={}), [])
        self.assertNotIn('GLPK 跳过策略', summary)
        self.assertIn('| lp_simplex | 0/1 |', summary)

    def test_competitor_only_statistics_do_not_compare_with_historical_lp(self):
        results = [dict(instance='afiro', solver='highs', repeat=1, status='PASS', wall_seconds=1.)]
        stats = benchmark.success_statistics(results, ['afiro'], 1, 10., ['highs'])
        self.assertEqual(stats['highs']['successes'], 1)
        self.assertEqual(stats['highs']['common_count'], 0)
        self.assertIsNone(stats['highs']['lp_common_mean'])

    def test_glpk_first_timeout_stops_subsequent_runs_in_large_datasets(self):
        config = dict(dataset='railway', glpk_stop_on_timeout=True, repeats=1)
        rows = [dict(instance='rail507', solver='glpk', status='PASS'),
                dict(instance='rail2536', solver='glpk', status='TIMEOUT')]
        for name in ['rail2536', 'rail2586', 'rail4284', 'rail4872']:
            self.assertIn('rail2536', benchmark.glpk_skip_reason(config, 'glpk', dict(name=name), rows))
            for solver in ['lp_simplex', 'scip', 'highs', 'clp']:
                self.assertIsNone(benchmark.glpk_skip_reason(config, solver, dict(name=name), rows))
        for status in ['PASS', 'AUDIT_FAILED', 'ERROR']:
            with self.subTest(status=status):
                self.assertIsNone(benchmark.glpk_skip_reason(config, 'glpk', dict(name='rail2586'),
                    [dict(instance='rail2536', solver='glpk', status=status)]))

    def test_netlib_never_skips_glpk_even_with_legacy_skip_settings(self):
        for dataset in ('netlib', 'netlib-infeasible'):
            with self.subTest(dataset=dataset):
                config = dict(dataset=dataset, glpk_stop_on_timeout=True,
                              glpk_gate_case='first', repeats=1)
                rows = [dict(instance='first', solver='glpk', status='TIMEOUT')]
                self.assertIsNone(benchmark.glpk_skip_reason(config, 'glpk', dict(name='next'), rows))

    def test_glpk_timeout_policy_is_selected_by_dataset(self):
        for dataset in benchmark.DATASETS:
            with self.subTest(dataset=dataset), TemporaryDirectory() as folder:
                root = Path(folder)
                argv = ['benchmark.py', '--dataset', dataset, '--solvers', 'glpk',
                        '--run-dir', str(root / 'run'), '--docs-dir', str(root / 'docs')]
                with patch.object(sys, 'argv', argv), \
                     patch.object(benchmark, 'select_models', return_value=[dict(name='case')]), \
                     patch.object(benchmark.shutil, 'which', return_value='/bin/glpsol'), \
                     patch.object(benchmark, 'environment', return_value={}), \
                     patch.object(benchmark, 'sources', return_value={}), \
                     patch.object(benchmark, 'input_metadata', return_value={}), \
                     patch.object(benchmark, 'checkpoint'), patch.object(benchmark, 'benchmark_case'):
                    benchmark.main()
                config = json.loads((root / 'run/config.json').read_text())
                self.assertEqual(config['glpk_stop_on_timeout'], dataset in ('pds', 'railway'))

    def test_netlib_rejects_glpk_gate_before_launching_workers(self):
        for dataset in ('netlib', 'netlib-infeasible'):
            with self.subTest(dataset=dataset):
                argv = ['benchmark.py', '--dataset', dataset, '--solvers', 'glpk',
                        '--glpk-gate-case', 'first']
                with patch.object(sys, 'argv', argv), patch('sys.stderr'), \
                     patch.object(benchmark, 'select_models', return_value=[dict(name='first')]), \
                     patch.object(benchmark.shutil, 'which') as executable, \
                     patch.object(benchmark, 'benchmark_case') as worker, \
                     self.assertRaises(SystemExit) as error:
                    benchmark.main()
                self.assertEqual(error.exception.code, 2)
                executable.assert_not_called()
                worker.assert_not_called()

    def test_reproduction_paths_work_from_repository_root(self):
        root = Path('/project/lp').resolve()
        with patch.object(benchmark, 'ROOT', root):
            config = dict(dataset='railway', models=[dict(name='rail507')],
                          binary=str(root / 'build with spaces/test_netlib'),
                          scip=str(root.parent / 'external/scip'), clp=str(root / 'clp'),
                          glpk='glpsol', solvers=['lp_simplex', 'scip'],
                          seconds=3600, glpk_gate_case='rail507')
            command = benchmark.reproduction_command(config)
            args = shlex.split(command)
            for option, expected in (('--binary', 'build with spaces/test_netlib'),
                                     ('--scip', 'scip'), ('--clp', './clp'),
                                     ('--glpk', 'glpsol'), ('--seconds', '3600'),
                                     ('--glpk-gate-case', 'rail507')):
                self.assertEqual(args[args.index(option) + 1], expected)
            self.assertEqual(config['binary'], str(root / 'build with spaces/test_netlib'))
            legacy = ('```sh\npython3 data/benchmark.py --dataset railway '
                      '--binary /project/lp/build/test_netlib --scip /opt/bin/scip\n```\n')
            for english in (False, True):
                report = benchmark.render_report('railway', {}, {'old': legacy}, english=english)
                self.assertIn('--binary build/benchmark/data/lp_solve_mps --scip scip', report)
                self.assertIn('BENCHMARKS=ON\n\ncmake --build', report)
                self.assertIn('--target lp_solve_mps -j 4\n\npython3', report)
                self.assertEqual(report.count('cmake --build'), 1)
                self.assertNotIn('/project/', report)
                self.assertNotIn('/opt/', report)

    def test_glpk_gate_only_skips_glpk_after_unsuccessful_gate(self):
        config = dict(dataset='railway', glpk_gate_case='rail507', repeats=1)
        row = dict(instance='rail507', solver='glpk', repeat=1, status='TIMEOUT', wall_seconds=3600.)
        self.assertIsNotNone(benchmark.glpk_skip_reason(config, 'glpk', dict(name='rail2586'), [row]))
        self.assertIsNone(benchmark.glpk_skip_reason(config, 'lp_simplex', dict(name='rail2586'), [row]))
        self.assertIsNone(benchmark.glpk_skip_reason(config, 'glpk', dict(name='rail507'), []))
        row['status'] = 'PASS'
        self.assertIsNone(benchmark.glpk_skip_reason(config, 'glpk', dict(name='rail2586'), [row]))
        with self.assertRaises(ValueError):
            benchmark.glpk_skip_reason(config, 'glpk', dict(name='rail2586'), [])

    def test_skips_are_separate_from_timeouts_failures_and_pending(self):
        results = [dict(instance='a', solver='lp_simplex', repeat=1, status='PASS', wall_seconds=1.),
                   dict(instance='a', solver='glpk', repeat=1, status='TIMEOUT', wall_seconds=3600.),
                   dict(instance='b', solver='lp_simplex', repeat=1, status='PASS', wall_seconds=2.),
                   dict(instance='b', solver='glpk', repeat=1, status='SKIPPED', wall_seconds=None)]
        stats = benchmark.success_statistics(results, ['a', 'b'], 1, 3600., ['lp_simplex', 'glpk'])
        self.assertEqual((stats['glpk']['timeouts'], stats['glpk']['failures'],
                          stats['glpk']['skipped'], stats['glpk']['pending']), (1, 0, 1, 0))
        self.assertEqual(stats['glpk']['common_count'], 0)
        self.assertIsNone(stats['glpk']['common_mean'])
        self.assertEqual(benchmark.aggregate(results, ['b'], 1, 3600., ['glpk'])[0]['glpk']['status'], 'SKIPPED')

    def test_overview_migrates_history_to_latest_run_in_both_languages(self):
        with TemporaryDirectory() as folder:
            docs = Path(folder)
            stored = {'afiro': '| afiro | 1 / 2 / 3 | 0 | PASS / 1 / 0 | PASS / 2 / 0 | PASS / 3 / 0 | NOT_RUN / — / — | NOT_RUN / — / — |'}
            latest = '| 测试信息 | 值 |\n| --- | --- |\n| 原始记录 | tmp/benchmarks/latest |\n'
            text = benchmark.render_report('netlib', stored, {'latest': latest})
            # Simulate an old report with a full historical block before the latest.
            text = text.replace('<!-- benchmark-run:latest -->',
                                '<!-- benchmark-run:old -->\nOLD RUN\n<!-- /benchmark-run -->\n\n'
                                '<!-- benchmark-run:latest -->')
            (docs / 'netlib.md').write_text(text)
            benchmark.overview(docs)
            generated = {p.name: p.read_text() for p in docs.iterdir()}
            for name, content in generated.items():
                self.assertNotIn('OLD RUN', content, name)
                self.assertNotIn('benchmark-run:old', content, name)
                self.assertEqual(content.count('<!-- benchmark-run:latest -->'), 1, name)
                self.assertIn('tmp/benchmarks/latest', content, name)
            self.assertEqual(benchmark.read_table(docs / 'netlib.md'), stored)
            self.assertIn(stored['afiro'], generated['netlib_en.md'])
            benchmark.overview(docs)
            self.assertEqual({p.name: p.read_text() for p in docs.iterdir()}, generated)

    def test_nested_driver_finds_cmake_build_root(self):
        with TemporaryDirectory() as folder:
            root = Path(folder).resolve()
            binary = root / 'data/lp_solve_mps'
            binary.parent.mkdir(parents=True)
            binary.touch()
            (root / 'CMakeCache.txt').write_text('CMAKE_BUILD_TYPE:STRING=Release\n')
            self.assertEqual(benchmark.build_root(binary), root)
            # Older build layouts must still record the same build tree.
            self.assertEqual(benchmark.build_root(root / 'data/test_netlib'), root)

    def setUp(self):
        self.lp = SimpleNamespace(
            num_col_=2, num_row_=1, col_cost_=[1., 2.], offset_=0.,
            col_lower_=[0., 0.], col_upper_=[2., math.inf],
            row_lower_=[3.], row_upper_=[math.inf],
            a_matrix_=SimpleNamespace(value_=[1., 1.], index_=[0, 0], start_=[0, 1, 2]))

    def test_original_model_audit(self):
        self.assertTrue(audit(self.lp, [2., 1.], 4., 4., 1e-7)['audit_pass'])
        for point, reported in [([3., 0.], 3.), ([0., 3.], 6.),
                                ([math.nan, 1.], 4.), ([2.], 4.),
                                ([2., 1.], 5.), ([2., 1.], math.nan)]:
            with self.subTest(point=point, reported=reported):
                self.assertFalse(audit(self.lp, point, reported, 4., 1e-7)['audit_pass'])

    def test_scip_omitted_zeros_and_unknown_names(self):
        with TemporaryDirectory() as folder:
            path = Path(folder) / 'solution'
            path.write_text('solution status: optimal solution found\nobjective value: 2\nx 2 (obj:1)\n')
            self.assertEqual(parse_scip_solution(path, ['x', 'y']), ([2., 0.], 2.))
            path.write_text('solution status: optimal solution found\nobjective value: 2\nz 2\n')
            with self.assertRaises(ValueError):
                parse_scip_solution(path, ['x', 'y'])

    def test_clp_binary_export_preserves_precision_and_rejects_truncation(self):
        with TemporaryDirectory() as folder:
            path = Path(folder) / 'solution'
            value = 1.2345678901234567
            data = struct.pack('=iid6d', 1, 2, -value, 0., 0., value, 0., 0., 0.)
            path.write_bytes(data)
            self.assertEqual(benchmark.parse_clp_solution(path, 1, 2), ([value, 0.], -value))
            for bad in (data[:-1], data + b'x', b''):
                path.write_bytes(bad)
                with self.assertRaises(ValueError):
                    benchmark.parse_clp_solution(path, 1, 2)
            path.write_bytes(data)
            with self.assertRaises(ValueError):
                benchmark.parse_clp_solution(path, 2, 1)

    def test_glpk_status_and_complete_coordinates_are_required(self):
        with TemporaryDirectory() as folder:
            path = Path(folder) / 'solution'
            good = 's bas 1 2 f f 4\nj 1 b 2 0\nj 2 b 1 0\ne o f\n'
            path.write_text(good)
            self.assertEqual(benchmark.parse_glpk_solution(path, 2), ([2., 1.], 4.))
            for bad in (good.replace('f f', 'f i'), good.replace('j 2 b 1 0\n', ''),
                        good.replace('j 2', 'j 1'), good.replace('j 2', 'j 3')):
                path.write_text(bad)
                with self.assertRaises(ValueError):
                    benchmark.parse_glpk_solution(path, 2)

    def test_glpk_compatible_input_preserves_original_and_model_data(self):
        with TemporaryDirectory() as folder:
            path = Path(folder) / 'model.mps'
            original = 'NAME MODEL\nOBJSENSE\n    MAX\nROWS\n N OBJ\n L LIMIT\nENDATA\n'
            path.write_text(original)
            compatible = benchmark.glpk_input(path)
            self.assertNotEqual(compatible, path)
            self.assertEqual(path.read_text(), original)
            self.assertEqual(compatible.read_text(), original.replace('OBJSENSE\n    MAX\n', ''))
            self.assertEqual(benchmark.glpk_input(compatible), compatible)

    def test_legacy_table_migration_preserves_three_solver_results(self):
        with TemporaryDirectory() as folder:
            path = Path(folder) / 'netlib.md'
            path.write_text('| 实例 | 行 | 列 | 非零元 | 参考目标 | lp | scip | highs | 重复 | 限时/秒 | 算法 | 运行与进度 |\n'
                            '| --- | --- |\n| m | 1 | 2 | 3 | 4 | PASS / 1 / 4 | TIMEOUT / — / — | PASS / 2 / 4 | 1 | 30 | auto/default/choose | oldrun |\n')
            row = benchmark.read_table(path)['m']
            self.assertIn('| m | 1 / 2 / 3 | 4 | PASS / 1 / 4 | TIMEOUT / — / — | PASS / 2 / 4 |', row)
            self.assertEqual(row.count('NOT_RUN'), 2)
            self.assertNotIn('oldrun', row)

    def test_timeout_and_failed_repeat_remain_visible(self):
        results = [dict(instance='m', solver=solver, status=status,
                        wall_seconds=elapsed, repeat=i + 1)
                   for solver in ('lp_simplex', 'scip', 'highs')
                   for i, (status, elapsed) in enumerate([('PASS', 1.), ('TIMEOUT', 10.)])]
        value = aggregate(results, ['m'], 2, 10.)[0]['lp_simplex']
        self.assertFalse(value['passed'])
        self.assertIsNone(value['median'])
        self.assertEqual(value['status'], 'PASS/TIMEOUT')
        self.assertEqual(value['penalized'], 50.5)

    def test_dotted_case_names_keep_distinct_artifacts(self):
        self.assertEqual(artifact(Path('pilot.ja-highs-1'), '.json').name,
                         'pilot.ja-highs-1.json')
        self.assertNotEqual(artifact(Path('pilot.ja-highs-1'), '.log'),
                            artifact(Path('pilot.ja-scip-1'), '.log'))

    def test_incomplete_report_cannot_publish(self):
        with TemporaryDirectory() as folder:
            run = Path(folder)
            (run / 'metadata.json').write_text(json.dumps(
                dict(config={}, environment={}, completed=False)))
            (run / 'results.json').write_text('[]')
            with self.assertRaisesRegex(ValueError, 'Incomplete'):
                report(run, run / 'docs')
            self.assertFalse((run / 'docs/README.md').exists())


class SummaryTests(unittest.TestCase):
    def test_success_geometric_means_use_explicit_common_set(self):
        data = {'a': [1., 4., 9.], 'b': [9., None, 1.], 'c': [None, 16., None]}
        results = [dict(instance=name, solver=solver, repeat=1,
                        status='PASS' if elapsed is not None else 'TIMEOUT',
                        wall_seconds=elapsed if elapsed is not None else 20.)
                   for name, times in data.items() for solver, elapsed in zip(benchmark.LEGACY_SOLVERS, times)]
        stats = benchmark.success_statistics(results, ['a', 'b', 'c', 'pending'], 1, 20., benchmark.LEGACY_SOLVERS)
        for solver, baseline, other, count in [
                ('lp_simplex', 3., 3., 2), ('scip', 1., 4., 1), ('highs', 3., 3., 2)]:
            self.assertAlmostEqual(stats[solver]['lp_common_mean'], baseline)
            self.assertAlmostEqual(stats[solver]['common_mean'], other)
            self.assertEqual(stats[solver]['successes'], 2)
            self.assertEqual(stats[solver]['common_count'], count)
            self.assertEqual(stats[solver]['pending'], 1)
        self.assertEqual(stats['scip']['timeouts'], 1)

    def test_pairwise_statistics_ignore_unrelated_solver_failures(self):
        results = [dict(instance=name, solver=solver, repeat=1,
                        status='PASS', wall_seconds=elapsed)
                   for name, solver, elapsed in [
                       ('a', 'lp_simplex', 1.), ('b', 'lp_simplex', 9.),
                       ('a', 'scip', 4.), ('b', 'scip', 16.)]]
        # HiGHS has no completed results, so the all-solver intersection
        # is empty. This must not erase the two valid lp/SCIP comparisons.
        stats = benchmark.success_statistics(results, ['a', 'b'], 1, 20., benchmark.LEGACY_SOLVERS)
        self.assertEqual(stats['scip']['common_count'], 2)
        self.assertAlmostEqual(stats['scip']['lp_common_mean'], 3.)
        self.assertAlmostEqual(stats['scip']['common_mean'], 8.)
        self.assertEqual(stats['highs']['common_count'], 0)
        self.assertIsNone(stats['highs']['lp_common_mean'])
        self.assertIsNone(stats['highs']['common_mean'])

    def test_partial_repeats_and_deadline_do_not_count_as_success(self):
        rows = [dict(instance='m', solver=solver, repeat=repeat, status=status, wall_seconds=elapsed)
                for solver, repeat, status, elapsed in [('lp_simplex', 1, 'PASS', 1.),
                    ('scip', 1, 'PASS', 1.), ('scip', 2, 'TIMEOUT', 10.),
                    ('highs', 1, 'PASS', 11.), ('highs', 2, 'PASS', 1.)]]
        stats = benchmark.success_statistics(rows, ['m'], 2, 10.)
        for solver in benchmark.SOLVERS:
            self.assertEqual(stats[solver]['successes'], 0)
            self.assertIsNone(stats[solver]['lp_common_mean'])
            self.assertIsNone(stats[solver]['common_mean'])

    def test_relative_time_is_lp_divided_by_solver_and_handles_zero(self):
        config = dict(models=[dict(name='m')], solvers=benchmark.LEGACY_SOLVERS,
                      repeats=1, seconds=10., lp_algorithm='auto',
                      scip_algorithm='default', highs_algorithm='choose',
                      tolerance=1e-7, iterations=100)
        metadata = dict(config=config, environment={})
        for lp_time, solver_time, expected in [(1., 4., '0.25'), (0., 2., '0'), (2., 0., '—')]:
            results = [dict(instance='m', solver=solver, repeat=1,
                            status='PASS', wall_seconds=elapsed)
                       for solver, elapsed in [('lp_simplex', lp_time), ('scip', solver_time)]]
            text = benchmark.run_summary(Path('run'), metadata, results)
            row = next(line for line in text.splitlines() if line.startswith('| SCIP / SoPlex |'))
            self.assertEqual(row.strip('|').split('|')[-1].strip(), expected)
            self.assertNotIn('待测次数', text)
            self.assertIn('相对耗时（lp_simplex / 各求解器）', text)

    def test_pds_summary_displays_numeric_case_order_without_changing_run_order(self):
        models = [dict(name=name) for name in ['pds-100', 'pds-80', 'pds-20']]
        config = dict(dataset='pds', suite_label='pds 全集', models=models,
                      solvers=['lp_simplex'], repeats=1, seconds=300.,
                      lp_algorithm='auto', scip_algorithm='default', highs_algorithm='choose',
                      tolerance=1e-7, iterations=300000)
        text = benchmark.run_summary(Path('run'), dict(config=config, environment={}), [])
        self.assertIn('| 案例范围 | pds 全集；pds-20, pds-80, pds-100 |', text)
        self.assertEqual([model['name'] for model in models], ['pds-100', 'pds-80', 'pds-20'])

    def test_geometric_mean_zero_and_empty(self):
        self.assertIsNone(benchmark.geometric_mean([]))
        self.assertEqual(benchmark.geometric_mean([0., 1.]), 0.)
        with self.assertRaises(ValueError):
            benchmark.geometric_mean([-1.])


class LifecycleTests(unittest.TestCase):
    def test_gate_skip_does_not_launch_worker(self):
        models = self.models(['other'])
        config = dict(dataset='railway', models=models, repeats=1, solvers=['lp_simplex', 'glpk'], glpk_gate_case='pilot.ja')
        self.results.append(dict(instance='pilot.ja', solver='glpk', repeat=1, status='TIMEOUT', wall_seconds=3600.))
        with patch.object(benchmark, 'run_one', return_value=dict(instance='other', solver='lp_simplex', repeat=1, status='PASS', wall_seconds=1.)) as run:
            benchmark.benchmark_case(self.config_path, config, models[0], 1, {}, self.results, self.metadata)
        self.assertEqual(run.call_count, 1)
        self.assertEqual(run.call_args.args[2], 'lp_simplex')
        skipped = next(r for r in self.results if r['instance'] == 'other' and r['solver'] == 'glpk')
        self.assertEqual(skipped['status'], 'SKIPPED')
        self.assertIsNone(skipped['wall_seconds'])

    def setUp(self):
        self.directory = TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        raw = self.root / 'data/netlib/raw/feasible'
        raw.mkdir(parents=True)
        self.raw = raw
        self.raw.joinpath('pilot.ja.mps.gz').write_bytes(gzip.compress(b'NAME test\nENDATA\n'))
        self.raw.joinpath('other.mps.gz').write_bytes(gzip.compress(b'NAME other\nENDATA\n'))
        source_root = Path(__file__).resolve().parents[1]
        shutil.copyfile(source_root / 'data/netlib/prepare.py', self.root / 'data/netlib/prepare.py')
        (self.root / 'data/netlib/feasible_gurobi_1e-8.csv').write_text(
            'name,status,pobj\npilot.ja,OPTIMAL,0\nother,OPTIMAL,0\n')
        self.run = self.root / 'run'
        (self.run / 'logs').mkdir(parents=True)
        self.config_path = self.run / 'config.json'
        self.results = []
        self.metadata = dict(model_changed_during_run=False)
        self.patch = patch.object(benchmark, 'ROOT', self.root)
        self.patch.start()
        self.addCleanup(self.patch.stop)

    def models(self, cases=None):
        return benchmark.select_models(argparse.ArgumentParser(), cases)

    def test_full_suite_is_selected_from_raw_without_generated_mps(self):
        models = self.models()
        self.assertEqual([v['name'] for v in models], ['other', 'pilot.ja'])
        self.assertFalse(list(self.root.rglob('*.mps')))
        self.assertEqual(self.models(['pilot.ja'])[0]['name'], 'pilot.ja')

    def execute(self, worker, repeats=2):
        models = self.models(['pilot.ja'])
        self.config = dict(models=models, repeats=repeats)
        benchmark.dump(self.config_path, self.config)
        with patch.object(benchmark, 'run_one', side_effect=worker):
            benchmark.benchmark_case(self.config_path, self.config, models[0], 0,
                                     {}, self.results, self.metadata)

    def test_same_model_for_three_solvers_and_repeats_then_cleanup(self):
        paths, solvers = [], []
        # A user-generated MPS outside our temporary directory is never touched.
        existing = self.root / 'data/netlib/feasible/pilot.ja.mps'
        existing.parent.mkdir()
        existing.write_bytes(b'user data')
        def worker(config_path, config, solver, model, repeat, env):
            path = Path(model['path'])
            self.assertTrue(path.is_file())
            self.assertEqual(len(list(self.run.rglob('*.mps'))), 1)
            paths.append(path)
            solvers.append((solver, repeat))
            return dict(instance=model['name'], solver=solver, repeat=repeat,
                        status='PASS', wall_seconds=1.)
        self.execute(worker)
        self.assertEqual(len(set(paths)), 1)
        self.assertEqual(set(solvers), {(s, r) for s in benchmark.SOLVERS for r in (1, 2)})
        self.assertEqual(len(self.results), 10)
        self.assertFalse(list(self.run.glob('mps-*')))
        self.assertEqual(existing.read_bytes(), b'user data')
        self.assertFalse(self.metadata['model_changed_during_run'])

    def test_interrupt_cleans_mps_and_keeps_completed_result(self):
        def worker(config_path, config, solver, model, repeat, env):
            if solver == 'scip':
                raise KeyboardInterrupt()
            return dict(instance=model['name'], solver=solver, repeat=repeat,
                        status='PASS', wall_seconds=1.)
        with self.assertRaises(KeyboardInterrupt):
            self.execute(worker)
        self.assertFalse(list(self.run.glob('mps-*')))
        self.assertEqual(len(json.loads((self.run / 'results.json').read_text())), 1)

    def test_solver_error_does_not_skip_other_solvers(self):
        def worker(config_path, config, solver, model, repeat, env):
            if solver == 'lp_simplex':
                raise ValueError('adapter failed')
            return dict(instance=model['name'], solver=solver, repeat=repeat,
                        status='PASS', wall_seconds=1.)
        self.execute(worker)
        self.assertEqual(len(self.results), 10)
        self.assertEqual(sum(row['status'] == 'ERROR' for row in self.results), 2)
        self.assertFalse(list(self.run.glob('mps-*')))

    def test_corrupt_raw_records_generation_failure_and_cleans_partial_mps(self):
        (self.raw / 'pilot.ja.mps.gz').write_bytes(b'bad gzip')
        with patch.object(benchmark, 'run_one') as worker:
            self.execute(worker)
            worker.assert_not_called()
        self.assertEqual(len(self.results), 10)
        self.assertTrue(all(row['status'] == 'PREPARE_FAILED' for row in self.results))
        self.assertFalse(list(self.run.glob('mps-*')))

    def test_pending_cases_have_no_timing_statistics(self):
        row = aggregate([], ['pending'], 3, 10.)[0]
        for solver in benchmark.SOLVERS:
            self.assertEqual(row[solver]['status'], 'PENDING')
            self.assertIsNone(row[solver]['penalized'])

    def test_progress_report_lists_pending_cases_and_remains_incomplete(self):
        models = self.models()
        config = dict(models=models, solvers=benchmark.SOLVERS, repeats=2, seconds=10., suite_label='test subset',
                      lp_algorithm='auto', scip_algorithm='default', highs_algorithm='choose',
                      tolerance=1e-7, iterations=100)
        environment = dict(os='test', cpu='test', memory_bytes='unknown', logical_cpus=1,
                           python='3.12', numpy='test', scipy='test', versions=dict(
            lp_simplex=dict(git_head='abc', git_status='', compiler='cc', cmake_cache={}),
            scip=dict(version_output='SCIP test'), highs=dict(version='test', binding='test')))
        metadata = dict(config=config, environment=environment, completed=False,
                        state='进行中', finished_cases=['other'], **{key: False for key in (
                            'source_changed_during_run', 'binary_changed_during_run',
                            'library_changed_during_run', 'model_changed_during_run',
                            'competitor_changed_during_run')})
        results = [dict(instance='other', solver=solver, repeat=repeat, status='PASS', wall_seconds=1.)
                   for solver in benchmark.SOLVERS for repeat in (1, 2)]
        benchmark.dump(self.run / 'metadata.json', metadata)
        benchmark.dump(self.run / 'results.json', results)
        docs = self.run / 'docs'
        report(self.run, docs, allow_partial=True)
        text = (docs / 'netlib.md').read_text()
        self.assertIn('进行中；1/2 个案例', text)
        self.assertIn('## 测试与统计口径', text)
        self.assertIn('test / test', text)
        self.assertIn('与 lp_simplex 共同成功数', text)
        self.assertNotIn('各自成功案例几何平均/秒', text)
        self.assertIn('共同案例几何平均/秒（lp_simplex / 本求解器）', text)
        self.assertTrue(next(line for line in text.splitlines()
                             if line.startswith('| lp_simplex |')).endswith('| — | — | — |'))
        self.assertIn('| SCIP / SoPlex | 1/2 | 0 | 0 | 1 | 1 / 1 | 1 |', text)
        self.assertIn('cc', text)
        index = (docs / 'README.md').read_text()
        self.assertIn('## Netlib feasible', index)
        self.assertIn('## PDS', index)
        self.assertIn('## Railway 连续 LP 松弛', index)
        self.assertIn('与 lp_simplex 共同成功数', index)
        self.assertIn('test / test', index)
        self.assertIn('--cases other pilot.ja', index)
        self.assertIn('## 复现说明', index)
        self.assertNotIn(benchmark.TABLE_HEADER.splitlines()[0], index)
        self.assertNotIn('## 逐例最新结果', index)
        self.assertEqual({p.name for p in docs.iterdir()}, {'README.md', 'README_en.md', 'netlib.md', 'netlib_en.md',
                                                                  'pds.md', 'pds_en.md', 'railway.md', 'railway_en.md',
                                                                  'netlib-infeasible.md', 'netlib-infeasible_en.md'})
        english = (docs / 'README_en.md').read_text()
        self.assertIn('[English](./README_en.md)', index)
        self.assertIn('[简体中文](./README.md)', english)
        self.assertIn('[Project home](../../README_en.md)', english)
        self.assertIn('[Documentation](../README_en.md)', english)
        english_report = (docs / 'netlib_en.md').read_text()
        self.assertIn('[English](./netlib_en.md)', text)
        self.assertIn('[简体中文](./netlib.md)', english_report)
        self.assertIn('[Benchmark summary](./README_en.md)', english_report)
        self.assertIn('## Latest per-instance results', english_report)
        self.assertIn('Running; 1/2 instances', english_report)
        self.assertIn('PENDING / — / —', english_report)
        self.assertIn('Not yet run;', (docs / 'pds_en.md').read_text())
        self.assertIn('(netlib_en.md)', english)
        self.assertNotRegex(english_report.replace('简体中文', ''), r'[\u4e00-\u9fff]')
        # Data cells and instance identifiers remain unchanged across translations.
        zh_rows = benchmark.read_table(docs / 'netlib.md')
        en_rows = {line.split('|')[1].strip(): line for line in english_report.splitlines()
                   if any(line.startswith('| ' + name + ' |') for name in zh_rows)}
        for name, row in zh_rows.items():
            self.assertEqual(row.split('|')[1:-2], en_rows[name].split('|')[1:-2])
        self.assertIn('Running; 1/2 instances', english)
        self.assertNotIn('### 运行 ', index)
        self.assertNotIn('### Run ', english)
        self.assertIn('| 时间（UTC） |', index)
        self.assertIn('| Time (UTC) |', english)
        self.assertIn('Not yet run;', english)
        self.assertIn('## Railway continuous LP relaxations', english)
        self.assertIn('unknown bytes / 1', english)
        self.assertIn('Geometric means on common set / s (lp_simplex / solver)', english)
        self.assertNotRegex(english.replace('简体中文', ''), r'[\u4e00-\u9fff]')
        # Localization preserves solver statistics and reproduction commands verbatim.
        solver_rows = lambda text: [line for line in text.splitlines()
                                   if line.count('|') > 4 and any(line.startswith('| ' + label + ' |')
                                          for label in benchmark.LABELS.values())]
        self.assertEqual(solver_rows(index), solver_rows(english))
        self.assertEqual(re.findall(r'```sh\n(.*?)```', index, re.S)[:1],
                         re.findall(r'```sh\n(.*?)```', english, re.S)[:1])
        self.assertIn('PENDING / — / —', benchmark.read_table(docs / 'netlib.md')['pilot.ja'])
        self.assertNotIn('（已完成 ', text)
        # A new report replaces older run summaries instead of accumulating them.
        target = docs / 'netlib.md'
        target.write_text(target.read_text().replace('<!-- benchmark-run:',
                          '<!-- benchmark-run:old -->\nOLD RUN\n<!-- /benchmark-run -->\n\n'
                          '<!-- benchmark-run:', 1))
        metadata['state'] = '已中断'
        benchmark.dump(self.run / 'metadata.json', metadata)
        report(self.run, docs, allow_partial=True)
        self.assertNotIn('OLD RUN', target.read_text())
        self.assertEqual(target.read_text().count('<!-- benchmark-run:'), 1)
        self.assertIn('已中断；1/2 个案例', (docs / 'netlib.md').read_text())
        self.assertIn('已中断；1/2 个案例', (docs / 'README.md').read_text())
        self.assertIn('Interrupted; 1/2 instances', (docs / 'README_en.md').read_text())
        self.assertIn('Interrupted; 1/2 instances', (docs / 'netlib_en.md').read_text())
        # Updating a subset preserves earlier rows and the other two datasets.
        other_tables = {name: (docs / (name + '.md')).read_bytes() for name in ('pds', 'railway')}
        previous = benchmark.read_table(docs / 'netlib.md')['pilot.ja']
        metadata['config']['models'] = [models[0]]
        benchmark.dump(self.run / 'metadata.json', metadata)
        report(self.run, docs, allow_partial=True)
        self.assertEqual(benchmark.read_table(docs / 'netlib.md')['pilot.ja'], previous)
        for name, content in other_tables.items():
            self.assertEqual((docs / (name + '.md')).read_bytes(), content)
        # The index can be regenerated from Markdown alone after raw logs are removed.
        original = (docs / 'README.md').read_bytes()
        original_english = (docs / 'README_en.md').read_bytes()
        translated_reports = {name: (docs / (name + '_en.md')).read_bytes()
                              for name in benchmark.DATASETS}
        (self.run / 'metadata.json').unlink()
        (self.run / 'results.json').unlink()
        (docs / 'README.md').unlink()
        (docs / 'README_en.md').unlink()
        for name in benchmark.DATASETS:
            (docs / (name + '_en.md')).unlink()
        benchmark.overview(docs)
        self.assertEqual((docs / 'README.md').read_bytes(), original)
        self.assertEqual((docs / 'README_en.md').read_bytes(), original_english)
        for name, content in translated_reports.items():
            self.assertEqual((docs / (name + '_en.md')).read_bytes(), content)

    def test_pds_selection_uses_its_reference_file_and_manifest(self):
        base = self.root / 'data/pds'
        (base / 'raw').mkdir(parents=True)
        (base / 'raw/pds-30.bz2').write_bytes(bz2.compress(b'raw'))
        (base / 'manifest.csv').write_text(
            'instance,constraints,variables,constraint_nonzeros,mps_sha256\n'
            'pds-30,2,3,4,abc\n')
        reference = base / 'reference_objectives.csv'
        reference.write_text('name,status,pobj,source_url\npds-30,OPTIMAL,21385445736,https://example.com\n')
        models = benchmark.select_models(argparse.ArgumentParser(), None, 'pds')
        self.assertEqual([m['name'] for m in models], ['pds-30'])
        self.assertEqual(models[0]['reference'], 21385445736.)
        self.assertEqual(models[0]['rows'], 2)
        self.assertEqual(models[0]['expected_mps_sha256'], 'abc')
        self.assertEqual(models[0]['dataset'], 'pds')
        metadata = benchmark.input_metadata('pds')
        self.assertEqual(metadata['data/pds/reference_objectives.csv'], benchmark.sha(reference))
        self.assertNotIn('data/netlib/feasible_gurobi_1e-8.csv', metadata)

    def test_pds_execution_and_reports_use_numeric_case_order(self):
        base = self.root / 'data/pds'
        (base / 'raw').mkdir(parents=True)
        for name in ('pds-100', 'pds-90', 'pds-20'):
            (base / 'raw' / (name + '.bz2')).write_bytes(bz2.compress(b'raw'))
        (base / 'reference_objectives.csv').write_text(
            'name,status,pobj\npds-100,OPTIMAL,1\npds-90,OPTIMAL,2\npds-20,OPTIMAL,3\n')
        expected = ['pds-20', 'pds-90', 'pds-100']
        for cases in (None, ['pds-100', 'pds-90', 'pds-20']):
            models = benchmark.select_models(argparse.ArgumentParser(), cases, 'pds')
            self.assertEqual([m['name'] for m in models], expected)
        # Display sorting preserves every existing result cell in both languages.
        stored = {name: f'| {name} | 1 / 2 / 3 | 1 | PASS / 2 / 1 |'
                  for name in reversed(expected)}
        for english in (False, True):
            report = benchmark.render_report('pds', stored, {}, english=english)
            rows = [line for line in report.splitlines() if line.startswith('| pds-')]
            self.assertEqual(rows, [stored[name] for name in expected])
        # Netlib retains its existing lexical ordering.
        self.assertEqual(benchmark.ordered_cases('netlib', ['pds-20', 'pds-100']),
                         ['pds-100', 'pds-20'])

    def test_railway_execution_and_reports_use_numeric_case_order(self):
        base = self.root / 'data/railway/raw'
        base.mkdir(parents=True)
        expected = ['rail507', 'rail516', 'rail582', 'rail2536',
                    'rail2586', 'rail4284', 'rail4872']
        for name in reversed(expected):
            (base / (name + '.gz')).write_bytes(gzip.compress(b'raw'))
        for cases in (None, list(reversed(expected)), sorted(expected)):
            models = benchmark.select_models(argparse.ArgumentParser(), cases, 'railway')
            self.assertEqual([m['name'] for m in models], expected)
        subset = ['rail4872', 'rail516', 'rail2536']
        models = benchmark.select_models(argparse.ArgumentParser(), subset, 'railway')
        self.assertEqual([m['name'] for m in models], ['rail516', 'rail2536', 'rail4872'])
        self.assertEqual(subset, ['rail4872', 'rail516', 'rail2536'])
        stored = {name: f'| {name} | 1 / 2 / 3 | — | PASS / 2 / 1 |'
                  for name in sorted(expected)}
        block = '| 案例范围 | railway 全集；' + ', '.join(sorted(expected)) + ' |\n'
        for english in (False, True):
            report = benchmark.render_report('railway', stored, {'old': block}, english=english)
            rows = [line for line in report.splitlines() if line.startswith('| rail')]
            self.assertEqual(rows, [stored[name] for name in expected])
            self.assertIn(', '.join(expected), report)

    def test_pds_rejects_missing_invalid_or_duplicate_references(self):
        base = self.root / 'data/pds'
        (base / 'raw').mkdir(parents=True)
        (base / 'raw/pds-30.bz2').write_bytes(bz2.compress(b'raw'))
        reference = base / 'reference_objectives.csv'
        invalid = [None, '', 'name,status,pobj\nother,OPTIMAL,1\n',
                   'name,status,pobj\npds-30,OPTIMAL,nan\n',
                   'name,status,pobj\npds-30,OPTIMAL,inf\n',
                   'name,status,pobj\npds-30,OPTIMAL,invalid\n',
                   'name,status,pobj\npds-30,UNKNOWN,1\n',
                   'name,status,pobj\npds-30,OPTIMAL,1\npds-30,OPTIMAL,2\n']
        for text in invalid:
            with self.subTest(csv=text):
                if text is not None:
                    reference.write_text(text)
                with patch('sys.stderr'), self.assertRaises(SystemExit) as error:
                    benchmark.select_models(argparse.ArgumentParser(), None, 'pds')
                self.assertEqual(error.exception.code, 2)

    def test_reference_metadata_keeps_netlib_and_railway_sources_separate(self):
        netlib = benchmark.input_metadata('netlib')
        self.assertIn('data/netlib/feasible_gurobi_1e-8.csv', netlib)
        self.assertIsNone(benchmark.reference_path('railway'))
        self.assertNotIn('data/netlib/feasible_gurobi_1e-8.csv', benchmark.input_metadata('railway'))

    def test_railway_dispatch_generates_continuous_mps_and_cleans_it(self):
        base = self.root / 'data/railway'
        (base / 'raw').mkdir(parents=True)
        (base / 'raw/railtest.gz').write_bytes(gzip.compress(b'2 1\n1 2 1 2\n'))
        source_root = Path(__file__).resolve().parents[1]
        shutil.copyfile(source_root / 'data/railway/prepare.py', base / 'prepare.py')
        models = benchmark.select_models(argparse.ArgumentParser(), None, 'railway')
        self.assertIsNone(models[0]['reference'])
        with benchmark.prepared_model(models[0], self.run) as path:
            self.assertTrue(path.is_file())
            self.assertIn(' UP BND1', path.read_text())
            self.assertEqual(path.parent.parent, self.run)
        self.assertFalse(list(self.run.rglob('*.mps')))
        self.assertFalse(list(self.run.glob('mps-*')))

    def test_optional_reference_still_requires_original_model_audit(self):
        lp = SimpleNamespace(
            num_col_=1, num_row_=1, col_cost_=[1.], offset_=0.,
            col_lower_=[0.], col_upper_=[1.], row_lower_=[1.], row_upper_=[math.inf],
            a_matrix_=SimpleNamespace(value_=[1.], index_=[0], start_=[0, 1]))
        result = audit(lp, [1.], 1., None, 1e-7)
        self.assertTrue(result['audit_pass'])
        self.assertFalse(result['reference_checked'])
        self.assertIsNone(result['reference_error'])
        self.assertFalse(audit(lp, [2.], 2., None, 1e-7)['audit_pass'])
        self.assertFalse(audit(lp, [1.], 2., None, 1e-7)['audit_pass'])

    def test_missing_run_cannot_be_marked_complete_even_with_partial_option(self):
        data = dict(config=dict(models=[dict(name='m')], repeats=1), environment={},
                    completed=True, **{key: False for key in (
                        'source_changed_during_run', 'binary_changed_during_run',
                        'library_changed_during_run', 'model_changed_during_run',
                        'competitor_changed_during_run')})
        benchmark.dump(self.run / 'metadata.json', data)
        benchmark.dump(self.run / 'results.json', [])
        with self.assertRaisesRegex(ValueError, 'Missing'):
            report(self.run, self.run / 'docs', allow_partial=True)


if __name__ == '__main__':
    unittest.main()
