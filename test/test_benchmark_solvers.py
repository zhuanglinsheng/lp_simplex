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
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'data'))
import benchmark
from benchmark import aggregate, artifact, audit, parse_scip_solution, report


class BenchmarkTests(unittest.TestCase):
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
            binary = root / 'data/netlib/test_netlib'
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
                                                                  'pds.md', 'pds_en.md', 'railway.md', 'railway_en.md'})
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
        # Other datasets retain their existing lexical ordering.
        self.assertEqual(benchmark.ordered_cases('netlib', ['pds-20', 'pds-100']),
                         ['pds-100', 'pds-20'])

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
