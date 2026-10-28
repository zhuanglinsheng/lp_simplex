"""Measurement reuse, invalidation, and report persistence."""
from contextlib import nullcontext
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'data'))
import benchmark


class MeasurementCacheTests(unittest.TestCase):
    def setUp(self):
        self.model = dict(name='afiro', dataset='netlib', raw_sha256='input', reference=0.)
        self.config = dict(dataset='netlib', models=[self.model], solvers=['lp_simplex'],
                           repeats=1, seconds=10., tolerance=1e-7, iterations=100,
                           lp_algorithm='dual-revised', scip_algorithm='default',
                           highs_algorithm='choose', glpk_stop_on_timeout=False)
        self.environment = dict(os='test', architecture='test', cpu='test',
                                versions={'lp_simplex': dict(binary_sha256='launcher',
                                    solver_libraries_sha256={'library': 'solver'}),
                                          'highs': dict(binding_sha256='highs')})
        self.row = dict(instance='afiro', solver='lp_simplex', repeat=1, status='PASS',
                        wall_seconds=1.23456789, objective=0., audit_pass=True)

    def entries(self, rows=None):
        return benchmark.update_result_cache([], self.config, self.environment, rows or [self.row])

    def test_cache_preserves_precision_and_survives_both_report_languages(self):
        cache = self.entries()
        with TemporaryDirectory() as folder:
            docs = Path(folder)
            text = benchmark.render_report('netlib', {}, {'run': 'summary'}, cache=cache)
            (docs / 'netlib.md').write_text(text)
            benchmark.overview(docs)
            for language in ('netlib.md', 'netlib_en.md'):
                saved = benchmark.read_result_cache(docs / language)
                self.assertEqual(saved, cache)
                row = benchmark.cached_measurement(saved, self.config, self.environment,
                                                   self.model, 'lp_simplex', 1)
                self.assertEqual(row['wall_seconds'], 1.23456789)

    def test_library_changes_invalidate_only_our_solver(self):
        highs = dict(self.row, solver='highs')
        cache = self.entries([self.row, highs])
        self.environment['versions']['lp_simplex']['solver_libraries_sha256']['library'] = 'new solver'
        self.assertIsNone(benchmark.cached_measurement(cache, self.config, self.environment,
                                                     self.model, 'lp_simplex', 1))
        self.assertEqual(benchmark.cached_measurement(cache, self.config, self.environment,
                                                     self.model, 'highs', 1), highs)

    def test_imported_failure_without_timing_keeps_competitor_statistics(self):
        failed = dict(self.row, solver='highs', status='SOLVER_FAILED', wall_seconds=None)
        cache = self.entries([self.row, failed])
        metadata = dict(config=self.config, environment=self.environment)
        view, compared = benchmark.comparison_results(metadata, [self.row], cache)
        self.assertIn('highs', view['config']['solvers'])
        self.assertIn(failed, compared)
        for status in ('PASS', 'TIMEOUT'):
            invalid = dict(failed, status=status)
            self.assertIsNone(benchmark.cached_measurement(self.entries([invalid]), self.config,
                                                          self.environment, self.model, 'highs', 1))

    def test_our_iteration_budget_does_not_invalidate_another_solver(self):
        highs = dict(self.row, solver='highs')
        cache = self.entries([self.row, highs])
        self.config['iterations'] += 1
        self.assertIsNone(benchmark.cached_measurement(cache, self.config, self.environment,
                                                     self.model, 'lp_simplex', 1))
        self.assertEqual(benchmark.cached_measurement(cache, self.config, self.environment,
                                                     self.model, 'highs', 1), highs)

    def test_inputs_settings_machine_and_repetitions_are_separate(self):
        cache = self.entries()
        for item, key, value in [(self.model, 'raw_sha256', 'changed input'),
                                 (self.config, 'seconds', 20.),
                                 (self.config, 'lp_algorithm', 'auto'),
                                 (self.environment, 'cpu', 'another CPU')]:
            with self.subTest(key=key):
                previous = item[key]
                item[key] = value
                self.assertIsNone(benchmark.cached_measurement(cache, self.config, self.environment,
                                                             self.model, 'lp_simplex', 1))
                item[key] = previous
        self.assertIsNone(benchmark.cached_measurement(cache, self.config, self.environment,
                                                     self.model, 'lp_simplex', 2))

    def test_cached_case_does_not_generate_mps_or_start_a_solver(self):
        with TemporaryDirectory() as folder:
            results = []
            metadata = dict(environment=self.environment, result_cache=self.entries())
            with patch.object(benchmark, 'prepared_model', side_effect=AssertionError('MPS generation')), \
                 patch.object(benchmark, 'run_one', side_effect=AssertionError('solver launch')):
                benchmark.benchmark_case(Path(folder) / 'config.json', self.config,
                                         self.model, 0, {}, results, metadata)
            self.assertEqual(results, [self.row])

    def test_changed_binary_retests_only_that_solver(self):
        self.config['solvers'] = ['lp_simplex', 'highs']
        highs = dict(self.row, solver='highs')
        metadata = dict(environment=self.environment, result_cache=self.entries([self.row, highs]))
        self.environment['versions']['lp_simplex']['binary_sha256'] = 'new launcher'
        with TemporaryDirectory() as folder:
            path = Path(folder) / 'model.mps'
            path.write_text('model')
            self.model['sha256'] = benchmark.sha(path)
            results = []
            with patch.object(benchmark, 'prepared_model', return_value=nullcontext(path)), \
                 patch.object(benchmark, 'run_one', return_value=self.row) as run:
                benchmark.benchmark_case(Path(folder) / 'config.json', self.config,
                                         self.model, 0, {}, results, metadata)
            self.assertEqual(run.call_count, 1)
            self.assertEqual(run.call_args.args[2], 'lp_simplex')
            self.assertEqual(results, [self.row, highs])

    def test_rerun_is_an_explicit_override(self):
        metadata = dict(environment=self.environment, result_cache=self.entries())
        self.config['rerun'] = ['lp_simplex']
        with TemporaryDirectory() as folder:
            path = Path(folder) / 'model.mps'
            path.write_text('model')
            self.model['sha256'] = benchmark.sha(path)
            with patch.object(benchmark, 'prepared_model', return_value=nullcontext(path)), \
                 patch.object(benchmark, 'run_one', return_value=self.row) as run:
                benchmark.benchmark_case(Path(folder) / 'config.json', self.config,
                                         self.model, 0, {}, [], metadata)
            run.assert_called_once()

    def test_summary_recomputes_comparisons_from_compatible_measurements(self):
        highs = dict(self.row, solver='highs', wall_seconds=2.)
        cache = self.entries([highs])
        with TemporaryDirectory() as folder:
            root = Path(folder)
            docs, run = root / 'docs', root / 'run'
            docs.mkdir()
            run.mkdir()
            (docs / 'netlib.md').write_text(benchmark.render_report('netlib', {}, {}, cache=cache))
            metadata = dict(config=self.config, environment=self.environment, completed=True,
                            **{key: False for key in ('source_changed_during_run', 'binary_changed_during_run',
                               'library_changed_during_run', 'model_changed_during_run', 'competitor_changed_during_run')})
            benchmark.dump(run / 'metadata.json', metadata)
            benchmark.dump(run / 'results.json', [dict(self.row, wall_seconds=1.)])
            benchmark.report(run, docs)
            text = (docs / 'netlib.md').read_text()
            self.assertIn('| HiGHS | 1/1 | 0 | 0 | 1 | 1 / 2 | 0.5 |', text)
            self.assertNotIn('历史结果', text)

    def test_explicit_import_checks_settings_and_never_adopts_our_old_build(self):
        self.config['solvers'] = list(benchmark.SOLVERS)
        versions = self.environment['versions']
        for solver in ('scip', 'clp', 'glpk'):
            versions[solver] = dict(version_output=solver + ' v1', binary_sha256=solver)
        versions['highs'].update(version='1.0', binding='test')
        results = [dict(self.row, solver=solver) for solver in benchmark.SOLVERS]
        block = benchmark.run_summary(Path('old'), dict(config=self.config, environment=self.environment), results)
        stored = {'afiro': '| afiro | 1 / 2 / 3 | 0 | ' + ' | '.join(['PASS / 1 / 0'] * 5) + ' |'}
        with TemporaryDirectory() as folder:
            path = Path(folder) / 'netlib.md'
            path.write_text(benchmark.render_report('netlib', stored, {'old': block}))
            cache = benchmark.adopt_report_results(path, self.config, self.environment)
            self.assertEqual({entry['result']['solver'] for entry in cache}, set(benchmark.SOLVERS[1:]))
            self.config['seconds'] = 20.
            with self.assertRaises(ValueError):
                benchmark.adopt_report_results(path, self.config, self.environment)
            self.config['seconds'] = 10.
            versions['clp']['version_output'] = 'clp v2'
            with self.assertRaises(ValueError):
                benchmark.adopt_report_results(path, self.config, self.environment)



if __name__ == '__main__':
    unittest.main()
