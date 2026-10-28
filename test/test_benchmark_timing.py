"""CLI timing boundaries and paired performance summaries."""
import json
import math
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'data'))
import benchmark
import compare_performance


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



if __name__ == '__main__':
    unittest.main()
