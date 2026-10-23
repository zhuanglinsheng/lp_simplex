"""Regression checks for audit failures and benchmark accounting."""
import json
import math
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest

from benchmark_solvers import aggregate, artifact, audit, parse_scip_solution, report


class BenchmarkTests(unittest.TestCase):
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


if __name__ == '__main__':
    unittest.main()
