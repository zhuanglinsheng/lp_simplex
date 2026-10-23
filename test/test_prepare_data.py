#!/usr/bin/env python3
"""Exercise data preparation failures without touching repository datasets."""
import bz2
import gzip
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


class PrepareDataTests(unittest.TestCase):
    def setUp(self):
        self.work = tempfile.TemporaryDirectory()
        self.addCleanup(self.work.cleanup)
        self.root = Path(self.work.name)
        (self.root / 'raw').mkdir()
        self.output = self.root / 'output'
        self.output.mkdir()

    def run_script(self, dataset, *args):
        script = self.root / 'prepare.py'
        shutil.copyfile(ROOT / 'data' / dataset / 'prepare.py', script)
        return subprocess.run([sys.executable, str(script), *args,
                               '--output-dir', str(self.output)],
                              capture_output=True, text=True)

    def test_gzip_crc_failure_preserves_existing_mps(self):
        raw = gzip.compress(b'NAME          EXAMPLE\nENDATA\n')
        (self.root / 'raw' / 'example.mps.gz').write_bytes(raw[:-4])
        existing = self.output / 'example.mps'
        existing.write_bytes(b'existing model')
        result = self.run_script('netlib', 'example')
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(existing.read_bytes(), b'existing model')
        self.assertEqual(list(self.output.iterdir()), [existing])

    def test_unknown_instance_rejected_before_creating_output(self):
        result = self.run_script('netlib', '../escape')
        self.assertEqual(result.returncode, 2)
        self.assertIn('Unknown instance', result.stderr)
        self.assertEqual(list(self.output.iterdir()), [])

    def test_selection_required(self):
        result = self.run_script('railway')
        self.assertEqual(result.returncode, 2)
        self.assertIn('Specify instance names or --all', result.stderr)

    def test_truncated_railway_preserves_existing_mps(self):
        (self.root / 'raw' / 'railtest.gz').write_bytes(gzip.compress(b'2 1\n1 2 1\n'))
        existing = self.output / 'railtest.mps'
        existing.write_bytes(b'existing model')
        result = self.run_script('railway', 'railtest')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('Truncated railway input', result.stderr)
        self.assertEqual(existing.read_bytes(), b'existing model')
        self.assertEqual(list(self.output.iterdir()), [existing])

    def test_pds_compile_failure_preserves_existing_mps(self):
        (self.root / 'raw' / 'pds-test.bz2').write_bytes(bz2.compress(b'input'))
        (self.root / 'raw' / 'emps.c').write_text('#error intentional compile failure\n')
        existing = self.output / 'pds-test.mps'
        existing.write_bytes(b'existing model')
        result = self.run_script('pds', 'pds-test')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('intentional compile failure', result.stderr)
        self.assertEqual(existing.read_bytes(), b'existing model')
        self.assertEqual(list(self.output.iterdir()), [existing])

    def test_railway_generates_continuous_bounded_model(self):
        (self.root / 'raw' / 'railtest.gz').write_bytes(gzip.compress(b'2 1\n2 2 1 2\n'))
        result = self.run_script('railway', 'railtest')
        self.assertEqual(result.returncode, 0, result.stderr)
        text = (self.output / 'railtest.mps').read_text()
        self.assertIn(' G  R0000001\n', text)
        self.assertIn(' G  R0000002\n', text)
        self.assertIn(' UP BND1      X0000001             1\n', text)
        self.assertNotIn('INTORG', text)
        self.assertTrue(text.endswith('ENDATA\n'))


if __name__ == '__main__':
    unittest.main()
