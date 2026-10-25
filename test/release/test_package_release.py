#!/usr/bin/env python3
"""Windows archive naming and transitive DLL closure, without a Windows host."""
import importlib.util
import json
import os
import shutil
from pathlib import Path
import tempfile
import unittest
from unittest import mock
import zipfile

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location('packager', ROOT / 'tools/package_release.py')
PACKAGER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(PACKAGER)


class WindowsArchive(unittest.TestCase):
    def test_name_and_dependency_closure(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            build = root / 'build'
            build.mkdir()
            (build / 'CMakeCache.txt').write_text(
                'CMAKE_BUILD_TYPE:STRING=Release\nCMAKE_PROJECT_VERSION:STATIC=1.0.0\n'
                'CMAKE_C_COMPILER:FILEPATH=cl.exe\n')
            (build / 'lp_simplexConfig.cmake').write_text(
                ''.join(f'set(lp_simplex_WITH_{key} OFF)\n' for key in ['KLU', 'SPQR', 'GMP']))
            for name in ['LICENSE', 'LICENSE.GPL-3.0', 'INSTALL.md', 'CHANGELOG.md',
                         'README.md', 'README_en.md']:
                (root / name).write_text('fixture\n')
            shutil.copytree(ROOT / 'tools/licenses/windows-runtime',
                            root / 'tools/licenses/windows-runtime')
            dependency = root / 'dependencies'
            (dependency / 'bin').mkdir(parents=True)
            for name in ['blas.dll', 'fortran.dll', 'unused.dll']:
                (dependency / 'bin' / name).write_bytes(b'fixture')
            system = root / 'system' / 'System32'
            system.mkdir(parents=True)
            (system / 'KERNEL32.dll').write_bytes(b'system')
            dependencies = {'lp_simplex.dll': 'blas.dll', 'blas.dll': 'fortran.dll',
                            'fortran.dll': 'KERNEL32.dll'}

            def output(*command):
                if command[:2] == ('git', 'status'):
                    return ''
                if command[:2] == ('git', 'rev-parse'):
                    return 'a' * 40
                if command[0] == 'cmake':
                    return 'cmake version fixture'
                if command[0] == 'dumpbin':
                    name = Path(command[-1]).name
                    # The last dumpbin line previously overwrote the archive name.
                    return f'Image has dependencies:\n    {dependencies[name]}\nSummary\n    181000 .text'
                raise AssertionError(command)

            def install(command, **kwargs):
                self.assertEqual(command[:2], ['cmake', '--install'])
                stage = Path(command[-1])
                (stage / 'bin').mkdir(parents=True)
                (stage / 'bin' / 'lp_simplex.dll').write_bytes(b'fixture')

            arguments = ['package_release.py', '--build-dir', str(build),
                         '--output-dir', str(root / 'dist'), '--platform', 'windows-x64',
                         '--dependency-prefix', str(dependency), '--require-clean']
            with mock.patch.object(PACKAGER, 'ROOT', root), \
                    mock.patch.object(PACKAGER, 'output', side_effect=output), \
                    mock.patch.object(PACKAGER.platform, 'system', return_value='Windows'), \
                    mock.patch.object(PACKAGER.platform, 'machine', return_value='AMD64'), \
                    mock.patch.object(PACKAGER.subprocess, 'run', side_effect=install), \
                    mock.patch.dict(os.environ, {'SystemRoot': str(system.parent)}), \
                    mock.patch('sys.argv', arguments):
                PACKAGER.main()
            name = 'lp_simplex-1.0.0-windows-x64'
            archive = root / 'dist' / (name + '.zip')
            self.assertTrue(archive.is_file())
            with zipfile.ZipFile(archive) as bundle:
                dlls = {Path(path).name for path in bundle.namelist() if path.endswith('.dll')}
                self.assertEqual(dlls, {'lp_simplex.dll', 'blas.dll', 'fortran.dll'})
                info = json.loads(bundle.read(name + '/BUILD-INFO.json'))
                self.assertEqual(info['platform'], 'windows-x64')
                self.assertFalse(info['dirty'])
                for notice in ['COPYING.RUNTIME', 'COPYING.LIB', 'NOTICE.txt']:
                    self.assertIn(name + '/share/licenses/third-party/gcc-runtime/' + notice,
                                  bundle.namelist())


if __name__ == '__main__':
    unittest.main()
