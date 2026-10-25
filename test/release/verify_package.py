#!/usr/bin/env python3
"""Verify the actual archive after extraction and relocation."""
import argparse
import json
import re
from pathlib import Path
import shutil
import subprocess
import tarfile
import tempfile

ROOT = Path(__file__).resolve().parents[2]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('archive', type=Path)
    parser.add_argument('--version', required=True)
    parser.add_argument('--require-full', action='store_true')
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix='lp-release-') as directory:
        work = Path(directory)
        original = work / 'extracted'
        original.mkdir()
        with tarfile.open(args.archive.resolve()) as bundle:
            bundle.extractall(original, filter='data')
        roots = list(original.iterdir())
        if len(roots) != 1 or not roots[0].is_dir():
            raise RuntimeError('Archive must contain exactly one package directory')
        package = work / 'relocated package'
        roots[0].rename(package)
        shutil.rmtree(original)
        manifest = json.loads((package / 'BUILD-INFO.json').read_text())
        if manifest['version'] != args.version:
            raise RuntimeError('Archive version mismatch')
        if args.require_full and manifest['capabilities'] != {
                'KLU': 'ON', 'SPQR': 'ON', 'GMP': 'ON'}:
            raise RuntimeError('Release omitted a required numerical backend')
        for name in ['LICENSE', 'LICENSE.GPL-3.0', 'INSTALL.md', 'CHANGELOG.md',
                     'include/lp_simplex/version.h']:
            if not (package / name).is_file():
                raise RuntimeError(f'Package omitted {name}')
        header = (package / 'include/lp_simplex/version.h').read_text()
        if f'#define LP_SIMPLEX_VERSION "{args.version}"' not in header:
            raise RuntimeError('Installed header version mismatch')
        configs = list(package.rglob('lp_simplexConfig.cmake'))
        if len(configs) != 1:
            raise RuntimeError('Expected exactly one installed package config')
        resolved = dict(re.findall(r'set\(lp_simplex_WITH_(\w+) (ON|OFF)\)',
                                   configs[0].read_text()))
        if resolved != manifest['capabilities']:
            raise RuntimeError('Manifest/installed capabilities mismatch')
        # Copy the consumer so its sources do not use the original checkout.
        source = work / 'consumer'
        shutil.copytree(ROOT / 'test/package_consumer', source)
        build = work / 'consumer-build'
        subprocess.run(['cmake', '-S', str(source), '-B', str(build),
                        f'-DPACKAGE_ROOT={package}',
                        f'-DEXPECTED_VERSION={args.version}', '-DTEST_CXX=ON',
                        '-DCMAKE_BUILD_TYPE=Release',
                        '-DCMAKE_FIND_USE_PACKAGE_REGISTRY=OFF',
                        '-DCMAKE_FIND_PACKAGE_NO_PACKAGE_REGISTRY=ON'], check=True,
                       cwd=work)
        subprocess.run(['cmake', '--build', str(build), '--parallel', '2'],
                       check=True, cwd=work)
        subprocess.run(['ctest', '--test-dir', str(build), '--output-on-failure'],
                       check=True, cwd=work)
    print(f'Extracted and relocated package {args.version}: C/C++ passed')


if __name__ == '__main__':
    main()
