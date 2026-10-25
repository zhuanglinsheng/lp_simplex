#!/usr/bin/env python3
"""Install a tested Release build and archive its resolved capabilities."""
import argparse
import json
import platform
from pathlib import Path
import re
import shutil
import subprocess
import tarfile

ROOT = Path(__file__).resolve().parents[1]


def output(*command):
    return subprocess.check_output(command, cwd=ROOT, text=True).strip()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build-dir', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, default=ROOT / 'build-dist')
    parser.add_argument('--require-clean', action='store_true')
    parser.add_argument('--platform', required=True,
                        choices=['linux-x64', 'macos-x64', 'macos-arm64'])
    args = parser.parse_args()
    native = {('Linux', 'x86_64'): 'linux-x64',
              ('Darwin', 'x86_64'): 'macos-x64',
              ('Darwin', 'arm64'): 'macos-arm64'}.get(
                  (platform.system(), platform.machine()))
    if args.platform != native:
        parser.error(f'Expected native platform {native}')
    dirty = bool(output('git', 'status', '--porcelain', '--untracked-files=normal'))
    if args.require_clean and dirty:
        parser.error('Official packages require a clean checkout')
    build = args.build_dir.resolve()
    cache = {}
    for line in (build / 'CMakeCache.txt').read_text().splitlines():
        match = re.match(r'([^#/:][^:]*):[^=]+=(.*)$', line)
        if match:
            cache[match[1]] = match[2]
    if cache.get('CMAKE_BUILD_TYPE') != 'Release':
        parser.error('Only Release builds may be packaged')
    version = cache['CMAKE_PROJECT_VERSION']
    name = f'lp_simplex-{version}-{args.platform}'
    destination = args.output_dir.resolve()
    destination.mkdir(parents=True, exist_ok=True)
    stage = destination / name
    if stage.exists():
        parser.error(f'Refusing to overwrite {stage}')
    subprocess.run(['cmake', '--install', str(build), '--prefix', str(stage)],
                   check=True)
    for filename in ['LICENSE', 'LICENSE.GPL-3.0', 'INSTALL.md', 'CHANGELOG.md',
                     'README.md', 'README_en.md']:
        shutil.copy2(ROOT / filename, stage / filename)
    config = (build / 'lp_simplexConfig.cmake').read_text()
    capabilities = dict(re.findall(r'set\(lp_simplex_WITH_(\w+) (ON|OFF)\)', config))
    if set(capabilities) != {'KLU', 'SPQR', 'GMP'}:
        raise RuntimeError('Incomplete installed capability metadata')
    dependencies = {key: value for key, value in cache.items()
                    if re.search(r'(BLAS|LAPACK|GMP|KLU|SPQR|CHOLMOD|SUITESPARSE)', key)
                    and (key.endswith('_LIBRARY') or key.startswith('LP_SIMPLEX_'))
                    and not value.endswith('-NOTFOUND')}
    libraries = sorted(stage.rglob('liblp_simplex.*'))
    library = next((path for path in libraries if not path.is_symlink()), None)
    if library is None:
        raise RuntimeError('Installed shared library is missing')
    linkage = output('otool', '-L', str(library)) if platform.system() == 'Darwin' else output('ldd', str(library))
    if str(ROOT) in '\n'.join(linkage.splitlines()[1:]):
        raise RuntimeError('Installed library retains a checkout/build path')
    if 'not found' in linkage:
        raise RuntimeError('Installed library has unresolved dynamic dependencies')
    manifest = {
        'version': version, 'platform': args.platform,
        'commit': output('git', 'rev-parse', 'HEAD'),
        'dirty': dirty,
        'cmake': output('cmake', '--version').splitlines()[0],
        'compiler': output(cache['CMAKE_C_COMPILER'], '--version').splitlines()[0],
        'capabilities': capabilities,
        'build_dependency_paths': dependencies,
        'runtime_dependencies': linkage.splitlines()[1:] if platform.system() == 'Darwin' else linkage.splitlines(),
    }
    (stage / 'BUILD-INFO.json').write_text(json.dumps(manifest, indent=2) + '\n')
    archive = destination / f'{name}.tar.gz'
    with tarfile.open(archive, 'w:gz') as bundle:
        bundle.add(stage, arcname=name)
    print(archive)


if __name__ == '__main__':
    main()
