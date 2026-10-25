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
import os

ROOT = Path(__file__).resolve().parents[1]


def output(*command):
    return subprocess.check_output(command, cwd=ROOT, text=True).strip()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build-dir', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, default=ROOT / 'build-dist')
    parser.add_argument('--require-clean', action='store_true')
    parser.add_argument('--dependency-prefix', type=Path)
    parser.add_argument('--platform', required=True,
                        choices=['linux-x64', 'macos-x64', 'macos-arm64', 'windows-x64'])
    args = parser.parse_args()
    native = {('Linux', 'x86_64'): 'linux-x64',
              ('Darwin', 'x86_64'): 'macos-x64',
              ('Darwin', 'arm64'): 'macos-arm64',
              ('Windows', 'AMD64'): 'windows-x64'}.get(
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
    if platform.system() == 'Windows':
        if not args.dependency_prefix:
            parser.error('Windows packaging requires --dependency-prefix')
        dependency = args.dependency_prefix.resolve()
        for notice in (dependency / 'share').glob('*/copyright'):
            target = stage / 'share/licenses/third-party' / notice.parent.name
            target.mkdir(parents=True, exist_ok=True)
            shutil.copy2(notice, target / 'copyright')
        libraries = [stage / 'bin/lp_simplex.dll']
    else:
        libraries = sorted(stage.rglob('liblp_simplex.*'))
    library = next((path for path in libraries if not path.is_symlink()), None)
    if library is None:
        raise RuntimeError('Installed shared library is missing')
    if platform.system() == 'Windows':
        linkage = '\n'.join(line.strip() for line in output('dumpbin', '/DEPENDENTS', str(library)).splitlines()
                            if re.fullmatch(r'[A-Za-z0-9_.+-]+\.dll', line.strip(), re.I))
        available = {dll.name.lower(): dll for dll in (dependency / 'bin').glob('*.dll')}
        bundled = {library.name.lower()}
        system = Path(os.environ['SystemRoot']) / 'System32'
        pending = [library]
        while pending:
            dll = pending.pop()
            for line in output('dumpbin', '/DEPENDENTS', str(dll)).splitlines():
                name = line.strip()
                if not re.fullmatch(r'[A-Za-z0-9_.+-]+\.dll', name, re.I):
                    continue
                key = name.lower()
                if key in bundled:
                    continue
                if key in available:
                    target = stage / 'bin' / available[key].name
                    shutil.copy2(available[key], target)
                    bundled.add(key)
                    pending.append(target)
                    continue
                if ((system / name).exists() or
                        key.startswith(('api-ms-', 'ext-ms-'))):
                    continue
                raise RuntimeError(f'Unbundled Windows dependency: {name} ({dll.name})')
    elif platform.system() == 'Darwin':
        linkage = output('otool', '-L', str(library))
    else:
        linkage = output('ldd', str(library))
    if str(ROOT) in '\n'.join(linkage.splitlines()[1:]):
        raise RuntimeError('Installed library retains a checkout/build path')
    if 'not found' in linkage:
        raise RuntimeError('Installed library has unresolved dynamic dependencies')
    compiler_files = sorted(build.glob('CMakeFiles/*/CMakeCCompiler.cmake'))
    compiler_info = compiler_files[-1].read_text() if compiler_files else ''
    compiler_id = re.search(r'set\(CMAKE_C_COMPILER_ID "([^"]+)"\)', compiler_info)
    compiler_version = re.search(r'set\(CMAKE_C_COMPILER_VERSION "([^"]+)"\)', compiler_info)
    compiler = ' '.join(match[1] for match in [compiler_id, compiler_version] if match)
    manifest = {
        'version': version, 'platform': args.platform,
        'commit': output('git', 'rev-parse', 'HEAD'),
        'dirty': dirty,
        'cmake': output('cmake', '--version').splitlines()[0],
        'compiler': compiler or cache['CMAKE_C_COMPILER'],
        'capabilities': capabilities,
        'build_dependency_paths': dependencies,
        'runtime_dependencies': linkage.splitlines()[1:] if platform.system() == 'Darwin' else linkage.splitlines(),
    }
    (stage / 'BUILD-INFO.json').write_text(json.dumps(manifest, indent=2) + '\n')
    if platform.system() == 'Windows':
        archive = Path(shutil.make_archive(str(destination / name), 'zip',
                                          root_dir=destination, base_dir=name))
    else:
        archive = destination / f'{name}.tar.gz'
        with tarfile.open(archive, 'w:gz') as bundle:
            bundle.add(stage, arcname=name)
    print(archive)


if __name__ == '__main__':
    main()
