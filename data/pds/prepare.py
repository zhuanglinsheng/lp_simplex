#!/usr/bin/env python3
"""Generate selected MPS files from local raw data using only Python's standard library."""
import argparse
import os
from pathlib import Path
import tempfile

ROOT = Path(__file__).resolve().parent


def select(parser, available):
    parser.add_argument('instances', nargs='*', help='Instance names (without compression/MPS suffix)')
    parser.add_argument('--all', action='store_true', help='Generate every instance')
    parser.add_argument('--output-dir', type=Path, default=DEFAULT_OUTPUT,
                        help='Output directory (default: %(default)s)')
    args = parser.parse_args()
    if args.all and args.instances:
        parser.error('Choose instance names or --all, not both')
    if not args.all and not args.instances:
        parser.error('Specify instance names or --all. Available: ' + ', '.join(available))
    selected = list(available) if args.all else args.instances
    for name in selected:
        if name not in available:
            parser.error('Unknown instance: ' + name + '. Available: ' + ', '.join(available))
    return args, list(dict.fromkeys(selected))


def atomic_write(destination, writer):
    destination.parent.mkdir(parents=True, exist_ok=True)
    fd, filename = tempfile.mkstemp(prefix='.' + destination.name + '.', dir=destination.parent)
    temporary = Path(filename)
    try:
        with os.fdopen(fd, 'wb') as output:
            writer(output)
        if temporary.stat().st_size == 0:
            raise ValueError('Empty MPS output: ' + str(destination))
        temporary.replace(destination)
    finally:
        temporary.unlink(missing_ok=True)

import bz2
import shlex
import shutil
import subprocess

DEFAULT_OUTPUT = ROOT / 'mps'


def main():
    available = {p.stem: p for p in sorted((ROOT / 'raw').glob('*.bz2'))}
    parser = argparse.ArgumentParser(description=__doc__)
    args, names = select(parser, available)
    with tempfile.TemporaryDirectory(prefix='lp-pds-') as directory:
        work = Path(directory)
        converter = work / 'emps'
        compiler = shlex.split(os.environ.get('CC', 'cc'))
        subprocess.run(compiler + ['-O2', '-Wno-parentheses', str(ROOT / 'raw' / 'emps.c'),
                                   '-o', str(converter)], check=True)
        for name in names:
            # Spooling avoids holding large MPC data in memory and ensures a bad
            # bzip2 stream cannot leave a converter process or partial MPS behind.
            with tempfile.TemporaryFile() as expanded:
                with bz2.open(available[name], 'rb') as source:
                    shutil.copyfileobj(source, expanded)
                expanded.seek(0)
                def write(output):
                    subprocess.run([str(converter)], stdin=expanded, stdout=output, check=True)
                destination = args.output_dir / (name + '.mps')
                atomic_write(destination, write)
                print(destination)


if __name__ == '__main__':
    main()
