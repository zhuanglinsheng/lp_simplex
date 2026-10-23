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

import gzip

DEFAULT_OUTPUT = ROOT / 'mps'


def convert(source, output):
    # Stream both compressed input and MPS output, including million-column cases.
    with gzip.open(source, 'rt', encoding='ascii') as stream:
        tokens = (int(value) for line in stream for value in line.split())
        def take():
            try:
                return next(tokens)
            except StopIteration:
                raise ValueError('Truncated railway input: ' + str(source)) from None
        m, n = take(), take()
        if not (0 < m < 10_000_000 and 0 < n < 10_000_000):
            raise ValueError('Dimensions exceed the eight-character MPS name convention')
        def emit(text):
            output.write(text.encode('ascii'))
        emit(f'NAME          {source.stem}\nROWS\n N  OBJ\n')
        for i in range(1, m + 1):
            emit(f' G  R{i:07d}\n')
        emit('COLUMNS\n')
        covered = set()
        for j in range(1, n + 1):
            cost, degree = take(), take()
            if cost not in (1, 2) or not 1 <= degree <= 12:
                raise ValueError(f'Invalid cost/degree in column {j}')
            rows = [take() for _ in range(degree)]
            if len(set(rows)) != degree or not all(1 <= i <= m for i in rows):
                raise ValueError(f'Invalid row indices in column {j}')
            covered.update(rows)
            column = f'X{j:07d}'
            emit(f'    {column:<8}  {"OBJ":<8}  {cost:>12}\n')
            for i in rows:
                emit(f'    {column:<8}  R{i:07d}  {1:>12}\n')
        if next(tokens, None) is not None:
            raise ValueError('Unexpected trailing data')
        if len(covered) != m:
            raise ValueError('Input has uncovered rows')
        emit('RHS\n')
        for i in range(1, m + 1):
            emit(f'    {"RHS1":<8}  R{i:07d}  {1:>12}\n')
        emit('BOUNDS\n')
        for j in range(1, n + 1):
            emit(f' UP {"BND1":<8}  X{j:07d}  {1:>12}\n')
        emit('ENDATA\n')


def main():
    available = {p.stem: p for p in sorted((ROOT / 'raw').glob('rail*.gz'))}
    parser = argparse.ArgumentParser(description=__doc__)
    args, names = select(parser, available)
    for name in names:
        destination = args.output_dir / (name + '.mps')
        atomic_write(destination, lambda output: convert(available[name], output))
        print(destination)


if __name__ == '__main__':
    main()
