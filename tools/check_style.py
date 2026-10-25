#!/usr/bin/env python3
"""Check first-party formatting, file headers and script whitespace."""
import argparse
import io
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tokenize

ROOT = Path(__file__).resolve().parents[1]
EXTENSIONS = {'.c', '.h', '.cpp', '.hpp'}


def source_files():
    """Exclude downloaded data sources, build products and local experiments."""
    paths = [path for folder in ('src', 'include', 'test')
             for path in (ROOT / folder).rglob('*')
             if path.is_file() and path.suffix in EXTENSIONS]
    paths.extend(path for path in (ROOT / 'data').iterdir()
                 if path.is_file() and path.suffix in EXTENSIONS)
    paths.extend(path for path in (ROOT / 'data/simple_examples').rglob('*')
                 if path.is_file() and path.suffix in EXTENSIONS)
    return sorted(paths)


def cmake_files():
    paths = {ROOT / 'CMakeLists.txt'}
    paths.update((ROOT / 'cmake').rglob('*.cmake*'))
    for folder in ('src', 'include', 'test', 'data'):
        for pattern in ('CMakeLists.txt', '*.cmake'):
            paths.update(path for path in (ROOT / folder).rglob(pattern)
                         if 'raw' not in path.relative_to(ROOT).parts)
    return sorted(path for path in paths if path.is_file())


def auxiliary_files():
    paths = {ROOT / '.cmake-format.json'}
    for folder in ('tools', 'test', 'data', '.github'):
        paths.update(path for path in (ROOT / folder).rglob('*')
                     if path.is_file() and path.suffix in
                     {'.py', '.yaml', '.yml', '.json'} and
                     'raw' not in path.relative_to(ROOT).parts)
    return sorted(paths)


def whitespace_issues(raw):
    problems = []
    if b'\r' in raw:
        problems.append('expected LF line endings')
    if not raw.endswith(b'\n'):
        problems.append('missing final newline')
    if any(line.rstrip(' \t') != line
           for line in raw.decode('utf-8').splitlines()):
        problems.append('trailing whitespace')
    return problems


def auxiliary_issues(path):
    raw = path.read_bytes()
    text = raw.decode('utf-8')
    problems = whitespace_issues(raw)
    if any('\t' in line[:len(line) - len(line.lstrip())]
           for line in text.splitlines()):
        problems.append('expected spaces for indentation')
    if path.suffix == '.py':
        try:
            # Continuation alignment and multiline literal contents are not
            # suite indentation. Let Python identify actual block prefixes.
            tokens = tokenize.generate_tokens(io.StringIO(text).readline)
            levels = [0]
            for token in tokens:
                if token.type == tokenize.INDENT:
                    width = len(token.string)
                    if set(token.string) != {' '} or width != levels[-1] + 4:
                        problems.append('expected four spaces per Python block level')
                    levels.append(width)
                elif token.type == tokenize.DEDENT:
                    levels.pop()
            compile(text, str(path), 'exec')
        except (SyntaxError, tokenize.TokenError) as error:
            problems.append('invalid Python: ' + str(error))
    elif path.suffix == '.json':
        try:
            expected = json.dumps(json.loads(text), ensure_ascii=False,
                                  indent=2) + '\n'
            if text != expected:
                problems.append('expected two-space JSON formatting')
        except ValueError as error:
            problems.append('invalid JSON: ' + str(error))
    return problems


def checked_formatter(parser, requested, version_file, package):
    executable = shutil.which(requested)
    expected = (ROOT / version_file).read_text().strip()
    if executable is None:
        parser.error(requested + ' is missing; install ' +
                     package + '==' + expected)
    result = subprocess.run([executable, '--version'], capture_output=True,
                            text=True, check=True)
    version = re.search(r'(\d+\.\d+\.\d+)', result.stdout)
    if not version or version.group(1) != expected:
        parser.error('expected ' + requested + ' ' + expected + '; got ' +
                     result.stdout.strip())
    return executable


def file_issues(path):
    raw = path.read_bytes()
    text = raw.decode('utf-8')
    problems = []
    leading = re.match(r'(?:\s*/\*.*?\*/)*', text, re.S).group()
    documentation = re.findall(r'/\*\*.*?\*/', leading, re.S)
    header_gap = re.match(r'/\*(?!\*).*?\*/(\s*)/\*\*', text, re.S)
    if header_gap and header_gap.group(1) != '\n\n':
        problems.append('expected one blank line before file documentation')
    if not any(re.search(r'@file\s+' + re.escape(path.name) + r'\s', block)
               and re.search(r'@brief\s+\S', block)
               for block in documentation):
        problems.append('missing leading @file basename / @brief block')
    problems.extend(whitespace_issues(raw))
    return problems


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '--fix', action='store_true',
        help='format files in place; header text needs human review')
    parser.add_argument('--clang-format', default='clang-format',
                        help='formatter executable (default: PATH)')
    parser.add_argument('--cmake-format', default='cmake-format',
                        help='CMake formatter executable (default: PATH)')
    args = parser.parse_args()
    formatter = checked_formatter(parser, args.clang_format,
                                  '.clang-format-version', 'clang-format')
    files = source_files()
    cmake = cmake_files()
    auxiliary = auxiliary_files()
    cmake_formatter = None
    if cmake:
        cmake_formatter = checked_formatter(
            parser, args.cmake_format, '.cmake-format-version', 'cmakelang')
    if args.fix:
        subprocess.run([formatter, '-i', *map(str, files)], check=True, cwd=ROOT)
        if cmake:
            subprocess.run([cmake_formatter, '-i', *map(str, cmake)],
                           check=True, cwd=ROOT)
    failed = False
    for path in files:
        for problem in file_issues(path):
            print(str(path.relative_to(ROOT)) + ': ' + problem, file=sys.stderr)
            failed = True
    for path in auxiliary:
        for problem in auxiliary_issues(path):
            print(str(path.relative_to(ROOT)) + ': ' + problem, file=sys.stderr)
            failed = True
    result = subprocess.run([formatter, '--dry-run', '--Werror',
                             *map(str, files)], cwd=ROOT)
    cmake_result = None
    if cmake:
        cmake_result = subprocess.run(
            [cmake_formatter, '--check', *map(str, cmake)], cwd=ROOT)
    if failed or result.returncode or (cmake_result and cmake_result.returncode):
        return 1
    print('Style check passed: ' + str(len(files)) + ' C/C++ files; clang-format ' +
          (ROOT / '.clang-format-version').read_text().strip() + '; ' +
          str(len(cmake)) + ' CMake files; ' + str(len(auxiliary)) +
          ' Python/YAML/JSON files; file headers and whitespace checked.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
