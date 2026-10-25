#!/usr/bin/env python3
"""Check source ownership boundaries without requiring a solver build."""
from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'src'
INCLUDE = re.compile(r'^\s*#\s*include\s*"([^"\n]+)"', re.MULTILINE)
KERNEL_ENTRIES = {
    'simplex_dual.h', 'simplex_dual_kernel.h', 'simplex_pan.h',
    'simplex_pan_bounded.h', 'simplex_tableau_solver.h',
}


def check(root=SOURCE):
    """Resolve local includes and reject dependencies on orchestration state."""
    paths = sorted(p for p in root.rglob('*') if p.suffix in {'.c', '.h'})
    headers = {}
    errors = []
    for path in paths:
        if path.suffix == '.h':
            headers.setdefault(path.name, []).append(path)
    for path in paths:
        owner = path.relative_to(root).parts[0]
        for name in INCLUDE.findall(path.read_text()):
            direct = (path.parent / name).resolve()
            matches = [direct] if direct.is_file() else headers.get(Path(name).name, [])
            if len(matches) != 1:
                errors.append(f'{path}: unresolved or ambiguous include {name}')
                continue
            target = matches[0]
            try:
                target_owner = target.relative_to(root).parts[0]
            except ValueError:
                errors.append(f'{path}: private include escapes src: {name}')
                continue
            if owner != 'core' and target_owner == 'core':
                errors.append(f'{path}: {owner} must not depend on core: {name}')
            if owner == 'runtime' and target.name in KERNEL_ENTRIES:
                errors.append(f'{path}: runtime must not invoke a solver: {name}')
            if target.name == 'simplex_auto_state.h' and path.name != 'simplex_auto.c':
                errors.append(f'{path}: AUTO state belongs to simplex_auto.c')
            if path.name == 'simplex_cover.c' and target_owner in {'bda', 'degeneracy'}:
                errors.append(f'{path}: cover must not borrow the BDA implementation')
    return paths, errors


def main():
    paths, errors = check()
    if errors:
        print('\n'.join(errors), file=sys.stderr)
        return 1
    print(f'Architecture check passed: {len(paths)} source files; '
          'core ownership, runtime/solver boundary and private AUTO state.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
