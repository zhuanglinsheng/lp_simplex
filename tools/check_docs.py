#!/usr/bin/env python3
"""Check repository documentation links and local Markdown fragments."""
from html import unescape
from pathlib import Path
import re
import sys
from urllib.parse import unquote, urlsplit

ROOT = Path(__file__).resolve().parents[1]
LINK = re.compile(r'!?\[[^\]\n]*\]\((<[^>]+>|[^\s)]+)(?:\s+"[^"]*")?\)')
HEADING = re.compile(r'^ {0,3}#{1,6}\s+(.+?)\s*#*\s*$')
ANCHOR = re.compile(r'<a\s+[^>]*(?:id|name)=["\']([^"\']+)["\']', re.I)


def prose_lines(text):
    """Ignore fenced examples, whose Markdown-looking text is not a link."""
    fence = None
    for number, line in enumerate(text.splitlines(), 1):
        marker = re.match(r'^\s{0,3}(`{3,}|~{3,})', line)
        if marker:
            value = marker.group(1)
            if fence is None:
                fence = value
            elif value[0] == fence[0] and len(value) >= len(fence):
                fence = None
            continue
        if fence is None:
            yield number, line


def anchors(text):
    result = set()
    counts = {}
    for _, line in prose_lines(text):
        result.update(ANCHOR.findall(line))
        heading = HEADING.match(line)
        if not heading:
            continue
        label = re.sub(r'\[([^\]]+)\]\([^)]*\)', r'\1', heading.group(1))
        label = re.sub(r'<[^>]*>', '', unescape(label)).lower()
        slug = re.sub(r'[^\w\- ]', '', label).replace(' ', '-')
        count = counts.get(slug, 0)
        candidate = slug if not count else f'{slug}-{count}'
        while candidate in result:
            count += 1
            candidate = f'{slug}-{count}'
        counts[slug] = count + 1
        result.add(candidate)
    return result


def documents(root):
    files = {root / 'README.md', root / 'README_en.md'}
    files.update((root / 'docs').rglob('*.md'))
    files.update(path for path in (root / 'data').rglob('README*.md')
                 if 'raw' not in path.relative_to(root).parts)
    return sorted(path for path in files if path.is_file())


def check(root=ROOT):
    root = root.resolve()
    files = documents(root)
    cache = {}
    errors = []
    checked = 0
    for path in files:
        for number, line in prose_lines(path.read_text()):
            for match in LINK.finditer(line):
                href = match.group(1).strip('<>')
                parts = urlsplit(href)
                if parts.scheme or parts.netloc:
                    continue
                checked += 1
                target = (path.parent / unquote(parts.path)).resolve() if parts.path else path
                location = f'{path.relative_to(root)}:{number}'
                try:
                    target.relative_to(root)
                except ValueError:
                    errors.append(f'{location}: target escapes repository {href}')
                    continue
                if not target.exists():
                    errors.append(f'{location}: missing target {href}')
                    continue
                fragment = unquote(parts.fragment)
                if fragment and target.suffix == '.md':
                    if target not in cache:
                        cache[target] = anchors(target.read_text())
                    if fragment not in cache[target]:
                        errors.append(f'{location}: missing fragment {href}')
    return files, checked, errors


def main():
    files, checked, errors = check()
    if errors:
        print('\n'.join(errors), file=sys.stderr)
        return 1
    print(f'Documentation check passed: {len(files)} pages; '
          f'{checked} local links and fragments.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
