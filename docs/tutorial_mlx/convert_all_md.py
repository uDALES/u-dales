#!/usr/bin/env python3
"""Merged tools script for docs/tutorial_mlx

This single script replaces the previous separate helpers and provides
the full workflow in one place:

- normalize_markdown(text) -> str
- fix fences/lists/anchors
- lint-check
- collapse duplicate .post.md files
- process all tutorial .md files

Behavior: by default the script WILL overwrite `.md` files in the
directory (backups of originals are written to `archive/`). Use
--dry-run to write `.post.md` files instead and avoid changing
originals.

Location: docs/tutorial_mlx/convert_all_md.py
"""

from __future__ import annotations

import argparse
import re
import shutil
from pathlib import Path
from typing import List, Optional


SCRIPT_DIR = Path(__file__).resolve().parent
ARCHIVE_DIR = SCRIPT_DIR / 'archive'


def ensure_archive() -> None:
    ARCHIVE_DIR.mkdir(exist_ok=True)


def slugify_heading(h: str) -> str:
    s = h.strip().lower()
    s = re.sub(r'[^a-z0-9\s\-]', '', s)
    s = re.sub(r'\s+', '-', s).strip('-')
    return s


def fix_fence_spacing_text(text: str) -> str:
    lines = text.splitlines()
    out = []
    i = 0
    n = len(lines)
    while i < n:
        line = lines[i]
        m = re.match(r'^(?P<indent>\s*)(?P<fence>`{3,}.*)$', line)
        if m:
            indent = m.group('indent')
            fence_line = m.group('fence')
            if len(out) > 0 and out[-1].strip() != '':
                out.append(indent)
            out.append(indent + fence_line.rstrip())
            i += 1
            while i < n:
                out.append(lines[i].rstrip())
                if re.match(r'^(?P<i>\s*)`{3,}', lines[i]):
                    break
                i += 1
            if i + 1 < n and lines[i + 1].strip() != '':
                out.append(indent)
            i += 1
            continue
        else:
            out.append(line.rstrip())
            i += 1
    while len(out) > 0 and out[-1].strip() == '':
        out.pop()
    return '\n'.join(out) + '\n'


def fix_list_spacing_text(text: str) -> str:
    lines = text.splitlines()
    out = []
    in_list = False
    for i, line in enumerate(lines):
        is_list = bool(re.match(r'^\s*[-+*]\s+', line))
        if is_list and not in_list:
            if len(out) > 0 and out[-1].strip() != '':
                out.append('')
            out.append(line.rstrip())
            in_list = True
            continue
        if not is_list and in_list:
            out.append('')
            in_list = False
        out.append(line.rstrip())
    while len(out) > 0 and out[-1].strip() == '':
        out.pop()
    return '\n'.join(out) + '\n'


def insert_explicit_anchors(text: str) -> str:
    fragments = set(re.findall(r'\]\(#([^\)]+)\)', text))
    if not fragments:
        return text
    headings = []
    lines = text.splitlines()
    for idx, line in enumerate(lines):
        m = re.match(r'^(#{1,6})\s*(.+?)(?:\s*\{#([^}]+)\})?\s*$', line)
        if m:
            headings.append((idx, m.group(1), m.group(2).strip(), m.group(3)))

    slug_to_idx = {}
    for idx, hashes, title, explicit in headings:
        slug = slugify_heading(title)
        slug_to_idx.setdefault(slug, []).append((idx, hashes, title, explicit))

    changed = False

    def normalize_frag_name(f: str) -> str:
        f2 = re.sub(r'^(m|h|mh)[_\-]?', '', f, flags=re.IGNORECASE)
        f2 = re.sub(r'^(section|sec)[_\-]?', '', f2, flags=re.IGNORECASE)
        f2 = re.sub(r'[^a-z0-9]+', '-', f2.lower()).strip('-')
        return f2

    for frag in sorted(fragments):
        if frag in re.findall(r'\{#([^}]+)\}', text):
            continue
        norm_frag = normalize_frag_name(frag)
        if frag in slug_to_idx:
            for idx, hashes, title, explicit in slug_to_idx[frag]:
                if explicit:
                    break
                lines[idx] = f"{hashes} {title} {{#{frag}}}"
                changed = True
                break
            continue
        if norm_frag in slug_to_idx:
            for idx, hashes, title, explicit in slug_to_idx[norm_frag]:
                if explicit:
                    break
                lines[idx] = f"{hashes} {title} {{#{frag}}}"
                changed = True
                break
            if changed:
                continue

        pattern = re.compile(rf'^\s*<a\s+id="{re.escape(frag)}"\s*>\s*</a>\s*$', flags=re.IGNORECASE)
        for idx, ln in enumerate(lines):
            if pattern.match(ln):
                def find_heading_near(k):
                    for j in range(k, min(len(lines), k + 13)):
                        if re.match(r'^(#{1,6})\s+', lines[j]):
                            return j
                    for j in range(k, max(-1, k - 13), -1):
                        if re.match(r'^(#{1,6})\s+', lines[j]):
                            return j
                    return None

                hidx = find_heading_near(idx + 1)
                if hidx is None:
                    hidx = find_heading_near(idx - 1)
                if hidx is not None:
                    m = re.match(r'^(#{1,6})\s*(.+?)(?:\s*\{#([^}]+)\})?\s*$', lines[hidx])
                    if m:
                        hashes = m.group(1)
                        title = m.group(2).strip()
                        explicit = m.group(3)
                        if not explicit:
                            lines[hidx] = f"{hashes} {title} {{#{frag}}}"
                            changed = True
                            lines[idx] = ''
                            break
    if changed:
        return '\n'.join(lines) + '\n'
    return text


def normalize_markdown(text: str) -> str:
    lines = text.splitlines()
    out_lines: List[str] = []
    in_fence = False
    pending_anchor: Optional[str] = None
    seen_first_h1 = False
    in_list_block = False

    def push_line(line: str):
        if line.strip() == "":
            if len(out_lines) >= 2 and out_lines[-1].strip() == "" and out_lines[-2].strip() == "":
                return
        out_lines.append(line)

    i = 0
    while i < len(lines):
        raw = lines[i]
        line = raw
        if re.match(r'^```{3,}', line):
            fence_start = re.match(r'^(?P<ticks>`{3,})(?P<lang>.*)$', line)
            if fence_start:
                lang = fence_start.group('lang').strip()
                if lang.lower().startswith('matlabtextoutput'):
                    lang = 'text'
                line = '```' + (('' if lang == '' else lang))

            if not in_fence and len(out_lines) > 0 and out_lines[-1].strip() != "":
                out_lines.append("")

            push_line(line.rstrip())
            if not in_fence:
                in_fence = True
            else:
                in_fence = False
                if i + 1 < len(lines) and lines[i+1].strip() != "":
                    out_lines.append("")
            i += 1
            in_list_block = False
            continue

        if in_fence:
            push_line(line.rstrip())
            i += 1
            continue

        if i == 0 and re.match(r'^```\s*markdown\s*$', line, flags=re.IGNORECASE):
            i += 1
            while i < len(lines) and lines[i].strip() == "":
                i += 1
            continue

        m_anchor = re.match(r'^\s*<a\s+id="([^\"]+)"\s*>\s*</a>\s*$', line)
        if m_anchor:
            pending_anchor = m_anchor.group(1)
            i += 1
            continue

        m_heading = re.match(r'^(?P<hashes>#{1,6})\s*(?P<title>.*?)(\s*\{#.*\})?\s*$', line)
        if m_heading:
            hashes = m_heading.group('hashes')
            title = m_heading.group('title').rstrip()
            if hashes == '#':
                if not seen_first_h1:
                    seen_first_h1 = True
                else:
                    hashes = '##'
            if len(out_lines) > 0 and out_lines[-1].strip() != "":
                out_lines.append("")
            if pending_anchor:
                line = f"{hashes} {title} {{#{pending_anchor}}}"
                pending_anchor = None
            else:
                line = f"{hashes} {title}"
            push_line(line.rstrip())
            if i + 1 < len(lines):
                nxt = lines[i+1].strip()
                if nxt != '' and not re.match(r'^[-+*]\s+', nxt) and not re.match(r'^```', nxt):
                    out_lines.append("")
            i += 1
            in_list_block = False
            continue

        m_list = re.match(r'^(?P<indent>\s*)(?P<marker>[-+*])\s+(?P<rest>.*)$', line)
        if m_list:
            indent = m_list.group('indent')
            marker = m_list.group('marker')
            rest = m_list.group('rest').rstrip()
            level = max(0, len(indent) // 2)
            new_indent = '  ' * level
            if not in_list_block:
                if len(out_lines) > 0 and out_lines[-1].strip() != "":
                    out_lines.append("")
                in_list_block = True
            new_line = f"{new_indent}{marker} {rest}"
            push_line(new_line.rstrip())
            i += 1
            if i < len(lines):
                m_next = re.match(r'^\s*[-+*]\s+', lines[i])
                if not m_next:
                    if len(out_lines) == 0 or out_lines[-1].strip() != "":
                        out_lines.append("")
                    in_list_block = False
            else:
                if len(out_lines) == 0 or out_lines[-1].strip() != "":
                    out_lines.append("")
                in_list_block = False
            continue

        line = re.sub(r'\\-', '-', line)
        line = re.sub(r'\\_', '_', line)
        line = re.sub(r'\\([\`*_{}\[\]()#+\-.!])', r'\1', line)
        line = line.rstrip()
        if in_list_block and line.strip() == "":
            if len(out_lines) == 0 or out_lines[-1].strip() != "":
                out_lines.append("")
            in_list_block = False
        push_line(line)
        i += 1

    while len(out_lines) > 0 and out_lines[-1].strip() == "":
        out_lines.pop()
    result = "\n".join(out_lines) + "\n"

    for _ in range(2):
        prev = result
        result = fix_fence_spacing_text(result)
        result = fix_list_spacing_text(result)
        result = insert_explicit_anchors(result)
        result = '\n'.join([ln.rstrip() for ln in result.splitlines()]) + '\n'
        if result == prev:
            break
    return result


def cmd_lint_check(file_path: Path) -> int:
    text = file_path.read_text(encoding='utf8')
    lines = text.splitlines()
    issues = []
    for i, l in enumerate(lines):
        if l.endswith(' '):
            issues.append((i + 1, 'Trailing whitespace'))
    h1s = [i + 1 for i, l in enumerate(lines) if re.match(r'^#\s+[^#]', l)]
    if len(h1s) > 1:
        issues.append((h1s[1], f'Multiple H1 headings (first at line {h1s[0]})'))
    fence_positions = [i for i, l in enumerate(lines) if l.strip().startswith('```')]
    for pos in fence_positions:
        if pos > 0 and lines[pos - 1].strip() != '':
            issues.append((pos + 1, 'No blank line before fenced code block'))
        j = pos + 1
        while j < len(lines) and not lines[j].strip().startswith('```'):
            j += 1
        if j < len(lines):
            if j + 1 < len(lines) and lines[j + 1].strip() != '':
                issues.append((j + 1, 'No blank line after fenced code block'))
    in_list = False
    for i, l in enumerate(lines):
        if re.match(r'^\s*[-+*]\s+', l):
            if not in_list:
                if i > 0 and lines[i - 1].strip() != '':
                    issues.append((i + 1, 'No blank line before list block'))
                in_list = True
        else:
            if in_list:
                if l.strip() != '':
                    issues.append((i + 1, 'No blank line after list block'))
                in_list = False
    fragments = re.findall(r'\]\(#([^\)]+)\)', text)
    anchors = set(re.findall(r'\{#([^}]+)\}', text))
    for m in re.finditer(r'^(#{1,6})\s*(.+?)(?:\s*\{#([^}]+)\})?\s*$', text, flags=re.M):
        title = m.group(2).strip()
        explicit = m.group(3)
        if explicit:
            anchors.add(explicit)
        else:
            anchors.add(slugify_heading(title))
    for f in sorted(set(fragments)):
        if f not in anchors:
            issues.append((0, f'Missing anchor for fragment: {f}'))
    if not issues:
        print('No issues found in', file_path)
        return 0
    print(f'Found {len(issues)} issues in {file_path}:')
    for lineno, msg in issues:
        if lineno:
            print(f'  Line {lineno}: {msg}')
        else:
            print(f'  {msg}')
    return 1


def collapse_duplicate_post_md() -> None:
    ensure_archive()
    changed = []
    for p in sorted(SCRIPT_DIR.iterdir()):
        if p.is_file() and '.post.md.post.md' in p.name:
            parts = p.name.split('.post.md')
            base = parts[0]
            target_name = base + '.post.md'
            target = SCRIPT_DIR / target_name
            if target.exists():
                bak = ARCHIVE_DIR / (target.name + '.dup.bak')
                if not bak.exists():
                    shutil.copy2(target, bak)
            shutil.copy2(p, target)
            dupbak = ARCHIVE_DIR / (p.name + '.bak')
            if not dupbak.exists():
                shutil.copy2(p, dupbak)
            try:
                p.unlink()
            except Exception:
                pass
            changed.append((p.name, target.name))
    if changed:
        print('Collapsed duplicates:')
        for s, t in changed:
            print(' -', s, '->', t)
    else:
        print('No duplicate .post.md files found')


def collect_md_files() -> list[Path]:
    files = []
    for p in sorted(SCRIPT_DIR.glob('*.md')):
        if p.name in ('convert_all_md.py', 'm2md.m'):
            continue
        files.append(p)
    return files


def backup(path: Path) -> Path:
    ensure_archive()
    target = ARCHIVE_DIR / (path.name + '.orig')
    if target.exists():
        i = 1
        while True:
            candidate = ARCHIVE_DIR / f"{path.name}.orig.{i}"
            if not candidate.exists():
                target = candidate
                break
            i += 1
    shutil.copy2(path, target)
    return target


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog='convert_all_md')
    parser.add_argument('--dry-run', action='store_true', help='Do not overwrite originals; write .post.md files instead')
    parser.add_argument('--collapse-duplicates', action='store_true', help='Collapse .post.md.post.md duplicates before processing')
    args = parser.parse_args(argv)

    if args.collapse_duplicates:
        collapse_duplicate_post_md()

    md_files = collect_md_files()
    if not md_files:
        print('No markdown files found in', SCRIPT_DIR)
        return 0

    total_changed = 0
    total_issues = 0

    for p in md_files:
        print('Processing', p.name)
        text = p.read_text(encoding='utf8')
        new = normalize_markdown(text)

        out_path = p if not args.dry_run else p.with_suffix(p.suffix + '.post.md')

        if new != text:
            total_changed += 1
            if not args.dry_run:
                b = backup(p)
                p.write_text(new, encoding='utf8')
                print('  Updated (backup at {})'.format(b.name))
            else:
                out_path.write_text(new, encoding='utf8')
                print('  Wrote dry-run:', out_path.name)
        else:
            print('  No changes')

        target_for_fence = out_path
        fixed = fix_fence_spacing_text(target_for_fence.read_text(encoding='utf8'))
        if fixed != target_for_fence.read_text(encoding='utf8'):
            target_for_fence.write_text(fixed, encoding='utf8')
            print('  Fixed fences in', target_for_fence.name)

        rc = cmd_lint_check(target_for_fence)
        if rc:
            total_issues += 1

    print('\nSummary:')
    print('  Files processed:', len(md_files))
    print('  Files changed:', total_changed)
    print('  Files with lint issues (after fixes):', total_issues)

    if total_issues:
        print('\nNote: this is a conservative automated pass. For a final authoritative run,')
        print('execute your project markdownlint configuration and fix any remaining')
        print('issues interactively or provide the config so I can adapt the script.')

    return 1 if total_issues else 0


if __name__ == '__main__':
    raise SystemExit(main())
