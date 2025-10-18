#!/usr/bin/env python3
"""Convert and clean markdown files in docs/tutorial_mlx in a single pass.

This script combines the lint passes from md_convert.py and the
anchor/header analysis from analyze_md_links.py. It operates on the
original `.md` files (not `.post.md` intermediates) and overwrites them by
default. Use --backup to keep a .bak copy of the original.
"""
from __future__ import annotations

import os
import re
import unicodedata
from pathlib import Path
from typing import List, Dict, Optional, Tuple

# ---- slug generation (same rules used in analyze_md_links.py) ----
def mkdocs_slug(text: str) -> str:
    s = re.sub(r'\\([\\`*_{}\[\]()#+.\-!~])', r"\1", text)
    s = re.sub(r'[`*~]+', '', s)
    s = unicodedata.normalize('NFKD', s)
    s = ''.join(ch for ch in s if not unicodedata.combining(ch))
    s = s.lower()
    # remove periods entirely (user requested behaviour)
    s = s.replace('.', '')
    # preserve underscores and hyphens; replace other non-alnum with hyphen
    s = re.sub(r'[^a-z0-9_-]+', '-', s)
    s = s.strip('-')
    return s

# ---- lint passes from md_convert.py (inlined) ----
_ESC_RE = re.compile(r'(?<!\\)\\([_\*\[\]\-])')

def final_trailing_space_sweep(text: str) -> str:
    out = []
    for ln in text.splitlines(True):
        if ln.endswith('\n'):
            core, nl = ln[:-1], '\n'
        else:
            core, nl = ln, ''
        m = re.search(r'([ \t]+)$', core)
        if m:
            trailing = m.group(1)
            if '\t' in trailing:
                out.append(re.sub(r'[ \t]+$', '', core) + nl)
            else:
                if len(trailing) == 2:
                    out.append(core + nl)
                else:
                    out.append(re.sub(r'[ \t]+$', '', core) + nl)
        else:
            out.append(core + nl)
    return ''.join(out)

def collapse_blank_runs_to_one(text: str) -> str:
    return re.sub(r"\n(?:[ \t]*\n){1,}", "\n\n", text)

def normalize_list_marker_space(text: str) -> str:
    out = []
    in_code = False
    for raw in text.splitlines(True):
        line = raw
        if line.strip().startswith('```'):
            in_code = not in_code
            out.append(line)
            continue
        if not in_code and '<a' not in line.lower():
            m = re.match(r'^(\s*(?:>\s*)*)(\s*)([-*+]|\d+\.)\s+(.*\S.*)$', line.rstrip('\n'))
            if m:
                line = f"{m.group(1)}{m.group(2)}{m.group(3)} {m.group(4)}\n"
        out.append(line)
    return ''.join(out)

def blanks_around_headings(text: str) -> str:
    lines = text.splitlines()
    i = 0
    while i < len(lines):
        if re.match(r'^(#{1,6})\s+\S', lines[i]):
            if i > 0 and lines[i-1].strip() != '':
                lines.insert(i, '')
                i += 1
            if i+1 < len(lines) and lines[i+1].strip() != '':
                lines.insert(i+1, '')
        i += 1
    return '\n'.join(lines) + ('\n' if text.endswith('\n') else '')

def blanks_around_fences(text: str) -> str:
    lines = text.splitlines()
    i = 0
    while i < len(lines):
        if lines[i].startswith('```'):
            if i > 0 and lines[i-1].strip() != '':
                lines.insert(i, '')
                i += 1
            j = i + 1
            while j < len(lines) and not lines[j].startswith('```'):
                j += 1
            if j < len(lines)-1 and lines[j+1].strip() != '':
                lines.insert(j+1, '')
            i = j
        i += 1
    return '\n'.join(lines) + ('\n' if text.endswith('\n') else '')

def demote_multiple_h1(text: str) -> str:
    count = 0
    out = []
    for line in text.splitlines(True):
        if re.match(r'^#\s+\S', line):
            count += 1
            if count >= 2:
                line = re.sub(r'^#(\s+)', r'##\1', line)
        out.append(line)
    return ''.join(out)

def unescape_selected_chars_except_matlab_text(text: str) -> str:
    out_lines = []
    in_fence = False
    fence_lang = ''
    for raw in text.splitlines(True):
        line = raw
        if line.lstrip().startswith('```'):
            if not in_fence:
                lang = line.strip()[3:].strip().split()[0] if len(line.strip()) > 3 else ''
                fence_lang = lang.lower()
                in_fence = True
            else:
                in_fence = False
                fence_lang = ''
            out_lines.append(line)
            continue
        if in_fence and fence_lang in ('matlab', 'text'):
            out_lines.append(line)
        else:
            out_lines.append(_ESC_RE.sub(r'\1', line))
    return ''.join(out_lines)

# ---- anchor/header analysis (inlined + simplified) ----
ANCHOR_RE = re.compile(r'<a\s+([^>]*)>', re.IGNORECASE)
ID_RE = re.compile(r'(?:id|name)\s*=\s*["\']([^"\']+)["\']', re.IGNORECASE)
HEADER_RE = re.compile(r'^(#{1,6})\s*(.+?)\s*(?:#*\s*)$')

def find_anchors(lines: List[str]) -> List[Dict]:
    anchors = []
    for i, line in enumerate(lines, start=1):
        for m in ANCHOR_RE.finditer(line):
            attrs = m.group(1)
            idm = ID_RE.search(attrs)
            if idm:
                anchors.append({'id': idm.group(1), 'line': i, 'tag': m.group(0), 'attrs': attrs.strip()})
    return anchors

def find_headers(lines: List[str]) -> List[Dict]:
    headers = []
    for i, line in enumerate(lines, start=1):
        m = HEADER_RE.match(line)
        if m:
            headers.append({'level': len(m.group(1)), 'text': m.group(2).strip(), 'line': i})
    return headers

def closest_header_for_anchor(anchor_line: int, headers: List[Dict]) -> Optional[Dict]:
    if not headers:
        return None
    return min(headers, key=lambda h: (abs(h['line'] - anchor_line), h['line']))

def dedupe_adjacent_lines_outside_fences(lines: List[str]) -> List[str]:
    out: List[str] = []
    in_fence = False
    last_line: Optional[str] = None
    fence_re = re.compile(r"^(```|~~~)")
    for line in lines:
        if fence_re.match(line):
            in_fence = not in_fence
            out.append(line)
            last_line = line
            continue
        if in_fence:
            out.append(line)
            last_line = line
            continue
        if last_line is not None and line == last_line:
            continue
        out.append(line)
        last_line = line
    return out

def apply_analysis_and_fixes(text: str) -> str:
    # Normalize certain fence languages: convert matlabTextOutput -> text
    # This keeps the content but standardizes the fence language for rendering.
    text = re.sub(r"(?mi)^```\s*matlabTextOutput\b", '```text', text)

    # operate on list-of-lines for fine-grained line numbers
    lines = text.splitlines(True)
    anchors = find_anchors(lines)
    headers = find_headers(lines)

    id_to_slug = {}
    for h in headers:
        id_to_slug[h['text']] = mkdocs_slug(h['text'])

    compact_to_slug = {}
    for h in headers:
        s = re.sub(r'\\([\\`*_{}\[\]()#+.\-!~])', r"\1", h['text'])
        s = re.sub(r'[`*~]+', '', s)
        s = unicodedata.normalize('NFKD', s)
        s = ''.join(ch for ch in s if not unicodedata.combining(ch))
        s = s.lower()
        compact = re.sub(r'[^a-z0-9_-]+', '', s)
        slug = mkdocs_slug(h['text'])
        if compact and compact != slug:
            compact_to_slug[compact] = slug
            compact_to_slug[compact.lower()] = slug

    new_lines: List[str] = []
    anchor_lines = {a['line']: a for a in anchors}

    for i, line in enumerate(lines, start=1):
        stripped = line.strip()
        if i in anchor_lines:
            if re.fullmatch(r'<a\s+[^>]*>\s*</a>', stripped, flags=re.IGNORECASE):
                continue
        if HEADER_RE.match(line):
            new_line = re.sub(r"\s*\{#[-_a-zA-Z0-9]+\}\s*$", "", line)
            if new_line != line:
                line = new_line

        # replace inline ](#ID) style links
        def _repl_all(mobj: re.Match) -> str:
            key = mobj.group(1)
            slug = None
            # try direct id map (case-sensitive then lower)
            if key in id_to_slug:
                slug = id_to_slug[key]
            elif key.lower() in id_to_slug:
                slug = id_to_slug[key.lower()]
            if not slug and headers:
                nearest = closest_header_for_anchor(i, headers)
                if nearest:
                    slug = mkdocs_slug(nearest['text'])
            if not slug:
                slug = compact_to_slug.get(key) or compact_to_slug.get(key.lower())
            if slug:
                return f'](#{slug})'
            return mobj.group(0)

        line = re.sub(r'\]\(#([^\)]+)\)', _repl_all, line)
        new_lines.append(line)

    # post-process: dedupe adjacent identical lines outside fences, collapse blanks
    post_lines = dedupe_adjacent_lines_outside_fences(new_lines)
    text = ''.join(post_lines)
    text = collapse_blank_runs_to_one(text)

    # run lint passes (three iterations as before)
    for _ in range(3):
        newer = normalize_list_marker_space(text)
        newer = blanks_around_headings(newer)
        newer = blanks_around_fences(newer)
        newer = collapse_blank_runs_to_one(newer)
        newer = demote_multiple_h1(newer)
        newer = final_trailing_space_sweep(newer)
        if newer == text:
            break
        text = unescape_selected_chars_except_matlab_text(newer)

    # ensure final newline
    if not text.endswith('\n'):
        text += '\n'
    return text

def process_all(root: Path, backup: bool = False):
    files = [p for p in root.rglob('*.md') if not p.name.endswith('.post.md')]
    for p in sorted(files):
        print(f'Processing: {p.relative_to(root)}')
        text = p.read_text(encoding='utf-8', errors='replace')
        new_text = apply_analysis_and_fixes(text)
        if backup:
            bak = p.with_suffix(p.suffix + '.bak')
            bak.write_text(text, encoding='utf-8')
        p.write_text(new_text, encoding='utf-8')

def main():
    """Run the conversion/cleanup pipeline on the tutorial_mlx folder only.

    This script intentionally takes no CLI arguments and does not call any
    external scripts. It overwrites the original `.md` files in the same
    directory as this script (no .post.md intermediates). Backups are not
    created by default to match your request.
    """
    script_dir = Path(__file__).parent
    process_all(script_dir, backup=False)

if __name__ == '__main__':
    main()
