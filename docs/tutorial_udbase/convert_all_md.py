#!/usr/bin/env python3
"""Single-file pipeline: run md_convert then analyze_md_links in sequence.

This script bundles the two-stage conversion pipeline used for the
`docs/tutorial_mlx` tutorials. It intentionally runs the original
md-convert (lint/cleanup) stage first and then runs the analyzer stage
that removes HTML anchor-only lines and rewrites inline anchor links to
computed header slugs.

    Important behavior and precautions
    - The md-convert stage overwrites the original `.md` files in-place.
        No `.bak` files are kept by default. This keeps the directory clean
        but means the change is destructive unless you have a backup
        (git, copy, etc.).
- The analyzer stage operates on the overwritten `.md` files and will
    remove anchor-only `<a...></a>` lines and rewrite `](#ID)` links to
    `](#computed-slug)` where possible.

Recommended precaution: run `git status` / commit or make a copy of the
`docs/tutorial_mlx` files before running this script if you want a
rollback option.
"""
from __future__ import annotations

import os
import re
import unicodedata
import sys
from pathlib import Path
from typing import List, Dict, Optional, Tuple

# ------------------------------------------------------------
# md_convert logic (copied & slightly adapted)
# ------------------------------------------------------------

_ESC_RE = re.compile(r'(?<!\\)\\([_\*\[\]\-])')

def final_trailing_space_sweep(text: str) -> str:
    out = []
    for ln in text.splitlines(True):
        if ln.endswith("\n"):
            core, nl = ln[:-1], "\n"
        else:
            core, nl = ln, ""
        m = re.search(r"([ \t]+)$", core)
        if m:
            trailing = m.group(1)
            if "\t" in trailing:
                out.append(re.sub(r"[ \t]+$", "", core) + nl)
            else:
                if len(trailing) == 2:
                    out.append(core + nl)
                else:
                    out.append(re.sub(r"[ \t]+$", "", core) + nl)
        else:
            out.append(core + nl)
    return "".join(out)

def collapse_blank_runs_to_one(text: str) -> str:
    return re.sub(r"\n(?:[ \t]*\n){1,}", "\n\n", text)

def normalize_list_marker_space(text: str) -> str:
    out = []
    in_code = False
    for raw in text.splitlines(True):
        line = raw
        if line.strip().startswith("```"):
            in_code = not in_code
            out.append(line)
            continue
        if not in_code and "<a" not in line.lower():
            m = re.match(r'^(\s*(?:>\s*)*)(\s*)([-*+]|\d+\.)\s+(.*\S.*)$', line.rstrip("\n"))
            if m:
                line = f"{m.group(1)}{m.group(2)}{m.group(3)} {m.group(4)}\n"
        out.append(line)
    return "".join(out)

def blanks_around_headings(text: str) -> str:
    lines = text.splitlines()
    i = 0
    while i < len(lines):
        if re.match(r'^(#{1,6})\s+\S', lines[i]):
            if i > 0 and lines[i-1].strip() != "":
                lines.insert(i, "")
                i += 1
            if i+1 < len(lines) and lines[i+1].strip() != "":
                lines.insert(i+1, "")
        i += 1
    return "\n".join(lines) + ("\n" if text.endswith("\n") else "")

def blanks_around_fences(text: str) -> str:
    lines = text.splitlines()
    i = 0
    while i < len(lines):
        if lines[i].startswith("```"):
            if i > 0 and lines[i-1].strip() != "":
                lines.insert(i, "")
                i += 1
            j = i + 1
            while j < len(lines) and not lines[j].startswith("```"):
                j += 1
            if j < len(lines)-1 and lines[j+1].strip() != "":
                lines.insert(j+1, "")
            i = j
        i += 1
    return "\n".join(lines) + ("\n" if text.endswith("\n") else "")

def demote_multiple_h1(text: str) -> str:
    count = 0
    out = []
    for line in text.splitlines(True):
        if re.match(r"^#\s+\S", line):
            count += 1
            if count >= 2:
                line = re.sub(r"^#(\s+)", r"##\1", line)
        out.append(line)
    return "".join(out)

def lint_fix_pass(text: str) -> str:
    t = text
    t = normalize_list_marker_space(t)
    t = blanks_around_headings(t)
    t = blanks_around_fences(t)
    t = collapse_blank_runs_to_one(t)
    t = demote_multiple_h1(t)
    t = final_trailing_space_sweep(t)
    return t

def ensure_final_newline(text: str) -> str:
    return text if text.endswith("\n") else text + "\n"

def trim_trailing_blank_lines(text: str) -> str:
    lines = text.splitlines(True)
    while lines and lines[-1].strip() == "":
        lines.pop()
    if not lines:
        return "\n"
    last = lines[-1]
    if last.endswith("\n"):
        return "".join(lines)
    return "".join(lines) + "\n"

def unescape_selected_chars_except_matlab_text(text: str) -> str:
    out_lines = []
    in_fence = False
    fence_lang = ""
    for raw in text.splitlines(True):
        line = raw
        if line.lstrip().startswith("```"):
            if not in_fence:
                lang = line.strip()[3:].strip().split()[0] if len(line.strip()) > 3 else ""
                fence_lang = lang.lower()
                in_fence = True
            else:
                in_fence = False
                fence_lang = ""
            out_lines.append(line)
            continue
        if in_fence and fence_lang in ("matlab", "text"):
            out_lines.append(line)
        else:
            out_lines.append(_ESC_RE.sub(r"\1", line))
    return "".join(out_lines)

def process_md_convert_all(root: Path) -> None:
    """Run the md_convert stage: lint/cleanup and overwrite original .md files.

    This pipeline intentionally overwrites the source `.md` files in-place
    (no backup files are produced) to keep the tutorial folder clean. Make
    a backup (git, copy, etc.) if you need a rollback option.
    """
    # process all .md files recursively under the tutorial folder and
    # overwrite them in-place (no `.bak` files are created).
    files = [p for p in root.rglob("*.md")]
    for p in sorted(files):
        text = p.read_text(encoding="utf-8", errors="replace")

        # Run lint passes (anchors/fragments untouched)
        for _ in range(3):
            newer = lint_fix_pass(text)
            if newer == text:
                break
            text = newer

            # Unescape requested characters, excluding matlab/text fences
            text = unescape_selected_chars_except_matlab_text(text)

            # Trim trailing blank-only lines at EOF (fix MD012 at document end)
            text = trim_trailing_blank_lines(text)
        # Overwrite the original .md file with the converted text
        p.write_text(text, encoding="utf-8")

# ------------------------------------------------------------
# analyze_md_links logic (copied & slightly adapted)
# ------------------------------------------------------------

ANCHOR_RE = re.compile(r"<a\s+([^>]*)>", re.IGNORECASE)
ID_RE = re.compile(r"(?:id|name)\s*=\s*[\"']([^\"']+)[\"']", re.IGNORECASE)
HEADER_RE = re.compile(r"^(#{1,6})\s*(.+?)\s*(?:#*\s*)$")

def mkdocs_slug(text: str) -> str:
    s = re.sub(r"\\([\\`*_{}\[\]()#+.\-!~])", r"\1", text)
    s = re.sub(r"[`*~]+", "", s)
    s = unicodedata.normalize('NFKD', s)
    s = ''.join(ch for ch in s if not unicodedata.combining(ch))
    s = s.lower()
    s = s.replace('.', '')
    # Replace any run of non-allowed chars with a hyphen, then collapse
    # consecutive hyphens to a single one. Finally strip leading/trailing
    # hyphens so the slug is clean.
    s = re.sub(r'[^a-z0-9_-]+', '-', s)
    s = re.sub(r'-{2,}', '-', s)
    s = s.strip('-')
    return s

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
    best = min(headers, key=lambda h: (abs(h['line'] - anchor_line), h['line']))
    return best

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

def apply_mappings_inplace(filepath: str, mappings: List[Dict], headers: List[Dict]) -> None:
    with open(filepath, 'r', encoding='utf-8') as f:
        lines = f.readlines()

    anchor_lines = {m['anchor_line']: m for m in mappings}
    header_by_line = {m['header_line']: m for m in mappings if m.get('header_line')}
    id_to_slug = {m['id']: m.get('header_id') for m in mappings if m.get('header_id')}
    id_to_slug_lower = {k.lower(): v for k, v in id_to_slug.items()}

    # compact variants
    compact_to_slug = {}
    for h in headers:
        text = h['text']
        s = re.sub(r"\\([\\`*_{}\[\]()#+.\-!~])", r"\1", text)
        s = re.sub(r'[`*~]+', '', s)
        s = unicodedata.normalize('NFKD', s)
        s = ''.join(ch for ch in s if not unicodedata.combining(ch))
        s = s.lower()
        compact = re.sub(r'[^a-z0-9_-]+', '', s)
        slug = mkdocs_slug(text)
        if compact and compact != slug:
            compact_to_slug[compact] = slug
            compact_to_slug[compact.lower()] = slug

    new_lines: List[str] = []

    for i, line in enumerate(lines, start=1):
        stripped = line.strip()
        if i in anchor_lines:
            if re.fullmatch(r'<a\s+[^>]*>\s*</a>', stripped, flags=re.IGNORECASE):
                continue
        if HEADER_RE.match(line):
            new_line = re.sub(r"\s*\{#[-_a-zA-Z0-9]+\}\s*$", "", line)
            if new_line != line:
                line = new_line

        if id_to_slug or headers:
            def _repl_link(mobj: re.Match) -> str:
                key = mobj.group(1)
                slug = None
                if id_to_slug:
                    slug = id_to_slug.get(key) or id_to_slug_lower.get(key.lower())
                if not slug and headers:
                    nearest = closest_header_for_anchor(i, headers)
                    if nearest:
                        slug = mkdocs_slug(nearest['text'])
                if slug:
                    return f'](#{slug})'
                return mobj.group(0)

            def _repl_all(mobj: re.Match) -> str:
                before = mobj.group(1)
                replaced = _repl_link(mobj)
                if replaced != mobj.group(0):
                    return replaced
                key = before
                slug2 = compact_to_slug.get(key) or compact_to_slug.get(key.lower())
                if slug2:
                    return f'](#{slug2})'
                return mobj.group(0)

            line = re.sub(r'\]\(#([^\)]+)\)', _repl_all, line)

        new_lines.append(line)

    # conservative dedupe and cleanup
    post_lines = dedupe_adjacent_lines_outside_fences(new_lines)

    def _collapse_blank_runs_to_one(text: str) -> str:
        return re.sub(r"\n(?:[ \t]*\n){1,}", "\n\n", text)

    def _trim_trailing_blank_lines(text: str) -> str:
        lines2 = text.splitlines(True)
        while lines2 and lines2[-1].strip() == "":
            lines2.pop()
        if not lines2:
            return "\n"
        last = lines2[-1]
        if last.endswith("\n"):
            return "".join(lines2)
        return "".join(lines2) + "\n"

    text = "".join(post_lines)
    text = _collapse_blank_runs_to_one(text)
    text = _trim_trailing_blank_lines(text)

    with open(filepath, 'w', encoding='utf-8') as f:
        f.write(text)

def analyze_one_postmd(fp: str) -> Tuple[List[Dict], List[Dict], List[Dict]]:
    with open(fp, 'r', encoding='utf-8') as f:
        lines = f.readlines()

    anchors = find_anchors(lines)
    headers = find_headers(lines)

    # The analyzer operates on `.md` files that were updated by the
    # conversion stage. If a file contains no anchors we leave it alone.

    mappings = []
    for a in anchors:
        h = closest_header_for_anchor(a['line'], headers)
        mappings.append({
            'id': a['id'],
            'anchor_line': a['line'],
            'header': h['text'] if h else None,
            'header_level': h['level'] if h else None,
            'header_line': h['line'] if h else None,
            'header_id': mkdocs_slug(h['text']) if h else None,
        })

    return anchors, headers, mappings

def process_analyze_all(root: Path) -> None:
    # operate on .md files (we overwrite originals in the md-convert stage)
    files = [os.path.join(root, name) for name in os.listdir(root) if name.endswith('.md')]
    files = sorted(files)
    if not files:
        print('No .md files found in the analysis directory.')
        return
    for fp in files:
        print(f"Processing: {os.path.basename(fp)}")
        anchors, headers, mappings = analyze_one_postmd(fp)
        apply_mappings_inplace(fp, mappings, headers)

        if anchors:
            print('Anchors found:')
            for a in anchors:
                print(f"  {a['id']} @ line {a['line']}")
            print("")

        if mappings:
            print('Anchor -> computed header_id (nearest header)')
            for m in mappings:
                old = m['id']
                new = m.get('header_id') or ''
                hdr = m.get('header') or ''
                hline = m.get('header_line') or ''
                print(f"  {old} -> {new}  (header: '{hdr}' @ line {hline})")
            print("")

        if headers:
            print('Header -> slug')
            for h in headers:
                slug = mkdocs_slug(h['text'])
                print(f"{h['text']} -> {slug}")
            print("")

# ------------------------------------------------------------
# Main: run md_convert stage then analyze stage
# ------------------------------------------------------------

def main():
    script_dir = Path(__file__).parent
    # Stage 1: md_convert
    print('Running md convert stage...')
    process_md_convert_all(script_dir)

    # Stage 2: analyze and apply mappings to `.md` files (we overwrite
    # originals in the conversion stage)
    print('Running analyzer stage...')
    process_analyze_all(script_dir)

if __name__ == '__main__':
    main()
