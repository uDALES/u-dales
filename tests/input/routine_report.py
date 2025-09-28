#!/usr/bin/env python3
"""
Simple routine_report: walk src/, filter lines containing any of the keywords
(subroutine, json%get, mpi_bcast) and write filtered files into tests/json/tmp.
Case-insensitive. Keeps original filename and writes one filtered file per source.
"""
import os
import re

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../'))
SRC_DIR = os.path.join(ROOT, 'src')
# write outputs into tmp2 as requested
OUT_DIR = os.path.join(os.path.dirname(__file__), 'tmp2')

# Patterns:
# - keep files only if they contain at least one of these (inclusion test)
INCLUDE_KEYWORDS = [r'json%get', r'mpi_bcast']
INCLUDE_PAT = re.compile('|'.join(INCLUDE_KEYWORDS), re.IGNORECASE)

# - but include lines that match either subroutine or the above patterns
LINE_KEYWORDS = [r'subroutine', r'json%get', r'mpi_bcast']
LINE_PAT = re.compile('|'.join(LINE_KEYWORDS), re.IGNORECASE)
os.makedirs(OUT_DIR, exist_ok=True)

count_files = 0
count_lines_in = 0
count_lines_out = 0

for root, dirs, files in os.walk(SRC_DIR):
    for fname in files:
        if not fname.endswith('.f90'):
            continue
        inpath = os.path.join(root, fname)
        rel = os.path.relpath(inpath, SRC_DIR)

        # read the file (we will analyse in-memory; don't write per-file summaries)
        try:
            with open(inpath, 'r', errors='ignore') as inf:
                lines = inf.readlines()
        except Exception:
            # skip unreadable files
            continue

        count_lines_in += len(lines)
        # Determine if file should be included at all (has json%get or MPI_BCAST)
        if not INCLUDE_PAT.search(''.join(lines)):
            continue

        count_files += 1

print(f'Filtered {count_files} files; scanned {count_lines_in} lines, wrote {count_lines_out} lines to {OUT_DIR}')

# --- Additional analysis: parse the generated .sum files to extract variables
# read by read_XXX_json and broadcast_XXX_parameters and produce a report
REPORT_PATH = os.path.join(OUT_DIR, 'namelist_variable_report.md')

# regexes to find subroutine headers and variable calls
RE_SUBROUTINE = re.compile(r"^\s*subroutine\s+(\w+)", re.IGNORECASE)
RE_JSONGET = re.compile(r"json%get\s*\(\s*'([A-Za-z0-9_]+)\.([A-Za-z0-9_]+)'", re.IGNORECASE)
RE_BCAST = re.compile(r"MPI_BCAST\s*\(\s*([^,\)\s]+)", re.IGNORECASE)

from collections import defaultdict

# mapping: namelist -> set(vars) for reads and broadcasts
reads = defaultdict(set)
bcasts = defaultdict(set)

def section_from_routine(rname):
    # Expect patterns like read_XXX_json or broadcast_XXX_parameters
    m = re.match(r"read_([a-z0-9_]+)_json", rname, re.IGNORECASE)
    if m:
        return m.group(1).lower()
    m = re.match(r"broadcast_([a-z0-9_]+)_parameters", rname, re.IGNORECASE)
    if m:
        return m.group(1).lower()
    return None

# Scan source files directly, analyze in-memory (no temporary .sum files)
for root, dirs, files in os.walk(SRC_DIR):
    for fname in files:
        if not fname.endswith('.f90'):
            continue
        path = os.path.join(root, fname)
        try:
            with open(path, 'r', errors='ignore') as f:
                lines = f.readlines()
        except Exception:
            continue

        # skip files without json%get or MPI_BCAST
        if not INCLUDE_PAT.search(''.join(lines)):
            continue

        cur_routine = None
        for ln in lines:
            sm = RE_SUBROUTINE.match(ln)
            if sm:
                cur_routine = sm.group(1)
                continue
            if not cur_routine:
                continue

            sec = section_from_routine(cur_routine)

            # JSON keys: if routine doesn't map to a section, fall back to the
            # SECTION token present in the json%get(...) call
            for jm in RE_JSONGET.finditer(ln):
                jsec = jm.group(1).lower()
                var = jm.group(2).lower()
                # prefer section from routine name; otherwise use the SECTION
                # token in the json%get call. If neither is present, skip the
                # entry — we don't want to collect variables under a None/unknown
                # key.
                use_sec = sec or jsec
                if not use_sec:
                    continue
                reads[use_sec].add(var)

            # MPI_BCAST entries - only record broadcasts when the current
            # routine maps to a named broadcast_XXX_parameters subroutine
            # (i.e. section_from_routine returned a section). This avoids
            # collecting unrelated MPI_BCAST uses into an 'unknown' bucket.
            if sec:
                for bm in RE_BCAST.finditer(ln):
                    var = bm.group(1)
                    var = re.sub(r"\(.*\)", '', var)
                    var = var.strip().lower()
                    # normalize common suffixes so broadcast variable names
                    # match canonical variable names where possible
                    var = re.sub(r"(_json)$", '', var)
                    bcasts[sec].add(var)

# write report
with open(REPORT_PATH, 'w') as rep:
    rep.write('# Namelist variable report\n\n')
    # remove any falsy keys (None, '') that might have crept in
    all_secs = sorted(s for s in set(list(reads.keys()) + list(bcasts.keys())) if s)
    # only show sections that actually have detected variables (omit empty ones
    # such as an empty 'unknown' bucket)
    visible_secs = [s for s in all_secs if reads.get(s) or bcasts.get(s)]
    if visible_secs:
        rep.write('Sections found: ' + ', '.join(visible_secs) + '\n\n')
    else:
        rep.write('Sections found: (none)\n\n')

    for sec in visible_secs:
        rlist = sorted(reads.get(sec, []))
        blist = sorted(bcasts.get(sec, []))
        # single long line: header plus JSON and BCAST lists
        json_str = ', '.join(rlist) if rlist else '(none)'
        bcast_str = ', '.join(blist) if blist else '(none)'
        # write Read from JSON line, blank line, then Broadcast line for readability
        rep.write(f'## {sec} \nRead from JSON: {json_str}\n\nBroadcast via MPI_BCAST: {bcast_str}\n')

print(f'Wrote namelist variable report to {REPORT_PATH}')

# --- Compare against canonical namelist variable lists from variable_index_full.md
IDX_PATH = os.path.join(os.path.dirname(__file__), 'variable_index_full.md')
canonical = defaultdict(set)
RE_SECTION_HDR = re.compile(r"^##\s+([A-Za-z0-9_\-]+)", re.IGNORECASE)
RE_TABLE_VAR = re.compile(r"^\|\s*`?([A-Za-z0-9_]+)`?\s*\|")
if os.path.exists(IDX_PATH):
    cur_sec = None
    with open(IDX_PATH, 'r') as idxf:
        for ln in idxf:
            m = RE_SECTION_HDR.match(ln)
            if m:
                cur_sec = m.group(1).split()[0].lower()
                continue
            if cur_sec:
                mv = RE_TABLE_VAR.match(ln)
                if mv:
                    name = mv.group(1).lower()
                    if name == 'variable':
                        continue
                    canonical[cur_sec].add(name)

    # Write a detailed per-namelist table: Variable | InIndex | JSON_Read | MPI_BCAST
    with open(REPORT_PATH, 'a') as rep:
        rep.write('\n# Comparison with canonical namelist index\n\n')
        secs = sorted(set(list(canonical.keys()) + list(reads.keys()) + list(bcasts.keys())))
        # Only include sections in the comparison if we detected any reads or
        # broadcasts for them. This hides canonical-only sections that have no
        # detected activity (e.g. an 'unknown' section that only appears in the
        # canonical index but has no JSON reads or MPI_BCASTs in the code).
        comp_secs = [s for s in secs if reads.get(s) or bcasts.get(s)]
        if not comp_secs:
            rep.write('(No detected variables to compare against canonical index.)\n')
        for sec in comp_secs:
            rep.write(f'## {sec}\n\n')
            can = canonical.get(sec, set())
            detected_reads = reads.get(sec, set())
            detected_bcasts = bcasts.get(sec, set())

            # Build a combined variable list (canonical first, then extras)
            combined = list(sorted(can))
            extras = sorted((detected_reads | detected_bcasts) - can)
            combined.extend(extras)

            rep.write('| Variable | In index | JSON read | MPI_BCAST |\n')
            rep.write('|---|---:|:---:|:---:|\n')
            for v in combined:
                in_index = '✓' if v in can else ''
                jread = '✓' if v in detected_reads else ''
                bcast = '✓' if v in detected_bcasts else ''
                rep.write(f'| `{v}` | {in_index} | {jread} | {bcast} |\n')

            # Summary lines
            rep.write('\n')
            rep.write(f'- Canonical variables in index: {len(can)}\n')
            rep.write(f'- JSON reads detected: {len(detected_reads)} (including {len(detected_reads & can)} in index)\n')
            rep.write(f'- MPI_BCAST detected: {len(detected_bcasts)} (including {len(detected_bcasts & can)} in index)\n')
            if extras:
                rep.write('\n**Detected variables not in index:**\n\n')
                for v in extras:
                    rep.write(f'- {v}\n')
            rep.write('\n')
    print(f'Appended detailed comparison to {REPORT_PATH}')
else:
    print(f'Canonical index {IDX_PATH} not found; skipping comparison')
