#!/usr/bin/env python3
"""Convert parameters_diff.txt into a simple table (CSV + Markdown).

Reads `tests/input/parameters_diff.txt` and writes:
 - `tests/input/parameters_diff_table.csv` (variable, namoptions_value, schema_value)
 - `tests/input/parameters_diff_table.md` (markdown table)

Run from repo root:
  python3 tests/input/diff_to_table.py
"""

from pathlib import Path
import csv
import re


repo = Path(__file__).resolve().parents[2]
diff_txt = repo / 'tests' / 'input' / 'parameters_diff.txt'
out_csv = repo / 'tests' / 'input' / 'parameters_diff_table.csv'
out_md = repo / 'tests' / 'input' / 'parameters_diff_table.md'

if not diff_txt.exists():
    print('Diff file not found:', diff_txt)
    raise SystemExit(2)

lines = diff_txt.read_text().splitlines()

rows = []

# Patterns
mismatch_re = re.compile(r"^([A-Za-z0-9_]+\.[A-Za-z0-9_]+) \((?:default|input)\) expected: (.+?)  found: (.+)$")
missing_re = re.compile(r"^([A-Za-z0-9_]+\.[A-Za-z0-9_]+) found: (.+)$")

for ln in lines:
    ln = ln.strip()
    if not ln:
        continue
    m = mismatch_re.match(ln)
    if m:
        var = m.group(1)
        expected = m.group(2).strip()
        found = m.group(3).strip()
        # normalize: section uppercase, variable lowercase
        section, varname = var.split('.')
        key = f"{section.upper()}.{varname.lower()}"
        rows.append((key, found, expected))
        continue
    m2 = missing_re.match(ln)
    if m2 and ln.startswith('Missing defaults') is False and not ln.startswith('Mismatches'):
        var = m2.group(1)
        found = m2.group(2).strip()
        section, varname = var.split('.')
        key = f"{section.upper()}.{varname.lower()}"
        rows.append((key, found, ''))

# Write CSV
with out_csv.open('w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['variable', 'namoptions_value', 'schema_value'])
    for r in rows:
        w.writerow(r)

# Write Markdown
with out_md.open('w') as f:
    f.write('| variable | namoptions_value | schema_value |\n')
    f.write('|---|---:|---:|\n')
    for var, found, expected in rows:
        f.write(f'| {var} | {found} | {expected} |\n')

print('Wrote', out_csv, 'and', out_md)
