#!/usr/bin/env python3
"""Compare generated JSON against input values and schema defaults.

Writes a diff file `parameters_diff.txt` in the `tests/input` directory by default.
Exit status: 0 if no differences found, 1 otherwise.

Usage: python3 tests/input/diff_generated.py [--generated PATH] [--input PATH] [--schema PATH] [--out PATH]
If arguments are omitted sensible defaults inside the repository are used:
 - generated: tests/input/parameters_json
 - input: first parameters.* in tests/input (prefer parameters.default)
 - schema: docs/schemas/udales_input_schema.json
 - out: tests/input/parameters_diff.txt
"""

import json
import sys
from pathlib import Path
import argparse


def normalize_value(v):
    """Normalize a value coming from the AWK-generated JSON.

    - strip trailing commas from strings
    - interpret Fortran-style floats and integers
    - convert 'T'/'F' and similar to booleans
    - leave other strings as-is
    """
    if isinstance(v, bool):
        return v
    if isinstance(v, (int, float)):
        return v
    if not isinstance(v, str):
        return v

    s = v.strip()
    # remove trailing comma if present
    if s.endswith(','):
        s = s[:-1].strip()

    if s == '':
        return ''

    low = s.lower()
    if low in ('.true.', 'true', 't'):
        return True
    if low in ('.false.', 'false', 'f'):
        return False

    # Try integer
    try:
        if all(ch.isdigit() or ch in '+-' for ch in s):
            return int(s)
    except Exception:
        pass

    # Try float (handles Fortran E-format)
    try:
        return float(s)
    except Exception:
        pass

    return s


def compare(generated, input_cfg, schema):
    mismatches = []
    missing_defaults = []

    schema_props = schema.get('properties', {})

    # Helper: find a section in schema_props case-insensitively
    def find_schema_section(name):
        # exact match
        if name in schema_props:
            return schema_props[name]
        # case-insensitive
        lname = name.lower()
        for k, v in schema_props.items():
            if k.lower() == lname:
                return v
        return None

    # Helper: find a property entry in a section schema case-insensitively
    def find_prop_entry(sec_props, pname):
        if not isinstance(sec_props, dict):
            return None
        if pname in sec_props:
            return sec_props[pname]
        lp = pname.lower()
        for k, v in sec_props.items():
            if k.lower() == lp:
                return v
        return None

    # Helper: find value in input section case-insensitively
    def find_input_value(input_section, pname):
        if not isinstance(input_section, dict):
            return None
        if pname in input_section:
            return input_section[pname]
        lp = pname.lower()
        for k, v in input_section.items():
            if k.lower() == lp:
                return v
        return None

    for section, params in generated.items():
        if not isinstance(params, dict):
            continue
        input_section = input_cfg.get(section, {})
        sec_schema = find_schema_section(section) or {}
        sec_props = sec_schema.get('properties', {}) if isinstance(sec_schema, dict) else {}

        for pname, pval in params.items():
            gen_val = normalize_value(pval)

            inp_raw = find_input_value(input_section, pname)
            if inp_raw is not None:
                inp_val = normalize_value(inp_raw)
                if inp_val != gen_val:
                    mismatches.append((section, pname, 'input', inp_val, gen_val))
            else:
                p_schema = find_prop_entry(sec_props, pname) or {}
                # Some parameters are required and intentionally have no default; skip them
                if section.upper() == 'RUN' and pname.lower() in ('nprocx', 'nprocy'):
                    continue
                if 'default' in p_schema:
                    default = p_schema['default']
                    # Normalize schema default if it's a string representing a number/boolean
                    if isinstance(default, str):
                        default_norm = normalize_value(default)
                    else:
                        default_norm = default

                    # For floats compare with tolerance
                    if isinstance(default_norm, float) and isinstance(gen_val, (int, float)):
                        if abs(float(gen_val) - float(default_norm)) > 1e-9:
                            mismatches.append((section, pname, 'default', default_norm, gen_val))
                    else:
                        if default_norm != gen_val:
                            mismatches.append((section, pname, 'default', default_norm, gen_val))
                else:
                    missing_defaults.append((section, pname, gen_val))

    return mismatches, missing_defaults


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('--generated', help='Generated JSON file', default=None)
    parser.add_argument('--input', help='Input JSON file (parameters.*)', default=None)
    parser.add_argument('--schema', help='Schema JSON file', default=None)
    parser.add_argument('--out', help='Output diff file', default=None)
    args = parser.parse_args(argv)

    repo = Path(__file__).resolve().parents[2]
    tests_input = repo / 'tests' / 'input'

    generated_path = Path(args.generated) if args.generated else tests_input / 'parameters_json'

    if args.input:
        input_path = Path(args.input)
    else:
        # prefer parameters.default then first parameters.* (excluding parameters_json)
        candidates = sorted([p for p in tests_input.glob('parameters.*') if p.name != 'parameters_json'])
        input_path = None
        for c in candidates:
            if c.name == 'parameters.default':
                input_path = c
                break
        if input_path is None and candidates:
            input_path = candidates[0]

    schema_path = Path(args.schema) if args.schema else repo / 'docs' / 'schemas' / 'udales_input_schema.json'
    out_path = Path(args.out) if args.out else tests_input / 'parameters_diff.txt'

    if not generated_path.exists():
        print(f'ERROR: generated file not found: {generated_path}', file=sys.stderr)
        return 2
    if input_path is None or not input_path.exists():
        print(f'ERROR: input file not found: {input_path}', file=sys.stderr)
        return 2
    if not schema_path.exists():
        print(f'ERROR: schema file not found: {schema_path}', file=sys.stderr)
        return 2

    with generated_path.open('r') as f:
        gen = json.load(f)
    with input_path.open('r') as f:
        inp = json.load(f)
    with schema_path.open('r') as f:
        schema = json.load(f)

    mismatches, missing_defaults = compare(gen, inp, schema)

    # Write report
    out_lines = []
    out_lines.append(f'Mismatches: {len(mismatches)}')
    for sec, name, source, expected, found in mismatches:
        out_lines.append(f'{sec}.{name} ({source}) expected: {expected!r}  found: {found!r}')

    out_lines.append('')
    out_lines.append(f'Missing defaults (no schema default to compare): {len(missing_defaults)}')
    for sec, name, found in missing_defaults:
        out_lines.append(f'{sec}.{name} found: {found!r}')

    out_text = '\n'.join(out_lines) + '\n'
    out_path.write_text(out_text)
    print(f'Wrote diff report to: {out_path}')

    if mismatches or missing_defaults:
        print('Differences found; exiting with failure (1)')
        return 1
    print('No differences found; exiting success (0)')
    return 0


if __name__ == '__main__':
    sys.exit(main())
