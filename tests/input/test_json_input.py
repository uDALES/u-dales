"""Standalone JSON validation test.

This script performs one test run and exits with 0 on success or skip,
and non-zero on failure. It is runnable with `python3 tests/input/test_json_input.py`.
"""
import json
import shutil
import subprocess
import sys
import os
import argparse
import tempfile
from pathlib import Path
import re
from types import SimpleNamespace


def _normalize_value(v):
    # Normalization removed: return the value as-is for strict comparisons
    return v


def _compare_generated(generated, input_cfg, schema):
    input_mismatches = []
    default_mismatches = []
    missing_defaults = []

    schema_props = schema.get('properties', {})

    def find_schema_section(name):
        # Case-sensitive lookup only
        return schema_props.get(name)

    def find_prop_entry(sec_props, pname):
        if not isinstance(sec_props, dict):
            return None
        # Case-sensitive lookup only
        return sec_props.get(pname)

    def find_input_value(input_section, pname):
        if not isinstance(input_section, dict):
            return None
        # Case-sensitive lookup only
        return input_section.get(pname)

    for section, params in generated.items():
        if not isinstance(params, dict):
            continue
        input_section = input_cfg.get(section, {})
        sec_schema = find_schema_section(section) or {}
        sec_props = sec_schema.get('properties', {}) if isinstance(sec_schema, dict) else {}

        for pname, pval in params.items():
            gen_val = _normalize_value(pval)

            inp_raw = find_input_value(input_section, pname)
            if inp_raw is not None:
                inp_val = _normalize_value(inp_raw)
                # Numeric comparison with relative tolerance
                try:
                    if isinstance(inp_val, (int, float)) and isinstance(gen_val, (int, float)):
                        denom = max(abs(float(inp_val)), 1e-12)
                        relerr = abs(float(gen_val) - float(inp_val)) / denom
                        if relerr > 1e-3:
                            input_mismatches.append((section, pname, 'input', inp_val, gen_val))
                    else:
                        if inp_val != gen_val:
                            input_mismatches.append((section, pname, 'input', inp_val, gen_val))
                except Exception:
                    if inp_val != gen_val:
                        input_mismatches.append((section, pname, 'input', inp_val, gen_val))
                # when present in input we do not also compare to defaults
            else:
                p_schema = find_prop_entry(sec_props, pname) or {}
                if section.upper() == 'RUN' and pname.lower() in ('nprocx', 'nprocy'):
                    continue
                if 'default' in p_schema:
                    default = p_schema['default']
                    if isinstance(default, str):
                        default_norm = _normalize_value(default)
                    else:
                        default_norm = default

                    # Numeric comparison with relative tolerance
                    try:
                        if isinstance(default_norm, (int, float)) and isinstance(gen_val, (int, float)):
                            denom = max(abs(float(default_norm)), 1e-12)
                            relerr = abs(float(gen_val) - float(default_norm)) / denom
                            if relerr > 1e-3:
                                default_mismatches.append((section, pname, 'default', default_norm, gen_val))
                        else:
                            if default_norm != gen_val:
                                default_mismatches.append((section, pname, 'default', default_norm, gen_val))
                    except Exception:
                        if default_norm != gen_val:
                            default_mismatches.append((section, pname, 'default', default_norm, gen_val))
                else:
                    missing_defaults.append((section, pname, gen_val))

    # Report input-defined mismatches first, then default mismatches
    mismatches = input_mismatches + default_mismatches
    return mismatches, missing_defaults


"""Standalone JSON validation test.

This script performs one test run and exits with 0 on success or skip,
and non-zero on failure. It is runnable with `python3 tests/input/test_json_input.py`.
"""

import json
import shutil
import subprocess
import sys
import os
from pathlib import Path


def main(argv=None):
    parser = argparse.ArgumentParser(description='Run u-DALES JSON test or comparator-only mode')
    parser.add_argument('--comparator-only', action='store_true', help='Skip running u-DALES and run comparator on existing generated JSON')
    parser.add_argument('--generated', help='Path to generated JSON to compare (overrides defaults)')
    parser.add_argument('--mode', choices=['default', 'random'], default='default', help='Which input to run: default (parameters.default) or random (parameters.random)')
    args = parser.parse_args(argv)

    repo_root = Path(__file__).resolve().parents[2]

    # Prefer executables from PATH, fall back to UD_TOPDIR/bin, then repo_root/bin
    udales_path = shutil.which('u-dales')
    ud_nam2json_path = shutil.which('ud_nam2json')

    ud_topdir = os.environ.get('UD_TOPDIR')

    if udales_path:
        udales_exe = Path(udales_path)
    elif ud_topdir:
        udales_exe = Path(ud_topdir) / 'bin' / 'u-dales'
    else:
        udales_exe = repo_root / 'bin' / 'u-dales'

    if ud_nam2json_path:
        ud_nam2json = Path(ud_nam2json_path)
    elif ud_topdir:
        ud_nam2json = Path(ud_topdir) / 'bin' / 'ud_nam2json'
    else:
        ud_nam2json = repo_root / 'bin' / 'ud_nam2json'

    # Prefer schema and input directory from UD_TOPDIR if set
    if ud_topdir:
        schema_path = Path(ud_topdir) / 'docs' / 'schemas' / 'udales_input_schema.json'
        input_dir = Path(ud_topdir) / 'tests' / 'input'
    else:
        schema_path = repo_root / 'docs' / 'schemas' / 'udales_input_schema.json'
        input_dir = repo_root / 'tests' / 'input'

    # Check schema exists always
    if not schema_path.exists():
        print(f'SKIP: schema not found: {schema_path}')
        return 0

    comparator_only = bool(args.comparator_only)

    # If not running comparator-only, ensure runtime prerequisites exist
    if not comparator_only:
        if shutil.which('mpirun') is None:
            print('SKIP: mpirun not available')
            return 0
        if not udales_exe.exists():
            print(f'SKIP: u-DALES executable not found: {udales_exe}')
            return 0
        if not ud_nam2json.exists():
            print(f'SKIP: ud_nam2json not found: {ud_nam2json}')
            return 0

    # Determine which inputs exist
    default_input = input_dir / 'parameters.default'
    random_input = input_dir / 'parameters.random'
    if not default_input.exists():
        print('SKIP: parameters.default not found in tests/input')
        return 0

    # Select inputs to run based on --mode. 'default' runs only parameters.default.
    # 'random' runs parameters.random if it exists, otherwise falls back to default.
    inputs_to_run = []
    if args.mode == 'default':
        inputs_to_run = [(default_input, 'default')]
    else:  # random
        if random_input.exists():
            inputs_to_run = [(random_input, 'random')]
        else:
            # Fallback to default if random not present
            inputs_to_run = [(default_input, 'default')]

    # Use an ephemeral TemporaryDirectory for the workdir so it is removed
    # automatically when this script exits. This avoids leaving stale debug
    # artifacts behind. If you need persistent tmp2, set UD_TOPDIR and the
    # script will continue to use tests/input/tmp2 as before.
    tempdir = tempfile.TemporaryDirectory(prefix='udales-tmp2-')
    workdir = Path(tempdir.name)
    # If comparator-only and args.generated points into workdir, remember to preserve it
    gen_spec = Path(args.generated).resolve() if comparator_only and args.generated else None

    # Load schema once
    with open(schema_path, 'r') as f:
        schema = json.load(f)

    # Iterate over requested inputs
    overall_ret = 0
    for inp_path, tag in inputs_to_run:
        # copy the input into tmp2 with its basename so comparator finds it
        work_input = workdir / inp_path.name
        try:
            shutil.copy2(inp_path, work_input)
        except Exception:
            pass

        # Load the specific input JSON for runtime parameters (nprocx/nprocy)
        with open(work_input, 'r') as f:
            input_config = json.load(f)

        # Determine generated JSON path. In comparator-only mode prefer args.generated or existing files.
        if comparator_only:
            if args.generated:
                generated_json_path = Path(args.generated)
            else:
                # prefer tmp2 then repo tests/input
                candidate = workdir / 'parameters_json'
                if candidate.exists():
                    generated_json_path = candidate
                else:
                    generated_json_path = input_dir / 'parameters_json'
            if not generated_json_path.exists():
                print(f'ERROR: generated JSON not found for comparator-only: {generated_json_path}')
                return 2
        else:
            nprocx = int(input_config.get('RUN', {}).get('nprocx', 1))
            nprocy = int(input_config.get('RUN', {}).get('nprocy', 1))
            total_procs = nprocx * nprocy

            # Run u-DALES with the simple mpirun invocation in the workdir
            try:
                cmd1 = f"mpirun -np {total_procs} u-dales {work_input.name}"
                print('Running:', cmd1)
                proc = subprocess.run(cmd1, cwd=str(workdir), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, timeout=300)
            except Exception as e:
                print('ERROR: running u-DALES failed to start:', e)
                return 2

            if proc.returncode != 0:
                print('u-DALES failed:')
                print(proc.stdout)
                print(proc.stderr)
                return 2

            # Now run ud_nam2json on the expected namelist filename (plain name)
            try:
                cmd2 = f"ud_nam2json namoptions_json parameters_json"
                print('Running:', cmd2, 'in', str(workdir))
                proc2 = subprocess.run(cmd2, cwd=str(workdir), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, timeout=30)
            except Exception as e:
                print('ERROR: running ud_nam2json failed to start:', e)
                return 2

            if proc2.returncode != 0:
                print('ud_nam2json failed:')
                print(proc2.stdout)
                print(proc2.stderr)
                return 2

            generated_json_path = workdir / 'parameters_json'
            if not generated_json_path.exists():
                print(f'ERROR: expected json file {generated_json_path} not found')
                return 2

    # Use the repository comparator tools to produce a canonical diff and markdown table.
    # 1) Run diff_generated.py pointing at the generated JSON in tmp2. Write textual diff to tmp2.
    # 2) Copy textual diff to repo tests/input so diff_to_table.py can read it.
    # 3) Run diff_to_table.py to produce the markdown table, then copy the markdown back to tmp2.

        diff_generated = repo_root / 'tests' / 'input' / 'diff_generated.py'
        # Fallback: some versions keep comparator in tests/input/attic
        if not diff_generated.exists():
            alt = repo_root / 'tests' / 'input' / 'attic' / 'diff_generated.py'
            if alt.exists():
                diff_generated = alt

        # Prepare containers in case comparator fails to set them
        mismatches = []
        missing_defaults = []

        # Run comparator in-process to keep this test self-contained
        # print('Running in-process comparator for', tag)
        proc_cmp = SimpleNamespace(returncode=0, stdout='', stderr='')
        try:
            if not generated_json_path.exists():
                proc_cmp.returncode = 2
                proc_cmp.stderr = f'ERROR: generated file not found: {generated_json_path}'
            else:
                gen = json.loads(generated_json_path.read_text())
                inp = json.loads(work_input.read_text())
                sch = json.loads(schema_path.read_text())
                mismatches, missing_defaults = _compare_generated(gen, inp, sch)

                out_lines = []
                out_lines.append(f'Mismatches: {len(mismatches)}')
                for sec, name, source, expected, found in mismatches:
                    out_lines.append(f'{sec}.{name} ({source}) expected: {expected!r}  found: {found!r}')

                out_lines.append('')
                out_lines.append(f'Missing defaults (no schema default to compare): {len(missing_defaults)}')
                for sec, name, found in missing_defaults:
                    out_lines.append(f'{sec}.{name} found: {found!r}')

                out_text = '\n'.join(out_lines) + '\n'
                # Do not write textual diff file; generate Markdown directly from mismatches/missing_defaults
                proc_cmp.stdout = 'Generated in-memory diff report'
                if mismatches or missing_defaults:
                    proc_cmp.returncode = 1
        except Exception as e:
            proc_cmp.returncode = 2
            proc_cmp.stderr = str(e)
        # Print comparator output for diagnostics
        if proc_cmp.stdout:
            print(proc_cmp.stdout)
        if proc_cmp.stderr:
            print('Comparator stderr:', proc_cmp.stderr)
        if proc_cmp.returncode != 0:
            print('Comparator found differences or failed for', tag, '; return code:', proc_cmp.returncode)

        # Generate Markdown directly from mismatches and missing_defaults (no textual diff file)
        out_md_repo = input_dir / f'test_json_{tag}.md'
        out_md_work = workdir / f'test_json_{tag}.md'
        try:
            rows = []
            for sec, name, source, expected, found in mismatches:
                key = f"{sec.upper()}.{name.lower()}"
                if source == 'input':
                    rows.append((key, found, expected, ''))
                else:
                    rows.append((key, found, '', expected))
            for sec, name, found in missing_defaults:
                key = f"{sec.upper()}.{name.lower()}"
                rows.append((key, found, '', ''))

            with out_md_repo.open('w') as f:
                    f.write('| variable | input_value | namoptions_value | schema_value |\n')
                    f.write('|---|---:|---:|---:|\n')
                    for var, found, input_val, schema_val in rows:
                        f.write(f'| {var} | {input_val} | {found} | {schema_val} |\n')
            try:
                shutil.copy2(out_md_repo, out_md_work)
            except Exception:
                with out_md_work.open('w') as f2:
                    f2.write(out_md_repo.read_text())
        except Exception as e:
            print('Warning: could not write Markdown:', e)

        # Accumulate return code (non-zero indicates a failure/diff)
        if proc_cmp.returncode != 0:
            overall_ret = proc_cmp.returncode

    # Exit with the most severe comparator return code
    if overall_ret != 0:
        return overall_ret

    print('PASS: generated JSON matches input values and schema defaults for all tests')
    return 0

if __name__ == '__main__':
    sys.exit(main())