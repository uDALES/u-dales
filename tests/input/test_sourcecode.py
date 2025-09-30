#!/usr/bin/env python3
"""
uDALES Variable Analysis - Test Source Code Runner

This script is a renamed variant of generate_reports.py used by CI/tests.
It produces two reports:
1. test_sourcecode_full.md - full variable index
2. test_sourcecode_status.md - color-coded status report

Exits with code 0 when no errors (suitable for unit tests), otherwise 1.
"""

import re
import sys
import time
import argparse
from pathlib import Path

# Import the main indexer class
sys.path.append(str(Path(__file__).parent))
from udales_variable_indexer import UdalesVariableIndexer


def analyze_variable_status(indexer: UdalesVariableIndexer):
    """Analyze variables and categorize them by status - only namelist variables."""
    full_support = []
    warnings = []
    errors = []

    all_namelist_vars = set()
    for vars_set in indexer.namelists.values():
        all_namelist_vars.update(vars_set)

    namelist_broadcasts = indexer.broadcasts & all_namelist_vars

    duplicate_json_reads = {var: count for var, count in indexer.json_read_counts.items() if count > 1}
    duplicate_broadcasts = {var: count for var, count in indexer.broadcast_counts.items() if count > 1}
    duplicate_namelists = {var: count for var, count in indexer.namelist_counts.items() if count > 1}

    stats = {
        'total_namelist_vars': len(all_namelist_vars),
        'total_json_reads': len(indexer.json_reads),
        'total_broadcasts': len(indexer.broadcasts),
        'namelist_broadcasts': len(namelist_broadcasts),
        'non_namelist_broadcasts': len(indexer.broadcasts) - len(namelist_broadcasts),
        'total_schema_vars': sum(len(vars_set) for vars_set in indexer.schema_vars.values()),
        'json_broadcast_mismatch': abs(len(indexer.json_reads) - len(namelist_broadcasts)),
        'duplicate_json_reads': duplicate_json_reads,
        'duplicate_broadcasts': duplicate_broadcasts
    }

    for namelist in sorted(indexer.namelists.keys()):
        namelist_vars = indexer.namelists[namelist]
        schema_vars_nl = indexer.schema_vars.get(namelist, set())

        for var in sorted(namelist_vars):
            var_info = {
                'namelist': namelist,
                'variable': var,
                'has_json_read': var in indexer.json_reads,
                'has_broadcast': var in indexer.broadcasts,
                'in_schema': var in schema_vars_nl
            }

            json_count = duplicate_json_reads.get(var, 1)
            broadcast_count = duplicate_broadcasts.get(var, 1)
            namelist_count = duplicate_namelists.get(var, 1)
            has_duplicates = json_count > 1 or broadcast_count > 1 or namelist_count > 1
            if has_duplicates:
                # If the only duplicate source is JSON reads, classify as a warning
                if json_count > 1 and broadcast_count == 1 and namelist_count == 1:
                    var_info['warning_type'] = 'duplicate_json_reads'
                    var_info['warning_detail'] = f"Duplicate JSON reads: {json_count}x"
                    warnings.append(var_info)
                    continue

                # Otherwise treat as an error (multiple duplicate sources or broadcast/namelist duplicates)
                dup_sources = []
                if namelist_count > 1:
                    dup_sources.append(f"namelist: {namelist_count}x")
                if json_count > 1:
                    dup_sources.append(f"json: {json_count}x")
                if broadcast_count > 1:
                    dup_sources.append(f"broadcast: {broadcast_count}x")
                var_info['error_type'] = 'duplicate_operations'
                var_info['error_detail'] = f"Duplicate operations in: {', '.join(dup_sources)}"
                errors.append(var_info)
                continue

            if (var_info['has_json_read'] and var_info['has_broadcast'] and var_info['in_schema']):
                full_support.append(var_info)
            elif not var_info['in_schema']:
                var_info['warning_type'] = 'schema_mismatch'
                var_info['warning_detail'] = 'In namelist but missing from schema'
                warnings.append(var_info)
            else:
                var_info['error_type'] = 'json_broadcast_mismatch'
                if not var_info['has_json_read']:
                    var_info['error_detail'] = 'Defined in namelist but no JSON reading support'
                else:
                    var_info['error_detail'] = 'Has JSON read but no MPI broadcast'
                errors.append(var_info)

    return full_support, warnings, errors, stats


def generate_status_report(indexer: UdalesVariableIndexer) -> str:
    full_support, warnings, errors, stats = analyze_variable_status(indexer)

    report = []
    report.append("# uDALES Variable Status Report")
    report.append("")
    report.append("This report shows the status of variables with color coding:")
    report.append("")
    report.append("- 🟢 GREEN: Full support (Namelist + JSON + Broadcast + Schema)")
    report.append("- 🟠 ORANGE: Warnings (Schema mismatches, duplicates)")
    report.append("- 🔴 RED: Errors (JSON/Broadcast issues)")
    report.append("")

    report.append("## Executive Summary")
    report.append("")
    report.append(f"- **Total Namelist Variables**: {stats['total_namelist_vars']}")
    report.append(f"- **JSON Support**: {stats['total_json_reads']} variables")
    report.append(f"- **Namelist Broadcasts**: {stats['namelist_broadcasts']} variables")
    report.append(f"- **Non-namelist Broadcasts**: {stats['non_namelist_broadcasts']} variables (internal/computed)")
    report.append(f"- **Schema Coverage**: {stats['total_schema_vars']} variables")
    report.append("")
    report.append(f"**Status Distribution:**")
    report.append(f"- 🟢 **{len(full_support)} variables** with full support")
    report.append(f"- 🟠 **{len(warnings)} variables** with warnings")
    report.append(f"- 🔴 **{len(errors)} variables** with errors")
    report.append("")

    report.append("## Variables by Namelist")
    report.append("*Only showing variables that are defined in Fortran namelists*")
    report.append("")

    for namelist in sorted(indexer.namelists.keys()):
        namelist_vars = indexer.namelists[namelist]
        schema_vars_nl = indexer.schema_vars.get(namelist, set())

        if not namelist_vars:
            continue

        green_vars = []
        orange_vars = []
        red_vars = []

        # Build maps of variables that have errors or warnings (duplicate operations etc.)
        error_map = {e['variable']: e for e in errors}
        warning_map = {w['variable']: w for w in warnings}

        for var in sorted(namelist_vars):
            var_info = {
                'variable': var,
                'has_json_read': var in indexer.json_reads,
                'has_broadcast': var in indexer.broadcasts,
                'in_schema': var in schema_vars_nl
            }


            # If this variable has a recorded error (duplicates etc.), show as red
            if var in error_map:
                err = error_map[var]
                red_vars.append((var, err.get('error_type', 'error'), err.get('error_detail', '')))
                continue

            # If this variable has a recorded warning, show as orange with detail
            if var in warning_map:
                warn = warning_map[var]
                orange_vars.append((var, warn.get('warning_type', 'warning'), warn.get('warning_detail', '')))
                continue

            if (var_info['has_json_read'] and var_info['has_broadcast'] and var_info['in_schema']):
                green_vars.append(var)
            elif not var_info['in_schema']:
                orange_vars.append((var, 'schema_missing'))
            else:
                red_vars.append((var, 'json_broadcast_mismatch', ''))

        report.append(f"### {namelist}")
        report.append("")

        if green_vars:
            green_list = ", ".join([f"`{var}`" for var in green_vars])
            report.append(f"🟢 **Full Support**: {green_list}")
            report.append("")

        if orange_vars:
            report.append("🟠 **Warnings**:")
            for item in orange_vars:
                # item can be (var, 'schema_missing') or (var, 'duplicate_json_reads', detail)
                if isinstance(item, tuple):
                    if len(item) == 2:
                        var, warning_type = item
                        if warning_type == 'schema_missing':
                            report.append(f"- `{var}` (missing from schema)")
                        else:
                            report.append(f"- `{var}` ({warning_type})")
                    elif len(item) == 3:
                        var, warning_type, detail = item
                        if warning_type == 'duplicate_json_reads':
                            report.append(f"- `{var}` (duplicate JSON reads: {detail})")
                        else:
                            report.append(f"- `{var}` ({warning_type}: {detail})")
                else:
                    report.append(f"- `{item}` (missing from schema)")
            report.append("")

        if red_vars:
            report.append("🔴 **Errors**:")
            for item in red_vars:
                if isinstance(item, tuple) and len(item) == 3:
                    var, etype, detail = item
                else:
                    var = item
                    etype = 'error'
                    detail = ''

                json_status = "✓" if var in indexer.json_reads else "✗"
                broadcast_status = "✓" if var in indexer.broadcasts else "✗"
                schema_status = "✓" if var in schema_vars_nl else "✗"

                if etype == 'duplicate_operations':
                    report.append(f"- `{var}` (JSON: {json_status}, Broadcast: {broadcast_status}, Schema: {schema_status}) - Duplicate operations: {detail}")
                elif etype == 'json_broadcast_mismatch':
                    if var not in indexer.json_reads:
                        report.append(f"- `{var}` (JSON: {json_status}, Broadcast: {broadcast_status}, Schema: {schema_status}) - No JSON support")
                    elif var not in indexer.broadcasts:
                        report.append(f"- `{var}` (JSON: {json_status}, Broadcast: {broadcast_status}, Schema: {schema_status}) - No MPI broadcast")
                    else:
                        report.append(f"- `{var}` (JSON: {json_status}, Broadcast: {broadcast_status}, Schema: {schema_status})")
                else:
                    report.append(f"- `{var}` (JSON: {json_status}, Broadcast: {broadcast_status}, Schema: {schema_status}) - {etype} {detail}")
            report.append("")

    # Recommendations section removed by request

    # Detailed Errors (if any)
    if errors:
        report.append("")
        report.append("## Errors (detailed)")
        report.append("")
        for e in errors:
            nl = e.get('namelist', '<unknown>')
            var = e.get('variable', '<unknown>')
            etype = e.get('error_type', 'error')
            detail = e.get('error_detail', '') or e.get('error_detail', '')
            report.append(f"- `{nl}.{var}`: {etype} - {detail}")

    # Report variables present in schema but not in namelists, and vice versa
    schema_only = set()
    for sec, vars_set in indexer.schema_vars.items():
        nl_vars = indexer.namelists.get(sec, set())
        for v in vars_set:
            if v not in nl_vars:
                schema_only.add(v)

    namelist_only = set()
    for sec, vars_set in indexer.namelists.items():
        schema_vars = indexer.schema_vars.get(sec, set())
        for v in vars_set:
            if v not in schema_vars:
                namelist_only.add(v)

    if schema_only:
        report.append("")
        report.append("## Variables in schema but not in namelists")
        report.append("")
        for v in sorted(schema_only):
            report.append(f"- `{v}`")

    if namelist_only:
        report.append("")
        report.append("## Variables in namelists but not in schema")
        report.append("")
        for v in sorted(namelist_only):
            report.append(f"- `{v}`")

    return "\n".join(report)


def main():
    parser = argparse.ArgumentParser(description="Generate uDALES variable reports")
    parser.add_argument(
        "--fortran-file",
        dest="fortran_file",
        help="Path to a Fortran file to parse (absolute or relative to src directory). If omitted, defaults to readparameters.f90",
        default=None,
    )
    args = parser.parse_args()

    script_dir = Path(__file__).parent
    src_dir = (script_dir / "../../src").resolve()
    schema_file = (script_dir / "../../docs/schemas/udales_input_schema.json").resolve()

    # Minimal output mode: don't print verbose progress to stdout
    start_time = time.time()
    indexer = UdalesVariableIndexer(str(src_dir), str(schema_file))

    # Determine which Fortran file to parse. Accept either an absolute path
    # or a path relative to the src directory (the latter is convenient for CI).
    if args.fortran_file:
        candidate = Path(args.fortran_file)
        if not candidate.is_absolute():
            candidate = src_dir / args.fortran_file
        readparams = candidate
    else:
        readparams = (indexer.src_dir / 'readparameters.f90')

    if readparams.exists():
        indexer.parse_fortran_files(str(readparams))
    else:
        print(f"  Error: {readparams} not found — aborting analysis.")
        sys.exit(1)

    indexer.parse_json_schema()

    # --- Apply ud_nam2json-style mappings directly to the indexer structures ---
    # Embedded default mapping (kept in ud_nam2json). Each entry is (src, dst)
    EMBEDDED_MAP_TEXT = '''
NAMSUBGRID               SUBGRID
NAMCHECKSIM.tcheck       OUTPUT.tcheck
SUBGRID.sg_cs            SUBGRID.cs
'''

    def load_mappings(script_dir_path):
        simple_map = {}
        qualified_map = {}
        lines = EMBEDDED_MAP_TEXT.strip().splitlines()
        try:
            mapfile = (script_dir_path / "../../docs/schemas/nam2jsonmap.txt").resolve()
            if mapfile.exists():
                with open(mapfile, 'r') as mf:
                    lines += [ln for ln in mf.read().splitlines()]
        except Exception:
            pass

        for raw in lines:
            line = raw.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            src = parts[0]
            dst = parts[-1]
            if '.' in src:
                ssec, skey = src.split('.', 1)
                if '.' in dst:
                    dsec, dkey = dst.split('.', 1)
                else:
                    dsec, dkey = dst, skey
                qualified_map[(ssec.upper(), skey.lower())] = (dsec.upper(), dkey.lower())
            else:
                simple_map[src.lower()] = dst
        return simple_map, qualified_map

    script_dir = Path(__file__).parent
    simple_map, qualified_map = load_mappings(script_dir)

    # Apply simple (section) mappings: move/rename sections when source exists
    # simple_map keys are lowercased source tokens
    for src_low, dst in list(simple_map.items()):
        src_sec = src_low.upper()
        # only treat as section move if source exists as a top-level namelist
        if src_sec in indexer.namelists and '.' not in dst:
            dst_sec = dst.upper()
            # move/merge namelist variables (destination takes precedence)
            src_vars = indexer.namelists.pop(src_sec, set())
            dst_vars = indexer.namelists.get(dst_sec, set())
            for v in src_vars:
                if v not in dst_vars:
                    dst_vars.add(v)
            indexer.namelists[dst_sec] = dst_vars

            # move file association if dest doesn't have one
            src_fn = indexer.namelist_files.pop(src_sec, None)
            if src_fn and dst_sec not in indexer.namelist_files:
                indexer.namelist_files[dst_sec] = src_fn

            # move schema vars similarly
            src_schema = indexer.schema_vars.pop(src_sec, set())
            dst_schema = indexer.schema_vars.get(dst_sec, set())
            for v in src_schema:
                if v not in dst_schema:
                    dst_schema.add(v)
            indexer.schema_vars[dst_sec] = dst_schema

    # Apply qualified mappings: SECTION.key -> DSTSECTION.newkey
    for (ssec, skey), (dsec, dkey) in list(qualified_map.items()):
        # move variable from ssec to dsec if present
        if ssec in indexer.namelists and skey in indexer.namelists[ssec]:
            indexer.namelists[ssec].discard(skey)
            indexer.namelists.setdefault(dsec, set()).add(dkey)
            # move file association if dest doesn't have one
            src_fn = indexer.namelist_files.get(ssec, None)
            if src_fn and dsec not in indexer.namelist_files:
                indexer.namelist_files[dsec] = src_fn
        # move schema var if present
        if skey in indexer.schema_vars.get(ssec, set()):
            indexer.schema_vars[ssec].discard(skey)
            indexer.schema_vars.setdefault(dsec, set()).add(dkey)

    # --- end mapping application ---

    parse_time = time.time() - start_time

    # Generate Full Report
    full_report = indexer.generate_markdown_report()
    full_output = script_dir / "test_sourcecode_full.md"
    with open(full_output, 'w') as f:
        f.write(full_report)
    # Generate Status Report
    status_report = generate_status_report(indexer)
    status_output = script_dir / "test_sourcecode_status.md"
    with open(status_output, 'w') as f:
        f.write(status_report)
    full_support, warnings, errors, stats = analyze_variable_status(indexer)

    # Minimal summary output
    print(f"Reports written: {full_output.name}, {status_output.name}")
    print(f"Summary: {len(full_support)} full, {len(warnings)} warnings, {len(errors)} errors")

    if len(errors) == 0:
        print("TEST PASSED")
    else:
        print(f"TEST FAILED: {len(errors)} error(s), {len(warnings)} warning(s)")

    exit_code = 0 if len(errors) == 0 else 1
    sys.exit(exit_code)


if __name__ == "__main__":
    main()
