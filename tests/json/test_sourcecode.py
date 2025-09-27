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

        for var in sorted(namelist_vars):
            var_info = {
                'variable': var,
                'has_json_read': var in indexer.json_reads,
                'has_broadcast': var in indexer.broadcasts,
                'in_schema': var in schema_vars_nl
            }

            if (var_info['has_json_read'] and var_info['has_broadcast'] and var_info['in_schema']):
                green_vars.append(var)
            elif not var_info['in_schema']:
                orange_vars.append((var, 'schema_missing'))
            else:
                red_vars.append(var)

        report.append(f"### {namelist}")
        report.append("")

        if green_vars:
            green_list = ", ".join([f"`{var}`" for var in green_vars])
            report.append(f"🟢 **Full Support**: {green_list}")
            report.append("")

        if orange_vars:
            report.append("🟠 **Warnings**:")
            for item in orange_vars:
                if isinstance(item, tuple) and len(item) == 2:
                    var, warning_type = item
                    report.append(f"- `{var}` (missing from schema)")
                else:
                    report.append(f"- `{item}` (missing from schema)")
            report.append("")

        if red_vars:
            report.append("🔴 **Errors**:")
            for var in red_vars:
                json_status = "✓" if var in indexer.json_reads else "✗"
                broadcast_status = "✓" if var in indexer.broadcasts else "✗"
                schema_status = "✓" if var in schema_vars_nl else "✗"

                if var not in indexer.json_reads:
                    report.append(f"- `{var}` (JSON: {json_status}, Broadcast: {broadcast_status}, Schema: {schema_status}) - No JSON support")
                elif var not in indexer.broadcasts:
                    report.append(f"- `{var}` (JSON: {json_status}, Broadcast: {broadcast_status}, Schema: {schema_status}) - No MPI broadcast")
                else:
                    report.append(f"- `{var}` (JSON: {json_status}, Broadcast: {broadcast_status}, Schema: {schema_status})")
            report.append("")

    report.append("## Recommendations")
    report.append("")
    if errors:
        report.append("### High Priority")
        report.append("1. **Add JSON reading support** for variables in namelists")
        report.append("2. **Add MPI broadcast calls** for JSON-read variables")
        report.append("")

    if warnings:
        report.append("### Medium Priority")
        report.append("3. **Update JSON schema** to include missing namelist variables")
        report.append("4. **Review schema variables** that don't correspond to namelists")
        report.append("")

    report.append("### General")
    report.append(f"5. **Improve JSON coverage**: Currently {len(indexer.json_reads)}/{stats['total_namelist_vars']} namelist variables support JSON")
    report.append(f"6. **Focus on namelist variables**: {stats['json_broadcast_mismatch']} namelist variables lack JSON support despite being broadcast")

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

    print("🔍 uDALES Variable Analysis - Test Source Code")
    print("=" * 50)
    print(f"Source directory: {src_dir}")
    print(f"Schema file: {schema_file}")
    print()

    print("🔍 Parsing source files and schema...")
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
    parse_time = time.time() - start_time

    print(f"✅ Analysis complete ({parse_time:.2f}s) - data loaded in memory for report generation")

    # Generate Full Report
    print("📊 Generating Full Report...")
    full_report = indexer.generate_markdown_report()
    full_output = script_dir / "test_sourcecode_full.md"
    with open(full_output, 'w') as f:
        f.write(full_report)
    print(f"✅ Full report generated: {full_output}")

    # Generate Status Report
    print("🎯 Generating Status Report...")
    status_report = generate_status_report(indexer)
    status_output = script_dir / "test_sourcecode_status.md"
    with open(status_output, 'w') as f:
        f.write(status_report)
    print(f"✅ Status report generated: {status_output}")

    full_support, warnings, errors, stats = analyze_variable_status(indexer)
    print(f"   🟢 {len(full_support)} variables with full support")
    print(f"   🟠 {len(warnings)} variables with warnings")
    print(f"   🔴 {len(errors)} variables with errors")

    if len(errors) == 0:
        print("\n✅ TEST PASSED: No errors found")
    else:
        print(f"\n❌ TEST FAILED: {len(errors)} error(s), {len(warnings)} warning(s)")

    print("\n" + "=" * 50)
    print("📋 Generated Reports Summary:")
    print("=" * 50)
    for name in [full_output.name, status_output.name]:
        size_kb = (script_dir / name).stat().st_size / 1024
        print(f"✅ {name:25} | {size_kb:6.1f} KB")

    exit_code = 0 if len(errors) == 0 else 1
    sys.exit(exit_code)


if __name__ == "__main__":
    main()
