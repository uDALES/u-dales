#!/usr/bin/env python3
"""
Static input-contract test.

Checks that variables defined in Fortran namelists are still:
1) present in the JSON schema used for tooling/editor support
2) broadcast where appropriate
"""

import argparse
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).parent))
from udales_variable_indexer import UdalesVariableIndexer


def analyze_variable_status(indexer: UdalesVariableIndexer):
    full_support = []
    warnings = []
    errors = []

    all_namelist_vars = set()
    for vars_set in indexer.namelists.values():
        all_namelist_vars.update(vars_set)

    namelist_broadcasts = indexer.broadcasts & all_namelist_vars
    duplicate_broadcasts = {var: count for var, count in indexer.broadcast_counts.items() if count > 1}
    duplicate_namelists = {var: count for var, count in indexer.namelist_counts.items() if count > 1}

    stats = {
        "total_namelist_vars": len(all_namelist_vars),
        "total_broadcasts": len(indexer.broadcasts),
        "namelist_broadcasts": len(namelist_broadcasts),
        "non_namelist_broadcasts": len(indexer.broadcasts) - len(namelist_broadcasts),
        "total_schema_vars": sum(len(vars_set) for vars_set in indexer.schema_vars.values()),
    }

    for namelist in sorted(indexer.namelists.keys()):
        schema_vars_nl = indexer.schema_vars.get(namelist, set())
        for var in sorted(indexer.namelists[namelist]):
            var_info = {
                "namelist": namelist,
                "variable": var,
                "has_broadcast": var in indexer.broadcasts,
                "in_schema": var in schema_vars_nl,
            }

            broadcast_count = duplicate_broadcasts.get(var, 1)
            namelist_count = duplicate_namelists.get(var, 1)
            if broadcast_count > 1 or namelist_count > 1:
                dup_sources = []
                if namelist_count > 1:
                    dup_sources.append(f"namelist: {namelist_count}x")
                if broadcast_count > 1:
                    dup_sources.append(f"broadcast: {broadcast_count}x")
                var_info["error_type"] = "duplicate_operations"
                var_info["error_detail"] = ", ".join(dup_sources)
                errors.append(var_info)
                continue

            if var_info["has_broadcast"] and var_info["in_schema"]:
                full_support.append(var_info)
            elif not var_info["in_schema"]:
                var_info["warning_type"] = "schema_mismatch"
                var_info["warning_detail"] = "In namelist but missing from schema"
                warnings.append(var_info)
            else:
                var_info["error_type"] = "broadcast_missing"
                var_info["error_detail"] = "Defined in namelist but no MPI broadcast"
                errors.append(var_info)

    return full_support, warnings, errors, stats


def write_reports(indexer: UdalesVariableIndexer, full_support, warnings, errors, stats):
    full_report = []
    full_report.append("# uDALES Variable Index Report")
    full_report.append("")
    full_report.append("| Namelist | Variable | Broadcast | Schema |")
    full_report.append("|----------|----------|-----------|--------|")
    for namelist in sorted(indexer.namelists.keys()):
        schema_vars_nl = indexer.schema_vars.get(namelist, set())
        for var in sorted(indexer.namelists[namelist]):
            broadcast_check = "✓" if var in indexer.broadcasts else "✗"
            schema_check = "✓" if var in schema_vars_nl else "✗"
            full_report.append(f"| `{namelist}` | `{var}` | {broadcast_check} | {schema_check} |")
    (Path(__file__).parent / "test_sourcecode_full.md").write_text("\n".join(full_report))

    status_report = []
    status_report.append("# uDALES Variable Status Report")
    status_report.append("")
    status_report.append(f"- Total namelist variables: {stats['total_namelist_vars']}")
    status_report.append(f"- Namelist broadcasts: {stats['namelist_broadcasts']}")
    status_report.append(f"- Non-namelist broadcasts: {stats['non_namelist_broadcasts']}")
    status_report.append(f"- Schema variables: {stats['total_schema_vars']}")
    status_report.append(f"- Full support: {len(full_support)}")
    status_report.append(f"- Warnings: {len(warnings)}")
    status_report.append(f"- Errors: {len(errors)}")
    status_report.append("")
    if warnings:
        status_report.append("## Warnings")
        status_report.append("")
        for warning in warnings:
            status_report.append(f"- `{warning['namelist']}.{warning['variable']}`: {warning['warning_detail']}")
        status_report.append("")
    if errors:
        status_report.append("## Errors")
        status_report.append("")
        for error in errors:
            status_report.append(f"- `{error['namelist']}.{error['variable']}`: {error['error_type']} - {error['error_detail']}")
        status_report.append("")
    (Path(__file__).parent / "test_sourcecode_status.md").write_text("\n".join(status_report))


def main(argv=None):
    parser = argparse.ArgumentParser(description="Run static source/schema checks for uDALES input variables")
    parser.add_argument("--src-dir", default="../../../src")
    parser.add_argument("--schema-file", default="../../../docs/schemas/udales_input_schema.json")
    parser.add_argument("--fortran-file", default=None)
    args = parser.parse_args(argv)

    script_dir = Path(__file__).parent
    src_dir = (script_dir / args.src_dir).resolve()
    schema_file = (script_dir / args.schema_file).resolve()

    indexer = UdalesVariableIndexer(str(src_dir), str(schema_file))
    if args.fortran_file:
        candidate = Path(args.fortran_file)
        if not candidate.is_absolute():
            candidate = src_dir / args.fortran_file
        readparams = candidate
    else:
        readparams = src_dir / "readparameters.f90"

    indexer.parse_fortran_files(str(readparams))
    indexer.parse_json_schema()

    embedded_map_text = """
NAMSUBGRID               SUBGRID
NAMCHECKSIM.tcheck       OUTPUT.tcheck
SUBGRID.sg_cs            SUBGRID.cs
"""

    simple_map = {}
    qualified_map = {}
    map_lines = embedded_map_text.strip().splitlines()
    map_file = (script_dir / "../../../docs/schemas/nml_mapping.txt").resolve()
    if map_file.exists():
        map_lines.extend(map_file.read_text().splitlines())

    for raw in map_lines:
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 2:
            continue
        src = parts[0]
        dst = parts[-1]
        if "." in src:
            ssec, skey = src.split(".", 1)
            if "." in dst:
                dsec, dkey = dst.split(".", 1)
            else:
                dsec, dkey = dst, skey
            qualified_map[(ssec.upper(), skey.lower())] = (dsec.upper(), dkey.lower())
        else:
            simple_map[src.lower()] = dst

    for src_low, dst in list(simple_map.items()):
        src_sec = src_low.upper()
        if src_sec in indexer.namelists and "." not in dst:
            dst_sec = dst.upper()
            src_vars = indexer.namelists.pop(src_sec, set())
            dst_vars = indexer.namelists.get(dst_sec, set())
            for var in src_vars:
                if var not in dst_vars:
                    dst_vars.add(var)
            indexer.namelists[dst_sec] = dst_vars

            src_fn = indexer.namelist_files.pop(src_sec, None)
            if src_fn and dst_sec not in indexer.namelist_files:
                indexer.namelist_files[dst_sec] = src_fn

            src_schema = indexer.schema_vars.pop(src_sec, set())
            dst_schema = indexer.schema_vars.get(dst_sec, set())
            for var in src_schema:
                if var not in dst_schema:
                    dst_schema.add(var)
            indexer.schema_vars[dst_sec] = dst_schema

    for (ssec, skey), (dsec, dkey) in list(qualified_map.items()):
        if ssec in indexer.namelists and skey in indexer.namelists[ssec]:
            indexer.namelists[ssec].discard(skey)
            indexer.namelists.setdefault(dsec, set()).add(dkey)
            src_fn = indexer.namelist_files.get(ssec, None)
            if src_fn and dsec not in indexer.namelist_files:
                indexer.namelist_files[dsec] = src_fn
        if skey in indexer.schema_vars.get(ssec, set()):
            indexer.schema_vars[ssec].discard(skey)
            indexer.schema_vars.setdefault(dsec, set()).add(dkey)

    full_support, warnings, errors, stats = analyze_variable_status(indexer)
    write_reports(indexer, full_support, warnings, errors, stats)

    if errors:
        print("test_sourcecode: FAILED")
        for error in errors:
            print(f"  {error['namelist']}.{error['variable']}: {error['error_detail']}")
        return 1

    print("test_sourcecode: PASSED")
    if warnings:
        print(f"  warnings: {len(warnings)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
