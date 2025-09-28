#!/usr/bin/env python3
"""
uDALES Variable Indexer - Compact Report Generator

This script generates a compact summary version of the variable analysis.
"""

import re
import sys
from pathlib import Path

# Import the main indexer class
sys.path.append(str(Path(__file__).parent))
from udales_variable_indexer import UdalesVariableIndexer


def generate_compact_report(indexer: UdalesVariableIndexer) -> str:
    """Generate a compact summary report."""
    report = []
    
    report.append("# uDALES Variable Index - Compact Summary")
    report.append("")
    report.append("## Overview")
    report.append("")
    
    # Get all namelists
    all_namelists = set(indexer.namelists.keys()) | set(indexer.schema_vars.keys())
    
    total_vars = sum(len(vars_set) for vars_set in indexer.namelists.values())
    json_vars = len(indexer.json_reads)
    broadcast_vars = len(indexer.broadcasts)
    schema_vars = sum(len(vars_set) for vars_set in indexer.schema_vars.values())
    
    report.append(f"- **Namelists found**: {len(all_namelists)}")
    report.append(f"- **Total namelist variables**: {total_vars}")
    report.append(f"- **Variables with JSON support**: {json_vars}")
    report.append(f"- **Variables broadcast via MPI**: {broadcast_vars}")
    report.append(f"- **Variables in JSON schema**: {schema_vars}")
    report.append("")
    
    # Summary table per namelist
    report.append("## Namelist Summary")
    report.append("")
    report.append("| Namelist | Total Vars | JSON Read | Broadcast | In Schema |")
    report.append("|----------|------------|-----------|-----------|-----------|")
    
    for namelist in sorted(all_namelists):
        namelist_vars = indexer.namelists.get(namelist, set())
        schema_vars_nl = indexer.schema_vars.get(namelist, set())
        all_vars = namelist_vars | schema_vars_nl
        
        if not all_vars:
            continue
            
        json_count = sum(1 for var in all_vars if var in indexer.json_reads)
        broadcast_count = sum(1 for var in all_vars if var in indexer.broadcasts)
        schema_count = len(schema_vars_nl)
        
        report.append(f"| {namelist} | {len(all_vars)} | {json_count} | {broadcast_count} | {schema_count} |")
    
    report.append("")
    
    # Variables with full support (all 4 aspects)
    report.append("## Variables with Full Support")
    report.append("*(Namelist + JSON + Broadcast + Schema)*")
    report.append("")
    
    full_support_vars = []
    for namelist in sorted(all_namelists):
        namelist_vars = indexer.namelists.get(namelist, set())
        schema_vars_nl = indexer.schema_vars.get(namelist, set())
        all_vars = namelist_vars | schema_vars_nl
        
        for var in sorted(all_vars):
            if (var in namelist_vars and 
                var in indexer.json_reads and 
                var in indexer.broadcasts and 
                var in schema_vars_nl):
                full_support_vars.append(f"{namelist}.{var}")
    
    if full_support_vars:
        for var in full_support_vars:
            report.append(f"- `{var}`")
    else:
        report.append("*No variables found with full support across all aspects.*")
    
    report.append("")
    
    # Variables missing from schema
    report.append("## Variables Missing from JSON Schema")
    report.append("*(In namelists but not in schema)*")
    report.append("")
    
    missing_from_schema = set()
    for namelist, vars_set in indexer.namelists.items():
        schema_vars_nl = indexer.schema_vars.get(namelist, set())
        missing_from_schema.update(vars_set - schema_vars_nl)
    
    if missing_from_schema:
        for var in sorted(missing_from_schema):
            report.append(f"- `{var}`")
    else:
        report.append("*All namelist variables are present in schema.*")
    
    report.append("")
    
    # Variables in schema but not in namelists
    report.append("## Schema Variables Missing from Namelists")
    report.append("*(In schema but not in namelists)*")
    report.append("")
    
    missing_from_namelists = set()
    for namelist, vars_set in indexer.schema_vars.items():
        namelist_vars = indexer.namelists.get(namelist, set())
        missing_from_namelists.update(vars_set - namelist_vars)
    
    if missing_from_namelists:
        for var in sorted(missing_from_namelists):
            report.append(f"- `{var}`")
    else:
        report.append("*All schema variables are present in namelists.*")
    
    return "\n".join(report)


def main():
    """Generate full, compact, and status reports."""
    script_dir = Path(__file__).parent
    src_dir = (script_dir / "../../src").resolve()
    schema_file = (script_dir / "../../docs/schemas/udales_input_schema.json").resolve()
    
    print("Generating uDALES Variable Analysis Reports...")
    print(f"Source directory: {src_dir}")
    print(f"Schema file: {schema_file}")
    print()
    
    # Run analysis
    indexer = UdalesVariableIndexer(src_dir, schema_file)
    indexer.parse_fortran_files()
    indexer.parse_json_schema()
    
    # Generate compact report
    compact_report = generate_compact_report(indexer)
    
    # Write compact report
    compact_output = script_dir / "variable_index_compact.md"
    with open(compact_output, 'w') as f:
        f.write(compact_report)
    
    print(f"Compact report generated: {compact_output}")
    
    # Generate the full report
    full_report = indexer.generate_markdown_report()
    full_output = script_dir / "variable_index_full.md"
    with open(full_output, 'w') as f:
        f.write(full_report)
    
    print(f"Full report generated: {full_output}")
    
    # Note: Status report (with color coding) is generated by generate_status_report.py
    print("For color-coded status report, run: python3 generate_status_report.py")


if __name__ == "__main__":
    main()