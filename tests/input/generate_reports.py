#!/usr/bin/env python3
"""
uDALES Variable Analysis - Unified Report Generator

This script generates v    report.append("- 🟢 **GREEN**: Full support (Namelist + JSON + Broadcast + Schema)")
    report.append("- 🟠 **ORANGE**: Warnings (Schema mismatches, 🟡 Duplicate operations)")  
    report.append("- 🔴 **RED**: Errors (JSON/Broadcast issues)")ble analysis reports:
1. Full detailed report with all variables per namelist
2. Color-coded status report with warnings and errors

Automatically generates both reports without requiring command-line arguments.
"""

import re
import sys
import time
from pathlib import Path

# Import the main indexer class
sys.path.append(str(Path(__file__).parent))
from udales_variable_indexer import UdalesVariableIndexer



def analyze_variable_status(indexer: UdalesVariableIndexer):
    """Analyze variables and categorize them by status - only namelist variables."""
    
    # Status categories - only for variables that exist in namelists
    full_support = []      # Green: All 4 aspects covered
    warnings = []          # Orange: Schema mismatches
    errors = []            # Red: JSON/broadcast issues
    
    # Get all namelist variables
    all_namelist_vars = set()
    for vars_set in indexer.namelists.values():
        all_namelist_vars.update(vars_set)
    
    # Count only broadcasts of namelist variables (exclude internal/computed variables)
    namelist_broadcasts = indexer.broadcasts & all_namelist_vars
    
    # Detect duplicates
    duplicate_json_reads = {var: count for var, count in indexer.json_read_counts.items() if count > 1}
    duplicate_broadcasts = {var: count for var, count in indexer.broadcast_counts.items() if count > 1}
    duplicate_namelists = {var: count for var, count in indexer.namelist_counts.items() if count > 1}
    
    # Track overall statistics - only for namelist variables
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
    
    # Analyze each namelist variable
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
            
            # Check for duplicates across all three lists
            json_count = duplicate_json_reads.get(var, 1)
            broadcast_count = duplicate_broadcasts.get(var, 1)
            namelist_count = duplicate_namelists.get(var, 1)
            has_duplicates = json_count > 1 or broadcast_count > 1 or namelist_count > 1

            # If duplicates exist in any list, treat as error with details of which lists
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
            
            # Determine status - all variables here are in namelists by definition
            if (var_info['has_json_read'] and var_info['has_broadcast'] and var_info['in_schema']):
                if has_duplicates:
                    # ORANGE: Has full support but with duplicates
                    var_info['warning_type'] = 'duplicate_operations'
                    duplicate_info = []
                    if json_count > 1:
                        duplicate_info.append(f"JSON: {json_count}x")
                    if broadcast_count > 1:
                        duplicate_info.append(f"Broadcast: {broadcast_count}x")
                    var_info['warning_detail'] = f"Duplicate operations: {', '.join(duplicate_info)}"
                    warnings.append(var_info)
                else:
                    # GREEN: Full support without duplicates
                    full_support.append(var_info)
            
            elif not var_info['in_schema']:
                # ORANGE: In namelist but missing from schema
                var_info['warning_type'] = 'schema_mismatch'
                var_info['warning_detail'] = 'In namelist but missing from schema'
                warnings.append(var_info)
            
            else:
                # RED: In namelist and schema but missing JSON/broadcast
                var_info['error_type'] = 'json_broadcast_mismatch'
                if not var_info['has_json_read']:
                    var_info['error_detail'] = 'Defined in namelist but no JSON reading support'
                else:
                    var_info['error_detail'] = 'Has JSON read but no MPI broadcast'
                errors.append(var_info)
    
    return full_support, warnings, errors, stats


def generate_status_report(indexer: UdalesVariableIndexer) -> str:
    """Generate a color-coded status report."""
    
    full_support, warnings, errors, stats = analyze_variable_status(indexer)
    
    report = []
    
    report.append("# uDALES Variable Status Report")
    report.append("")
    report.append("This report shows the status of variables with color coding:")
    report.append("")
    report.append("- 🟢 **GREEN**: Full support (Namelist + JSON + Broadcast + Schema)")
    report.append("- 🟠 **ORANGE**: Warnings (Schema mismatches, duplicates)")
    report.append("- 🔴 **RED**: Errors (JSON/Broadcast issues)")
    report.append("")
    
    # Executive Summary
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
    

    
    # Critical Issues Section
    if stats['json_broadcast_mismatch'] > 10:
        report.append("## 🚨 Critical Configuration Issue")
        report.append("")
        report.append(f"**JSON Coverage Gap**: {stats['total_json_reads']} JSON reads vs {stats['namelist_broadcasts']} namelist broadcasts")
        report.append(f"**Missing JSON Support**: {stats['json_broadcast_mismatch']} namelist variables")
        report.append("")
        report.append("This means many namelist variables are broadcast to all MPI processes but cannot be configured via JSON.")
        report.append(f"Additionally, {stats['non_namelist_broadcasts']} internal/computed variables are also broadcast (this is normal).")
        report.append("")
    

    
    # Group variables by namelist - only show variables that exist in namelists
    report.append("## Variables by Namelist")
    report.append("*Only showing variables that are defined in Fortran namelists*")
    report.append("")
    
    for namelist in sorted(indexer.namelists.keys()):
        namelist_vars = indexer.namelists[namelist]
        schema_vars_nl = indexer.schema_vars.get(namelist, set())
        
        if not namelist_vars:
            continue
            
        # Categorize variables for this namelist - only process namelist variables
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
            
            # Check for duplicates (handled globally in analyze_variable_status)
            json_count = stats['duplicate_json_reads'].get(var, 1)
            broadcast_count = stats['duplicate_broadcasts'].get(var, 1)
            
            # All variables here are in namelist by definition
            if (var_info['has_json_read'] and var_info['has_broadcast'] and var_info['in_schema']):
                green_vars.append(var)
            elif not var_info['in_schema']:
                orange_vars.append((var, 'schema_missing'))  # In namelist but missing from schema
            else:
                red_vars.append(var)  # In namelist and schema but missing JSON/broadcast
        
        # Output section for this namelist
        file_location = indexer.namelist_files.get(namelist, "unknown file")
        report.append(f"### {namelist} Namelist *(defined in {file_location})*")
        report.append("")
        
        # Green variables - compact comma-separated format
        if green_vars:
            green_list = ", ".join([f"`{var}`" for var in green_vars])
            report.append(f"🟢 **Full Support**: {green_list}")
            report.append("")
        
        # Orange variables - warnings (namelist variables missing from schema or with duplicates)
        if orange_vars:
            report.append("🟠 **Warnings**:")
            for item in orange_vars:
                if isinstance(item, tuple) and len(item) == 4:
                    # Duplicate variable: (var, 'duplicate', json_count, broadcast_count)
                    var, warning_type, json_count, broadcast_count = item
                    duplicate_info = []
                    if json_count > 1:
                        duplicate_info.append(f"JSON: {json_count}x")
                    if broadcast_count > 1:
                        duplicate_info.append(f"Broadcast: {broadcast_count}x")
                    duplicate_text = ", ".join(duplicate_info)
                    report.append(f"- `{var}` (duplicate operations: {duplicate_text})")
                elif isinstance(item, tuple) and len(item) == 2:
                    # Schema missing variable: (var, 'schema_missing')
                    var, warning_type = item
                    report.append(f"- `{var}` (missing from schema)")
                else:
                    # Legacy single variable case
                    report.append(f"- `{item}` (missing from schema)")
            report.append("")
        
        # Red variables - errors (namelist variables with JSON/broadcast issues)
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
    

    
    # Recommendations
    report.append("## 📋 Recommendations")
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
    """Generate variable analysis reports."""
    
    script_dir = Path(__file__).parent
    src_dir = (script_dir / "../../src").resolve()
    schema_file = (script_dir / "../../docs/schemas/udales_input_schema.json").resolve()
    
    print("🔍 uDALES Variable Analysis - Report Generator")
    print("=" * 50)
    print(f"Source directory: {src_dir}")
    print(f"Schema file: {schema_file}")
    print()
    
    # Run analysis once - parse all Fortran files and JSON schema
    # This data will be reused in memory for all subsequent report generation
    print("🔍 Parsing source files and schema...")
    
    start_time = time.time()
    indexer = UdalesVariableIndexer(str(src_dir), str(schema_file))

    # Only process the canonical readparameters.f90 file for all extractions
    # (namelists, JSON reads and MPI broadcasts). Do not scan other Fortran files.
    readparams = (indexer.src_dir / 'readparameters.f90')
    if readparams.exists():
        indexer.parse_fortran_files(str(readparams))
    else:
        print(f"  Error: {readparams} not found — aborting analysis.")
        sys.exit(1)

    # Parse JSON schema as before
    indexer.parse_json_schema()
    parse_time = time.time() - start_time
    
    print(f"✅ Analysis complete ({parse_time:.2f}s) - data loaded in memory for report generation")
    
    # Generate reports
    reports_generated = []
    
    # Generate Full Report
    print("📊 Generating Full Report...")
    full_report = indexer.generate_markdown_report()
    full_output = script_dir / "variable_index_full.md"
    with open(full_output, 'w') as f:
        f.write(full_report)
    reports_generated.append(("Full Report", full_output))
    print(f"✅ Full report generated: {full_output}")
    
    # Generate Status Report
    print("🎯 Generating Status Report...")
    status_report = generate_status_report(indexer)
    status_output = script_dir / "variable_status_report.md"
    with open(status_output, 'w') as f:
        f.write(status_report)
    reports_generated.append(("Status Report", status_output))
    
    # Show quick summary for status report
    full_support, warnings, errors, stats = analyze_variable_status(indexer)
    print(f"✅ Status report generated: {status_output}")
    print(f"   🟢 {len(full_support)} variables with full support")
    print(f"   🟠 {len(warnings)} variables with warnings")
    print(f"   🔴 {len(errors)} variables with errors")

    # Final test result: pass if there are no errors
    if len(errors) == 0:
        print("\n✅ TEST PASSED: No errors found")
    else:
        print(f"\n❌ TEST FAILED: {len(errors)} error(s), {len(warnings)} warning(s)")

    # Summary
    print("\n" + "=" * 50)
    print("📋 Generated Reports Summary:")
    print("=" * 50)
    
    for title, filepath in reports_generated:
        size_kb = filepath.stat().st_size / 1024
        print(f"✅ {title:15} | {filepath.name:25} | {size_kb:6.1f} KB")
    
    print("\n🎯 Recommended: Start with 'variable_status_report.md' for color-coded overview")
    print("📖 Full documentation: VARIABLE_INDEXER_README.md")
    # Determine exit code for unit tests: pass if no errors, fail otherwise
    exit_code = 0 if len(errors) == 0 else 1
    # Exit with appropriate code so unit tests can assert pass/fail
    sys.exit(exit_code)

if __name__ == "__main__":
    main()