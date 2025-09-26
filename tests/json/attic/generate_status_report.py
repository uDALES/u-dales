#!/usr/bin/env python3
"""
uDALES Variable Indexer - Status Report Generator

This script generates a color-coded status report showing variables with:
- Green: Full support (all 4 categories)
- Orange: Warnings (schema mismatches)
- Red: Errors (JSON/broadcast mismatches, missing JSON support)
"""

import re
import sys
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
    
    # Track overall statistics - only for namelist variables
    stats = {
        'total_namelist_vars': len(all_namelist_vars),
        'total_json_reads': len(indexer.json_reads),
        'total_broadcasts': len(indexer.broadcasts),
        'namelist_broadcasts': len(namelist_broadcasts),
        'non_namelist_broadcasts': len(indexer.broadcasts) - len(namelist_broadcasts),
        'total_schema_vars': sum(len(vars_set) for vars_set in indexer.schema_vars.values()),
        'json_broadcast_mismatch': abs(len(indexer.json_reads) - len(namelist_broadcasts))
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
            
            # Determine status - all variables here are in namelists by definition
            if (var_info['has_json_read'] and var_info['has_broadcast'] and var_info['in_schema']):
                # GREEN: Full support
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
    report.append("- 🟠 **ORANGE**: Warnings (Schema mismatches)")  
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
            
            # All variables here are in namelist by definition
            if (var_info['has_json_read'] and var_info['has_broadcast'] and var_info['in_schema']):
                green_vars.append(var)
            elif not var_info['in_schema']:
                orange_vars.append(var)  # In namelist but missing from schema
            else:
                red_vars.append(var)  # In namelist and schema but missing JSON/broadcast
        
        # Output section for this namelist
        report.append(f"### {namelist} Namelist")
        report.append("")
        
        # Green variables - compact comma-separated format
        if green_vars:
            green_list = ", ".join([f"`{var}`" for var in green_vars])
            report.append(f"🟢 **Full Support**: {green_list}")
            report.append("")
        
        # Orange variables - warnings (namelist variables missing from schema)
        if orange_vars:
            report.append("🟠 **Warnings**:")
            for var in orange_vars:
                report.append(f"- `{var}` (missing from schema)")
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
    """Generate the status report."""
    script_dir = Path(__file__).parent
    src_dir = (script_dir / "../../src").resolve()
    schema_file = (script_dir / "../../docs/schemas/udales_input_schema.json").resolve()
    
    print("Generating uDALES Variable Status Report...")
    print(f"Source directory: {src_dir}")
    print(f"Schema file: {schema_file}")
    print()
    
    # Run analysis
    indexer = UdalesVariableIndexer(src_dir, schema_file)
    indexer.parse_fortran_files()
    indexer.parse_json_schema()
    
    # Generate status report
    status_report = generate_status_report(indexer)
    
    # Write report
    output_file = script_dir / "variable_status_report.md"
    with open(output_file, 'w') as f:
        f.write(status_report)
    
    print(f"Status report generated: {output_file}")
    
    # Quick summary
    full_support, warnings, errors, stats = analyze_variable_status(indexer)
    print(f"\nQuick Summary:")
    print(f"  🟢 {len(full_support)} variables with full support")
    print(f"  🟠 {len(warnings)} variables with warnings")
    print(f"  🔴 {len(errors)} variables with errors")


if __name__ == "__main__":
    main()