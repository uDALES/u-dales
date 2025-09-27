#!/usr/bin/env python3
"""
uDALES Variable Indexer

This script analyzes the uDALES Fortran source code and generates a comprehensive
index of variables that are:
1) Defined in namelists (and in which namelist)
2) Being read from JSON configuration files
3) Being broadcast via MPI
4) Present in the JSON schema

The output is a markdown table per namelist showing the status of each variable.

Usage: python udales_variable_indexer.py [--src-dir SRC_DIR] [--schema-file SCHEMA_FILE] [--output OUTPUT_FILE]
"""

import argparse
import json
import os
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Set, Tuple, Optional

class UdalesVariableIndexer:
    """Main class for analyzing uDALES variables."""
    
    def __init__(self, src_dir: str, schema_file: str):
        self.src_dir = Path(src_dir)
        self.schema_file = Path(schema_file)
        
        # Data structures to store findings
        self.namelists: Dict[str, Set[str]] = defaultdict(set)  # namelist -> {variables}
        self.namelist_files: Dict[str, str] = {}  # namelist -> filename where defined
        self.json_reads: Set[str] = set()  # variables read from JSON
        self.broadcasts: Set[str] = set()  # variables that are broadcast
        self.schema_vars: Dict[str, Set[str]] = defaultdict(set)  # namelist -> {variables in schema}
        
        # Duplicate tracking
        from collections import Counter
        self.json_read_counts: Counter = Counter()  # count occurrences of JSON reads
        self.broadcast_counts: Counter = Counter()  # count occurrences of broadcasts
        
    def parse_fortran_files(self):
        """Parse all Fortran files to extract namelist definitions, JSON reads, and broadcasts."""
        print("Parsing Fortran source files...")
        
        for f90_file in self.src_dir.glob("*.f90"):
            # Skip backup files
            if f90_file.name.endswith('.backup') or '.backup.' in f90_file.name:
                continue
            print(f"  Processing {f90_file.name}")
            try:
                with open(f90_file, 'r', encoding='utf-8', errors='ignore') as f:
                    content = f.read()
                
                self._extract_namelists(content, f90_file.name)
                self._extract_json_reads(content)
                self._extract_broadcasts(content)
                
            except Exception as e:
                print(f"    Warning: Could not process {f90_file}: {e}")
    
    def _extract_namelists(self, content: str, filename: str):
        """Extract namelist definitions from Fortran source."""
        # Find namelist declarations more precisely
        # Pattern: namelist/NAMELISTNAME/ &
        #          var1, var2, var3, &
        #          var4, var5
        
        # Common Fortran keywords and intrinsics to filter out
        fortran_keywords = {
            'abs', 'acos', 'asin', 'atan', 'atan2', 'ceiling', 'cos', 'cosh', 'exp', 
            'floor', 'log', 'log10', 'max', 'min', 'mod', 'nint', 'real', 'sin', 
            'sinh', 'sqrt', 'tan', 'tanh', 'allocated', 'allocatable', 'allocate', 
            'deallocate', 'associated', 'present', 'size', 'shape', 'lbound', 'ubound',
            'integer', 'real', 'logical', 'character', 'complex', 'double', 'precision',
            'parameter', 'intent', 'optional', 'pointer', 'target', 'save', 'public', 
            'private', 'protected', 'volatile', 'asynchronous', 'bind', 'abstract',
            'dimension', 'kind', 'len', 'if', 'then', 'else', 'elseif', 'endif', 
            'select', 'case', 'default', 'endselect', 'do', 'enddo', 'while', 'exit',
            'cycle', 'where', 'elsewhere', 'endwhere', 'forall', 'endforall', 'goto',
            'continue', 'stop', 'pause', 'return', 'call', 'function', 'subroutine',
            'program', 'module', 'contains', 'use', 'only', 'implicit', 'none',
            'include', 'data', 'common', 'equivalence', 'external', 'intrinsic',
            'sequence', 'type', 'endtype', 'interface', 'endinterface', 'procedure',
            'true', 'false', 'and', 'or', 'not', 'eqv', 'neqv', 'lt', 'le', 'eq',
            'ne', 'gt', 'ge', 'read', 'write', 'print', 'open', 'close', 'inquire',
            'backspace', 'endfile', 'rewind', 'flush', 'wait', 'format', 'iostat',
            'err', 'end', 'file', 'unit', 'rec', 'fmt', 'nml', 'advance', 'size',
            'eor', 'pos', 'stream', 'access', 'form', 'recl', 'blank', 'position',
            'action', 'delim', 'pad', 'sign', 'round', 'decimal', 'encoding', 'pending',
            'exists', 'opened', 'number', 'named', 'name', 'sequential', 'direct',
            'formatted', 'unformatted', 'nextrec', 'iolength', 'readwrite', 'readonly',
            'j', 'k', 'i', 'n', 'm', 'l', 'ii', 'jj', 'kk', 'nn', 'mm', 'll',
            'equal', 'equals', 'error', 'exist', 'comm3d', 'mpierr', 'mpi_real',
            'my_real', 'info', 'does', 'be', 'at', 'in', 'to', 'by', 'for', 'with',
            'on', 'of', 'from', 'into', 'within', 'through', 'during', 'before',
            'after', 'above', 'below', 'up', 'down', 'out', 'off', 'over', 'under'
        }
        
        # Look for namelist declarations with better pattern matching
        lines = content.split('\n')
        current_namelist = None
        collecting_vars = False
        
        debug = bool(os.environ.get('UD_INDEXER_DEBUG'))
        for idx, line in enumerate(lines, start=1):
            raw_line = line.rstrip()
            # Remove comments for parsing
            proc_line = re.sub(r'!.*$', '', line).strip()
            if not proc_line:
                continue
                
            # Check for namelist start
            namelist_match = re.match(r'\s*namelist\s*/\s*(\w+)\s*/\s*&?\s*(.*)', line, re.IGNORECASE)
            if namelist_match:
                current_namelist = namelist_match.group(1).upper()
                # Track which file this namelist is defined in
                self.namelist_files[current_namelist] = filename
                var_part = namelist_match.group(2)
                if debug:
                    print(f"DEBUG [_extract_namelists] {filename}:{idx}: NAMELIST START -> {raw_line}")
                if var_part:
                    # show the part on the same line if present
                    if debug:
                        print(f"DEBUG [_extract_namelists] {filename}:{idx}: NAMELIST VARS ON START -> {var_part}")
                    self._extract_vars_from_line(var_part, current_namelist, fortran_keywords, debug=debug)

                # Check if this namelist line has continuation
                if raw_line.endswith('&'):
                    collecting_vars = True
                else:
                    # Namelist is complete on this line, don't collect from subsequent lines
                    collecting_vars = False
                    current_namelist = None
                continue
                
            # If we're collecting variables and find another namelist, stop
            if collecting_vars and re.match(r'\s*namelist\s*/', line, re.IGNORECASE):
                collecting_vars = False
                current_namelist = None
                continue
                
            # If we're collecting variables for a namelist
            if collecting_vars and current_namelist:
                # Always process the current continuation line (strip trailing & if present)
                var_line = re.sub(r'&\s*$', '', raw_line).strip()
                if var_line:
                    if debug:
                        print(f"DEBUG [_extract_namelists] {filename}:{idx}: NAMELIST CONT -> {var_line}")
                    self._extract_vars_from_line(var_line, current_namelist, fortran_keywords, debug=debug)

                # If this line does not end with an ampersand, the namelist ended here
                if not raw_line.endswith('&'):
                    collecting_vars = False
                    current_namelist = None
    
    def _extract_vars_from_line(self, line: str, namelist: str, keywords: set, debug: bool = False):
        """Extract variable names from a line of namelist declaration."""
        # Split by comma and whitespace, filter valid identifiers
        potential_vars = re.findall(r'\b([a-zA-Z_][a-zA-Z0-9_]*)\b', line)
        
        for var in potential_vars:
            var_lower = var.lower()
            # Filter out Fortran keywords and single letters (usually loop indices)
            if (var_lower not in keywords and 
                len(var) > 1 and 
                not var.isdigit() and
                not re.match(r'^[a-z]$', var_lower)):  # Skip single letter variables
                self.namelists[namelist].add(var_lower)
                if debug:
                    print(f"DEBUG [_extract_vars_from_line] added {var_lower} to {namelist}")
    
    def _extract_json_reads(self, content: str):
        """Extract variables that are read from JSON configuration."""
        # Look for patterns like: json%get('NAMELIST.variable', temp_var, found)
        # Use single pattern that covers both direct and call forms
        patterns = [
            r"(?:call\s+)?json%get\s*\(\s*['\"](\w+)\.(\w+)['\"]"
        ]
        
        for pattern in patterns:
            matches = re.findall(pattern, content, re.IGNORECASE)
            for namelist, var in matches:
                var_lower = var.lower()
                self.json_reads.add(var_lower)
                self.json_read_counts[var_lower] += 1
    
    def _extract_broadcasts(self, content: str):
        """Extract variables that are broadcast via MPI."""
        # Look for MPI_BCAST calls
        # Pattern: call MPI_BCAST(variable_name, ...)
        broadcast_patterns = [
            r'(?:call\s+)?MPI_BCAST\s*\(\s*([a-zA-Z_]\w*)'
        ]
        
        for pattern in broadcast_patterns:
            matches = re.findall(pattern, content, re.IGNORECASE)
            for var in matches:
                # Filter out obvious non-variables
                if var.lower() not in ['comm3d', 'mpierr', 'mpi_integer', 'mpi_real', 'my_real', 'mpi_logical']:
                    var_lower = var.lower()
                    self.broadcasts.add(var_lower)
                    self.broadcast_counts[var_lower] += 1
    
    def parse_json_schema(self):
        """Parse the JSON schema to extract variable definitions."""
        print(f"Parsing JSON schema: {self.schema_file}")
        
        try:
            with open(self.schema_file, 'r') as f:
                schema = json.load(f)
            
            # Navigate the schema structure
            if 'properties' in schema:
                for namelist_name, namelist_def in schema['properties'].items():
                    namelist_name = namelist_name.upper()
                    
                    if isinstance(namelist_def, dict) and 'properties' in namelist_def:
                        for var_name in namelist_def['properties'].keys():
                            self.schema_vars[namelist_name].add(var_name.lower())
                            
        except Exception as e:
            print(f"Warning: Could not parse schema file: {e}")
    
    def generate_markdown_report(self) -> str:
        """Generate a comprehensive markdown report."""
        report = []
        
        report.append("# uDALES Variable Index Report")
        report.append("")
        report.append("This report shows the status of variables across different aspects of the uDALES code:")
        report.append("")
        report.append("- **Namelist**: Variable is defined in a Fortran namelist (✓/✗)")
        report.append("- **JSON Read**: Variable is read from JSON configuration (count or ✗)")
        report.append("- **Broadcast**: Variable is broadcast via MPI (count or ✗)")
        report.append("- **Schema**: Variable is present in the JSON schema (✓/✗)")
        report.append("")
        
        # Get all namelists (only from Fortran code, not schema-only namelists)
        all_namelists = set(self.namelists.keys())
        
        for namelist in sorted(all_namelists):
            # Add file information if available
            file_info = ""
            if namelist in self.namelist_files:
                file_info = f" *(defined in {self.namelist_files[namelist]})*"
            
            report.append(f"## {namelist} Namelist{file_info}")
            report.append("")
            
            # Get variables for this namelist (only include variables actually in Fortran namelists)
            namelist_vars = self.namelists.get(namelist, set())
            schema_vars = self.schema_vars.get(namelist, set())
            # Only report variables that are actually in Fortran namelists
            all_vars = namelist_vars
            
            if not all_vars:
                report.append("*No variables found for this namelist.*")
                report.append("")
                continue
            
            # Create table header
            report.append("| Variable | Namelist | JSON Read | Broadcast | Schema |")
            report.append("|----------|----------|-----------|-----------|--------|")
            
            # Add variables to table
            for var in sorted(all_vars):
                namelist_check = "✓" if var in namelist_vars else "✗"
                
                # Show count for JSON reads
                if var in self.json_reads:
                    json_count = self.json_read_counts.get(var, 1)
                    json_check = f"✓ ({json_count}x)" if json_count > 1 else "✓"
                else:
                    json_check = "✗"
                
                # Show count for broadcasts
                if var in self.broadcasts:
                    broadcast_count = self.broadcast_counts.get(var, 1)
                    broadcast_check = f"✓ ({broadcast_count}x)" if broadcast_count > 1 else "✓"
                else:
                    broadcast_check = "✗"
                
                schema_check = "✓" if var in schema_vars else "✗"
                
                report.append(f"| `{var}` | {namelist_check} | {json_check} | {broadcast_check} | {schema_check} |")
            
            report.append("")
        
        # Add summary statistics
        report.append("## Summary Statistics")
        report.append("")
        
        total_namelist_vars = sum(len(vars) for vars in self.namelists.values())
        total_json_vars = len(self.json_reads)
        total_broadcast_vars = len(self.broadcasts)
        total_schema_vars = sum(len(vars) for vars in self.schema_vars.values())
        
        report.append(f"- **Total namelist variables**: {total_namelist_vars}")
        report.append(f"- **Total JSON-readable variables**: {total_json_vars}")
        report.append(f"- **Total broadcast variables**: {total_broadcast_vars}")
        report.append(f"- **Total schema variables**: {total_schema_vars}")
        report.append(f"- **Total namelists found**: {len(all_namelists)}")
        report.append("")
        
        # Add potential issues
        report.append("## Potential Issues")
        report.append("")
        
        # Variables in schema but not in namelists
        schema_only_vars = set()
        for namelist, vars_set in self.schema_vars.items():
            namelist_vars = self.namelists.get(namelist, set())
            schema_only_vars.update(vars_set - namelist_vars)
        
        if schema_only_vars:
            report.append("### Variables in schema but not in namelists:")
            for var in sorted(schema_only_vars):
                report.append(f"- `{var}`")
            report.append("")
        
        # Variables in namelists but not in schema
        namelist_only_vars = set()
        for namelist, vars_set in self.namelists.items():
            schema_vars = self.schema_vars.get(namelist, set())
            namelist_only_vars.update(vars_set - schema_vars)
        
        if namelist_only_vars:
            report.append("### Variables in namelists but not in schema:")
            for var in sorted(namelist_only_vars):
                report.append(f"- `{var}`")
            report.append("")
        
        return "\n".join(report)
    
    def run_analysis(self) -> str:
        """Run the complete analysis and return the markdown report."""
        self.parse_fortran_files()
        self.parse_json_schema()
        return self.generate_markdown_report()


def main():
    parser = argparse.ArgumentParser(
        description="Analyze uDALES variables across namelists, JSON, broadcasts, and schema"
    )
    parser.add_argument(
        "--src-dir", 
        default="../../src", 
        help="Path to the src directory containing Fortran files"
    )
    parser.add_argument(
        "--schema-file", 
        default="../../docs/schemas/udales_input_schema.json",
        help="Path to the JSON schema file"
    )
    parser.add_argument(
        "--output", 
        default="variable_index_report.md",
        help="Output file for the markdown report"
    )
    
    args = parser.parse_args()
    
    # Resolve paths relative to script location
    script_dir = Path(__file__).parent
    src_dir = (script_dir / args.src_dir).resolve()
    schema_file = (script_dir / args.schema_file).resolve()
    
    if not src_dir.exists():
        print(f"Error: Source directory not found: {src_dir}")
        sys.exit(1)
    
    if not schema_file.exists():
        print(f"Error: Schema file not found: {schema_file}")
        sys.exit(1)
    
    print(f"Source directory: {src_dir}")
    print(f"Schema file: {schema_file}")
    print(f"Output file: {args.output}")
    print()
    
    # Run analysis
    indexer = UdalesVariableIndexer(src_dir, schema_file)
    report = indexer.run_analysis()
    
    # Write report
    with open(args.output, 'w') as f:
        f.write(report)
    
    print(f"\nReport generated: {args.output}")


if __name__ == "__main__":
    main()