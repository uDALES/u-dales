#!/usr/bin/env python3
"""
Schema completeness checker for u-DALES.

This script runs u-DALES with the comprehensive configuration and analyzes 
the output to find parameters that might be missing from the schema.
"""

import re
import subprocess
import sys
from pathlib import Path
from typing import Set, List, Dict, Any
import json

class SchemaCompletenessChecker:
    """Check for parameters missing from the schema by analyzing u-DALES output."""
    
    def __init__(self):
        self.project_dir = Path(__file__).resolve().parents[2]
        self.test_dir = Path(__file__).parent / "comprehensive_test_output"
        self.schema_path = self.project_dir / "docs" / "schemas" / "udales_input_schema.json"
        
        # Load schema parameters
        with open(self.schema_path, 'r') as f:
            schema = json.load(f)
        
        self.schema_params = set()
        for section_data in schema.get('properties', {}).values():
            if isinstance(section_data, dict) and 'properties' in section_data:
                self.schema_params.update(section_data['properties'].keys())
        
        print(f"Loaded {len(self.schema_params)} parameters from schema")
    
    def find_udales_executable(self) -> Path:
        """Find the u-DALES executable."""
        possible_paths = [
            self.project_dir / "build" / "release" / "u-dales",
            self.project_dir / "build" / "debug" / "u-dales",
            self.project_dir / "u-dales",
            Path("./u-dales")
        ]
        
        for exe_path in possible_paths:
            if exe_path.exists():
                return exe_path
        
        raise FileNotFoundError(f"u-DALES executable not found. Looked in: {[str(p) for p in possible_paths]}")
    
    def run_udales_with_config(self, config_type: str) -> tuple:
        """Run u-DALES with comprehensive configuration."""
        original_cwd = Path.cwd()
        
        try:
            # Change to test directory
            self.test_dir.absolute()
            if not self.test_dir.exists():
                raise FileNotFoundError(f"Test directory not found: {self.test_dir}")
            
            import os
            os.chdir(self.test_dir)
            
            # Find executable
            udales_exe = self.find_udales_executable()
            
            # Set up configuration
            if config_type == 'json':
                config_file = self.test_dir / "config_comprehensive.json"
                if config_file.exists():
                    # Create symlink to config.json
                    config_link = self.test_dir / "config.json"
                    if config_link.exists():
                        config_link.unlink()
                    config_link.symlink_to(config_file.name)
                else:
                    raise FileNotFoundError(f"Comprehensive JSON config not found: {config_file}")
            
            elif config_type == 'namelist':
                namelist_file = self.test_dir / "namoptions.comprehensive"
                if namelist_file.exists():
                    # Create symlink to standard namelist name
                    namelist_link = self.test_dir / "namoptions.001"
                    if namelist_link.exists():
                        namelist_link.unlink()
                    namelist_link.symlink_to(namelist_file.name)
                else:
                    raise FileNotFoundError(f"Comprehensive namelist not found: {namelist_file}")
            
            print(f"Running u-DALES with {config_type} configuration...")
            print(f"Working directory: {self.test_dir}")
            print(f"Executable: {udales_exe}")
            
            # Run u-DALES with timeout
            result = subprocess.run(
                [str(udales_exe)],
                capture_output=True,
                text=True,
                timeout=60  # 1 minute timeout
            )
            
            return result.returncode, result.stdout, result.stderr
            
        except subprocess.TimeoutExpired:
            return -1, "", "u-DALES execution timed out (60s)"
        except Exception as e:
            return -2, "", f"Error running u-DALES: {str(e)}"
        finally:
            os.chdir(original_cwd)
    
    def extract_parameters_from_output(self, stdout: str, stderr: str) -> Set[str]:
        """Extract parameter names from u-DALES output."""
        all_output = stdout + "\n" + stderr
        found_params = set()
        
        # Common patterns for parameter output in Fortran programs
        patterns = [
            # Direct parameter assignments
            r'(\w+)\s*=\s*[^\s]+',
            # Reading parameter messages
            r'(?:Reading|Read|Setting|Set)\s+(\w+)(?:\s*[:=]\s*[^\s]+)?',
            # Parameter initialization messages  
            r'(?:Initialize|Init|Parameter)\s+(\w+)',
            # Namelist parameter patterns
            r'&\w+[^/]*\s+(\w+)\s*=',
            # Error messages about parameters
            r'(?:Unknown|Invalid|Missing)\s+(?:parameter|variable)\s+[\'"]?(\w+)[\'"]?',
            # JSON parameter patterns
            r'[\'"](\w+)[\'"]\s*:\s*[^\s,}]+',
        ]
        
        for pattern in patterns:
            matches = re.finditer(pattern, all_output, re.IGNORECASE | re.MULTILINE)
            for match in matches:
                param_name = match.group(1).strip()
                # Filter out obvious non-parameters
                if (len(param_name) > 2 and 
                    param_name.isalnum() and 
                    not param_name.isupper() and  # Avoid constants
                    not param_name.isdigit()):    # Avoid numbers
                    found_params.add(param_name)
        
        # Also look for specific u-DALES parameter patterns
        udales_patterns = [
            r'namelist\s+(\w+)',
            r'variable\s+(\w+)',
            r'option\s+(\w+)',
        ]
        
        for pattern in udales_patterns:
            matches = re.finditer(pattern, all_output, re.IGNORECASE)
            for match in matches:
                found_params.add(match.group(1))
        
        return found_params
    
    def analyze_missing_parameters(self, found_params: Set[str]) -> Dict[str, Any]:
        """Analyze parameters found in output vs schema."""
        missing_from_schema = found_params - self.schema_params
        in_schema_not_found = self.schema_params - found_params
        
        # Filter out likely false positives
        filtered_missing = set()
        for param in missing_from_schema:
            # Skip if it looks like a module name, function, or other non-parameter
            if (not param.startswith('mod') and 
                not param.startswith('sub') and
                not param.endswith('_mod') and
                len(param) > 2 and
                param.islower()):
                filtered_missing.add(param)
        
        return {
            'schema_params': len(self.schema_params),
            'found_in_output': len(found_params),
            'missing_from_schema': sorted(filtered_missing),
            'in_schema_not_found': sorted(in_schema_not_found),
            'overlap': len(self.schema_params & found_params)
        }
    
    def run_comprehensive_check(self) -> Dict[str, Any]:
        """Run comprehensive parameter completeness check."""
        print("=" * 80)
        print("SCHEMA COMPLETENESS CHECK")
        print("=" * 80)
        
        results = {}
        
        # Test both input methods
        for config_type in ['json', 'namelist']:
            print(f"\n{config_type.upper()} INPUT TEST:")
            print("-" * 40)
            
            try:
                returncode, stdout, stderr = self.run_udales_with_config(config_type)
                
                print(f"Return code: {returncode}")
                if returncode == 0:
                    print("✓ u-DALES ran successfully")
                else:
                    print(f"⚠ u-DALES exited with code {returncode}")
                
                # Extract parameters from output
                found_params = self.extract_parameters_from_output(stdout, stderr)
                print(f"Found {len(found_params)} parameter references in output")
                
                # Analyze
                analysis = self.analyze_missing_parameters(found_params)
                
                print(f"Schema parameters: {analysis['schema_params']}")
                print(f"Output parameters: {analysis['found_in_output']}")
                print(f"Overlap: {analysis['overlap']}")
                print(f"Missing from schema: {len(analysis['missing_from_schema'])}")
                
                if analysis['missing_from_schema']:
                    print("Potentially missing parameters:")
                    for param in analysis['missing_from_schema'][:10]:  # Show first 10
                        print(f"  - {param}")
                    if len(analysis['missing_from_schema']) > 10:
                        print(f"  ... and {len(analysis['missing_from_schema']) - 10} more")
                
                results[config_type] = {
                    'returncode': returncode,
                    'stdout': stdout[:1000],  # First 1000 chars
                    'stderr': stderr[:1000],  # First 1000 chars  
                    'found_params': found_params,
                    'analysis': analysis
                }
                
            except Exception as e:
                print(f"✗ Error testing {config_type}: {e}")
                results[config_type] = {'error': str(e)}
        
        return results
    
    def generate_report(self, results: Dict[str, Any]):
        """Generate a comprehensive report."""
        print("\n" + "=" * 80)
        print("COMPLETENESS ANALYSIS REPORT")
        print("=" * 80)
        
        all_missing = set()
        
        for config_type, result in results.items():
            if 'analysis' in result:
                analysis = result['analysis']
                missing = set(analysis.get('missing_from_schema', []))
                all_missing.update(missing)
                
                print(f"\n{config_type.upper()} Results:")
                print(f"  Schema coverage: {analysis['overlap']}/{analysis['schema_params']} ({analysis['overlap']/analysis['schema_params']*100:.1f}%)")
                
                if missing:
                    print(f"  Potentially missing: {len(missing)} parameters")
        
        if all_missing:
            print(f"\nAll potentially missing parameters ({len(all_missing)}):")
            for param in sorted(all_missing):
                print(f"  - {param}")
            
            # Save to file
            report_file = self.test_dir / "missing_parameters_report.txt"
            with open(report_file, 'w') as f:
                f.write("Parameters potentially missing from schema:\n")
                f.write("=" * 50 + "\n")
                for param in sorted(all_missing):
                    f.write(f"{param}\n")
            
            print(f"\nDetailed report saved to: {report_file}")
        else:
            print("\n✓ No obviously missing parameters detected!")


def main():
    """Run the schema completeness check."""
    try:
        checker = SchemaCompletenessChecker()
        results = checker.run_comprehensive_check()
        checker.generate_report(results)
        
    except Exception as e:
        print(f"Error running completeness check: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()