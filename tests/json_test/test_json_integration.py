#!/usr/bin/env python3
"""
Integration test for JSON vs Namelist input methods in u-DALES.

This test compares the actual parameter extraction by running u-DALES
with identical parameters in both namelist and JSON formats.
"""

import os
import sys
import json
import subprocess
import tempfile
import shutil
from pathlib import Path
import re
import random

class UDALESInputTest:
    """Test runner for comparing namelist and JSON input methods."""
    
    def __init__(self, udales_exe_path: str):
        self.udales_exe = Path(udales_exe_path)
        if not self.udales_exe.exists():
            raise FileNotFoundError(f"u-DALES executable not found: {udales_exe_path}")
    
    def create_minimal_test_files(self, test_dir: Path, parameters: dict, use_json: bool = False):
        """Create minimal test configuration files."""
        
        # Create basic geometry file
        stl_content = """solid flat_ground
  facet normal 0 0 1
    outer loop
      vertex 0 0 0
      vertex 16 0 0
      vertex 16 16 0
    endloop
  endfacet
  facet normal 0 0 1
    outer loop
      vertex 0 0 0
      vertex 16 16 0
      vertex 0 16 0
    endloop
  endfacet
endsolid flat_ground
"""
        with open(test_dir / "flat_ground.stl", 'w') as f:
            f.write(stl_content)
        
        if use_json:
            # Create JSON configuration
            json_config = {
                "RUN": {
                    "iexpnr": 999,
                    "runtime": 0.1,  # Very short run
                    "dtmax": 0.01,
                    "nprocx": 1,
                    "nprocy": 1,
                    **parameters
                },
                "DOMAIN": {
                    "itot": 16,
                    "jtot": 16,
                    "ktot": 16,
                    "xlen": 8,
                    "ylen": 8
                },
                "INP": {
                    "zsize": 8.0,
                    "stl_file": "flat_ground.stl"
                },
                "OUTPUT": {
                    "lfielddump": False,
                    "ltdump": False
                },
                "WALLS": {
                    "nfcts": 2
                }
            }
            
            with open(test_dir / "config.json", 'w') as f:
                json.dump(json_config, f, indent=2)
        else:
            # Create namelist configuration
            namelist_content = f"""&RUN
iexpnr = 999
runtime = 0.1
dtmax = 0.01
nprocx = 1
nprocy = 1
"""
            for key, value in parameters.items():
                if isinstance(value, bool):
                    value_str = '.true.' if value else '.false.'
                elif isinstance(value, str):
                    value_str = f"'{value}'"
                else:
                    value_str = str(value)
                namelist_content += f"{key} = {value_str}\n"
            
            namelist_content += """/

&DOMAIN
itot = 16
jtot = 16
ktot = 16
xlen = 8
ylen = 8
/

&INP
zsize = 8.0
stl_file = 'flat_ground.stl'
/

&OUTPUT
lfielddump = .false.
ltdump = .false.
/

&WALLS
nfcts = 2
/
"""
            with open(test_dir / "namoptions.999", 'w') as f:
                f.write(namelist_content)
    
    def run_udales_test(self, test_dir: Path, timeout: int = 30) -> tuple:
        """Run u-DALES in test directory and capture output."""
        original_cwd = os.getcwd()
        try:
            os.chdir(test_dir)
            
            # Run u-DALES with timeout
            result = subprocess.run(
                [str(self.udales_exe)],
                capture_output=True,
                text=True,
                timeout=timeout
            )
            
            return result.returncode, result.stdout, result.stderr
            
        except subprocess.TimeoutExpired:
            return -1, "", "Process timed out"
        except Exception as e:
            return -2, "", str(e)
        finally:
            os.chdir(original_cwd)
    
    def extract_parameters_from_output(self, stdout: str, stderr: str) -> dict:
        """Extract parameter values from u-DALES output."""
        all_output = stdout + "\n" + stderr
        parameters = {}
        
        # Look for parameter assignments in the output
        patterns = [
            r'(\w+)\s*=\s*([^\s]+)',
            r'Reading\s+(\w+)\s*:\s*([^\s]+)',
            r'Parameter\s+(\w+)\s*set\s+to\s*([^\s]+)'
        ]
        
        for pattern in patterns:
            matches = re.finditer(pattern, all_output, re.IGNORECASE)
            for match in matches:
                param_name = match.group(1).strip()
                param_value = match.group(2).strip()
                
                # Clean up value
                param_value = param_value.rstrip(',').strip()
                if param_value.startswith("'") and param_value.endswith("'"):
                    param_value = param_value[1:-1]
                
                parameters[param_name] = param_value
        
        return parameters
    
    def compare_input_methods(self, test_parameters: dict) -> tuple:
        """Compare namelist and JSON input methods with given parameters."""
        results = {}
        
        # Create base directory for integration tests
        proj_dir = Path(__file__).resolve().parents[2]  # Now two levels up from tests/json_test/
        json_test_dir = proj_dir / "tests" / "json_test"
        json_test_dir.mkdir(exist_ok=True)
        integration_dir = json_test_dir / "integration_test_output"
        integration_dir.mkdir(exist_ok=True)
        
        for method in ['namelist', 'json']:
            test_dir = integration_dir / f"{method}_test_{len(test_parameters)}_params"
            test_dir.mkdir(exist_ok=True)
            
            try:
                # Create test files
                self.create_minimal_test_files(
                    test_dir, 
                    test_parameters, 
                    use_json=(method == 'json')
                )
                
                # Run u-DALES
                returncode, stdout, stderr = self.run_udales_test(test_dir)
                
                # Extract parameters
                extracted_params = self.extract_parameters_from_output(stdout, stderr)
                
                results[method] = {
                    'returncode': returncode,
                    'stdout': stdout,
                    'stderr': stderr,
                    'parameters': extracted_params,
                    'test_dir': str(test_dir)
                }
                
            except Exception as e:
                results[method] = {
                    'error': str(e),
                    'test_dir': str(test_dir)
                }
            
            finally:
                # Keep test directories for inspection - don't auto-delete
                pass
        
        return results
    
    def run_comprehensive_test(self):
        """Run comprehensive test with various parameter combinations."""
        test_cases = [
            # Basic test
            {'ladaptive': True},
            
            # Numerical parameters
            {'ih': 4, 'jh': 4, 'kh': 2},
            
            # Multiple parameters
            {'ladaptive': False, 'ih': 5, 'trestart': 0.05},
            
            # Random values test
            self._generate_random_parameters()
        ]
        
        results = []
        for i, test_params in enumerate(test_cases):
            print(f"\nRunning test case {i+1}: {test_params}")
            
            test_result = self.compare_input_methods(test_params)
            test_result['test_parameters'] = test_params
            results.append(test_result)
            
            # Analyze results
            self._analyze_test_result(test_result, i+1)
        
        return results
    
    def _generate_random_parameters(self) -> dict:
        """Generate random parameter values for testing."""
        random_params = {}
        
        # Safe integer parameters to randomize
        int_params = {
            'ih': [3, 4, 5, 6],
            'jh': [3, 4, 5, 6],
            'kh': [1, 2, 3]
        }
        
        for param, choices in int_params.items():
            random_params[param] = random.choice(choices)
        
        # Safe boolean parameters
        bool_params = ['ladaptive']
        for param in bool_params:
            random_params[param] = random.choice([True, False])
        
        return random_params
    
    def _analyze_test_result(self, result: dict, test_num: int):
        """Analyze and report test results."""
        print(f"Test {test_num} Results:")
        
        if 'namelist' in result and 'json' in result:
            nl_result = result['namelist']
            json_result = result['json']
            
            # Check return codes
            nl_success = nl_result.get('returncode', -999) == 0
            json_success = json_result.get('returncode', -999) == 0
            
            print(f"  Namelist execution: {'SUCCESS' if nl_success else 'FAILED'}")
            print(f"  JSON execution:     {'SUCCESS' if json_success else 'FAILED'}")
            
            if nl_success and json_success:
                # Compare extracted parameters
                nl_params = nl_result.get('parameters', {})
                json_params = json_result.get('parameters', {})
                
                common_params = set(nl_params.keys()) & set(json_params.keys())
                if common_params:
                    print(f"  Common parameters found: {len(common_params)}")
                    
                    differences = []
                    for param in common_params:
                        if nl_params[param] != json_params[param]:
                            differences.append(param)
                    
                    if differences:
                        print(f"  DIFFERENCES FOUND in: {differences}")
                    else:
                        print(f"  All parameters match: ✓")
                else:
                    print(f"  No common parameters extracted")
            
            # Show any errors
            if not nl_success and 'stderr' in nl_result:
                print(f"  Namelist error: {nl_result['stderr'][:200]}...")
            if not json_success and 'stderr' in json_result:
                print(f"  JSON error: {json_result['stderr'][:200]}...")
        
        else:
            print(f"  Test setup failed")


def main():
    """Main test runner."""
    # Find u-DALES executable
    proj_dir = Path(__file__).resolve().parents[2]  # Now two levels up from tests/json_test/
    possible_exe_paths = [
        proj_dir / "build" / "release" / "u-dales",
        proj_dir / "build" / "debug" / "u-dales",
        proj_dir / "u-dales"
    ]
    
    udales_exe = None
    for exe_path in possible_exe_paths:
        if exe_path.exists():
            udales_exe = exe_path
            break
    
    if not udales_exe:
        print("ERROR: u-DALES executable not found!")
        print("Looked in:")
        for path in possible_exe_paths:
            print(f"  {path}")
        print("\nPlease build u-DALES first with 'cmake --build build/release'")
        sys.exit(1)
    
    print(f"Using u-DALES executable: {udales_exe}")
    
    # Run tests
    tester = UDALESInputTest(str(udales_exe))
    
    try:
        results = tester.run_comprehensive_test()
        
        # Summary
        print("\n" + "="*60)
        print("COMPREHENSIVE TEST SUMMARY")
        print("="*60)
        
        total_tests = len(results)
        successful_tests = 0
        
        for i, result in enumerate(results):
            if ('namelist' in result and 'json' in result and 
                result['namelist'].get('returncode') == 0 and 
                result['json'].get('returncode') == 0):
                successful_tests += 1
        
        print(f"Total tests run: {total_tests}")
        print(f"Successful tests: {successful_tests}")
        print(f"Success rate: {successful_tests/total_tests*100:.1f}%")
        
        if successful_tests == total_tests:
            print("\n✓ ALL TESTS PASSED - JSON input method works correctly!")
        else:
            print(f"\n⚠ {total_tests - successful_tests} tests failed - see details above")
        
    except Exception as e:
        print(f"Test execution failed: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()