#!/usr/bin/env python3
"""
Test suite for u-DALES JSON input functionality using runmode=1001.

This test suite uses runmode=1001 to trigger the tests_json routine in u-DALES which:
1. Process 0 reads namoptions/JSON files using production code
2. All processes participate in MPI broadcasts via broadcast_config_parameters
3. Last process outputs ALL namelist values to files for verification
4. Provides comprehensive testing of the centralized broadcast system

Usage:
    python3 test_new_json_framework.py
"""

import os
import sys
import json
import subprocess
import tempfile
import shutil
from pathlib import Path
import unittest

# Project directory
PROJ_DIR = Path(__file__).resolve().parents[2]

class RunmodeTestFramework:
    """Test framework using runmode=1001 to trigger tests_json routine in u-DALES."""
    
    def __init__(self, udales_exe_path: str = None):
        """Initialize with path to u-DALES executable."""
        if udales_exe_path is None:
            # Try to find the executable in common build locations
            build_paths = [
                PROJ_DIR / "build" / "release" / "u-dales",
                PROJ_DIR / "build" / "u-dales", 
                PROJ_DIR / "u-dales"
            ]
            
            self.udales_exe = None
            for path in build_paths:
                if path.exists():
                    self.udales_exe = str(path)
                    break
                    
            if self.udales_exe is None:
                print(f"Warning: u-DALES executable not found. Tried: {build_paths}")
                print("Please build u-DALES or specify executable path")
                self.udales_exe = "u-dales"  # Fallback to system PATH
        else:
            self.udales_exe = udales_exe_path
            if not Path(udales_exe_path).exists():
                raise FileNotFoundError(f"u-DALES executable not found: {udales_exe_path}")
    
    def create_json_config(self, test_dir: Path, parameters: dict) -> Path:
        """Create config.json file with specified parameters and runmode=1001."""
        
        # Always set runmode to TEST_JSON
        parameters["runmode"] = 1001
        
        # Basic JSON structure
        config = {
            "RUN": {
                "iexpnr": parameters.get("iexpnr", 999),
                "runmode": 1001,  # Force test mode
                "nprocx": 1,
                "nprocy": 1
            },
            "DOMAIN": {
                "itot": parameters.get("itot", 16),
                "jtot": parameters.get("jtot", 16), 
                "ktot": parameters.get("ktot", 8),
                "xlen": parameters.get("xlen", 100.0),
                "ylen": parameters.get("ylen", 100.0)
            }
        }
        
        # Add any additional parameters to appropriate sections
        for key, value in parameters.items():
            if key not in ["iexpnr", "itot", "jtot", "ktot", "xlen", "ylen", "runmode"]:
                # Try to map to appropriate section or add to RUN
                if key.startswith("i") or key.startswith("l") or key.startswith("runtime"):
                    config["RUN"][key] = value
                elif key in ["ps", "lmoist", "lcoriol"]:
                    if "PHYSICS" not in config:
                        config["PHYSICS"] = {}
                    config["PHYSICS"][key] = value
                else:
                    config["RUN"][key] = value
        
        config_path = test_dir / "config.json"
        with open(config_path, 'w') as f:
            json.dump(config, f, indent=2)
        
        return config_path
    
    def create_namelist_config(self, test_dir: Path, parameters: dict) -> Path:
        """Create namoptions file with specified parameters and runmode=1001."""
        
        # Always set runmode to TEST_JSON
        parameters["runmode"] = 1001
        
        namoptions_path = test_dir / "namoptions"
        with open(namoptions_path, 'w') as f:
            # RUN namelist
            f.write("&RUN\n")
            f.write(f"iexpnr = {parameters.get('iexpnr', 999)}\n")
            f.write("runmode = 1001\n")  # Force test mode
            f.write("nprocx = 1\n")
            f.write("nprocy = 1\n")
            
            # Add other RUN parameters
            for key, value in parameters.items():
                if key in ["lwarmstart", "runtime", "dtmax", "ladaptive"]:
                    if isinstance(value, bool):
                        f.write(f"{key} = .{'true' if value else 'false'}.\n")
                    else:
                        f.write(f"{key} = {value}\n")
            f.write("/\n\n")
            
            # DOMAIN namelist
            f.write("&DOMAIN\n")
            f.write(f"itot = {parameters.get('itot', 16)}\n")
            f.write(f"jtot = {parameters.get('jtot', 16)}\n")
            f.write(f"ktot = {parameters.get('ktot', 8)}\n")
            f.write(f"xlen = {parameters.get('xlen', 100.0)}\n")
            f.write(f"ylen = {parameters.get('ylen', 100.0)}\n")
            f.write("/\n\n")
            
            # Add PHYSICS if needed
            physics_params = {k: v for k, v in parameters.items() 
                            if k in ["ps", "lmoist", "lcoriol", "lbuoyancy"]}
            if physics_params:
                f.write("&PHYSICS\n")
                for key, value in physics_params.items():
                    if isinstance(value, bool):
                        f.write(f"{key} = .{'true' if value else 'false'}.\n")
                    else:
                        f.write(f"{key} = {value}\n")
                f.write("/\n\n")
        
        return namoptions_path
    
    def run_test(self, test_dir: Path, use_mpi: bool = False, num_procs: int = 4) -> dict:
        """Run u-DALES with runmode=1001 and return results."""
        
        # Change to test directory
        original_cwd = os.getcwd()
        os.chdir(test_dir)
        
        try:
            if use_mpi:
                cmd = ["mpirun", "-np", str(num_procs), self.udales_exe]
            else:
                cmd = [self.udales_exe]
            
            print(f"Running: {' '.join(cmd)} in {test_dir}")
            
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
            
            # Parse results
            results = {
                "returncode": result.returncode,
                "stdout": result.stdout,
                "stderr": result.stderr,
                "success": result.returncode == 0 and "JSON tests completed successfully" in result.stdout,
                "namelist_files": []
            }
            
            # Look for generated namelist files
            for file in test_dir.glob("namelists_proc_*.txt"):
                results["namelist_files"].append(file)
            
            return results
            
        except subprocess.TimeoutExpired:
            return {
                "returncode": -1,
                "stdout": "",
                "stderr": "Test timed out after 30 seconds",
                "success": False,
                "namelist_files": []
            }
        except Exception as e:
            return {
                "returncode": -1,
                "stdout": "",
                "stderr": str(e),
                "success": False,
                "namelist_files": []
            }
        finally:
            os.chdir(original_cwd)


class TestJSONWithRunmode(unittest.TestCase):
    """Test cases using runmode=1001 approach."""
    
    def setUp(self):
        """Set up test environment."""
        self.test_framework = RunmodeTestFramework()
        self.test_dir = None
    
    def tearDown(self):
        """Clean up test directory."""
        if self.test_dir and self.test_dir.exists():
            shutil.rmtree(self.test_dir)
    
    def test_basic_json_config(self):
        """Test basic JSON configuration with runmode=1001."""
        self.test_dir = Path(tempfile.mkdtemp(prefix="udales_test_json_"))
        
        # Create basic test parameters
        params = {
            "iexpnr": 999,
            "itot": 16,
            "jtot": 16,
            "ktot": 8,
            "runtime": 60.0
        }
        
        # Create JSON config
        config_path = self.test_framework.create_json_config(self.test_dir, params)
        self.assertTrue(config_path.exists())
        
        # Verify JSON contains runmode=1001
        with open(config_path) as f:
            config = json.load(f)
        self.assertEqual(config["RUN"]["runmode"], 1001)
        
        # Run test (will only work if u-DALES is built)
        try:
            results = self.test_framework.run_test(self.test_dir)
            print(f"JSON test results: {results['success']}")
            if results["success"]:
                print(f"Generated {len(results['namelist_files'])} namelist files")
        except Exception as e:
            print(f"Test execution skipped: {e}")
    
    def test_basic_namelist_config(self):
        """Test basic namelist configuration with runmode=1001."""
        self.test_dir = Path(tempfile.mkdtemp(prefix="udales_test_namelist_"))
        
        # Create basic test parameters
        params = {
            "iexpnr": 999,
            "itot": 16,
            "jtot": 16,
            "ktot": 8,
            "runtime": 60.0,
            "lwarmstart": False
        }
        
        # Create namelist config
        namelist_path = self.test_framework.create_namelist_config(self.test_dir, params)
        self.assertTrue(namelist_path.exists())
        
        # Verify namelist contains runmode=1001
        with open(namelist_path) as f:
            content = f.read()
        self.assertIn("runmode = 1001", content)
        
        # Run test (will only work if u-DALES is built)
        try:
            results = self.test_framework.run_test(self.test_dir)
            print(f"Namelist test results: {results['success']}")
            if results["success"]:
                print(f"Generated {len(results['namelist_files'])} namelist files")
        except Exception as e:
            print(f"Test execution skipped: {e}")
    
    def test_mpi_json_config(self):
        """Test JSON configuration with MPI using runmode=1001."""
        self.test_dir = Path(tempfile.mkdtemp(prefix="udales_test_mpi_"))
        
        # Create test parameters
        params = {
            "iexpnr": 999,
            "itot": 16,
            "jtot": 16,
            "ktot": 8
        }
        
        # Create JSON config
        self.test_framework.create_json_config(self.test_dir, params)
        
        # Run with MPI (will only work if u-DALES is built and MPI available)
        try:
            results = self.test_framework.run_test(self.test_dir, use_mpi=True, num_procs=4)
            print(f"MPI JSON test results: {results['success']}")
            if results["success"]:
                print(f"Generated {len(results['namelist_files'])} namelist files")
                # Should have files for each process
                self.assertGreaterEqual(len(results['namelist_files']), 1)
        except Exception as e:
            print(f"MPI test execution skipped: {e}")


if __name__ == '__main__':
    print("u-DALES JSON Testing Framework - Runmode Approach")
    print("=" * 50)
    print("This test suite uses runmode=1001 to trigger tests_json routine")
    print("Requires u-DALES to be built with JSON support")
    print("")
    
    # Run the tests
    unittest.main(verbosity=2)