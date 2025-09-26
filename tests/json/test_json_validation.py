#!/usr/bin/env python3
"""
Comprehensive unit test for u-DALES JSON input functionality.

This test validates that:
1. JSON input files are correctly read by u-DALES with runmode=1001
2. Parameters are properly parsed and broadcast via MPI
3. Output values match input values for specified parameters
4. Unspecified parameters get proper default values from the schema

The test uses runmode=1001 which triggers the tests_json routine in u-DALES.
"""

import os
import json
import tempfile
import subprocess
import shutil
from pathlib import Path
import unittest

class UDalesJSONTest(unittest.TestCase):
    """Test class for u-DALES JSON input validation."""
    
    @classmethod
    def setUpClass(cls):
        """Set up test environment once for all tests."""
        cls.project_dir = Path(__file__).resolve().parents[2]
        cls.test_dir = Path(__file__).resolve().parent
        cls.udales_exe = cls.project_dir / "bin" / "u-dales"
        cls.ud_nam2json_exe = cls.project_dir / "bin" / "ud_nam2json"
        cls.schema_path = cls.project_dir / "docs" / "schemas" / "udales_input_schema.json"
        
        # Load schema for default value comparison
        with open(cls.schema_path, 'r') as f:
            cls.schema = json.load(f)
        
        # Verify executables exist
        if not cls.udales_exe.exists():
            raise FileNotFoundError(f"u-DALES executable not found: {cls.udales_exe}")
        if not cls.ud_nam2json_exe.exists():
            raise FileNotFoundError(f"ud_nam2json executable not found: {cls.ud_nam2json_exe}")
        
        print(f"Using u-DALES executable: {cls.udales_exe}")
        print(f"Using ud_nam2json: {cls.ud_nam2json_exe}")
    
    def setUp(self):
        """Set up for each individual test."""
        self.temp_dir = Path(tempfile.mkdtemp(prefix="udales_test_"))
        
    def tearDown(self):
        """Clean up after each test."""
        if self.temp_dir.exists():
            shutil.rmtree(self.temp_dir)
    
    def get_default_value(self, section_name, param_name):
        """Get default value for a parameter from the schema."""
        try:
            section_schema = self.schema["properties"][section_name]
            param_schema = section_schema["properties"][param_name]
            return param_schema.get("default")
        except KeyError:
            return None
    
    def run_udales_with_json(self, json_config_path):
        """
        Run u-DALES with JSON configuration using appropriate number of MPI processes.
        
        Returns:
            tuple: (return_code, stdout, stderr, output_files)
        """
        # Read config to determine MPI process count
        with open(json_config_path, 'r') as f:
            config = json.load(f)
        
        nprocx = config["RUN"]["nprocx"]
        nprocy = config["RUN"]["nprocy"]
        total_procs = nprocx * nprocy
        
        print(f"Running u-DALES with {total_procs} processes ({nprocx}x{nprocy})")
        
        # Copy JSON config to temp directory
        config_name = json_config_path.name
        temp_config = self.temp_dir / config_name
        shutil.copy2(json_config_path, temp_config)
        
        # Set up environment and run u-DALES with mpirun
        # Load necessary modules and source environment
        env_script = "/rds/general/user/mvr/home/udales/udales.env"
        
        cmd = f"""
        module load intel/2025a netCDF/4.9.2-iimpi-2023a netCDF-Fortran/4.6.1-iimpi-2023a FFTW/3.3.9-intel-2021a CMake/3.29.3-GCCcore-13.3.0 git/2.45.1-GCCcore-13.3.0 && \\
        source {env_script} && \\
        mpirun -n {total_procs} {self.udales_exe} {temp_config}
        """
        
        print(f"Command: {cmd.strip()}")
        
        result = subprocess.run(
            cmd,
            shell=True,  # Use shell to load modules and source environment
            cwd=self.temp_dir,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
            timeout=300  # 5 minute timeout
        )
        
        # Find output files
        output_files = list(self.temp_dir.glob("*"))
        
        return result.returncode, result.stdout, result.stderr, output_files
    
    def extract_output_namelists(self):
        """
        Extract namelist files from u-DALES output.
        
        The runmode=1001 should generate namelist files that we can convert back to JSON.
        """
        namelist_files = []
        for pattern in ["namoptions.*", "namelist.*", "*.namelist"]:
            namelist_files.extend(self.temp_dir.glob(pattern))
        
        return namelist_files
    
    def convert_namelist_to_json(self, namelist_path):
        """Convert a namelist file to JSON using ud_nam2json."""
        
        # Run ud_nam2json on the namelist
        result = subprocess.run(
            [str(self.ud_nam2json_exe), str(namelist_path)],
            cwd=namelist_path.parent,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True
        )
        
        if result.returncode != 0:
            raise RuntimeError(f"ud_nam2json failed: {result.stderr}")
        
        # Find the generated JSON file
        json_files = list(namelist_path.parent.glob("parameters.*"))
        if not json_files:
            raise RuntimeError("ud_nam2json did not generate expected JSON file")
        
        # Return the content of the generated JSON
        with open(json_files[0], 'r') as f:
            return json.load(f)
    
    def compare_configs(self, input_config, output_config, test_name):
        """
        Compare input and output configurations.
        
        Args:
            input_config: Original JSON configuration
            output_config: Configuration extracted from u-DALES output
            test_name: Name of the test for reporting
        """
        print(f"\n=== Comparing configurations for {test_name} ===")
        
        errors = []
        matches = []
        defaults_checked = []
        
        # Compare each section
        for section_name, section_data in input_config.items():
            if not isinstance(section_data, dict):
                continue
                
            print(f"\nSection: {section_name}")
            
            if section_name not in output_config:
                errors.append(f"Section {section_name} missing from output")
                continue
            
            output_section = output_config[section_name]
            
            # Compare each parameter in the input
            for param_name, input_value in section_data.items():
                if param_name in output_section:
                    output_value = output_section[param_name]
                    if input_value == output_value:
                        matches.append(f"{section_name}.{param_name}: {input_value} ✓")
                        print(f"  ✓ {param_name}: {input_value}")
                    else:
                        error_msg = f"{section_name}.{param_name}: expected {input_value}, got {output_value}"
                        errors.append(error_msg)
                        print(f"  ✗ {param_name}: expected {input_value}, got {output_value}")
                else:
                    error_msg = f"{section_name}.{param_name} missing from output"
                    errors.append(error_msg)
                    print(f"  ✗ {param_name}: missing from output")
        
        # Check that parameters not in input have default values
        for section_name, section_data in output_config.items():
            if not isinstance(section_data, dict):
                continue
                
            input_section = input_config.get(section_name, {})
            
            for param_name, output_value in section_data.items():
                if param_name not in input_section:
                    # This parameter wasn't in input, check if it matches default
                    default_value = self.get_default_value(section_name, param_name)
                    if default_value is not None:
                        if output_value == default_value:
                            defaults_checked.append(f"{section_name}.{param_name}: {output_value} (default) ✓")
                            print(f"  ✓ {param_name}: {output_value} (default)")
                        else:
                            error_msg = f"{section_name}.{param_name}: expected default {default_value}, got {output_value}"
                            errors.append(error_msg)
                            print(f"  ✗ {param_name}: expected default {default_value}, got {output_value}")
        
        # Summary
        print(f"\n=== Summary for {test_name} ===")
        print(f"Matches: {len(matches)}")
        print(f"Defaults verified: {len(defaults_checked)}")
        print(f"Errors: {len(errors)}")
        
        if errors:
            print("\nErrors found:")
            for error in errors[:10]:  # Show first 10 errors
                print(f"  - {error}")
            if len(errors) > 10:
                print(f"  ... and {len(errors) - 10} more errors")
        
        return len(errors) == 0, matches, defaults_checked, errors
    
    def run_json_test(self, config_filename):
        """Run a complete JSON test for a given configuration file."""
        
        config_path = self.test_dir / config_filename
        if not config_path.exists():
            self.fail(f"Configuration file not found: {config_path}")
        
        print(f"\n{'='*60}")
        print(f"Testing: {config_filename}")
        print(f"{'='*60}")
        
        # Load input configuration
        with open(config_path, 'r') as f:
            input_config = json.load(f)
        
        # Run u-DALES
        try:
            returncode, stdout, stderr, output_files = self.run_udales_with_json(config_path)
            
            print(f"u-DALES return code: {returncode}")
            if returncode != 0:
                print(f"STDOUT:\n{stdout}")
                print(f"STDERR:\n{stderr}")
                self.fail(f"u-DALES failed with return code {returncode}")
            
            # Extract namelist files from output
            namelist_files = self.extract_output_namelists()
            print(f"Found {len(namelist_files)} namelist files: {[f.name for f in namelist_files]}")
            
            if not namelist_files:
                # Look for any output files to debug
                print(f"Output files: {[f.name for f in output_files]}")
                self.fail("No namelist files found in u-DALES output")
            
            # Convert the first namelist back to JSON
            output_config = self.convert_namelist_to_json(namelist_files[0])
            
            # Compare configurations
            success, matches, defaults, errors = self.compare_configs(
                input_config, output_config, config_filename
            )
            
            self.assertTrue(success, f"Configuration comparison failed for {config_filename}")
            
        except subprocess.TimeoutExpired:
            self.fail(f"u-DALES timed out for {config_filename}")
        except Exception as e:
            self.fail(f"Test failed for {config_filename}: {e}")
    
    def test_baseline_config(self):
        """Test parameters.500 (baseline configuration)."""
        self.run_json_test("parameters.500")
    
    def test_enhanced_config(self):
        """Test parameters.501 (enhanced configuration)."""
        self.run_json_test("parameters.501")
    
    def test_random_config(self):
        """Test parameters.800 (fully random configuration)."""
        self.run_json_test("parameters.800")
    
    def test_example_configs(self):
        """Test all parameters.* files in the directory, using experiment number logic."""
        param_files = sorted(self.test_dir.glob("parameters.*"))
        for param_path in param_files:
            config_name = param_path.name
            # Extract experiment number (e.g., 001, 500, etc.)
            try:
                exp_num = int(config_name.split(".")[1])
            except (IndexError, ValueError):
                exp_num = None
            with self.subTest(config=config_name):
                # Case logic for different experiment numbers
                if exp_num == 500:
                    print("Detected baseline configuration (500)")
                elif exp_num == 501:
                    print("Detected enhanced configuration (501)")
                elif exp_num == 800:
                    print("Detected random configuration (800)")
                elif exp_num in {1,2,24,101,102,201,949,950,999}:
                    print(f"Detected example configuration ({exp_num})")
                else:
                    print(f"Detected unknown experiment number: {exp_num}")
                self.run_json_test(config_name)

if __name__ == "__main__":
    # Run tests with verbose output
    unittest.main(verbosity=2)