#!/usr/bin/env python3
"""
Run all JSON input tests for u-DALES.

This script runs all JSON-related unit tests and integration tests
in the proper order with clear output.
"""

import sys
import subprocess
import unittest
from pathlib import Path

def run_command(cmd, description):
    """Run a command and report success/failure."""
    print(f"\n{'='*60}")
    print(f"Running: {description}")
    print(f"Command: {' '.join(cmd)}")
    print('='*60)
    
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        print("✓ PASSED")
        if result.stdout:
            print("Output:")
            print(result.stdout)
        return True
    except subprocess.CalledProcessError as e:
        print("✗ FAILED")
        print("Error output:")
        print(e.stderr)
        if e.stdout:
            print("Standard output:")
            print(e.stdout)
        return False

def main():
    """Run all JSON tests."""
    print("u-DALES JSON Input Test Suite")
    print("=" * 60)
    
    # Change to the json_test directory
    test_dir = Path(__file__).parent
    original_dir = Path.cwd()
    
    try:
        print(f"Test directory: {test_dir}")
        
        # List of tests to run
        tests = [
            (["python3", "-m", "unittest", "test_json_input.TestJSONInput.test_schema_extraction", "-v"], 
             "JSON Schema Extraction Test"),
            (["python3", "-m", "unittest", "test_json_input.TestJSONInput.test_random_value_generation", "-v"], 
             "Random Value Generation Test"),
            (["python3", "test_comprehensive_schema.py"], 
             "Comprehensive Schema Test (all 224 parameters)"),
            (["python3", "demo_json_testing.py"], 
             "JSON Testing Demonstration"),
        ]
        
        results = []
        for cmd, description in tests:
            # Run test from the json_test directory
            full_cmd = ["bash", "-c", f"cd {test_dir} && {' '.join(cmd)}"]
            success = run_command(full_cmd, description)
            results.append((description, success))
        
        # Summary
        print(f"\n{'='*60}")
        print("TEST SUMMARY")
        print('='*60)
        
        passed = 0
        total = len(results)
        
        for description, success in results:
            status = "✓ PASSED" if success else "✗ FAILED"
            print(f"{status:10} {description}")
            if success:
                passed += 1
        
        print(f"\nOverall: {passed}/{total} tests passed")
        
        if passed == total:
            print("🎉 All JSON tests passed!")
            return 0
        else:
            print("⚠️  Some tests failed. Check output above for details.")
            return 1
            
    except Exception as e:
        print(f"Error running tests: {e}")
        return 1

if __name__ == "__main__":
    sys.exit(main())