#!/usr/bin/env python3
"""
Comprehensive test runner using runmode=1001 for u-DALES JSON functionality.

This script runs multiple test scenarios using the runmode=1001 approach:
1. Creates various parameter combinations 
2. Tests both JSON and namelist input methods
3. Runs with single process and MPI
4. Verifies that tests_json routine works correctly

Usage:
    python3 run_comprehensive_tests.py
"""

import sys
import json
import random
import tempfile
import shutil
from pathlib import Path
from test_runmode_framework import RunmodeTestFramework

def generate_test_scenarios():
    """Generate various test parameter scenarios."""
    
    scenarios = [
        {
            "name": "basic_cold_start",
            "params": {
                "iexpnr": 999,
                "itot": 16,
                "jtot": 16,
                "ktot": 8,
                "runtime": 60.0,
                "lwarmstart": False
            }
        },
        {
            "name": "larger_domain",
            "params": {
                "iexpnr": 888,
                "itot": 32,
                "jtot": 32,
                "ktot": 16,
                "xlen": 200.0,
                "ylen": 200.0,
                "runtime": 120.0
            }
        },
        {
            "name": "moist_physics",
            "params": {
                "iexpnr": 777,
                "itot": 16,
                "jtot": 16,
                "ktot": 8,
                "lmoist": True,
                "lcoriol": True,
                "ps": 101325.0
            }
        },
        {
            "name": "random_params",
            "params": {
                "iexpnr": random.randint(100, 999),
                "itot": random.choice([16, 32, 64]),
                "jtot": random.choice([16, 32, 64]),
                "ktot": random.choice([8, 16, 32]),
                "xlen": random.uniform(100.0, 500.0),
                "ylen": random.uniform(100.0, 500.0),
                "runtime": random.uniform(60.0, 300.0),
                "lmoist": random.choice([True, False]),
                "lcoriol": random.choice([True, False])
            }
        }
    ]
    
    return scenarios

def run_scenario(framework: RunmodeTestFramework, scenario: dict, test_type: str, use_mpi: bool = False):
    """Run a single test scenario."""
    
    print(f"\\n--- Running {scenario['name']} ({test_type}) {'with MPI' if use_mpi else 'single process'} ---")
    
    # Create temporary test directory
    test_dir = Path(tempfile.mkdtemp(prefix=f"udales_{scenario['name']}_{test_type}_"))
    
    try:
        # Create configuration
        if test_type == "json":
            config_path = framework.create_json_config(test_dir, scenario['params'])
            print(f"Created JSON config: {config_path}")
        else:
            config_path = framework.create_namelist_config(test_dir, scenario['params'])
            print(f"Created namelist config: {config_path}")
        
        # Run test
        results = framework.run_test(test_dir, use_mpi=use_mpi, num_procs=4 if use_mpi else 1)
        
        # Report results
        if results['success']:
            print(f"✅ SUCCESS - Generated {len(results['namelist_files'])} namelist files")
            
            # Show sample of generated files
            for i, namelist_file in enumerate(results['namelist_files'][:2]):  # Show first 2
                print(f"   Generated: {namelist_file.name}")
                if i == 0:  # Show content of first file
                    try:
                        with open(namelist_file, 'r') as f:
                            content = f.read()[:300]  # First 300 chars
                        print(f"   Sample content: {content}...")
                    except Exception:
                        pass
        else:
            print(f"❌ FAILED - Return code: {results['returncode']}")
            if results['stderr']:
                print(f"   Error: {results['stderr'][:200]}...")
            if results['stdout']:
                print(f"   Output: {results['stdout'][:200]}...")
        
        return results['success']
        
    except Exception as e:
        print(f"❌ EXCEPTION - {str(e)}")
        return False
        
    finally:
        # Clean up
        if test_dir.exists():
            shutil.rmtree(test_dir)

def main():
    """Run comprehensive test suite."""
    
    print("u-DALES Comprehensive JSON Testing - Runmode Approach")
    print("=" * 60)
    print("Testing JSON vs namelist input using runmode=1001")
    print("This triggers the tests_json routine for comprehensive validation")
    print("")
    
    # Initialize test framework
    try:
        framework = RunmodeTestFramework()
        print(f"Using u-DALES executable: {framework.udales_exe}")
    except Exception as e:
        print(f"Error initializing test framework: {e}")
        print("Please build u-DALES first:")
        print("  cd build/release && make")
        return 1
    
    # Generate test scenarios
    scenarios = generate_test_scenarios()
    print(f"Generated {len(scenarios)} test scenarios")
    
    # Track results
    total_tests = 0
    passed_tests = 0
    
    # Run each scenario with both JSON and namelist
    for scenario in scenarios:
        
        # Test with JSON
        total_tests += 1
        if run_scenario(framework, scenario, "json"):
            passed_tests += 1
        
        # Test with namelist
        total_tests += 1
        if run_scenario(framework, scenario, "namelist"):
            passed_tests += 1
        
        # Test with MPI (only for first scenario to save time)
        if scenario['name'] == 'basic_cold_start':
            total_tests += 1
            if run_scenario(framework, scenario, "json", use_mpi=True):
                passed_tests += 1
    
    # Final results
    print(f"\\n{'=' * 60}")
    print(f"COMPREHENSIVE TEST RESULTS")
    print(f"{'=' * 60}")
    print(f"Total tests run: {total_tests}")
    print(f"Tests passed:    {passed_tests}")
    print(f"Tests failed:    {total_tests - passed_tests}")
    print(f"Success rate:    {100 * passed_tests / total_tests:.1f}%")
    
    if passed_tests == total_tests:
        print("\\n🎉 ALL TESTS PASSED!")
        print("The runmode=1001 approach is working correctly.")
        return 0
    else:
        print(f"\\n⚠️  {total_tests - passed_tests} tests failed")
        print("Check u-DALES build and MPI availability")
        return 1

if __name__ == '__main__':
    sys.exit(main())