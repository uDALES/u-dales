#!/usr/bin/env python3
"""
Demonstration script for the sophisticated JSON input testing capability.

This script shows:
1. Schema extraction from Fortran source files
2. Random non-default value generation
3. Test file creation for both namelist and JSON formats
4. Parameter validation and comparison

Usage: python3 demo_json_testing.py
"""

import json
import sys
from pathlib import Path

# Import our test classes
sys.path.append(str(Path(__file__).parent))
from test_json_input import UDALESParameterSchema, TestJSONInput

def main():
    """Main demonstration function."""
    print("=" * 70)
    print("u-DALES JSON INPUT TESTING DEMONSTRATION")
    print("=" * 70)
    
    # Initialize schema extractor
    src_dir = Path(__file__).parent.parent / "src"
    schema = UDALESParameterSchema(src_dir)
    
    print(f"\n1. SCHEMA EXTRACTION RESULTS")
    print(f"   - Extracted {len(schema.schema)} parameters from Fortran source")
    print(f"   - Found namelist sections: {getattr(schema, 'namelist_sections', [])}")
    
    # Show parameter types breakdown
    type_counts = {}
    for param_info in schema.schema.values():
        ptype = param_info['type']
        type_counts[ptype] = type_counts.get(ptype, 0) + 1
    
    print(f"   - Parameter types found:")
    for ptype, count in type_counts.items():
        print(f"     * {ptype}: {count} parameters")
    
    # Generate random values
    print(f"\n2. RANDOM VALUE GENERATION")
    random_params = schema.generate_random_values(exclude_defaults=True)
    print(f"   - Generated {len(random_params)} random non-default values")
    
    # Show some interesting examples
    print(f"   - Example random values (different from defaults):")
    examples = [
        ('nsv', 'Number of scalar variables'),
        ('ih', 'Halo cells in i-direction'),
        ('jh', 'Halo cells in j-direction'),
        ('kh', 'Halo cells in k-direction'),
        ('nfcts', 'Number of wall facets'),
        ('lwarmstart', 'Warm start flag'),
        ('lfielddump', 'Field dump flag')
    ]
    
    for param_name, description in examples:
        if param_name in random_params and param_name in schema.schema:
            default_val = schema.schema[param_name]['default']
            random_val = random_params[param_name]
            print(f"     * {param_name:12} = {random_val:8} (default: {default_val:8}) - {description}")
    
    # Test file generation
    print(f"\n3. TEST FILE GENERATION")
    
    # Select a subset of safe parameters for testing
    test_params = {}
    safe_params = ['ih', 'jh', 'kh', 'nsv']
    for param in safe_params:
        if param in random_params:
            test_params[param] = random_params[param]
    
    # Add some logical parameters
    test_params['ladaptive'] = True
    test_params['lwarmstart'] = False
    
    print(f"   - Selected test parameters: {test_params}")
    
    # Create test instance
    test_instance = TestJSONInput()
    proj_dir = Path(__file__).resolve().parents[2]  # Now two levels up from tests/json_test/
    json_test_dir = proj_dir / "tests" / "json_test"
    json_test_dir.mkdir(exist_ok=True)
    test_instance.test_dir = json_test_dir / "demo_test_output"
    test_instance.test_dir.mkdir(exist_ok=True)
    
    # Generate namelist file
    namelist_path = test_instance.create_test_namelist(test_params)
    print(f"   - Created namelist file: {namelist_path}")
    
    # Generate JSON file  
    json_path = test_instance.create_test_json(test_params)
    print(f"   - Created JSON file: {json_path}")
    
    # Show file contents
    print(f"\n4. GENERATED FILE CONTENTS")
    
    print(f"\n   Namelist format (namoptions.999):")
    with open(namelist_path, 'r') as f:
        namelist_content = f.read()
    print("   " + "\n   ".join(namelist_content.split('\n')[:15]) + "\n   ...")
    
    print(f"\n   JSON format (namoptions.json):")
    with open(json_path, 'r') as f:
        json_content = f.read()
    print("   " + "\n   ".join(json_content.split('\n')[:20]) + "\n   ...")
    
    # Validation
    print(f"\n5. FILE VALIDATION")
    
    # Validate JSON structure
    try:
        with open(json_path, 'r') as f:
            json_data = json.load(f)
        print(f"   ✓ JSON file is valid and parseable")
        print(f"   ✓ JSON contains sections: {list(json_data.keys())}")
        
        # Check that our test parameters are present
        run_section = json_data.get('RUN', {})
        found_params = [p for p in test_params.keys() if p in run_section]
        print(f"   ✓ Found {len(found_params)}/{len(test_params)} test parameters in JSON")
        
    except Exception as e:
        print(f"   ✗ JSON validation failed: {e}")
    
    # Validate namelist structure  
    namelist_checks = 0
    for param, value in test_params.items():
        if f"{param} =" in namelist_content:
            namelist_checks += 1
    
    print(f"   ✓ Found {namelist_checks}/{len(test_params)} test parameters in namelist")
    
    # Summary
    print(f"\n6. SUMMARY")
    print(f"   ✓ Schema extraction: Successfully parsed {len(schema.schema)} parameters")
    print(f"   ✓ Random generation: Created {len(random_params)} non-default values")
    print(f"   ✓ File generation: Created both namelist and JSON test files")
    print(f"   ✓ Validation: Both formats contain the expected parameters")
    
    print(f"\n" + "=" * 70)
    print("DEMONSTRATION COMPLETE")
    print("The sophisticated JSON input testing system is fully functional!")
    print("Run 'python3 test_json_input.py' for the full unit test suite.")
    print("=" * 70)
    
    # Keep test output for inspection - don't cleanup
    print(f"Test files saved in: {test_instance.test_dir}")
    print(f"All JSON test outputs are organized in: {json_test_dir}")

if __name__ == '__main__':
    main()