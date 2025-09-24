#!/usr/bin/env python3
"""
JSON Input Validation Tool for u-DALES

This script validates u-DALES JSON input files against the official schema.

Usage:
    python3 validate_json_input.py config.json
    python3 validate_json_input.py --schema-info
"""

import json
import sys
from pathlib import Path

def load_schema():
    """Load the u-DALES JSON schema."""
    # Find schema relative to this script
    script_dir = Path(__file__).parent
    schema_path = script_dir / "schemas" / "udales_input_schema.json"
    
    if not schema_path.exists():
        # Try alternative locations
        alternatives = [
            script_dir.parent / "docs" / "schemas" / "udales_input_schema.json",
            script_dir / "udales_input_schema.json"
        ]
        
        for alt_path in alternatives:
            if alt_path.exists():
                schema_path = alt_path
                break
        else:
            raise FileNotFoundError(f"Could not find schema file. Looked in: {[str(p) for p in [schema_path] + alternatives]}")
    
    with open(schema_path, 'r') as f:
        return json.load(f)

def validate_json_file(json_file_path, schema):
    """Validate a JSON file against the schema."""
    try:
        import jsonschema
    except ImportError:
        print("Warning: jsonschema library not available. Install with: pip install jsonschema")
        print("Performing basic structure validation only...")
        return basic_validation(json_file_path, schema)
    
    # Load and validate
    with open(json_file_path, 'r') as f:
        config = json.load(f)
    
    try:
        jsonschema.validate(config, schema)
        return True, "JSON file is valid according to schema"
    except jsonschema.ValidationError as e:
        return False, f"Validation error: {e.message}"
    except jsonschema.SchemaError as e:
        return False, f"Schema error: {e.message}"

def basic_validation(json_file_path, schema):
    """Basic validation without jsonschema library."""
    with open(json_file_path, 'r') as f:
        config = json.load(f)
    
    # Check structure
    schema_props = schema.get('properties', {})
    issues = []
    
    # Check that config sections exist in schema
    for section in config.keys():
        if section not in schema_props:
            issues.append(f"Unknown section: {section}")
    
    # Check parameters in each section
    for section, section_data in config.items():
        if section in schema_props:
            schema_section = schema_props[section].get('properties', {})
            if isinstance(section_data, dict):
                for param in section_data.keys():
                    if param not in schema_section:
                        issues.append(f"Unknown parameter in {section}: {param}")
    
    if issues:
        return False, "Basic validation issues: " + "; ".join(issues)
    else:
        return True, "Basic structure validation passed"

def show_schema_info(schema):
    """Display information about the schema."""
    print("u-DALES JSON Schema Information")
    print("=" * 40)
    print(f"Title: {schema.get('title', 'N/A')}")
    print(f"Description: {schema.get('description', 'N/A')}")
    print(f"Schema version: {schema.get('$schema', 'N/A')}")
    print()
    
    properties = schema.get('properties', {})
    print(f"Available sections ({len(properties)}):")
    
    total_params = 0
    for section_name, section_data in properties.items():
        if isinstance(section_data, dict):
            section_props = section_data.get('properties', {})
            param_count = len(section_props)
            total_params += param_count
            
            print(f"  {section_name:12} - {param_count:3} parameters")
            if section_data.get('description'):
                print(f"               {section_data['description']}")
    
    print(f"\nTotal parameters: {total_params}")

def main():
    """Main function."""
    if len(sys.argv) < 2:
        print("Usage: python3 validate_json_input.py <json_file> | --schema-info")
        sys.exit(1)
    
    try:
        schema = load_schema()
    except Exception as e:
        print(f"Error loading schema: {e}")
        sys.exit(1)
    
    if sys.argv[1] == "--schema-info":
        show_schema_info(schema)
        return
    
    json_file = sys.argv[1]
    
    if not Path(json_file).exists():
        print(f"Error: File not found: {json_file}")
        sys.exit(1)
    
    try:
        is_valid, message = validate_json_file(json_file, schema)
        
        if is_valid:
            print(f"✓ {message}")
        else:
            print(f"✗ {message}")
            sys.exit(1)
            
    except Exception as e:
        print(f"Error during validation: {e}")
        sys.exit(1)

if __name__ == '__main__':
    main()