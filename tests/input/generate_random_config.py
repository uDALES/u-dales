#!/usr/bin/env python3
"""
Generate a random u-DALES configuration file based on the schema.
This script reads the udales_input_schema.json and generates random values
for each parameter according to their type and constraints.
"""

import json
import random
import os
import string
from pathlib import Path

def generate_random_string(length=10):
    """Generate a random string of given length."""
    return ''.join(random.choices(string.ascii_letters + string.digits, k=length))

def generate_random_value(param_def):
    """Generate a random value based on parameter definition."""
    param_type = param_def.get("type")
    
    if param_type == "integer":
        min_val = param_def.get("minimum", 1)
        max_val = param_def.get("maximum", 1000)
        
        # Handle special cases for enum values
        if "enum" in param_def:
            return random.choice(param_def["enum"])
        
        # Generate random integer in range
        return random.randint(min_val, max_val)
    
    elif param_type == "number":
        min_val = param_def.get("minimum", 0.0)
        max_val = param_def.get("maximum", 1000.0)
        
        # If no maximum specified, use a reasonable upper bound
        if "maximum" not in param_def:
            if min_val >= 0:
                max_val = max(100.0, min_val * 10)
            else:
                max_val = 100.0
        
        # Generate random float in range
        return round(random.uniform(min_val, max_val), 6)
    
    elif param_type == "boolean":
        return random.choice([True, False])
    
    elif param_type == "string":
        # Generate a reasonable random string
        if "stl_file" in param_def.get("description", "").lower():
            return f"geometry_{generate_random_string(6)}.stl"
        elif "file" in param_def.get("description", "").lower():
            return f"file_{generate_random_string(6)}.dat"
        else:
            return generate_random_string(8)
    
    else:
        # Fallback for unknown types
        return None

def generate_random_config(schema_path):
    """Generate a random configuration based on the schema."""
    
    # Load the schema
    with open(schema_path, 'r') as f:
        schema = json.load(f)
    
    config = {}
    
    # Process each section in the schema
    for section_name, section_def in schema.get("properties", {}).items():
        if section_def.get("type") == "object":
            config[section_name] = {}
            
            # Process each parameter in the section
            for param_name, param_def in section_def.get("properties", {}).items():
                random_value = generate_random_value(param_def)
                if random_value is not None:
                    config[section_name][param_name] = random_value
    
    # Ensure runmode is set to 1001 for testing
    if "RUN" in config:
        config["RUN"]["runmode"] = 1001
    
    # Ensure required parameters are present and make sense
    if "RUN" in config and "DOMAIN" in config:
        # Set fixed domain dimensions
        config["DOMAIN"]["itot"] = 16
        config["DOMAIN"]["jtot"] = 16  
        config["DOMAIN"]["ktot"] = 32
        
        # Set domain size - make sure it's reasonable and positive
        config["DOMAIN"]["xlen"] = 100.0  # 100 meters
        config["DOMAIN"]["ylen"] = 100.0  # 100 meters

        # Set nprocx and nprocy 
        proc_choices = [1, 1]
        nprocx = random.choice(proc_choices)
        # choose a different value for nprocy; fallback to 1 if none
        remaining = [p for p in proc_choices if p != nprocx]
        nprocy = random.choice(remaining) if remaining else 1

        config["RUN"]["nprocx"] = nprocx
        config["RUN"]["nprocy"] = nprocy
    
    return config

def add_random_params_to_config(base_config_path, schema_path, output_path):
    """Create a new config by adding random parameters to each section of an existing config."""
    
    # Load the base configuration
    with open(base_config_path, 'r') as f:
        base_config = json.load(f)
    
    # Load the schema to get available parameters
    with open(schema_path, 'r') as f:
        schema = json.load(f)
    
    modified_config = json.loads(json.dumps(base_config))  # Deep copy
    
    # For each section in the base config
    for section_name, section_data in base_config.items():
        if isinstance(section_data, dict):
            # Get schema definition for this section
            section_schema = schema.get("properties", {}).get(section_name, {})
            section_properties = section_schema.get("properties", {})
            
            # Find parameters that exist in schema but not in base config
            available_params = set(section_properties.keys())
            existing_params = set(section_data.keys())
            unused_params = available_params - existing_params
            
            # Add one random parameter if available
            if unused_params:
                random_param = random.choice(list(unused_params))
                param_def = section_properties[random_param]
                random_value = generate_random_value(param_def)
                
                if random_value is not None:
                    modified_config[section_name][random_param] = random_value
                    print(f"Added to {section_name}: {random_param} = {random_value}")
    
    return modified_config

def main():
    """Main function to generate random config files."""
    # Paths (computed relative to this script)
    # script is at <repo>/tests/input/generate_random_config.py -> repo root is parents[2]
    repo_root = Path(__file__).resolve().parents[2]
    schema_path = str(repo_root / "docs" / "schemas" / "udales_input_schema.json")
    # Write output to the same directory as this script (tests/input)
    script_dir = str(Path(__file__).resolve().parent)
    output_dir = script_dir
    output_path = os.path.join(output_dir, "parameters.random")

    print("Generating a fully-random u-DALES configuration...")
    print(f"Schema: {schema_path}")

    try:
        # Ensure output directory exists
        os.makedirs(output_dir, exist_ok=True)

        # Generate fully random configuration
        config = generate_random_config(schema_path)

        with open(output_path, 'w') as f:
            json.dump(config, f, indent=2)

        print(f"Created: {output_path} (sections: {len(config)})")
    except Exception as e:
        print(f"Error generating configuration file: {e}")
        return 1

    return 0

if __name__ == "__main__":
    exit(main())