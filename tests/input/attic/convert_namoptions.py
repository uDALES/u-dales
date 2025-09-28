#!/usr/bin/env python3
"""
Script to convert all namoption files from examples directory to JSON format
and add runmode=1001 parameter for testing.
"""

import os
import json
import subprocess
import shutil
from pathlib import Path

def convert_namoption_to_json(namoption_path, output_dir):
    """Convert a single namoption file to JSON with runmode=1001 added."""
    
    # Get the experiment number from the filename
    filename = Path(namoption_path).name
    if filename.startswith('namoptions.'):
        exp_num = filename.split('.')[-1]
    else:
        exp_num = 'unknown'
    
    print(f"Converting {filename} (experiment {exp_num})...")
    
    # Create temporary directory for conversion
    temp_dir = Path("/tmp") / f"convert_{exp_num}"
    temp_dir.mkdir(exist_ok=True)
    
    try:
        # Copy namoption file to temp directory
        temp_namoption = temp_dir / filename
        shutil.copy2(namoption_path, temp_namoption)
        
        # Run ud_nam2json
        ud_nam2json_path = "/rds/general/user/mvr/home/udales/u-dales/bin/ud_nam2json"
        result = subprocess.run(
            [ud_nam2json_path, str(temp_namoption)],
            cwd=temp_dir,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True
        )
        
        if result.returncode != 0:
            print(f"Error converting {filename}: {result.stderr}")
            return None
            
        # Read the generated JSON file
        json_file = temp_dir / f"parameters.{exp_num}"
        if not json_file.exists():
            print(f"Error: Expected output file {json_file} not found")
            return None
            
        with open(json_file, 'r') as f:
            config_data = json.load(f)
        
        # Add runmode=1001 to the RUN section
        if 'RUN' not in config_data:
            config_data['RUN'] = {}
        config_data['RUN']['runmode'] = 1001
        
        # Write to output directory
        output_file = Path(output_dir) / f"config_{exp_num}.json"
        with open(output_file, 'w') as f:
            json.dump(config_data, f, indent=2)
        
        print(f"Successfully created {output_file}")
        return output_file
        
    finally:
        # Clean up temporary directory
        shutil.rmtree(temp_dir, ignore_errors=True)

def main():
    """Main function to convert all namoption files."""
    
    examples_dir = Path("/rds/general/user/mvr/home/udales/u-dales/examples")
    output_dir = Path("/rds/general/user/mvr/home/udales/u-dales/tests/json_test")
    
    # Find all namoption files
    namoption_files = []
    for subdir in examples_dir.iterdir():
        if subdir.is_dir():
            for namoption_file in subdir.glob("namoptions.*"):
                # Skip the old file in 102
                if "old" not in namoption_file.name:
                    namoption_files.append(namoption_file)
    
    namoption_files.sort()
    
    print(f"Found {len(namoption_files)} namoption files to convert:")
    for nf in namoption_files:
        print(f"  {nf}")
    
    print(f"\nConverting to JSON with runmode=1001...")
    
    successful_conversions = []
    for namoption_file in namoption_files:
        result = convert_namoption_to_json(namoption_file, output_dir)
        if result:
            successful_conversions.append(result)
    
    print(f"\nConversion complete! Created {len(successful_conversions)} JSON config files:")
    for config_file in successful_conversions:
        print(f"  {config_file}")

if __name__ == "__main__":
    main()