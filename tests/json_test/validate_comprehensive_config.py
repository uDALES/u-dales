#!/usr/bin/env python3
"""
Validate that all 224 parameters in the comprehensive test have non-default values.
"""

import json
from pathlib import Path

def validate_comprehensive_test():
    """Validate the comprehensive test configuration."""
    
    # Load schema and comprehensive config
    project_dir = Path(__file__).resolve().parents[2]
    schema_path = project_dir / "docs" / "schemas" / "udales_input_schema.json"
    config_path = Path(__file__).parent / "comprehensive_test_output" / "config_comprehensive.json"
    
    with open(schema_path, 'r') as f:
        schema = json.load(f)
    
    with open(config_path, 'r') as f:
        config = json.load(f)
    
    print("COMPREHENSIVE CONFIGURATION VALIDATION")
    print("=" * 60)
    
    # Extract all schema parameters with defaults
    schema_params = {}
    for section_name, section_data in schema.get('properties', {}).items():
        if isinstance(section_data, dict) and 'properties' in section_data:
            for param_name, param_data in section_data['properties'].items():
                if isinstance(param_data, dict):
                    schema_params[param_name] = {
                        'section': section_name,
                        'default': param_data.get('default'),
                        'type': param_data.get('type')
                    }
    
    # Check each parameter in the generated config
    total_params = 0
    non_default_count = 0
    same_as_default = []
    missing_from_config = []
    
    for param_name, param_info in schema_params.items():
        total_params += 1
        section = param_info['section']
        default_val = param_info['default']
        
        if section in config and param_name in config[section]:
            config_val = config[section][param_name]
            
            # Check if value differs from default
            if config_val != default_val:
                non_default_count += 1
            else:
                same_as_default.append((param_name, section, default_val))
        else:
            missing_from_config.append((param_name, section))
    
    # Report results
    print(f"Total schema parameters: {total_params}")
    print(f"Parameters in config: {total_params - len(missing_from_config)}")
    print(f"Non-default values: {non_default_count}")
    print(f"Same as default: {len(same_as_default)}")
    print(f"Missing from config: {len(missing_from_config)}")
    
    success_rate = (non_default_count / total_params) * 100
    print(f"Success rate: {success_rate:.1f}%")
    
    if same_as_default:
        print(f"\nParameters still using default values ({len(same_as_default)}):")
        for param, section, default in same_as_default[:10]:
            print(f"  {param} [{section}] = {default}")
        if len(same_as_default) > 10:
            print(f"  ... and {len(same_as_default) - 10} more")
    
    if missing_from_config:
        print(f"\nParameters missing from config ({len(missing_from_config)}):")
        for param, section in missing_from_config[:10]:
            print(f"  {param} [{section}]")
        if len(missing_from_config) > 10:
            print(f"  ... and {len(missing_from_config) - 10} more")
    
    # Show some examples of successful non-default assignments
    print(f"\nExamples of successful non-default values:")
    examples_shown = 0
    for section_name, section_params in config.items():
        for param_name, config_val in section_params.items():
            if param_name in schema_params:
                default_val = schema_params[param_name]['default']
                param_type = schema_params[param_name]['type']
                
                if config_val != default_val and examples_shown < 10:
                    print(f"  {param_name:20} = {config_val:15} (default: {default_val}, type: {param_type})")
                    examples_shown += 1
                
                if examples_shown >= 10:
                    break
        if examples_shown >= 10:
            break
    
    print(f"\n{'✓' if success_rate > 95 else '⚠'} Comprehensive test validation complete!")
    
    return {
        'total_params': total_params,
        'non_default_count': non_default_count,
        'success_rate': success_rate,
        'same_as_default': same_as_default,
        'missing_from_config': missing_from_config
    }

if __name__ == '__main__':
    validate_comprehensive_test()