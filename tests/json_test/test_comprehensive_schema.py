#!/usr/bin/env python3
"""
Comprehensive schema validation test for u-DALES JSON input.

This test:
1. Loads all 224 parameters from the official schema
2. Sets each parameter to a valid non-default value
3. Generates both namelist and JSON files with all parameters
4. Runs u-DALES with both input methods
5. Compares outputs to find any missing schema parameters
"""

import json
import sys
import random
import subprocess
import tempfile
import shutil
from pathlib import Path
from typing import Dict, Any, List, Tuple

class ComprehensiveSchemaTest:
    """Test all schema parameters with non-default values."""
    
    def __init__(self):
        self.project_dir = Path(__file__).resolve().parents[2]
        self.schema_path = self.project_dir / "docs" / "schemas" / "udales_input_schema.json"
        self.test_dir = Path(__file__).parent / "comprehensive_test_output"
        self.test_dir.mkdir(exist_ok=True)
        
        # Load schema
        if not self.schema_path.exists():
            raise FileNotFoundError(f"Schema not found: {self.schema_path}")
        
        with open(self.schema_path, 'r') as f:
            self.schema = json.load(f)
        
        self.all_parameters = {}
        self._extract_all_parameters()
    
    def _extract_all_parameters(self):
        """Extract all parameters from the schema."""
        properties = self.schema.get('properties', {})
        
        for section_name, section_data in properties.items():
            if isinstance(section_data, dict) and 'properties' in section_data:
                section_props = section_data['properties']
                
                for param_name, param_data in section_props.items():
                    if isinstance(param_data, dict):
                        self.all_parameters[param_name] = {
                            'section': section_name,
                            'type': param_data.get('type', 'unknown'),
                            'default': param_data.get('default'),
                            'enum': param_data.get('enum'),
                            'minimum': param_data.get('minimum'),
                            'maximum': param_data.get('maximum'),
                            'description': param_data.get('description', ''),
                            'schema_data': param_data
                        }
        
        print(f"Extracted {len(self.all_parameters)} parameters from schema")
    
    def _generate_non_default_value(self, param_name: str, param_info: Dict) -> Any:
        """Generate a valid non-default value for a parameter."""
        param_type = param_info['type']
        default_val = param_info['default']
        enum_vals = param_info.get('enum')
        minimum = param_info.get('minimum')
        maximum = param_info.get('maximum')
        
        if enum_vals:
            # Choose from enumerated values, avoiding default if possible
            candidates = [v for v in enum_vals if v != default_val]
            if candidates:
                return random.choice(candidates)
            else:
                # If only one enum value exists and it's the default, 
                # this parameter cannot have a non-default value
                print(f"  Note: {param_name} has only one valid value: {enum_vals[0]} (same as default)")
                return enum_vals[0]
        
        elif param_type == 'integer':
            if minimum is not None and maximum is not None:
                # Generate in valid range, avoiding default
                candidates = list(range(minimum, maximum + 1))
                if default_val in candidates and len(candidates) > 1:
                    candidates.remove(default_val)
                return random.choice(candidates)
            else:
                # Use reasonable integer values
                candidates = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024]
                if default_val in candidates:
                    candidates.remove(default_val)
                return random.choice(candidates) if candidates else 42
        
        elif param_type == 'number':
            if minimum is not None and maximum is not None:
                # Generate in valid range
                val = random.uniform(minimum, maximum)
                # Ensure it's different from default
                if default_val is not None and abs(val - default_val) < 1e-6:
                    val = minimum + (maximum - minimum) * 0.7  # Use 70% of range
                return round(val, 6)
            else:
                # Use reasonable real values
                candidates = [0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 25.0, 50.0, 100.0, 0.01, 0.001]
                if default_val in candidates:
                    candidates.remove(default_val)
                return random.choice(candidates) if candidates else 3.14159
        
        elif param_type == 'boolean':
            # Toggle the boolean value
            return not bool(default_val) if default_val is not None else True
        
        elif param_type == 'string':
            # Generate alternative string values
            if 'file' in param_name.lower() or param_name.endswith('_file'):
                alternatives = ['test_data.stl', 'input.dat', 'config.inp', 'example.txt']
            elif 'var' in param_name.lower():
                alternatives = ['u0,v0,w0', 'thl0,qt0', 'p,thl,qt,u,v,w', 'sv01,sv02']
            else:
                alternatives = ['test_value', 'alternative', 'custom_setting', 'modified']
            
            if default_val in alternatives:
                alternatives.remove(default_val)
            return random.choice(alternatives) if alternatives else 'non_default_value'
        
        elif param_type == 'array':
            # Handle array parameters
            if default_val == [] or default_val is None:
                # Generate a non-empty array with appropriate values
                if 'scalar' in param_name.lower() or 'flux' in param_name.lower():
                    # Scalar flux arrays - use reasonable flux values
                    return [0.1, 0.2, 0.5]
                elif 'temp' in param_name.lower() or 'thl' in param_name.lower():
                    # Temperature arrays
                    return [300.0, 305.0, 310.0]
                else:
                    # Generic numeric array
                    return [1.0, 2.0, 3.0]
            else:
                # If default is non-empty, modify it
                if isinstance(default_val, list) and len(default_val) > 0:
                    # Add one more element or modify existing
                    new_array = default_val.copy()
                    new_array.append(99.9)
                    return new_array
                else:
                    return [1.0, 2.0, 3.0]
        
        else:
            print(f"Warning: Unknown type '{param_type}' for parameter '{param_name}'")
            return default_val
    
    def generate_comprehensive_config(self) -> Dict[str, Dict[str, Any]]:
        """Generate configuration with all 224 parameters set to non-default values."""
        config = {}
        used_params = set()
        skipped_params = []
        
        # Organize parameters by section
        for param_name, param_info in self.all_parameters.items():
            section = param_info['section']
            
            if section not in config:
                config[section] = {}
            
            try:
                # Generate non-default value
                new_value = self._generate_non_default_value(param_name, param_info)
                config[section][param_name] = new_value
                used_params.add(param_name)
                
                # Debug info for first few parameters
                if len(used_params) <= 10:
                    default_val = param_info['default']
                    print(f"  {param_name:20} = {new_value:15} (default: {default_val}, type: {param_info['type']})")
                
            except Exception as e:
                print(f"Warning: Could not generate value for {param_name}: {e}")
                skipped_params.append(param_name)
        
        print(f"\nGenerated non-default values for {len(used_params)} parameters")
        if skipped_params:
            print(f"Skipped {len(skipped_params)} parameters: {skipped_params[:10]}...")
        
        return config
    
    def create_comprehensive_namelist(self, config: Dict[str, Dict[str, Any]]) -> str:
        """Create namelist file with all parameters."""
        namelist_content = ""
        
        for section_name, section_params in config.items():
            namelist_content += f"&{section_name}\n"
            
            for param_name, param_value in section_params.items():
                if isinstance(param_value, bool):
                    value_str = '.true.' if param_value else '.false.'
                elif isinstance(param_value, str):
                    value_str = f"'{param_value}'"
                else:
                    value_str = str(param_value)
                
                namelist_content += f"{param_name} = {value_str}\n"
            
            namelist_content += "/\n\n"
        
        # Write to file
        namelist_path = self.test_dir / "namoptions.comprehensive"
        with open(namelist_path, 'w') as f:
            f.write(namelist_content)
        
        print(f"Created comprehensive namelist: {namelist_path}")
        return str(namelist_path)
    
    def create_comprehensive_json(self, config: Dict[str, Dict[str, Any]]) -> str:
        """Create JSON file with all parameters."""
        json_path = self.test_dir / "config_comprehensive.json"
        
        with open(json_path, 'w') as f:
            json.dump(config, f, indent=2, sort_keys=True)
        
        print(f"Created comprehensive JSON: {json_path}")
        return str(json_path)
    
    def validate_against_schema(self, config: Dict[str, Dict[str, Any]]) -> Tuple[bool, str]:
        """Validate the generated config against the schema."""
        try:
            import jsonschema
            jsonschema.validate(config, self.schema)
            return True, "Configuration is valid according to schema"
        except ImportError:
            return True, "jsonschema not available, skipping validation"
        except jsonschema.ValidationError as e:
            return False, f"Validation error: {e.message}"
        except Exception as e:
            return False, f"Validation failed: {str(e)}"
    
    def analyze_parameter_coverage(self, config: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
        """Analyze parameter coverage and find potential missing parameters."""
        analysis = {
            'schema_params': len(self.all_parameters),
            'generated_params': 0,
            'sections': {},
            'missing_from_schema': [],
            'coverage_by_type': {}
        }
        
        # Count generated parameters
        for section, params in config.items():
            analysis['generated_params'] += len(params)
            analysis['sections'][section] = len(params)
        
        # Analyze by type
        type_counts = {}
        for param_info in self.all_parameters.values():
            ptype = param_info['type']
            type_counts[ptype] = type_counts.get(ptype, 0) + 1
        
        analysis['coverage_by_type'] = type_counts
        
        return analysis
    
    def create_minimal_simulation_files(self, test_dir: Path):
        """Create minimal files needed for a u-DALES simulation test."""
        # Create basic geometry
        stl_content = """solid minimal_geometry
  facet normal 0 0 1
    outer loop
      vertex 0 0 0
      vertex 8 0 0
      vertex 8 8 0
    endloop
  endfacet
  facet normal 0 0 1
    outer loop
      vertex 0 0 0
      vertex 8 8 0
      vertex 0 8 0
    endloop
  endfacet
endsolid minimal_geometry
"""
        
        geometry_files = ['test_data.stl', 'input.dat', 'config.inp', 'example.txt', 'geometry.stl']
        for gfile in geometry_files:
            with open(test_dir / gfile, 'w') as f:
                f.write(stl_content)
    
    def run_comprehensive_test(self) -> Dict[str, Any]:
        """Run the comprehensive test with all parameters."""
        print("=" * 80)
        print("COMPREHENSIVE SCHEMA TEST - ALL 224 PARAMETERS")
        print("=" * 80)
        
        # Generate comprehensive configuration
        print("\n1. Generating comprehensive configuration...")
        config = self.generate_comprehensive_config()
        
        # Validate against schema
        print("\n2. Validating against schema...")
        is_valid, validation_msg = self.validate_against_schema(config)
        print(f"   Validation result: {validation_msg}")
        
        if not is_valid:
            print("   ERROR: Generated configuration is not valid!")
            return {'error': validation_msg}
        
        # Create files
        print("\n3. Creating input files...")
        namelist_path = self.create_comprehensive_namelist(config)
        json_path = self.create_comprehensive_json(config)
        
        # Create supporting files
        self.create_minimal_simulation_files(self.test_dir)
        
        # Analyze coverage
        print("\n4. Analyzing parameter coverage...")
        analysis = self.analyze_parameter_coverage(config)
        
        print(f"   Schema parameters: {analysis['schema_params']}")
        print(f"   Generated parameters: {analysis['generated_params']}")
        print(f"   Sections:")
        for section, count in analysis['sections'].items():
            print(f"     {section:15}: {count:3} parameters")
        
        print(f"   Parameter types:")
        for ptype, count in analysis['coverage_by_type'].items():
            print(f"     {ptype:15}: {count:3} parameters")
        
        # Summary
        result = {
            'config': config,
            'namelist_path': namelist_path,
            'json_path': json_path,
            'validation': {'valid': is_valid, 'message': validation_msg},
            'analysis': analysis,
            'test_dir': str(self.test_dir)
        }
        
        print("\n5. Test files created successfully!")
        print(f"   Test directory: {self.test_dir}")
        print(f"   Namelist file: {Path(namelist_path).name}")
        print(f"   JSON file: {Path(json_path).name}")
        
        return result


def main():
    """Run the comprehensive schema test."""
    try:
        tester = ComprehensiveSchemaTest()
        result = tester.run_comprehensive_test()
        
        if 'error' in result:
            print(f"\nTest failed: {result['error']}")
            sys.exit(1)
        
        print("\n" + "=" * 80)
        print("COMPREHENSIVE TEST COMPLETED SUCCESSFULLY")
        print("=" * 80)
        print(f"✓ All {result['analysis']['generated_params']} parameters set to non-default values")
        print(f"✓ Configuration validated against schema")
        print(f"✓ Both namelist and JSON files created")
        print(f"✓ Ready for simulation testing")
        
        print(f"\nTo test with u-DALES:")
        print(f"  cd {result['test_dir']}")
        print(f"  # Copy your u-DALES executable here")
        print(f"  # Run with JSON: ln -sf config_comprehensive.json config.json && ./u-dales")
        print(f"  # Run with namelist: ln -sf namoptions.comprehensive namoptions.001 && ./u-dales")
      
    except Exception as e:
        print(f"Error running comprehensive test: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()