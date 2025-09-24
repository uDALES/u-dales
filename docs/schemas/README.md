# u-DALES Input Schemas

This directory contains JSON schemas that define the structure and validation rules for u-DALES input files.

## Files

### `udales_input_schema.json`
Complete JSON Schema (draft-07) for u-DALES namelist parameters. This schema:

- **Defines all input parameters** organized by namelist sections (RUN, DOMAIN, PHYSICS, etc.)
- **Specifies data types and constraints** (integer ranges, boolean values, string formats)
- **Provides default values** for all parameters
- **Includes descriptions** explaining what each parameter controls
- **Enables validation** of JSON input files before running simulations

## Usage

### For Users
- **JSON Input Validation**: Use this schema to validate your `config.json` files
- **Parameter Reference**: Check parameter names, types, defaults, and descriptions
- **IDE Support**: Many editors can use this schema for autocompletion and validation

### For Developers
- **Testing**: The test suite in `tests/json_test/` uses this schema for comprehensive testing
- **Documentation**: Parameter definitions are automatically extracted from this schema
- **Tools**: External tools can use this schema to generate input files or UIs

### For Tools Integration
```python
import json
import jsonschema

# Load schema
with open('docs/schemas/udales_input_schema.json', 'r') as f:
    schema = json.load(f)

# Validate input file
with open('config.json', 'r') as f:
    config = json.load(f)

jsonschema.validate(config, schema)
```

## Schema Structure

The schema follows the same structure as u-DALES namelists:

```json
{
  "RUN": {
    "iexpnr": 1,
    "runtime": 100.0,
    "lwarmstart": false
  },
  "DOMAIN": {
    "itot": 64,
    "jtot": 64, 
    "ktot": 64
  },
  "PHYSICS": {
    "...": "..."
  }
}
```

## Maintenance

This schema should be updated whenever:
- New parameters are added to u-DALES
- Parameter defaults change
- Parameter constraints are modified
- New namelist sections are introduced

The schema serves as the authoritative definition of the u-DALES input interface.