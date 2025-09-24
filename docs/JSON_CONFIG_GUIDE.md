# u-DALES JSON Configuration Guide

u-DALES supports JSON input as a modern alternative to traditional Fortran namelists.

## Quick Start

1. **Create a configuration file**: `config.json`
2. **Use the same structure as namelists** but in JSON format
3. **Run u-DALES normally** - it will automatically detect and use the JSON file

## Example Configuration

```json
{
  "RUN": {
    "iexpnr": 1,
    "runtime": 100.0,
    "lwarmstart": false,
    "ladaptive": true
  },
  "DOMAIN": {
    "itot": 64,
    "jtot": 64,
    "ktot": 64,
    "xlen": 32.0,
    "ylen": 32.0
  },
  "OUTPUT": {
    "lfielddump": true,
    "tfielddump": 10.0,
    "fieldvars": "u0,v0,w0,thl0"
  }
}
```

## Validation

Validate your configuration before running:

```bash
python3 docs/validate_json_input.py config.json
```

## Schema Reference

Complete parameter definitions with types, defaults, and descriptions are available in:
- [`docs/schemas/udales_input_schema.json`](docs/schemas/udales_input_schema.json)

## Fallback Behavior

- If `config.json` exists, u-DALES uses JSON input
- If `config.json` doesn't exist, u-DALES falls back to traditional namelists
- Both methods produce identical results

## Benefits

✅ **Modern format**: Industry-standard JSON  
✅ **Validation**: Schema-based input validation  
✅ **Tool integration**: Easy integration with external tools  
✅ **Backward compatible**: Existing namelists continue to work

## Testing

Comprehensive test suite available in `tests/json_test/`:

```bash
cd tests/json_test
python3 test_json_input.py        # Unit tests
python3 test_json_integration.py  # Integration tests
python3 demo_json_testing.py      # Demonstration
```