# JSON Test Output Directory

This directory contains all outputs from the JSON input functionality tests for u-DALES.

## Directory Structure

- `unit_test_output/` - Output files from `test_json_input.py` unit tests
- `integration_test_output/` - Output files from `test_json_integration.py` integration tests  
- `demo_test_output/` - Output files from `demo_json_testing.py` demonstration

## Contents

Each test run creates subdirectories containing:
- Generated namelist files (`namoptions.*`)
- Generated JSON files (`config.json`)
- Test geometry files (`*.stl`)
- Simulation output files (when integration tests run successfully)
- Log files and error outputs

## Purpose

This organized output structure allows for:
- Easy inspection of generated test files
- Debugging of test failures
- Comparison between namelist and JSON input methods
- Verification that both input methods produce identical configurations

## Cleanup

Test outputs are preserved for inspection. To clean up:
```bash
rm -rf /path/to/udales/tests/json_test/*_output/
```

The `.gitignore` file ensures that temporary test outputs are not committed to the repository while preserving the directory structure and documentation.