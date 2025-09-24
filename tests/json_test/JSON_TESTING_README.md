# JSON Input Testing for u-DALES

This directory contains a comprehensive test suite for the JSON input functionality implemented in u-DALES.

## Test Files

### Core Test Suite

1. **`test_json_input.py`** - Main unit test suite
   - Extracts parameter schema from Fortran source files
   - Generates random non-default values for testing
   - Creates equivalent namelist and JSON input files
   - Validates file structure and parameter consistency
   - Runs 6 comprehensive unit tests

2. **`test_json_integration.py`** - Integration tests
   - Runs actual u-DALES simulations with both input methods
   - Compares parameter extraction and execution results
   - Tests with multiple parameter combinations including random values
   - Validates that both methods produce identical outputs

3. **`demo_json_testing.py`** - Demonstration script
   - Shows the complete functionality of the test system
   - Demonstrates schema extraction capabilities
   - Shows random value generation and file creation
   - Provides a quick overview of all testing features

## Key Features

### Schema Extraction
- Automatically parses Fortran source files to extract parameter definitions
- Identifies parameter types (integer, real, logical, string)
- Extracts default values from variable declarations
- Supports 22+ well-known u-DALES parameters

### Intelligent Random Testing
- Generates random values that differ from defaults
- Respects parameter type constraints
- Uses realistic value ranges for each parameter type
- Excludes default values to ensure meaningful testing

### Comprehensive Validation
- Creates both namelist and JSON formats with identical parameters
- Validates JSON structure and parseability  
- Ensures parameter consistency between formats
- Supports complex parameter combinations

## Usage

### Run Unit Tests
```bash
python3 test_json_input.py
```

### Run Integration Tests (requires compiled u-DALES)
```bash
python3 test_json_integration.py
```

### See Demonstration
```bash
python3 demo_json_testing.py
```

## Test Results

All tests pass successfully, demonstrating that:
- ✅ Schema extraction works correctly
- ✅ Random value generation excludes defaults
- ✅ Both input formats are created with correct structure
- ✅ Parameters are consistently represented in both formats
- ✅ JSON input method is equivalent to namelist method

## Example Output

The test suite successfully extracts parameters like:
- `nsv` (integer): Number of scalar variables
- `ih`, `jh`, `kh` (integer): Halo cells in each direction  
- `lwarmstart`, `lfielddump` (logical): Boolean configuration flags
- `runtime`, `dtmax` (real): Timing parameters
- `fname_options`, `stl_file` (string): File path parameters

Random non-default values are generated for all parameters, ensuring comprehensive testing coverage.

## Integration with u-DALES Build

These tests validate the JSON input functionality added to u-DALES through:
- Git submodule integration of json-fortran library
- CMake build system modifications with USE_JSON_INPUT flag
- Conditional compilation in `src/modstartup.f90`
- Backward compatibility with existing namelist method

The test suite confirms that the JSON input implementation works correctly and produces identical results to the traditional namelist approach.