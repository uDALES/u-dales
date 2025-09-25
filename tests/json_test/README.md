# u-DALES JSON Testing Framework

This directory contains the test suite for u-DALES JSON input functionality using the **runmode approach**.

## Overview

The testing framework leverages the integrated `runmode=1001` functionality in u-DALES to trigger comprehensive testing of JSON vs namelist input methods. When `runmode=1001` is set in configuration files, u-DALES executes the `tests_json` routine which:

1. **Process 0 reads configuration** - Uses actual production code (`readnamelists`/`readjsonconfig`)
2. **All processes participate in MPI broadcasts** - Tests the centralized `broadcast_config_parameters` routine  
3. **Last process outputs ALL parameters** - Creates files with complete parameter verification
4. **Program exits cleanly** - Indicates successful test completion

## Key Benefits

- ✅ **Direct testing of production code** - No mock implementations or separate test harnesses
- ✅ **Complete MPI broadcast verification** - Tests actual parallel communication patterns
- ✅ **Comprehensive parameter coverage** - All namelist sections automatically tested
- ✅ **Simple usage** - Just set `runmode=1001` in any configuration file
- ✅ **File-based output** - Clean namelist files written for inspection

## Test Files

### Main Test Framework
- `test_runmode_framework.py` - Core test framework using runmode=1001 approach
- `run_comprehensive_tests.py` - Comprehensive test runner with multiple scenarios

### Usage Examples

#### Quick Test
```bash
# Run basic unit tests
python3 test_runmode_framework.py

# Run comprehensive test suite  
python3 run_comprehensive_tests.py
```

#### Manual Test with u-DALES
```bash
# Create JSON config with runmode=1001
echo '{"RUN": {"iexpnr": 999, "runmode": 1001, "nprocx": 1, "nprocy": 1}, "DOMAIN": {"itot": 16, "jtot": 16, "ktot": 8}}' > config.json

# Run u-DALES (builds must include JSON support)
/path/to/u-dales

# Check generated namelist files
ls namelists_proc_*.txt
```

#### MPI Testing
```bash
# Create config and run with MPI
echo '{"RUN": {"runmode": 1001, "nprocx": 2, "nprocy": 2}, "DOMAIN": {"itot": 32, "jtot": 32, "ktot": 16}}' > config.json
mpirun -np 4 /path/to/u-dales

# Verify all processes generated files
ls namelists_proc_*.txt
```

## Expected Output

When using `runmode=1001`, u-DALES will:

1. **Read configuration** (JSON or namelist) on process 0
2. **Broadcast parameters** to all processes via `broadcast_config_parameters`
3. **Execute tests_json routine** which:
   - Tests namelist reading + MPI broadcast
   - Tests JSON reading + MPI broadcast  
   - Each process writes its parameter values to `namelists_proc_<id>.txt`
4. **Exit with success message**: "JSON tests completed successfully"

### Sample Output Files

Each process creates a file like `namelists_proc_3.txt` containing:
```fortran
&RUN
 IEXPNR=999,
 RUNMODE=1001,
 LWARMSTART=F,
 RUNTIME=60.0,
 /
&DOMAIN
 ITOT=16,
 JTOT=16,
 KTOT=8,
 XLEN=100.0,
 YLEN=100.0,
 /
```

## Requirements

- **u-DALES built with JSON support** - CMake automatically includes json-fortran
- **MPI available** (for multi-process testing) - Single process also works
- **Python 3** (for automated test scripts)

## Configuration Methods

### JSON Configuration
```json
{
  "RUN": {
    "iexpnr": 999,
    "runmode": 1001,
    "nprocx": 1,
    "nprocy": 1
  },
  "DOMAIN": {
    "itot": 16,
    "jtot": 16,
    "ktot": 8
  }
}
```

### Namelist Configuration
```fortran
&RUN
iexpnr = 999
runmode = 1001
nprocx = 1
nprocy = 1
/

&DOMAIN
itot = 16
jtot = 16
ktot = 8
/
```

## Troubleshooting

### Test Not Running
- **Check runmode**: Verify `runmode = 1001` is set in configuration
- **Check executable**: Ensure u-DALES is built with JSON support
- **Check files**: Ensure `config.json` or `namoptions` exists in working directory

### MPI Issues
- **Check MPI installation**: Verify `mpirun` is available
- **Check process count**: Ensure `-np` matches or exceeds `nprocx * nprocy`
- **Check permissions**: Ensure executable permissions on u-DALES binary

### No Output Files
- **Check success message**: Look for "JSON tests completed successfully"
- **Check return code**: Should be 0 for successful completion
- **Check directory**: Files written to current working directory

## Integration with CI/CD

The runmode approach is ideal for automated testing:

```bash
#!/bin/bash
# Simple CI test script
cd /path/to/test/dir
echo '{"RUN": {"runmode": 1001}, "DOMAIN": {"itot": 16, "jtot": 16, "ktot": 8}}' > config.json
timeout 60 /path/to/u-dales
if [ $? -eq 0 ]; then
    echo "✅ JSON tests passed"
else
    echo "❌ JSON tests failed"
    exit 1
fi
```

This approach provides comprehensive, reliable testing of the u-DALES JSON input system using the production codebase itself.