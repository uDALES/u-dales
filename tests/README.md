# Tests

The following integration tests are used to verify that any change to the source code does not alter simulation results between the `master` and working branch. **Tests should be run locally before pushing changes** and are run automatically at every commit through CI.

Tests are run and can be modified in two different ways by:

- specifying a whole case: run using configuration files from the `tests/cases` folder. Tests are run automatically for all case folders included in the `tests/cases` folder and can be easily expanded by adding a new case folders. Currently we run tests for case `103`.

- patching an existing example: run using configuration files from the `examples` folder but with a reduced simulation time to cut CPU costs. This is done by patching namelists with patch files included in the  `tests/patches` folder. Tests are run automatically for all patch files included in the `tests/patches` folder is detected and can be easily expanded by adding a new patch files corresponding to the additional example case to run. Currently we run tests for example case `001`, `102`, `201`, `501` and `502`.

## How to run tests

Make sure you have the required Python libraries installed on your system (see [DEVELOP.md](../DEVELOP.md)) and activate the `udales` conda environment. Then, to run the tests run the following commands from the `tests` directory:

```
python run_tests.py <branch_a> <branch_b> <build_type>
```

Where `<branch_a>` and `<branch_b>` are the two branches you want to compare and `<build_type>` is either `Debug` or `Release`.

E.g.
```
python run_tests.py master dmey/patch-1 Release
```

Tests outputs are saved under `tests/outputs` and currently include a boxplots of approximate errors for all four dimensional (3D space + time) quantities included in the uDALES output netCDF files.

## JSON Input Tests

Additional tests for the JSON input functionality are located in the `json_test/` subdirectory:

- **`json_test/test_json_input.py`**: Unit tests for JSON input schema extraction and parameter validation
- **`json_test/test_json_integration.py`**: Integration tests comparing namelist vs JSON input methods  
- **`json_test/demo_json_testing.py`**: Demonstration script showing all JSON testing capabilities
- **`json_test/test_comprehensive_schema.py`**: Comprehensive test of all 224 schema parameters

To run the JSON input tests:

```bash
cd json_test/

# Run all JSON tests with a single command
python3 run_all_json_tests.py

# Or run individual tests:

# Unit tests using Python unittest framework
python3 -m unittest test_json_input.py
python3 -m unittest test_json_integration.py

# Or run tests directly
python3 test_json_input.py        # Unit tests for JSON schema extraction
python3 test_json_integration.py  # Integration tests comparing namelist vs JSON
python3 demo_json_testing.py      # Demonstration of all features
python3 test_comprehensive_schema.py  # Test all 224 parameters with non-default values

# Validate a JSON input file
python3 ../../docs/validate_json_input.py config.json

# Show schema information  
python3 ../../docs/validate_json_input.py --schema-info
```

All JSON test outputs are organized in the `json_test/` directory structure. See `json_test/README.md` for detailed documentation.

