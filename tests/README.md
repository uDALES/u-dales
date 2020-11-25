# Tests

The following integration tests are used to verify that changes in the code do not alter the simulation results. These tests and automatically run on CI at every commit but should be run locally before pushing changes to the code.

Tests are run directly using configuration files from the `examples` folder but with a reduced simulation time to cut CPU costs. This is done by patching namelists with patch files included in the  `tests/patches` folder. Tests are run automatically for all patch files included in the `tests/patches` folder is detected and can be easily expanded by adding a new patch files corresponding to the additional example case to run. Currently we run tests for example case `001`, `102`, `201`, `501` and `502` 

## How to run tests

Make sure you have the required Python libraries installed on your system (see [DEVELOP.md](../DEVELOP.md)) and activate the `udales` conda environment. Then, to run the tests run the following commands:

```
python run_tests.py <branch_a> <branch_b> <build_type>
```

Where `<branch_a>` and `<branch_b>` are the two branches you want to compare and `<build_type>` is either `Debug` or `Release`.

E.g.
```
python run_tests.py master dmey/patch-1 Release
```

Tests outputs are saved under `tests/outputs` and currently include a boxplots of approximate errors for all four dimensional (3D space + time) quantities included in the uDALES output netCDF files.

