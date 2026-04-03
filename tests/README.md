# Tests

The `tests` directory is organized by test scope:

- `integration/`: end-to-end checks of interacting components
- `regression/`: case comparisons and branch-to-branch output checks
- `system/`: whole-code validation cases, including resource-heavy runs
- `unit/`: isolated tests when routines can be tested independently

## Integration

`tests/integration/input` contains checks for the input interface:

- `test_roundtrip_input.py` runs a round-trip integration test. It executes `u-dales`, writes the interpreted namelist state from the last rank, compares input and output namelists with `f90nml`, and optionally checks schema defaults.
- `test_sourcecode.py` performs a static consistency check between Fortran namelist definitions, MPI broadcasts, and the schema used for tooling/editor support.

## Regression

`tests/regression` contains branch-comparison and case-based regression assets:

- `cases/`: complete test cases
- `patches/`: patches applied to example cases to create cheaper regression runs
- `scripts/`: helper scripts for building models and comparing outputs
- `run_tests.py`: regression test entry point

To run regression tests:

```bash
cd tests/regression
python run_tests.py <branch_a> <branch_b> <build_type>
```

Where `<branch_a>` and `<branch_b>` are the two branches you want to compare and `<build_type>` is either `Debug` or `Release`.

Example:

```bash
cd tests/regression
python run_tests.py master dmey/patch-1 Release
```

## Notes

- Large HPC-scale validation runs should live under `tests/system/`.
- Python unit tests should live alongside the Python code, for example under `tools/python/tests`.
- `tests_sparse_ijk` and `tests_tree_input` currently contain old log files rather than active tests, so they have been left untouched.
