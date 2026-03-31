# Tests

The top-level `tests` directory is for repo-level tests and shared case assets.
It is primarily intended for checks that exercise compiled uDALES behaviour,
MPI runs, case inputs, or branch-to-branch outputs.

It is not the only test location in the repo:

- `tests/`: repo-level regression, integration, and system tests
- `tools/python/tests/`: Python package tests for `udbase`, `udprep`, `udgeom`, and related utilities
- `tools/python/examples/` or a dedicated dev area: exploratory or visual scripts that are not stable automated tests

Use the location that matches the thing being validated, not just the language
used to run the test.

The `tests` directory itself is organized by test scope:

- `integration/`: end-to-end checks of interacting components
- `regression/`: case comparisons and branch-to-branch output checks
- `system/`: whole-code validation cases, including resource-heavy runs
- `unit/`: isolated tests when routines can be tested independently

## Integration

`tests/integration` contains executable-driven end-to-end checks for this branch:

- `sparse_ijk/`: MPI validation for `read_sparse_ijk()` using `runmode = 1004`
- `tree_input/`: vegetation/tree input case assets used by solver and Python tooling tests

To run the sparse IBM test:

```bash
cd tests/integration/sparse_ijk
./run_test.sh
```

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

- Put new tests in `tests/` when they validate whole-code behaviour, MPI runs, executable outputs, or shared case fixtures.
- Put new tests in `tools/python/tests/` when they validate Python modules directly.
- Do not put exploratory or plotting-heavy development scripts in either automated test tree; keep those in `tools/python/examples/` or a dedicated dev location.
