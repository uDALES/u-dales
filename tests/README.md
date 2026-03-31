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

- `cases/`: shared committed case fixtures used by multiple repo-level tests
- `integration/`: end-to-end checks of interacting components
- `regression/`: case comparisons and branch-to-branch output checks
- `system/`: whole-code validation cases, including resource-heavy runs
- `unit/`: isolated tests when routines can be tested independently

## Integration

`tests/integration` contains executable-driven end-to-end checks for this branch:

- `directshortwave/`: Python-driven integration tests for direct shortwave on committed cases `100` and `525`
- `ibm_sparse_input/`: MPI validation for `read_sparse_ijk()` using `runmode = 1004`
- `tree_sparse_compare/`: MPI validation for sparse tree forcing using `runmode = 1005`

`tests/cases/` holds shared committed fixtures used by these tests. At present:

- `101/`: IBM sparse-input case used by `integration/ibm_sparse_input/`
- `100/`: no-tree direct shortwave reference case used by `integration/directshortwave/`
- `525/`: flat-terrain tree case used by `integration/directshortwave/` and `integration/tree_sparse_compare/`

To run the direct shortwave reference test:

```bash
source .venv/bin/activate
python tests/integration/directshortwave/test_directshortwave.py
```

To run the sparse IBM input test:

```bash
cd tests/integration/ibm_sparse_input
./run_test.sh
```

To run the sparse tree comparison test:

```bash
cd tests/integration/tree_sparse_compare
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
- If a test is driven by Python but relies on committed case fixtures under `tests/` or compares multiple implementations, prefer `tests/integration/`.
- Do not put exploratory or plotting-heavy development scripts in either automated test tree; keep those in `tools/python/examples/` or a dedicated dev location.
