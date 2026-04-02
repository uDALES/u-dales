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

## Test Layers

As this tree grows, tests should be organized by both scope and execution
contract. Each test suite should have a clear answer to:

- what it validates
- what it requires to run
- whether it is expected to run in GitHub Actions

The intended layers are:

- `unit/`: isolated routines with minimal dependencies and fast runtime
- `integration/`: end-to-end or cross-component checks that use committed
  fixtures or compare multiple implementations
- `system/`: heavier whole-code validation, often solver-driven and potentially
  resource intensive
- `regression/`: branch-to-branch or reference-output comparisons

These layers describe scope, not whether a test is required in every context.
That is controlled separately by the CI/runtime contract below.

## CI And Runtime Contract

The repository should distinguish between test layers that are expected to run
in GitHub Actions and heavier suites that are not.

Current intended categories are:

- `supported`: stable merge-gating coverage expected to run in standard CI
- `experimental`: important automated coverage for developing functionality
  that is not yet part of the primary merge gate
- `heavy`: resource-heavy or cluster-oriented coverage not expected in routine
  GitHub Actions

Today, GitHub Actions runs the curated `supported` selection:

- `.github/workflows/ci.yml` installs system dependencies including `gfortran`,
  MPI, NetCDF, Graphviz, and FFTW
- the workflow then runs `tests/run_tests.py supported`
- this supported selection currently includes lightweight Python unit tests,
  the branch-comparison regression harness, and the IBM sparse-input solver
  integration test
- experimental and heavy tests are not part of the required merge gate

For a curated top-level entry point, use `tests/run_tests.py`. The group
membership lives in `tests/test_suites.yml`:

```bash
python tests/run_tests.py supported --branch-a master --branch-b HEAD --build-type Release
python tests/run_tests.py experimental
python tests/run_tests.py all --branch-a master --branch-b HEAD --build-type Release
```

This means a test being placed under `tests/integration/` does not by itself
mean it runs in the supported merge gate. The suite should document whether it is:

- `supported`
- `experimental`
- `heavy`

The `tests` directory itself is organized by test scope:

- `cases/`: shared committed case fixtures used by multiple repo-level tests
- `integration/`: end-to-end checks of interacting components
- `regression/`: case comparisons and branch-to-branch output checks
- `system/`: whole-code validation cases, including resource-heavy runs
- `unit/`: isolated tests when routines can be tested independently

## Integration

`tests/integration` contains end-to-end checks for this branch. These tests may
be Python-driven or executable-driven, but they should validate interactions
between multiple components rather than one isolated API.

- `directshortwave/`: Python-driven preprocessing integration tests for direct
  shortwave on committed cases `100` and `525`
- `ibm_sparse_input/`: MPI validation for `read_sparse_ijk()` using `runmode = 1004`

`tests/cases/` holds shared committed fixtures used by these tests. At present:

- `101/`: IBM sparse-input case used by `integration/ibm_sparse_input/`
- `100/`: no-tree direct shortwave reference case used by `integration/directshortwave/`
- `525/`: flat-terrain tree case used by `integration/directshortwave/`
- `526/`: reduced tree case used by `regression/new_vegetation_module_against_v2.2/`

To run the direct shortwave reference test:

```bash
source ../.venv/bin/activate
python tests/integration/directshortwave/test_directshortwave.py
```

To run the sparse IBM input test:

```bash
cd tests/integration/ibm_sparse_input
./run_test.sh
```

## Regression

`tests/regression` contains branch-comparison and case-based regression assets:

- `david_tests/`: the older branch-comparison regression harness and its helper assets
- `new_vegetation_module_against_v2.2/`: the `526` legacy vegetation regression against release `v2.2.0`

At present there are two regression paths:

- `david_tests/`: an older branch-comparison build harness used by the supported suite
- `new_vegetation_module_against_v2.2/`: a dedicated solver-output regression for the new vegetation module against the `v2.2.0` release

To run regression tests:

```bash
cd tests/regression/david_tests
python run_tests.py <branch_a> <branch_b> <build_type>
```

Or via the top-level dispatcher:

```bash
python tests/run_tests.py supported --branch-a <branch_a> --branch-b <branch_b> --build-type <build_type>
```

Where `<branch_a>` and `<branch_b>` are the two branches you want to compare and `<build_type>` is either `Debug` or `Release`.

Example:

```bash
cd tests/regression/david_tests
python run_tests.py master dmey/patch-1 Release
```

## Notes

- Put new tests in `tests/` when they validate whole-code behaviour, MPI runs, executable outputs, or shared case fixtures.
- Put new tests in `tools/python/tests/` when they validate Python modules directly.
- If a test is driven by Python but relies on committed case fixtures under `tests/` or compares multiple implementations, prefer `tests/integration/`.
- Do not assume every test under `tests/` belongs in the default supported GitHub Actions gate. Document whether it is `supported`, `experimental`, or `heavy`.
- For Python tooling and preprocessing paths that are still under development, prefer documenting them as `experimental` rather than making them part of the required merge gate too early.
- For tests with nontrivial prerequisites, make those requirements explicit in the suite README and in the test code where possible.
- Do not put exploratory or plotting-heavy development scripts in either automated test tree; keep those in `tools/python/examples/` or a dedicated dev location.
- See `tests/ROADMAP.md` for the current phased testing roadmap and project status against it.
