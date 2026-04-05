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
  solver-facing integration checks, preprocessing/tooling integration checks,
  and the older branch-comparison regression harness
- on macOS, CI currently runs the lighter `supported-macos` selection instead,
  which excludes the branch-comparison regression while `master` remains
  incompatible with the newer Homebrew CMake helper-project path
- experimental and heavy tests are not part of the required merge gate

For a curated top-level entry point, use `tests/run_tests.py`. The group
membership lives in `tests/test_suites.yml`:

```bash
python tests/run_tests.py supported --branch-a master --branch-b HEAD --build-type Release
python tests/run_tests.py experimental
python tests/run_tests.py all --branch-a master --branch-b HEAD --build-type Release
```

## Test Manifest Schema

`tests/test_suites.yml` is the source of truth for curated automated test
selections. The current schema is intentionally small and stable.

Top-level structure:

- `description`: short human-readable summary of the manifest
- `groups`: mapping from selection name to a group definition

Group definition:

- `suites`: optional list of concrete suite entries to run
- `includes`: optional list of other group names to include

Suite definition:

- `label`: human-readable display name shown by `tests/run_tests.py`
- `class`: support policy for the suite
  - `supported`: part of the required, stable merge-gating path
  - `experimental`: useful automated coverage that is not yet part of the
    primary merge gate
- `kind`: scope of validation
  - `unit`: isolated API or routine checks with minimal dependencies
  - `integration`: multi-component or real-case checks, often using committed
    fixtures or external tools
  - `reference`: compare against committed reference outputs or another
    implementation
  - `system`: heavier whole-code validation when used
  - `regression`: branch/version/reference comparison workflows
- `component`: subsystem or ownership area, for example `udbase`, `udprep`,
  `udgeom`, `ibm`, or `solver`
- `platform`: expected runtime target such as `linux`, `macos`, `hpc`, or `any`
- `cost`: rough runtime tier such as `fast`, `medium`, or `slow`
- `command`: ordered command list passed directly to `subprocess.run`
- `env_<NAME>`: optional per-suite environment variable override injected into
  the child process as `<NAME>`

Current conventions:

- `class` describes execution policy, not technical scope
- `kind` describes scope, not whether the suite is required in CI
- tests that depend on committed real case directories under `tests/cases/` or
  compare Python against another implementation should usually be
  `kind: integration` or `kind: reference` even if they are Python-driven
- `supported-macos` is a temporary compatibility group and should be treated as
  an exception rather than the long-term schema shape

Planned evolution:

- the schema is expected to keep `class` and `kind` as separate concepts
- the runner currently tolerates legacy `purpose` entries for compatibility, but
  new manifest edits should use `kind`

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

When a suite depends on a specific `runmode`, keep that explicit in the
test-local namelist template rather than relying on whatever the case fixture
happens to contain.

## Integration

`tests/integration` contains end-to-end checks for this branch. These tests may
be Python-driven or executable-driven, but they should validate interactions
between multiple components rather than one isolated API.

- `directshortwave/`: Python-driven preprocessing integration tests for direct
  shortwave on committed cases `100` and `525`
- `ibm_sparse_input/`: MPI validation for `read_sparse_ijk()` using `runmode = 1004`
- `mpi_operators/`: direct MPI operator validation for `runmode = 1005` on
  the Xie/Castro case `100` across `1x1`, `2x1`, `1x2`, and `2x2`
- `processor_boundaries/`: MPI decomposition parity checks on the Xie/Castro
  no-tree case `100` and the vegetation case `526`
- `python_preproc_against_matlab/`: preprocessing parity test between the
  MATLAB and Python entry points on no-tree case `100`
- `udbase_against_matlab/`: Python-vs-MATLAB parity checks on committed cases
- `udprep/`: preprocessing integration checks on committed cases and binaries

`tests/cases/` holds shared committed fixtures used by these tests. At present:

- `101/`: IBM sparse-input case used by `integration/ibm_sparse_input/`
- `100/`: no-tree direct shortwave reference case used by `integration/directshortwave/`
- `100/`: also used by `integration/processor_boundaries/` for the Xie/Castro
  no-tree decomposition check
- `100/`: also used by `integration/mpi_operators/` for the direct MPI
  operator check
- `100/`: also used by `integration/python_preproc_against_matlab/`
- `525/`: flat-terrain tree case used by `integration/directshortwave/`
- `526/`: reduced tree case used by `regression/new_vegetation_module_against_v2.2/`
- `526/`: also used by `integration/processor_boundaries/` for the vegetation
  decomposition check

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

To run the direct MPI operator test:

```bash
cd tests/integration/mpi_operators
./run_test.sh
```

## Regression

`tests/regression` contains branch-comparison and case-based regression assets:

- `david_tests/`: the older branch-comparison regression harness and its helper assets
- `new_vegetation_module_against_v2.2/`: the `526` legacy vegetation regression against release `v2.2.0`
- `mpi_averaging_regression/`: branch/commit regression for MPI-sensitive dumped fields and diagnostics

At present there are three regression paths:

- `david_tests/`: an older branch-comparison build harness used by the supported suite
- `new_vegetation_module_against_v2.2/`: a dedicated solver-output regression for the new vegetation module against the `v2.2.0` release
- `mpi_averaging_regression/`: a compact branch-comparison regression for decomposition-sensitive dumped fields on cases `100` and `526`

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
