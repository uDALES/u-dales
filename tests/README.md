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

GPU coverage has dedicated selections because it requires NVHPC and CUDA
hardware:

- `gpu-build`: isolated CPU and GPU builds for the requested build type
- `gpu-smoke`: Debug build, CUDA initializer self-test, and two tiny parity cases
- `gpu-nightly`: Release build and representative IBM, vegetation, and SEB parity
- `gpu-mpi`: Debug two-rank parity on two GPUs
- `gpu-full`: the full single- and multi-GPU matrix

For Python-library-only validation, use the dedicated dispatcher selection:

```bash
python tests/run_tests.py python-library
```

This runs Python-library unit tests plus Python-driven library integration/
reference suites (including directshortwave, udprep integration, udbase
parity, and Python preprocessing parity). It skips solver/MPI and regression
suites.

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
python tests/run_tests.py gpu-smoke
python tests/run_tests.py gpu-nightly
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
- `gpu/`: strict same-commit CPU/GPU output parity on deterministic committed
  fixtures; see `tests/integration/gpu/README.md` for the case matrix, runner
  requirements, and CI contract

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
source tools/python/.venv/bin/activate
python tests/integration/directshortwave/test_directshortwave.py
```

To run the sparse IBM input test:

```bash
cd tests/integration/ibm_sparse_input
./run_test.sh
```

To run the Debug GPU smoke suite after loading NVHPC and activating the
canonical Python environment:

```bash
python tests/run_tests.py gpu-smoke
```

This builds CPU and GPU executables in separate CMake trees, runs the Debug CUDA
device self-test suite, executes both solvers from the same commit,
and compares every variable in the required NetCDF outputs. GPU suites are not
included in `supported` or `all`; they are selected explicitly on GPU hardware.

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

## Validation Before Committing, Pushing, or Merging

### Before committing

Inspect the working tree and check the patch for whitespace errors:

```bash
git status --short
git diff
git diff --check
```

Build the Debug solver and run the Python test stream and GPU matrix
configuration checks:

```bash
cmake -S . -B build/debug -DCMAKE_BUILD_TYPE=Debug
cmake --build build/debug -j 4

source tools/python/.venv/bin/activate
python tests/run_tests.py python-library
python tests/integration/gpu/run_gpu_tests.py full --validate-config
```

For solver or GPU-related changes, also validate every committed GPU-suite
fixture with the CPU executable:

```bash
python tests/integration/gpu/run_gpu_tests.py full \
  --cpu-only \
  --cpu-executable build/debug/u-dales
```

### Before pushing

After committing the changes, run the supported Debug stream and the complete
Release stream. Branch-comparison tests use committed revisions, so `HEAD`
must contain the changes being tested.

```bash
python tests/run_tests.py supported \
  --branch-a master \
  --branch-b HEAD \
  --build-type Debug

cmake -S . -B build/release -DCMAKE_BUILD_TYPE=Release
cmake --build build/release -j 4

python tests/run_tests.py all \
  --branch-a master \
  --branch-b HEAD \
  --build-type Release
```

For GPU-related changes, run at least the Debug smoke suite on a machine with
a visible NVIDIA GPU and the required NVHPC environment:

```bash
python tests/run_tests.py gpu-smoke
```

### Before merging branches

Start from a clean working tree, update the remote references, and inspect the
incoming change summary:

```bash
git status --short
git fetch origin
git diff --stat HEAD..origin/master
```

Perform the merge without creating its commit immediately:

```bash
git merge --no-commit --no-ff origin/master
```

After resolving conflicts, confirm that none remain and that the resolved
patch has no whitespace errors:

```bash
git diff --name-only --diff-filter=U
git diff --check
```

Repeat the **Before committing** checks before creating the merge commit, then
run the **Before pushing** checks before publishing it. When suitable GPU
hardware is available, finish with:

```bash
python tests/run_tests.py gpu-smoke
python tests/run_tests.py gpu-full
```

`gpu-full` requires four visible GPUs. With one GPU, use `gpu-nightly`; with
two GPUs, additionally use `gpu-mpi`.

### Coverage of this validation workflow

Running the complete sequence above covers:

- CPU Debug and Release compilation
- supported and experimental CPU tests
- Python tooling and preprocessing tests
- CPU regression against `master`
- GPU Debug and Release compilation
- Debug CUDA device self-tests
- all 20 CPU/GPU parity cases
- serial and multi-GPU MPI paths

If preprocessing source code changed, rebuild the preprocessing tools before
running their tests so that stale executables or extension modules cannot hide
a problem:

```bash
PREPROCESSING_PYTHON_EXECUTABLE="$(command -v python)" \
  bash tools/build_preprocessing.sh common preprocessing_tools

python tests/run_tests.py python-library
```

Branch-comparison tests operate on committed Git revisions. Consequently,
`HEAD` must contain the changes being validated before running the pre-push
regression commands with `--branch-b HEAD`.

Subject to the documented dependencies and four visible GPUs for `gpu-full`,
this sequence covers all tests registered in the curated uDALES test system.
It does not represent exhaustive testing of every possible namelist
configuration.
