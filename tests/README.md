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
- `input/`: namelist input interface checks. `test_roundtrip_input.py` runs
  `u-dales`, writes the interpreted namelist state from the last rank, and
  compares input and output namelists with `f90nml`; `test_sourcecode.py`
  performs a static consistency check between Fortran namelist definitions,
  MPI broadcasts, and the schema used for tooling/editor support
- `zgrid_metrics/`: vertical grid metric identities for a uniform and a
  stretched profile using `runmode = 1006`, including the rejection of a
  stretched grid with an equidistant-only Poisson configuration
- `mpi_operators/`: direct MPI operator validation for `runmode = 1005` on
  the Xie/Castro case `100` across `1x1`, `2x1`, `1x2`, and `2x2`
- `processor_boundaries/`: MPI decomposition parity checks on the Xie/Castro
  no-tree case `100` and the vegetation case `526`
- `python_preproc_against_matlab/`: preprocessing parity test between the
  MATLAB and Python entry points on no-tree case `100`
- `warmstart_facets/`: warm-start facet-state continuity regression for issue
  #269 on the energy-balance case `064`
- `stretched_grid/`: preprocessing and a short solver run of the energy-balance
  case `064` on an exponentially stretched vertical grid
- `sgs_statistics/`: SGS statistics validation for `runmode = 1007` on the
  stretched-grid case `300` across `1x1`, `2x1`, `1x2` and `2x2`
- `udbase_against_matlab/`: Python-vs-MATLAB parity checks on committed cases
- `udprep/`: preprocessing integration checks on committed cases and binaries

### MPI launcher and rank counts

The MPI runners share one launch helper, `tools/mpi_env.sh` for the shell
scripts and `tests/integration/mpi_env.py` for the Python tests. It resolves the
launcher (`MPIEXEC`, else Intel MPI's `mpiexec` next to `mpiifort`, else
`mpiexec` from `PATH`), passes `MPI_LAUNCH_EXTRA_ARGS` through unchanged, and
places Open MPI's session directory under `TMPDIR`.

uDALES runs one rank per CPU and never oversubscribes. Before every launch the
helper compares the requested ranks with the CPUs granted to the process and
fails with both numbers when the request cannot be met, for example

```
ERROR: requested 4 MPI ranks but only 2 CPUs are available (Slurm allocation).
```

The granted CPUs are, in this order:

- in a Slurm job: `SLURM_NTASKS x SLURM_CPUS_PER_TASK` (else
  `SLURM_CPUS_ON_NODE x SLURM_JOB_NUM_NODES`)
- in a PBS Pro job: `NCPUS x` the distinct hosts in `PBS_NODEFILE`, that is the
  allocated cores; the nodefile line count (one per `mpiprocs`) is only a fallback
- otherwise: the CPUs the operating system makes available to the process
  (the affinity mask on Linux, `hw.logicalcpu` on macOS)

The `supported` selection needs four. On nodes with hyperthreading ask the
scheduler for physical cores (Slurm: `--hint=nomultithread`): the launcher
binds one rank per core and refuses an allocation of four hardware threads on
two cores. Outside a scheduler Open MPI counts one slot per core; on a
hyperthreaded laptop with fewer cores than ranks it stops with "not enough
slots", and `MPI_LAUNCH_EXTRA_ARGS=--use-hwthread-cpus` lets it run one rank
per hardware thread instead (GitHub Actions does the same on its 4-vCPU runner
through the equivalent MCA variables in `.github/workflows/ci.yml`). The
hosted macOS CI runner has three CPUs, so the
`supported-macos` selection runs the MPI tests on at most two ranks
(`ZGRID_MODES`, `UDALES_TEST_NPROCX/NPROCY`) and leaves the fixed 2x2 ibm sparse
input test to Linux and HPC.

### Running the tests on your platform

Every platform needs the same three things: a built executable, the Python
environment, and (for the MPI suites) an allocation with as many CPUs as ranks.

- **Executable.** Build with `tools/build_executable.sh <system> release`
  (or `debug`). `tests/run_tests.py --build-type Release|Debug` passes the
  matching `build/<release|debug>/u-dales` to every MPI suite; when running a
  suite on its own set `UDALES_BUILD=<path to u-dales>` (the shell runners
  default to `build/debug/u-dales`, the Python runners to `build/release/u-dales`).
  Tests that preprocess a case also need the preprocessing tools
  (`tools/build_preprocessing.sh <system> preprocessing_tools`, and `VIEW3D_EXE`
  when View3D is not on `PATH`).
- **Python.** Activate `tools/python/.venv` (see `tools/python/README_VENV.md`);
  `tests/run_tests.py` uses whatever interpreter is active.
- **Runtime modules.** On module-based clusters the solver runs inside a shell
  that loads `UDALES_RUNTIME_MODULES` (space-separated list, the same modules the
  executable was built with). Set it in the job script.
- **Launcher.** `MPIEXEC` when the launcher is not `mpiexec` (Cray systems:
  `MPIEXEC=srun`); `MPI_LAUNCH_EXTRA_ARGS` for site-specific flags.

Platform notes:

- **Laptop (Linux or macOS).** No scheduler: the granted CPUs are what the OS
  reports. Four-rank tests need four cores (see the Open MPI slot note above).
  On macOS keep `TMPDIR` short; the helper sets it for Open MPI.
- **Slurm clusters (ECMWF Atos, ARCHER2).** Run the selection inside a batch or
  interactive job with `--ntasks=4 --cpus-per-task=1 --hint=nomultithread`
  (plus the account and QoS of the site). On ARCHER2 set `MPIEXEC=srun`.
- **PBS Pro clusters (Imperial CX3, HX1).** Request
  `-l select=1:ncpus=4:mpiprocs=4` so that both the allocation and the
  launcher see four slots; the helper reads `NCPUS` and `PBS_NODEFILE`.
- **GitHub Actions.** `.github/workflows/ci.yml` runs `supported` on Linux
  (4 vCPUs, hardware threads counted as CPUs) and `supported-macos` on macOS
  (3 CPUs). Nothing is oversubscribed there either.

A typical cluster job then runs

```bash
export UDALES_BUILD=$PWD/build/release/u-dales
export UDALES_RUNTIME_MODULES="<modules used for the build>"
source tools/python/.venv/bin/activate
python tests/run_tests.py supported --branch-a origin/master --build-type Release
```

The `david_tests` regression harness in the `supported` selection reports PASS
and skips when the working tree has uncommitted changes; commit first to run it.

`tests/cases/` holds shared committed fixtures used by these tests. At present:

- `101/`: IBM sparse-input case used by `integration/ibm_sparse_input/`
- `100/`: no-tree direct shortwave reference case used by `integration/directshortwave/`
- `100/`: also used by `integration/processor_boundaries/` for the Xie/Castro
  no-tree decomposition check
- `100/`: also used by `integration/mpi_operators/` for the direct MPI
  operator check
- `100/`: also used by `integration/zgrid_metrics/` for the vertical grid
  metric checks (uniform profile and a test-local stretched profile)
- `100/`: also used by `integration/python_preproc_against_matlab/`
- `525/`: flat-terrain tree case used by `integration/directshortwave/`
- `526/`: reduced tree case used by `regression/new_vegetation_module_against_v2.2/`
- `526/`: also used by `integration/processor_boundaries/` for the vegetation
  decomposition check
- `064/`: energy-balance case used by `integration/warmstart_facets/` for the
  warm-start facet continuity regression
- `064/`: also used by `integration/stretched_grid/`, which preprocesses and
  runs it again on a stretched vertical grid
- `300/`: small 16x16x32 case with a **stretched** vertical grid and
  `libm = .false.`, used by `integration/sgs_statistics/`. The non-uniform `zf`
  is the point of the fixture: on a uniform grid a `dzf`/`dzh` confusion in the
  SGS stencils is invisible. It carries only `namoptions.300`, `prof.inp.300`
  and `lscale.inp.300` — no IBM or facet input is needed.

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

To run the direct MPI operator test:

```bash
cd tests/integration/mpi_operators
./run_test.sh
```

To run the SGS statistics test (experimental, not part of the merge gate):

```bash
cd tests/integration/sgs_statistics
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
