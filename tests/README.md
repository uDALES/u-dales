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
- `gpu-smoke`: Debug build, CUDA device self-tests, four tiny parity cases, the
  surface-energy-balance parity case (the one that checks the facet flux handover)
  and two stretched-vertical-grid parity cases (one-equation SGS, kappa advection),
  which the uniform-grid cases cannot distinguish from an equidistant formula
- `gpu-nightly`: Release build and representative IBM, vegetation, and SEB parity
- `gpu-mpi`: Debug two-rank parity on two GPUs, the two-rank Poisson checks and
  the warm-start facet-state chains on the Debug GPU build with the device
  self-tests on
- `gpu-full`: the full single- and multi-GPU matrix, the four-rank Poisson checks
  and the warm-start facet-state chains on the Release GPU build

For Python-library-only validation, use the dedicated dispatcher selection:

```bash
python tests/run_tests.py python-library
```

This runs Python-library unit tests plus Python-driven library integration/
reference suites (including directshortwave, udprep integration, udbase
parity, Python preprocessing parity, the GPU harness unit tests, and the build
wrapper and runtime module-loading contracts). It skips solver/MPI and
regression suites.

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
- `build_tools/`: compiler-free contract tests for the CPU/GPU build-directory
  selection in `tools/build_executable.sh`
- `runtime_environment/`: contract tests ensuring that MPI integration tests
  load environment modules only when explicitly requested
- `ibm_sparse_input/`: MPI validation for `read_sparse_ijk()` using `runmode = 1004`
- `zgrid_metrics/`: vertical grid metric identities for a uniform and a
  stretched profile using `runmode = 1020`, including the rejection of a
  stretched grid with an equidistant-only Poisson configuration
- `mpi_operators/`: direct MPI operator validation for `runmode = 1005` on
  the Xie/Castro case `100` across `1x1`, `2x1`, `1x2`, and `2x2`
- `poisson/`: Poisson solver validation for `runmode = 1019` on case `103`:
  both solvers, four decompositions, CPU or GPU build, exact-solution check on
  every rank, plus `run_solver_equivalence.py`, which compares `ipoiss=3` on CPU
  and GPU builds against `ipoiss=0` on the CPU build within a tolerance
- `processor_boundaries/`: MPI decomposition parity checks on the Xie/Castro
  no-tree case `100` and the vegetation case `526`
- `python_preproc_against_matlab/`: preprocessing parity test between the
  MATLAB and Python entry points on no-tree case `100`
- `warmstart_facets/`: warm-start facet-state continuity regression for issue
  #269 on the energy-balance case `064`; with `UDALES_GPU_BUILD` it also
  warm-starts the GPU build from CPU and GPU checkpoints and compares the
  facet temperatures across builds
- `stretched_grid/`: preprocessing and a short solver run of the energy-balance
  case `064` on an exponentially stretched vertical grid
- `sgs_statistics/`: SGS statistics validation for `runmode = 1021` on the
  stretched-grid case `300` across `1x1`, `2x1`, `1x2` and `2x2`
- `lchem/`: the NO/NO2/O3 chemistry step of `lchem` must be carried into the
  next time step (issue #427); uniform scalars on the energy-balance case
  `064` against the implicit reaction map, single rank
- `udbase_against_matlab/`: Python-vs-MATLAB parity checks on committed cases
- `udprep/`: preprocessing integration checks on committed cases and binaries
- `gpu/`: strict same-commit CPU/GPU output parity on deterministic committed
  fixtures; see `tests/integration/gpu/README.md` for the case matrix, runner
  requirements, and CI contract

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
(`ZGRID_MODES`, `UDALES_TEST_NPROCX/NPROCY`, `UDALES_MAX_RANKS=2` and
`--max-ranks 2` for the Poisson suites) and leaves the fixed 2x2 ibm sparse
input and mpi operators tests to Linux and HPC.

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

The MPI integration tests inherit the caller's active compiler/MPI environment
and do not infer a cluster from the presence of the `module` command. Cluster
automation may explicitly request a module stack by exporting a space-separated
`UDALES_RUNTIME_MODULES` value before invoking the test. If the required modules
are already loaded, leave `UDALES_RUNTIME_MODULES` unset.

To run the Debug GPU smoke suite after loading NVHPC and activating the
canonical Python environment:

```bash
python tests/run_tests.py gpu-smoke
```

This builds CPU and GPU executables in separate CMake trees, runs the Debug CUDA
device self-test suite, executes both solvers from the same commit,
and compares every variable in the required NetCDF outputs. GPU suites are not
included in `supported` or `all`; they are selected explicitly on GPU hardware.
The Debug GPU smoke suite is currently local-only and is not run by GitHub
Actions.

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

## Validation Before Committing, Pushing, or Merging

### On the Imperial HX1 cluster

The whole sequence below runs with one command on HX1:

```bash
tests/hpc_run_all_tests.sh          # fast: about 30 min of run time
tests/hpc_run_all_tests.sh full     # everything
tests/hpc_run_all_tests.sh gpu      # the four builds, then gpu-smoke and gpu-nightly only
tests/hpc_run_all_tests.sh full --login   # orchestrate from the login node
tests/hpc_run_all_tests.sh full --login --cpu-on-gpu   # CPU-only jobs on the a100 queue
```

It submits an orchestrator job, so it carries on after the terminal is
closed. With `--login` the orchestration runs on the login node under `nohup`
instead: the builds happen there and the wait on the test jobs costs nothing,
so no PBS walltime has to be guessed for it (an orchestrator job must outlive
the queue waits of every job it submits, and a long request on the small-CPU
queue can itself sit unplaced for days). Either way the test jobs ask for a
little more than they have taken on HX1, from 30 min for `python` to 3 h for
`all`. With `--cpu-on-gpu` the CPU-only jobs (`python`, `supported`, `fixtures`,
`all`) ask the a100 queue for one GPU each instead of the small-CPU queue; use it
when that queue cannot place anything (in November 2026 all of its nodes were
reserved by another group until 15 December, and jobs of any walltime sat with
"Placement set is too small"). The orchestrator builds all four executables from clean build directories,
submits the test jobs that can run at the same time, then runs the ones that
must run alone one after the other, and writes a final report, stating the
mode and every selection run, with the path of every log under
`logdir/run_all_tests_<timestamp>/` (the same place
`tools/ud_compare_multiple_inputs.py` logs to).

`fast`, the default, covers every build, `python-library`, `gpu-smoke`,
`gpu-mpi`, `gpu-nightly` (so both GPU build types, CPU/GPU parity and the
restart round trip) and `supported` in Debug. `full` adds `gpu-full`, which
needs four free A100s and can queue for hours for five minutes of work;
`gpu-2node`, the Poisson self-test and the mpi parity selection with each
pair of ranks on different nodes (2 nodes x 2 GPUs, `place=scatter`); the
CPU-only run of the GPU fixture matrix, submitted with `place=excl` because
its 128³ Debug case slows 3.5× when it shares a node, so it waits for an idle
node and then takes about 40 minutes; and `all` in Release. Queue waits are
outside the timings above and are reported per job.

`--dry-run` writes the job scripts without submitting anything. Commit first
if you want the supported regression to run: it checks out `master` and
`HEAD` in place and skips itself on a dirty tree.

The rest of this section is the environment that script sets up, for running
selections by hand. On HX1 everything comes from Lmod modules, and three
things about the cluster interact with the test harness. Set the environment
once per shell or job, then run the commands in the following subsections
unchanged.

**1. Run the compiled-solver streams on a compute node.** `python-library` is
fine on a login node. Everything that launches `mpiexec` or builds the solver
(`supported`, `all`, every `gpu-*` selection) belongs in a job. For the CPU
streams a batch job is simplest; for the GPU selections either submit a batch
job or take an interactive node:

```bash
# two GPUs: enough for gpu-smoke, gpu-nightly and gpu-mpi
qsub -I -l select=1:ncpus=2:mpiprocs=2:mem=256gb:ngpus=2:gpu_type=A100 -l walltime=01:00:00
# gpu-full needs four
qsub -I -l select=1:ncpus=8:mpiprocs=8:mem=128gb:ngpus=4:gpu_type=A100 -l walltime=04:00:00
```

The `a100` queue is one node per job with a fair-use ceiling of 12 GPUs per
user, so `gpu-full` (four GPUs) is the largest single request.

**2. Environment for the CPU streams** (`python-library`, `supported`, `all`):

```bash
cd u-dales
module load intel/2023a netCDF/4.9.2-iimpi-2023a netCDF-Fortran/4.6.1-iimpi-2023a \
            FFTW/3.3.10-intel-compilers-2023.1.0 CMake/3.26.3-GCCcore-12.3.0
# The regression harnesses need git; the GPU compute nodes have none in
# /usr/bin, so load it and point the harnesses that default to /usr/bin/git at it.
module load git/2.41.0-GCCcore-12.3.0-nodocs
export UDALES_GIT="$(command -v git)"

# The solver integration suites load their runtime themselves, through
# tests/integration/common/runtime_modules.sh, from this variable.
export UDALES_RUNTIME_MODULES="intel/2023a netCDF/4.9.2-iimpi-2023a netCDF-Fortran/4.6.1-iimpi-2023a FFTW/3.3.10-intel-compilers-2023.1.0"

# The supported regression builds master and HEAD with plain cmake. Hand it
# the compiler and netCDF the hx1 target uses, or cmake picks ifx from PATH.
export UDALES_CMAKE_ARGS="-DCMAKE_Fortran_COMPILER=$(command -v mpiifort) -DNETCDF_DIR=$EBROOTNETCDF -DNETCDF_FORTRAN_DIR=$EBROOTNETCDFMINFORTRAN"

# The experimental regressions build older revisions through their own
# tools/build_executable.sh, which predates the hx1 target; "common" is the
# one target every revision has, and it relies on the modules loaded above.
export UDALES_BUILD_SYSTEM=common
export MPIEXEC="$(command -v mpiexec)"

# HX1 sets BASH_ENV to Lmod's init script, so every child bash re-defines the
# real `module` command - which defeats the stub that the build-wrapper and
# runtime-module contract tests rely on, and makes them fail against real
# module names. Dropping it is safe: the exported function survives.
unset BASH_ENV

source tools/python/.venv/bin/activate
```

**3. Environment for the GPU selections** (`gpu-smoke`, `gpu-mpi`,
`gpu-nightly`, `gpu-full`), on top of the block above:

```bash
NVHPC=/gpfs/easybuild/prod/software/NVHPC/23.7-CUDA-12.2.0
IMPI=/gpfs/easybuild/prod/software/impi/2021.9.0-intel-compilers-2023.1.0/mpi/2021.9.0/bin

# tests/integration/gpu/build_test_binaries.sh: which build targets, which
# MPI wrappers. It needs absolute paths, and the two must differ.
export UDALES_CPU_SYSTEM=hx1
export UDALES_GPU_SYSTEM=gpuhx1
export UDALES_CPU_FORTRAN_COMPILER=$IMPI/mpiifort
# The GPU MPI is "stack 3": OpenMPI 4.1.6 built for NVHPC 23.7 with IPv6,
# PBS (tm) launch and the site UCX, by tools/GPU_dependency/install_dependencies.sh.
# The OpenMPIs bundled with NVHPC only work on a single node here (see
# .github/skills/udales-exec/references/clusters.md, HX1 GPU).
OMPI=$PWD/tools/GPU_dependency/stack3-OpenMPI-4.1.6-tm/build/openmpi/openmpi
export UDALES_GPU_FORTRAN_COMPILER=$OMPI/bin/mpif90

# run_gpu_tests.py: one launcher per binary. Appended, so Intel's tools keep
# precedence for everything else.
export UDALES_CPU_MPIEXEC=$IMPI/mpiexec
export UDALES_GPU_MPIEXEC=$OMPI/bin/mpiexec
export PATH="$PATH:$OMPI/bin"

# UCX finds its CUDA transports only through UCX_MODULE_DIR plus the EB_UCX_*
# module lists (the site's UCX-CUDA module sets them); without them UCX takes
# GPU buffers for host memory and multi-rank runs segfault on small messages.
export UCX_MODULE_DIR=/gpfs/easybuild/prod/software/UCX-CUDA/1.14.1-GCCcore-12.3.0-CUDA-12.1.1/ucx
export EB_UCX_uct_MODULES=":ib:rdmacm:cma:cuda"
export EB_UCX_ucm_MODULES=":cuda"
export EB_UCX_uct_cuda_MODULES=":gdrcopy"
export LD_LIBRARY_PATH="$NVHPC/Linux_x86_64/23.7/cuda/12.2/lib64:/gpfs/easybuild/prod/software/GDRCopy/2.3.1-GCCcore-12.3.0/lib:$LD_LIBRARY_PATH"

# Parity compares CPU and GPU output to 1e-6; keep host threads out of it.
export OMP_NUM_THREADS=1
```

Do not load the NVHPC module in this shell: `build_executable.sh gpuhx1`
loads it in its own subshell, and the GPU executable reaches its libraries
through its RPATH plus the GCCcore that the Intel set already provides.

**4. Keep the working tree clean while `supported` or `all` runs.** The
branch-comparison regression checks out `master` and `HEAD` in place and
restores the checkout afterwards. With uncommitted changes it skips itself
with a warning rather than failing, and an edit made while it is running can
block the checkout. Commit first, then run.

**5. The GPU and CPU selections share `build/cpu/debug` and
`build/cpu/release`**: `build_test_binaries.sh` rebuilds both. Do not run a
GPU selection and a CPU selection at the same time. Chaining the batch jobs
with `qsub -W depend=afterany:<cpu-job-id>` serialises them.

**6. What to expect.** Run this way on HX1 (Intel 2023a for the CPU
builds, NVHPC 23.7 for the GPU ones), `python-library`, `supported` in both
build types, the CPU-only fixture run and every GPU selection pass. What
remains is in the `experimental` stream:

- the two `experimental` regressions: they copy a `2decomp-fft` directory
  that this branch no longer carries as a submodule, and fail before building.
  A property of the branch, not of the cluster.

Worth knowing even though no test fails on it any more: with the Intel module
set loaded, Intel's TBB is on `LD_LIBRARY_PATH` and numba silently prefers it
over OpenMP as its threading layer. TBB assigns `prange` chunks dynamically, so
results that come from per-thread partial sums (the directshortwave budget
totals, for instance) move by a few ULPs between runs. Set
`NUMBA_THREADING_LAYER=omp` if you need the preprocessing bit-reproducible.

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
cmake -S . -B build/cpu/debug -DCMAKE_BUILD_TYPE=Debug
cmake --build build/cpu/debug -j 4

source tools/python/.venv/bin/activate
python tests/run_tests.py python-library
python tests/integration/gpu/run_gpu_tests.py full --validate-config
```

Run the supported stream against that Debug build to exercise the compiled
solver, including the in-solver runmode tests that the Python stream above does
not reach. Keep the virtual environment activated: the dispatcher launches each
suite with the interpreter that started it, so running this outside the
environment fails on missing dependencies rather than on the code being tested.

```bash
python tests/run_tests.py supported
```

The branch-comparison regression inside this selection compares committed
revisions, so it is only meaningful once the changes are committed. The
`Before pushing` step below repeats the selection with explicit `--branch-a`
and `--branch-b` for that reason.

For solver or GPU-related changes, also validate every committed GPU-suite
fixture with the CPU executable:

```bash
python tests/integration/gpu/run_gpu_tests.py full \
  --cpu-only \
  --cpu-executable build/cpu/debug/u-dales
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

cmake -S . -B build/cpu/release -DCMAKE_BUILD_TYPE=Release
cmake --build build/cpu/release -j 4

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
