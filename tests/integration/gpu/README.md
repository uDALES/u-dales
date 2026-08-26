# GPU Porting Tests

This suite checks the GPU port against the CPU implementation built from the
same source revision. CPU output is the parity reference; it is deliberately
not read from `master`, because an intentional physics change on the current
branch must be exercised by both execution paths.

The suite is initially classified as `experimental`/`heavy`. Promote the smoke
test to a required check only after its tolerances and runner availability have
been stable for a suitable period.

## What is implemented

`case_matrix.json` defines 31 deterministic cases and five selections:

| Selection | Cases | Intended use |
| --- | --- | --- |
| `smoke` | four 8 x 8 x 8 cases | Debug development and trusted pull requests |
| `scalar-sources` | two serial and two two-rank cases | scalar-source parity and global positioning |
| `nightly` | 26 serial cases | single-GPU scheduled regression |
| `mpi` | dry and scalar-source two-rank X/Y cases | manually dispatched two-GPU check |
| `full` | nightly plus X, Y, and 2 x 2 MPI | manual four-GPU validation |

The current coverage is:

- dry momentum, second-order advection, Vreman closure, and neutral bottom wall
- both `ipoiss=0` and `ipoiss=3`
- temperature and moisture with second-order advection
- Kappa temperature and scalar advection, including the temperature tendency
  conversion kernels
- one-equation SGS TKE advection, production, diffusion, and integration
- fixed bottom temperature, moisture, and scalar flux kernels
- non-neutral flat-bottom Uno momentum/heat wall functions and their atomic
  heat-flux accumulation
- Coriolis and profile forcing, nonzero large-scale subsidence/radiative/
  moisture tendencies, shifted periodic-boundary forcing, `ifixuinf=1`, and
  all three gravity-damping modes
- X and Y profile inflow/outflow, pressure-top, no-slip top, and nonzero top
  temperature/moisture/scalar fluxes
- three-species NO/NO2/O3 chemistry
- point and finite-line passive-scalar sources, including global peak-location
  checks on the non-root partition of two-rank X and Y decompositions
- IBM wall functions, including the paths that update `iomomflux`, `iotflux`,
  and `obcTfluxA`
- the extended scalar halos involved in the former `sca1` mismatch
- vegetation forcing and tree output
- surface energy balance, sparse view factors, and periodic EB correction
- two-rank X and Y decompositions and a four-rank 2 x 2 decomposition
- the `fielddump` variables no other case requests: the wall stresses `tau_x`,
  `tau_y` and `tau_z`, the heat flux `thl_flux`, and the post-correction
  divergence `div`. Without this case those four transfers could be dropped
  entirely and every other case would still pass. `div` is there for the code
  path only - after the pressure correction it sits at round-off, so it cannot
  detect a stale velocity field
- a second fielddump case asking for `ty` and `tz` but not `tx` or `hf`. Those
  four transfers are conditional on fielddump naming them, and with only a case
  that asks for all four, a transfer keyed to the wrong code is invisible -
  every code is then either wanted everywhere or nowhere. The partial selection
  is what tells them apart
- the `statsdump` modes no other case enables: `lydump`, `lytdump`, `lxydump`
  and `lmintdump` in one case, `lkslicedump`, `lislicedump` and `ljslicedump`
  in another. Both carry moisture and the averages case carries three scalars,
  because several `ytdump` statistics are written whether or not the physics
  that assigns them is enabled

Every Debug GPU smoke/MPI run must also emit a CUDA device-suite marker from
every MPI rank. The opt-in `src/tests_cuda.f90` suite checks the extended-halo
initializer, the Kappa limiter, directional Kappa scalar advection, first-order
upwind scalar advection, both Kappa temperature tendency-copy directions, and
the driver-inlet boundary kernel. It also checks all eight X/Y periodic device
routines using coordinate-dependent interior values, poisoned halo cells, and
independently calculated expected values. These checks cover momentum and the
optional one-equation SGS fields, temperature and its Kappa-only extended
field, moisture, and multiple extended-grid scalars without relying only on
downstream NetCDF output differences.
The Kappa/upwind and temperature-copy checks run when the selected case has
allocated their required extended-grid arrays; the smoke selection deliberately
contains cases that do so. The runner enables the suite with
`UDALES_RUN_CUDA_SELFTEST=1`; ordinary Debug GPU simulations do not run it.

The `ported_routines` contract in `case_matrix.json` is checked against the
current Fortran source. At present it maps all 48 production CUDA global/device
routines (including device functions) and all 45 subroutines containing
executable OpenACC regions to one or more parity cases or to the Debug device
suite. Matrix validation fails when a new ported routine has no declared test,
a deleted routine leaves a stale entry, a mapped test does not exist, or a
namelist no longer selects the option asserted by its case. It also rejects the
explicitly unsupported GPU configurations `ipoiss=2` and `lpurif=.true.`.

## Restart round trip

`run_restart_roundtrip.py` is the one test here that is not a single run per
build. It runs `namoptions.103.restart-roundtrip` cold on the CPU and on the
GPU, then warm on each build from the restart that build's own cold run wrote
at `trestart`, and compares four ways: cold CPU against cold GPU (restart files
included), warm CPU against warm GPU, and cold against warm on each build at
the end time - the last field-dump record and the restart written there. On
the CPU that comparison is bit for bit: the case has no atomics, so a
difference of one bit is a lost or mis-restored piece of state. On the GPU it
is to 1e-11, because a warm start forms the slab averages with the host
reduction while the continuing run formed them on the device, and the two sum
in a different order. A restart that lost state would lose it identically on
both builds, which is why warm-versus-warm parity alone is not enough. No
restart fixture is committed; the files live and die in the run directory.

This test found two defects on its first runs: a `do j = j, je` typo that
left `thv0` uninitialised on a warm start, and a missing halo exchange on the
warm-start path that left the scalars' outer halo unfilled for the first step.

```bash
python tests/integration/gpu/run_restart_roundtrip.py --require-debug-selftest
```

It takes the same executable, work-root and summary options as
`run_gpu_tests.py`, and reads only the tolerances from the matrix.

## Comparison contract

`compare_outputs.py` compares every variable in every required output, not a
hand-selected list. Outputs are NetCDF, except the restart files: the
`initd` and `inits` prefixes select a reader for the Fortran unformatted
format `modsave::writerestartfiles` produces, which names each record after
the array it holds and compares it under the same rules - integer records
exactly, real ones under the tolerance a NetCDF variable of that name would
get. It fails when:

- a required output is missing
- the CPU and GPU file sets differ
- dimensions, variable names, shapes, or fill-value masks differ
- an unmasked value is NaN or infinite
- a numerical value violates `abs(gpu-cpu) <= atol + rtol*abs(cpu)`

The default tolerances and any variable-specific overrides are versioned in
`case_matrix.json`. Failure reports include the file, variable, maximum error,
and array index. The runner can write JSON, JUnit XML, logs, and the small
NetCDF outputs to an artifact directory.

Cases may additionally declare semantic `output_assertions`. The scalar-source
MPI cases require a nonzero final `sca1` peak at the configured global cell
centre and in the expected non-root rank-local fielddump file. This catches
shared CPU and GPU failures such as silently skipping the source or treating
its global coordinates as rank-local coordinates, which ordinary
same-decomposition parity would miss.

The hardware-independent comparator and matrix-validation tests are included
in the normal `python-library` stream, so ordinary CPU GitHub Actions still
protect the GPU testing infrastructure.

## Build isolation

CPU and GPU executables must use separate CMake caches. The build wrapper now
accepts these optional variables:

- `UDALES_BUILD_DIR`
- `UDALES_BUILD_JOBS`
- `UDALES_FORTRAN_COMPILER`

Without an override, `tools/build_executable.sh` writes CPU builds to
`build/cpu/<build-type>` and GPU builds to `build/gpu/<build-type>`.

`build_test_binaries.sh` uses them to produce:

```text
build/cpu/debug/u-dales
build/gpu/debug/u-dales
build/cpu/release/u-dales
build/gpu/release/u-dales
```

The default CPU compiler is `/usr/bin/mpif90`. Override
`UDALES_CPU_FORTRAN_COMPILER` when the CPU reference needs another compatible
MPI wrapper. The GPU compiler defaults to the `mpif90` in the active NVHPC
environment and can be overridden with `UDALES_GPU_FORTRAN_COMPILER`.

## Local execution

Use the repository's canonical Python environment and load NVHPC before
calling the top-level dispatcher:

```bash
source tools/python/.venv/bin/activate
module purge
module use /opt/nvidia/hpc_sdk/modulefiles
module load nvhpc/24.11
export UDALES_GPU_FORTRAN_COMPILER="$(command -v mpif90)"
export UDALES_CPU_MPIEXEC=/usr/bin/mpiexec
export UDALES_GPU_MPIEXEC="$(command -v mpiexec)"

python tests/run_tests.py gpu-smoke
```

The broader suites are:

```bash
python tests/run_tests.py gpu-nightly
python tests/run_tests.py gpu-mpi
python tests/run_tests.py gpu-full
```

For focused scalar-source development with existing executables, run:

```bash
python tests/integration/gpu/run_gpu_tests.py scalar-sources \
  --cpu-executable build/cpu/debug/u-dales \
  --gpu-executable build/gpu/debug/u-dales
```

That selection includes two-rank X and Y positioning cases and therefore uses
two GPUs with the default rank-to-device binding. For a correctness-only run
on one GPU, allow both MPI ranks to share it:

```bash
UDALES_GPU_BIND=0 python tests/integration/gpu/run_gpu_tests.py scalar-sources \
  --cpu-executable build/cpu/debug/u-dales \
  --gpu-executable build/gpu/debug/u-dales
```

`gpu-mpi` requires at least two visible GPUs and `gpu-full` requires four,
because the existing `tools/bind.sh` assigns one local MPI rank per device. Set
`UDALES_GPU_BIND=0` for correctness testing with shared-device ranks or when
device placement is managed externally.

Do not set `UDALES_RUN_CUDA_SELFTEST` manually when using this runner. The
`--require-debug-selftest` option owns the flag and verifies one pass marker
for every expected MPI rank.

The runner accepts distinct launcher environments for the executables:
`UDALES_CPU_MPIEXEC`, `UDALES_GPU_MPIEXEC`, `UDALES_CPU_MPI_ARGS`, and
`UDALES_GPU_MPI_ARGS`. They fall back to the older common `MPIEXEC` and
`MPI_LAUNCH_EXTRA_ARGS` variables.

To reuse executables without rebuilding, call the runner directly:

```bash
python tests/integration/gpu/run_gpu_tests.py smoke \
  --cpu-executable build/cpu/debug/u-dales \
  --gpu-executable build/gpu/debug/u-dales \
  --require-debug-selftest \
  --artifacts-dir /tmp/udales-gpu-artifacts
```

Validate the matrix without NetCDF dependencies, MPI, or a GPU:

```bash
python tests/integration/gpu/run_gpu_tests.py smoke --validate-config
```

Validate that the fixtures themselves run and produce all expected finite
outputs when only a CPU executable is available:

```bash
python tests/integration/gpu/run_gpu_tests.py smoke \
  --cpu-only --cpu-executable build/cpu/debug/u-dales
```

The runner stages fresh copies under `/tmp` by default. Set `--work-root` to a
node-local scratch directory when `/tmp` is unsuitable. It sets
`HDF5_USE_FILE_LOCKING=FALSE` unless the caller has already provided a value.

## GitHub Actions

The Debug smoke suite is local-only and is not run by GitHub Actions. Run it
before pushing GPU-related changes with:

```bash
python tests/run_tests.py gpu-smoke
```

`.github/workflows/gpu-ci.yml` contains two hardware jobs:

- Release nightly on the schedule or manual dispatch
- two-rank MPI only by manual dispatch on a runner labelled `multi-gpu`

The expected runner labels are:

```text
self-hosted, linux, x64, gpu, nvhpc
self-hosted, linux, x64, multi-gpu, nvhpc
```

The runner image/account must provide NVHPC 24.11, the GPU NetCDF/FFTW paths
used by `tools/build_executable.sh`, `/usr/bin/mpif90` for the CPU build, and
`tools/python/.venv` with NumPy and netCDF4. Set `UDALES_GPU_PYTHON` in the
runner service environment if the canonical venv is elsewhere.

The standard Ubuntu job validates the committed smoke matrix without running
the solver. Configure protected branches and runner-group access before
restoring `gpu-smoke` as a GitHub-hosted check.

Logs and comparison artifacts are retained for 14 days. Scheduled workflows
only run from GitHub's default branch, so the schedule becomes active there
after this workflow is merged to that branch.

## Adding a case

1. Prefer an existing fixture under `tests/cases/`; keep small test-specific
   namelists and overrides in this directory.
2. Use fixed timesteps and deterministic initialization for direct parity.
3. Give the case a focused `coverage` list and explicit required outputs.
4. Start in `nightly`; move it into `smoke` only if it is fast and reliable.
5. Add variable-specific tolerance only with a documented numerical reason.
6. Run the matrix validator and the hardware-independent unit tests.

Generated outputs and binaries must not be committed.
