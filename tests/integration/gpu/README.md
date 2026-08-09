# GPU Porting Tests

This suite checks the GPU port against the CPU implementation built from the
same source revision. CPU output is the parity reference; it is deliberately
not read from `master`, because an intentional physics change on the current
branch must be exercised by both execution paths.

The suite is initially classified as `experimental`/`heavy`. Promote the smoke
test to a required check only after its tolerances and runner availability have
been stable for a suitable period.

## What is implemented

`case_matrix.json` defines six deterministic cases and four selections:

| Selection | Cases | Intended use |
| --- | --- | --- |
| `smoke` | two 8 x 8 x 8 cases | Debug development and trusted pull requests |
| `nightly` | smoke plus IBM, vegetation, and SEB | single-GPU scheduled regression |
| `mpi` | two-rank X-decomposed case | manually dispatched two-GPU check |
| `full` | nightly plus MPI | manual validation on a multi-GPU machine |

The current coverage is:

- dry momentum, advection, subgrid closure, and neutral bottom-wall handling
- both `ipoiss=0` and `ipoiss=3`
- temperature, moisture, second-order thermodynamic advection, and Kappa
  scalar advection
- non-neutral bottom-wall heat and momentum fluxes
- IBM wall functions, including the paths that update `iomomflux`, `iotflux`,
  and `obcTfluxA`
- the extended scalar halos involved in the former `sca1` mismatch
- vegetation forcing and tree output
- surface energy balance, sparse view factors, and periodic EB correction
- a two-rank/multi-GPU X decomposition

Every Debug GPU smoke/MPI run must also emit a CUDA extended-halo initializer
self-test marker from every MPI rank. The Python runner opts into the dedicated
`src/tests_cuda.f90` check with `UDALES_RUN_CUDA_SELFTEST=1`; ordinary Debug GPU
simulations do not run it. A missing marker is a failure, preventing a Release
executable from accidentally being used for the Debug device check.

## Comparison contract

`compare_outputs.py` compares every variable in every required NetCDF output,
not a hand-selected list. It fails when:

- a required output is missing
- the CPU and GPU file sets differ
- dimensions, variable names, shapes, or fill-value masks differ
- an unmasked value is NaN or infinite
- a numerical value violates `abs(gpu-cpu) <= atol + rtol*abs(cpu)`

The default tolerances and any variable-specific overrides are versioned in
`case_matrix.json`. Failure reports include the file, variable, maximum error,
and array index. The runner can write JSON, JUnit XML, logs, and the small
NetCDF outputs to an artifact directory.

The hardware-independent comparator and matrix-validation tests are included
in the normal `python-library` stream, so ordinary CPU GitHub Actions still
protect the GPU testing infrastructure.

## Build isolation

CPU and GPU executables must use separate CMake caches. The build wrapper now
accepts these optional variables:

- `UDALES_BUILD_DIR`
- `UDALES_BUILD_JOBS`
- `UDALES_FORTRAN_COMPILER`

`build_test_binaries.sh` uses them to produce:

```text
build/common/debug/u-dales
build/gpu/debug/u-dales
build/common/release/u-dales
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

`gpu-mpi` and `gpu-full` require at least two visible GPUs because the existing
`tools/bind.sh` assigns one local MPI rank per device. Set
`UDALES_GPU_BIND=0` only when device placement is managed externally.

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
  --cpu-executable build/common/debug/u-dales \
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
  --cpu-only --cpu-executable build/common/debug/u-dales
```

The runner stages fresh copies under `/tmp` by default. Set `--work-root` to a
node-local scratch directory when `/tmp` is unsuitable. It sets
`HDF5_USE_FILE_LOCKING=FALSE` unless the caller has already provided a value.

## GitHub Actions

`.github/workflows/gpu-ci.yml` contains three hardware jobs:

- Debug smoke on pushes to `gpu` and trusted same-repository pull requests
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

The workflow intentionally refuses to execute a fork pull request on the
self-hosted GPU runner. The standard Ubuntu job still validates the committed
matrix for such pull requests. Configure protected branches and runner-group
access before making `gpu-smoke` a required check.

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
