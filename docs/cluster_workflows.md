# Cluster Workflows

This note records the cluster-side build, preprocessing, execution, and
analysis workflow that is already encoded in the repository scripts under
`tools/`. Most of it is written for the Imperial clusters (`icl` = CX3,
`hx1` = HX1); the [ECMWF HPC2020](#ecmwf-hpc2020) section covers the ECMWF
Atos system.

## Executable Build

Use the project wrapper instead of assembling a module stack by hand:

```bash
./tools/build_executable.sh icl debug
./tools/build_executable.sh icl release
```

Replace `icl` with `hx1`, `archer` or `ecmwf` on those machines. That script is
the source of truth for the solver build environment on each cluster. Each
branch loads the compiler, MPI, NetCDF and FFTW modules that the build expects,
plus CMake and (on the Imperial clusters) Git.

## Preprocessing Build

Set up the Python environment and build the full preprocessing toolchain with:

```bash
./tools/python/setup_venv.sh icl
```

This creates `tools/python/.venv`, installs the Python dependencies, and builds
View3D plus the f2py extension modules required by the Python preprocessing
route.

If you only need to rebuild the preprocessing binaries after the environment is
already available, use:

```bash
./tools/build_preprocessing.sh icl preprocessing_tools
```

Use `./tools/build_preprocessing.sh icl view3d` only when you deliberately want
to rebuild View3D without the f2py extension modules.

## Preprocessing Runs

`tools/write_inputs.sh` scans `namoptions.###` for `nompthreads` and uses it
for the preprocessing CPU request. If `nompthreads` is omitted, the
preprocessing default is `8`; if it appears more than once anywhere in the file,
the wrapper stops and asks for a single value. The wrapper exports the derived
value internally as `PREPROC_NCPU`; the default View3D configuration also uses
that value to choose the View3D OpenMP thread count unless
`VIEW3D_NUM_THREADS` is set explicitly. For preprocessing, `DA_TOOLSDIR`
defaults to the directory containing `write_inputs.sh` and `DA_EXPDIR` defaults
to the parent directory of the case directory unless these are set in
`config.sh` or the calling environment.

When submitting preprocessing to an Imperial HPC compute node with
`tools/write_inputs.sh <route> <case-directory> c`, the wrapper uses
`PREPROC_WALLTIME="24:00:00"` and `PREPROC_MEM="128gb"` unless these are set in
`config.sh` or the calling environment. These control the preprocessing PBS job
only and are separate from the solver job `WALLTIME` and `MEM` settings used by
`tools/hpc_execute.sh`. `PREPROC_MEM` must be written as a number followed by
lowercase `gb`, such as `128gb`; a unitless value such as `128` is rejected
before PBS submission. Unless `VIEW3D_MAX_DENSE_MATRIX_GIB` is set explicitly,
`tools/view3d_config.sh` derives the View3D dense-matrix guard from the
preprocessing memory request: requests above `16gb` leave 16 GiB for overhead,
while smaller requests use the requested GiB value.

## Batch Execution

Use `ud_run` rather than writing a new launcher command from scratch:

```bash
./bin/ud_run icl sim <case-directory>      # or hx1 on HX1
```

On `icl` and `hx1` this runs `tools/hpc_execute.sh` with the cluster set
explicitly (`UDALES_SYSTEM=cx3` or `hx1`), so the script does not have to detect
it. The case directory must provide `config.sh`, and the wrapper writes and
submits the PBS job script with the module stack and `mpirun` invocation that
the project expects.

Use the gather step to collect outputs after the run:

```bash
./bin/ud_run icl gather <case-directory>   # runs tools/hpc_gather.sh
```

On other machines use `archer` or `ecmwf` as the first argument; see
[Running uDALES](udales-simulation-setup.md) for the full list of machines and
the scripts `ud_run` calls.

## Python Environment

Create the project virtual environment with the setup script, which also
builds the preprocessing tools (View3D and the f2py extension modules):

```bash
./tools/python/setup_venv.sh icl
```

When activating the environment on the cluster, load the matching Python
module first so the runtime libraries are available:

```bash
module load Python/3.13.1-GCCcore-14.2.0
source tools/python/.venv/bin/activate
```

Use the same Python module for repo Python workflows on the cluster. In
particular, the `f2py`-based extensions in this repository are expected to be
built and run with the same Python runtime environment above rather than whichever `python3`
happens to be first on `PATH`.

## Interactive Analysis

On the login nodes, `HOME`, `$EPHEMERAL`, and login-node `$TMPDIR` may all be
backed by the shared RDS/GPFS filesystem. For large NetCDF/HDF5 reads, that can
lead to very slow or hanging bulk variable reads even when metadata access
works.

For interactive debugging of NetCDF outputs:

- copy the files to local `/tmp` first
- set `HDF5_USE_FILE_LOCKING=FALSE`
- then read them with `ncdump`, `ncks`, or Python `netCDF4`

This is especially relevant for regression comparison of `treedump.*.nc`
outputs.

## MPI Launcher Notes

Interactive MPI launches on the login nodes can behave differently from the PBS
job environment. In particular:

- the working `mpiexec` is not always the same launcher that appears first on
  `PATH` after loading modules
- batch-style output redirection patterns may fail interactively even when they
  work inside submitted jobs
- sandboxed agent sessions (Codex, Claude Code) can add another layer of
  difference: a launcher failure seen inside the sandbox may be a sandbox
  socket restriction rather than a real cluster-side problem

So for interactive debugging:

- prefer reproducing the environment from the repo wrappers
- if a login-node MPI launch fails inside an agent sandbox, retry it outside the sandbox
  before treating it as a solver or cluster configuration issue
- keep the launcher invocation minimal
- avoid changing MPI launcher behavior and output handling unless you have
  confirmed it works on the current node

## ECMWF HPC2020

The ECMWF Atos HPC2020 (hostnames such as `ac6-101`)
uses Slurm and ECMWF's own Lmod module tree under `/usr/local/apps`.

Build the solver on a login node:

```bash
./tools/build_executable.sh ecmwf release
```

This loads
`prgenv/intel intel/2021.4.0 intel-mpi/2021.4.0 netcdf4/4.10.0 fftw/3.3.10 cmake/4.2.4`
and compiles with `mpiifort`. The netCDF C and Fortran libraries share one
prefix (`/usr/local/apps/netcdf4/4.10.0/INTEL/2021.4`).

Run and gather:

```bash
./bin/ud_run ecmwf sim <case-directory>     # tools/ecmwf_execute.sh
./bin/ud_run ecmwf gather <case-directory>  # tools/ecmwf_gather.sh
```

- `ecmwf_execute.sh` copies the inputs and the executable to
  `$DA_WORKDIR/<exp>`, then submits `NNODE` nodes with `NCPU` ranks per node
  (at most 128, physical cores only) launched with `srun`. The job loads the
  runtime half of the build stack — `prgenv/intel`, `intel`, `intel-mpi`,
  `netcdf4` and `fftw`, but not CMake. It refuses to submit unless `DA_WORKDIR`
  is under `$SCRATCH`.
- `ecmwf_gather.sh` submits a one-task job that loads `prgenv/intel`,
  `intel/2021.4.0`, `netcdf4/4.10.0` and `nco/5.3.7` and runs
  `gather_outputs.sh`. It requires `MEM` (`16G` or `16gb`).
- `QOS` must be set in `config.sh` and applies to the simulation only; there is
  no default. Use `np` (exclusive nodes, 240 GB each) for real runs, `nf`
  (shared, 1 node, at most 128 cores / 128 GB) for small tests.
- The gather job is fixed to `nf` and ignores `QOS`: it is one NCO process, and
  on `np` it would be allocated — and billed — a whole 128-core node.
- `MEM` is requested by the gather job always, and by the simulation only when
  `QOS` is not `np` (shared nodes default to 8 GB). `16G` and `16gb` are both
  accepted and normalised for Slurm; a unit is required.
- The Slurm account is fixed by `ACCOUNT=` at the top of both scripts.

Filesystems: use `$SCRATCH` (50 TB) for run directories. It is purged 30 days
after last access, so move results to keep elsewhere (e.g. ECFS). `$HOME`
(10 GB), `$PERM` (500 GB) and `$HPCPERM` (1 TB) are too small for multi-TB
output, and none of the three is meant for parallel I/O. See the
[ECMWF HPC2020 user guide](https://confluence.ecmwf.int/spaces/UDOC/pages/240851027/HPC2020+User+Guide).

To check a job script by hand without queuing it, use
`sbatch --test-only <job-file>`. To dry-run the wrapper scripts themselves,
put a stub `sbatch` first on `PATH` that calls
`/usr/local/bin/sbatch --test-only "$@"` — ECMWF's `sbatch` is a site wrapper,
so the scripts must go through it.
