# Cluster Workflows

This note records the cluster-side build, preprocessing, execution, and
analysis workflow that is already encoded in the repository scripts under
`tools/`.

## Executable Build

Use the project wrapper instead of assembling a module stack by hand:

```bash
./tools/build_executable.sh icl debug
./tools/build_executable.sh icl release
```

That script is the source of truth for the solver build environment on the
cluster. It loads the compiler, MPI, NetCDF, FFTW, CMake, and Git modules that
the solver build expects.

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

Use the execution wrapper rather than writing a new launcher command from
scratch:

```bash
./tools/hpc_execute.sh <case-directory>
```

The case directory must provide `config.sh`, and the wrapper writes and submits
the PBS job script with the module stack and `mpiexec` invocation that the
project expects.

Use the gather wrapper to collect outputs after the run:

```bash
./tools/hpc_gather.sh <case-directory>
```

## Python Environment

Use the repo-local virtual environment created by `setup_venv.sh`:

```bash
tools/python/.venv
```

When activating that environment on the cluster, load the matching Python module
first so the runtime libraries are available:

```bash
module load Python/3.13.1-GCCcore-14.2.0
source tools/python/.venv/bin/activate
```

Use the same Python module for repo Python workflows on the cluster. In
particular, the `f2py`-based extensions in this repository are expected to be
built and run with the same Python runtime rather than whichever `python3`
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
- Codex sandboxed sessions can add another layer of difference: a launcher
  failure seen inside the sandbox may be a sandbox socket restriction rather
  than a real cluster-side problem

So for interactive debugging:

- prefer reproducing the environment from the repo wrappers
- if a login-node MPI launch fails inside Codex, retry it outside the sandbox
  before treating it as a solver or cluster configuration issue
- keep the launcher invocation minimal
- avoid changing MPI launcher behavior and output handling unless you have
  confirmed it works on the current node
