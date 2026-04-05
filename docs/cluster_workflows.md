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

Build View3D and related preprocessing binaries with:

```bash
./tools/build_preprocessing.sh icl
```

This is the cluster-side source of truth for preprocessing build setup.

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

Use the shared virtual environment at:

```bash
/rds/general/user/mvr/home/udales/.venv
```

When activating that environment on the cluster, load the matching Python
module first so the runtime libraries are available:

```bash
module load Python/3.9.6-GCCcore-11.2.0
source /rds/general/user/mvr/home/udales/.venv/bin/activate
```

Use that Python 3.9 module for repo Python workflows on the cluster. In
particular, the `f2py`-based extensions in this repository are expected to be
built and run with the Python 3.9 environment above rather than whichever
`python3` happens to be first on `PATH`.

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
