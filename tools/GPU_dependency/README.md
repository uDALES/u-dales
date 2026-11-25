# GPU dependency stacks for uDALES on Imperial HX1

The GPU build (`tools/build_executable.sh gpuhx1`) needs HDF5, netCDF-C,
netCDF-Fortran and FFTW compiled with the NVIDIA HPC SDK against one specific
MPI, and the job must launch with that MPI's `mpirun`. Three such stacks have
been built; each subdirectory holds the install scripts and a README with the
run order. Only the scripts are kept here, not the sources or installs.

| Directory | MPI | Use |
|---|---|---|
| `stack1-NVHPC-OpenMPI-3.1.5` | OpenMPI 3.1.5 bundled with NVHPC 23.7 | single node; the original stack |
| `stack2-NVHPC-OpenMPI-4.1.5` | OpenMPI 4.1.5 bundled with NVHPC 23.7 (UCX, CUDA) | single node; ~20 % faster than stack 1 |
| `stack3-OpenMPI-4.1.6-tm` | OpenMPI 4.1.6 built here with IPv6, PBS launch and the site UCX-CUDA | single and multi-node; the current default |

The bundled OpenMPIs cannot run across HX1 nodes: the node networks are
IPv6-only and those builds lack IPv6 support, and their bundled UCX cannot use
the InfiniBand cards. Stack 3 exists for that reason.

Where the stack is selected: the `gpuhx1` block of `tools/build_executable.sh`
(MPI wrapper and library prefixes), the `hx1:gpu` block of
`tools/hpc_execute.sh` (runtime `mpirun` and UCX environment) and the GPU
launcher lines of `tests/hpc_run_all_tests.sh`. All three default to stack 3
in `tools/GPU_dependency/stack3-OpenMPI-4.1.6-tm/build`, found relative to
the checkout, so after `install_dependencies.sh` has built stack 3 nothing
needs editing. `UDALES_GPU_STACK_DIR` points them at a build directory
elsewhere. All three keep the lines of the other stacks as comments so a
stack can be switched by moving the comment markers, always in all three
files together; the install log of a stack 1 or 2 build ends with the exact
lines to change.

## Building a stack with one command

**Build stack 3.** It is the recommended stack and the one the uDALES build
and launch scripts use: it is the only one that runs across nodes, and it is
also the fastest on a single node. Try stack 2, then stack 1, only if stack 3
does not build on your system (they need no OpenMPI build of their own, since
they use the OpenMPI shipped inside the NVHPC module).

On Imperial HX1 (submits a PBS job, loads the modules itself):

```bash
cd u-dales
tools/GPU_dependency/install_dependencies.sh                  # stack 3, 8 cores, 2 h walltime
tools/GPU_dependency/install_dependencies.sh --jobs 16        # more cores for the compilers and the job
tools/GPU_dependency/install_dependencies.sh --on-gpu-queue   # while the CPU queues are reserved (until 2026-12-15)
tools/GPU_dependency/install_dependencies.sh --dry-run        # print the PBS script, submit nothing
```

`hx1` is the default system, so it need not be typed. The job asks for
`select=1:ncpus=<jobs>:mem=32gb` and a walltime of two hours for stack 3, one
hour for stacks 1 and 2 (`--mem`, `--walltime`, or `--select` for the whole
string, change that); the short request is placed sooner. A build takes
about 45 minutes on 8 cores, 31 of them for OpenMPI.

On another machine with the NVHPC module (runs in the background with
`nohup`, so it survives closing the terminal):

```bash
cd u-dales
tools/GPU_dependency/install_dependencies.sh common               # stack 3, 8 cores
tools/GPU_dependency/install_dependencies.sh common --jobs 16
```

Only if stack 3 fails to build:

```bash
tools/GPU_dependency/install_dependencies.sh --stack 2            # add common for a non-HX1 machine
tools/GPU_dependency/install_dependencies.sh --stack 1
```

All module names are in one file, `tools/GPU_dependency/modules.sh`, chosen
by the `hx1`/`common` argument. With `common` the script loads the NVIDIA HPC
SDK the way its installer sets it up,

```bash
module use /opt/nvidia/hpc_sdk/modulefiles
module load nvhpc/24.11
```

and takes `cmake` and the zlib headers from the system (`cmake` and
`zlib1g-dev`/`zlib-devel` packages). If your SDK is another version or lives
elsewhere, change `NVHPC_MODULE` / `NVHPC_MODULEFILES` at the top of
`modules.sh`, or export them before running. The modules are loaded once
before anything is submitted or started: if the module is missing the script
stops right there and prints which `nvhpc` modules exist on the machine, the
file and line to change if one of them will do, and where to get the SDK if
none exists. The stack 3 OpenMPI is configured on a `common` machine with
IPv6 and CUDA, PBS launch only if PBS is installed, UCX only if `UCX_ROOT` is
set, and OpenMPI's internal PMIx, hwloc and libevent.

Either way the script returns at once and prints the command to follow the
log:

```
Submitted job 736884.pbs-6  (select=1:ncpus=8:mem=32gb:ngpus=1:gpu_type=A100, walltime 02:00:00)
PBS script: .../stack3-OpenMPI-4.1.6-tm/build/logs/install_20261125_105612.pbs

The log appears once the job starts. Follow it with:
    tail -f .../stack3-OpenMPI-4.1.6-tm/build/logs/install_20261125_105612.log
Job state:  qstat 736884.pbs-6
It ends with INSTALL COMPLETE or INSTALL FAILED.
```

### What the script does

`install_dependencies.sh` runs every stage of a stack in order (for stack 3
the OpenMPI first) as one job that outlives the terminal: on `hx1` a PBS job,
on `common` a `nohup` process. The stage scripts are copied into
`<stack>/build/` and run there, so sources, installed libraries and logs all
end up in that directory, which is git-ignored. Components already installed
there are skipped, so after a failure the same command resumes the build;
`--force` rebuilds everything. The log (`build/logs/install_<time>.log`) ends
with `INSTALL COMPLETE` or `INSTALL FAILED`. After a stack 3 build it says that
nothing needs editing; after a stack 1 or 2 build it lists, file by file, the
variables to change in `tools/build_executable.sh`, `tools/hpc_execute.sh` and
`tests/hpc_run_all_tests.sh`.

On `hx1` the tarballs are downloaded on the submitting node before the job is
queued, so the compute node needs no internet. Downloads check certificates;
the RHEL 8 CA store on HX1 lacks the root of `downloads.unidata.ucar.edu`'s
chain, so `download_helper.sh` retries with a newer bundle (`UDALES_CA_BUNDLE`
if set, otherwise the certifi bundle of `python3` or of `tools/python/.venv`).
