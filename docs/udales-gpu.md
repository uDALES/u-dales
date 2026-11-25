# uDALES on GPUs

uDALES has a GPU port written with OpenACC and CUDA Fortran, compiled with the
NVIDIA HPC SDK (NVHPC). It reads the same inputs and writes the same outputs
as the CPU code; the test suite compares the two to a tolerance of 1e-6 (see
[Testing and CI](./udales-testing.md)). Each MPI rank drives one GPU, so a GPU
run with `NNODE` nodes and `NGPU` GPUs per node needs
`nprocx * nprocy = NNODE * NGPU` in `namoptions`.

Supported systems:

| `<system>` of `tools/build_executable.sh` | Machine | Notes |
|---|---|---|
| `gpuhx1` | Imperial HX1, A100 nodes (`v1_a100` queue) | the reference GPU system; single and multi-node |
| `gpu` | a Linux workstation with an NVIDIA GPU and NVHPC | single node; library paths set in the script |

Getting a GPU build takes three steps: build the dependency libraries with
NVHPC, build the executable against them, and run with `NGPU` set in
`config.sh`.

## 1. Build the dependencies with NVHPC

The GPU executable needs HDF5, netCDF-C, netCDF-Fortran and FFTW compiled with
the NVHPC compilers, all against one and the same MPI, and the job has to
launch with that MPI's `mpirun`. The Intel- or GCC-built netCDF and FFTW
modules of a cluster cannot be mixed in, so these libraries are built from
source once per machine. The scripts for this live in
`tools/GPU_dependency/` (details in its
[README](https://github.com/uDALES/u-dales/blob/master/tools/GPU_dependency/README.md)).

Three MPI "stacks" are defined there:

| Stack | MPI | Use |
|---|---|---|
| 1 | OpenMPI 3.1.5 bundled with NVHPC | single node only |
| 2 | OpenMPI 4.1.5 bundled with NVHPC (UCX, CUDA) | single node only; faster than 1 |
| 3 | OpenMPI 4.1.6 built by the script with IPv6, PBS launch and UCX-CUDA | single **and multi-node**; the default |

**Build stack 3.** It is the only one that runs across HX1 nodes (the bundled
OpenMPIs lack IPv6 and HX1's node networks are IPv6-only) and the fastest on
one node. Try stack 2, then 1, only if stack 3 does not build on your machine.

One command builds a whole stack. On HX1 it is submitted as a PBS job that
survives closing the terminal:

```sh
# From the u-dales directory
tools/GPU_dependency/install_dependencies.sh                  # stack 3, 8 cores, 2 h walltime
tools/GPU_dependency/install_dependencies.sh --jobs 16        # more cores
tools/GPU_dependency/install_dependencies.sh --dry-run        # show the PBS script, submit nothing
```

On another Linux machine with the NVIDIA HPC SDK, the same script runs the
build in the background with `nohup`:

```sh
tools/GPU_dependency/install_dependencies.sh common            # stack 3, 8 cores
tools/GPU_dependency/install_dependencies.sh common --jobs 16
```

Only if stack 3 fails to build:

```sh
tools/GPU_dependency/install_dependencies.sh --stack 2         # add "common" off HX1
tools/GPU_dependency/install_dependencies.sh --stack 1
```

The script prints the `tail -f` command for its log and returns at once. The
log ends with `INSTALL COMPLETE` or `INSTALL FAILED`; after a failure the same
command resumes where it stopped (finished components are skipped; `--force`
rebuilds them). A stack 3 build takes about 45 minutes on 8 cores, two
thirds of it for OpenMPI.

Everything is installed under `tools/GPU_dependency/<stack>/build/`, which is
ignored by git. `tools/build_executable.sh`, `tools/hpc_execute.sh` and the
test scripts look for stack 3 in exactly that place, relative to the checkout,
so after a stack 3 build nothing needs editing. A build somewhere else is
reached through the environment variable `UDALES_GPU_STACK_DIR`. After a
stack 1 or 2 build the log lists, file by file, the lines to change.

Module names are kept in one file, `tools/GPU_dependency/modules.sh`. On HX1
it loads `NVHPC/23.7-CUDA-12.2.0` with CMake and zlib (and for the stack 3
OpenMPI the site UCX-CUDA, PMIx, hwloc and libevent). With `common` it loads
the SDK the way NVIDIA's installer sets it up,

```sh
module use /opt/nvidia/hpc_sdk/modulefiles
module load nvhpc/24.11
```

and takes `cmake` and the zlib headers from the system (install the `cmake`
and `zlib1g-dev` / `zlib-devel` packages). For another SDK version or
location change `NVHPC_MODULE` / `NVHPC_MODULEFILES` at the top of
`modules.sh`, or export them before running. The modules are loaded once
before anything is submitted; if the module is missing the script stops right
there, lists the `nvhpc` modules that do exist, and names the file and line to
change, or where to download the SDK if there is none.

## 2. Build the executable

On HX1:

```sh
# From the u-dales directory
tools/build_executable.sh gpuhx1 release      # -> build/gpu/release/u-dales
tools/build_executable.sh gpuhx1 debug        # -> build/gpu/debug/u-dales
```

The script loads NVHPC, takes MPI and the libraries from the stack 3 build
directory, and stops with a message pointing at `install_dependencies.sh` if
that directory is missing. GPU and CPU builds keep separate CMake caches
(`build/gpu/` and `build/cpu/`), so both can exist side by side.

The HX1 executable is built for compute capability 8.0 (A100) without PTX;
for other GPUs pass `-DUDALES_CUDA_ARCH=<cc>` as described in
[Installation](./udales-installation.md#build-on-hpcs).

On a workstation use the `gpu` target. Its block in
`tools/build_executable.sh` names the MPI wrapper and the netCDF and FFTW
prefixes directly; point them at the `netcdf-c-*/netcdfc`,
`netcdf-fortran-*/netcdff` and `fftw-*/fftw3` directories under the stack's
`build/` and set `UDALES_CUDA_ARCH` to your card (86 for an RTX 30 series,
89 for 40 series, 80 for A100):

```sh
tools/build_executable.sh gpu release
```

## 3. Run

A GPU run is a normal `ud_run` with a GPU executable and `NGPU` in
`config.sh`. On HX1:

```sh
# Contents of the config.sh file
# We assume you are running the following commands from your
# top-level project directory.

export DA_EXPDIR=$(pwd)/experiments                      # Experiments top-level directory
export DA_TOOLSDIR=$(pwd)/u-dales/tools                  # Directory of scripts
export DA_BUILD=$(pwd)/u-dales/build/gpu/release/u-dales # GPU build
export DA_WORKDIR=$EPHEMERAL                             # Output top-level directory
export NNODE=1                                           # Nodes
export NGPU=4                                            # GPUs per node: makes this a GPU run
export NCPU=4                                            # CPU cores per node, one per GPU is enough
export WALLTIME="01:00:00"                               # hours:minutes:seconds
export MEM="256gb"                                       # Memory request per node
# export PLACE=scatter                                   # multi-node: one chunk per node
# export GPU_TYPE=A100                                   # default; "" drops the constraint
```

with `nprocx * nprocy = NNODE * NGPU` in `namoptions` (here 4, for example
`nprocx = 1, nprocy = 4`), then

```sh
./u-dales/bin/ud_run hx1 sim experiments/009
./u-dales/bin/ud_run hx1 gather experiments/009
```

`tools/hpc_execute.sh` detects the GPU build from `DA_BUILD`, requests
`ngpus=NGPU` A100s per node, loads NVHPC and puts the stack 3 `mpirun` and
the UCX CUDA transports in the job environment, and launches
`mpirun --map-by ppr:NGPU:node` through `tools/bind.sh`, which gives each rank
its own GPU. Nothing has to be loaded by hand.

Multi-node runs (`NNODE > 1`) need stack 3; the bundled OpenMPIs cannot
launch across HX1 nodes. HX1 allows a user 12 running GPUs at a time, so the
largest job is two nodes with four A100s each. Measured on HX1: a 512^3 grid
runs on four GPUs of one node and gets slower, not faster, when spread over
two nodes (the 2DECOMP transposes cross the InfiniBand link), so use more
than one node only for cases that do not fit on one. A 1024^3 grid needs the
eight GPUs of two nodes (75 GB of the 80 GB per card), with
`nprocx = 1, nprocy = 8`, `NNODE=2`, `NGPU=4` and `PLACE=scatter`. For any
layout prefer `nprocx = 1`.

On a workstation, `ud_run common sim` (`tools/local_execute.sh`) recognises a
GPU build and launches it through `bind.sh` as well; set `NCPU` to the number
of GPUs and make sure the `mpiexec` on `PATH` is the one the stack was built
with (`MPIEXEC` overrides it).

## Checking a GPU build

`tests/hpc_run_all_tests.sh gpu` builds the CPU and GPU executables from clean
and runs the `gpu-smoke` and `gpu-nightly` parity selections, which compare
GPU and CPU results case by case. See [Testing and CI](./udales-testing.md)
and `tests/integration/gpu/README.md` for the case matrix and for running the
harness by hand.
