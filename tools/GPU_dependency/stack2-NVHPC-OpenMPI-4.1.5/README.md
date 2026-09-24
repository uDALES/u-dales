# Stack 2: NVHPC 23.7 with its bundled OpenMPI 4.1.5

Same libraries and scripts as stack 1, compiled against the OpenMPI 4.1.5 that
NVHPC 23.7 also ships (`comm_libs/12.2/openmpi4/openmpi-4.1.5`), which has
UCX and CUDA support. About 20 % faster per step than stack 1 on one node and
free of the openib warnings.

**Scope: single-node GPU runs only.** Like stack 1 it is built without IPv6,
so it cannot launch across HX1 nodes, and its bundled UCX cannot load its
InfiniBand module against the system rdma-core. For multi-node runs use
stack 3.

The only difference from the stack 1 scripts is the `NVHPC_MPI` line of each
`s1` script (the stack 1 line is kept above it as a comment).

## What gets built

| Library | Version | Notes |
|---|---:|---|
| HDF5 | 2.1.1 | CMake build, parallel, Fortran, zlib/DEFLATE |
| netCDF-C | 4.10.1 | netCDF-4, parallel I/O, no DAP |
| netCDF-Fortran | 4.6.4 | against the netCDF-C above |
| FFTW | 3.3.11 | MPI and OpenMP, double precision |

Order: HDF5, netCDF-C, netCDF-Fortran; FFTW is independent. Each library has
four stages, `s0` download, `s1` environment, `s2` configure, `s3` build and
install. The `s2`/`s3` scripts of netCDF-C, netCDF-Fortran and FFTW source
their `s1` themselves; HDF5's `s3` does not, so source `hdf5_s1` first.

## How to run

Copy this directory to where the libraries should live (the scripts place
sources and installs beside themselves), then:

```bash
./hdf5_s0_download.sh
source ./hdf5_s1_set_env_for_install.sh
./hdf5_s2_configure_install.sh
./hdf5_s3_install.sh

./netcdf_c_s0_download.sh
./netcdf_c_s2_configure_install.sh
./netcdf_c_s3_build_install.sh

./netcdf_fortran_s0_download.sh
./netcdf_fortran_s2_configure_install.sh
./netcdf_fortran_s3_build_install.sh

./fftw_s0_download.sh
./fftw_s2_configure_install.sh
./fftw_s3_build_install.sh
```

`BUILD_JOBS=16` raises the parallelism of the `s3` steps. If `wget` is
missing, fetch the archives with `curl -fL -o <archive> <url>` and untar them
by hand. Installed prefixes: `hdf5-2.1.1/hdf5`, `netcdf-c-4.10.1/netcdfc`,
`netcdf-fortran-4.6.4/netcdff`, `fftw-3.3.11/fftw3`; verify with
`nc-config --has-parallel4` (yes) and `grep "Parallel HDF5" hdf5-2.1.1/hdf5/lib/libhdf5.settings`.

## Using the result

Point the `gpuhx1` block of `tools/build_executable.sh` at the installed
prefixes (`NETCDF_DIR`, `NETCDF_FORTRAN_DIR`, `FFTWDIR`) and the matching MPI
wrapper directory, and the `hx1:gpu` block of `tools/hpc_execute.sh` at the
matching `mpirun`. The three must belong to the same stack: the executable,
2DECOMP&FFT, these libraries and the launcher all share one `libmpi`.
