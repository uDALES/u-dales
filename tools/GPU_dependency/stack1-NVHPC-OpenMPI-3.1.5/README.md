# Stack 1: NVHPC 23.7 with its bundled OpenMPI 3.1.5

Dependency libraries for the GPU build of uDALES on Imperial HX1, compiled
with the NVIDIA HPC SDK 23.7 (CUDA 12.2) and the OpenMPI 3.1.5 it ships in
`comm_libs/mpi`.

**Scope: single-node GPU runs only.** This OpenMPI is built without IPv6 (HX1's
node networks are IPv6-only, so it cannot launch across nodes) and without UCX.
It also packs non-contiguous GPU datatypes one element at a time, which is why
uDALES exchanges its device halos through packed buffers (issue #377).

Modules used by the `s1` scripts: `NVHPC/23.7-CUDA-12.2.0`,
`zlib/1.2.13-GCCcore-12.3.0`, and `CMake/3.26.3-GCCcore-12.3.0` for HDF5.

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
