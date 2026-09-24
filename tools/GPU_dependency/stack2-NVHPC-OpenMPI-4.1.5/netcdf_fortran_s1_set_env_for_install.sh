#!/usr/bin/env bash
set -euo pipefail

workdir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
export HDF5_ROOT="$workdir/hdf5-2.1.1/hdf5"
export NETCDF_ROOT="$workdir/netcdf-c-4.10.1/netcdfc"

module purge
module load NVHPC/23.7-CUDA-12.2.0
module load zlib/1.2.13-GCCcore-12.3.0

# OpenMPI 3.1.5 (NVHPC default, no UCX, no working InfiniBand path on HX1):
# export NVHPC_MPI="$EBROOTNVHPC/Linux_x86_64/23.7/comm_libs/mpi/bin"
# OpenMPI 4.1.5 with UCX and CUDA support, also shipped with NVHPC 23.7:
export NVHPC_MPI="$EBROOTNVHPC/Linux_x86_64/23.7/comm_libs/12.2/openmpi4/openmpi-4.1.5/bin"
export CC="$NVHPC_MPI/mpicc"
export CXX="$NVHPC_MPI/mpicxx"
export FC="$NVHPC_MPI/mpif90"
export F77="$FC"

export CFLAGS="-fPIC"
export CXXFLAGS="-fPIC"
export FFLAGS="-fPIC"
export FCFLAGS="-fPIC"

export CPPFLAGS="-I$NETCDF_ROOT/include -I$HDF5_ROOT/include -I$EBROOTZLIB/include"
export LDFLAGS="-L$NETCDF_ROOT/lib -L$HDF5_ROOT/lib -L$EBROOTZLIB/lib -Wl,-rpath,$NETCDF_ROOT/lib -Wl,-rpath,$HDF5_ROOT/lib -Wl,-rpath,$EBROOTZLIB/lib"
export LD_LIBRARY_PATH="$NETCDF_ROOT/lib:$HDF5_ROOT/lib:$EBROOTZLIB/lib:${LD_LIBRARY_PATH:-}"
export PATH="$NETCDF_ROOT/bin:$HDF5_ROOT/bin:$PATH"
export PKG_CONFIG_PATH="$NETCDF_ROOT/lib/pkgconfig:$HDF5_ROOT/lib/pkgconfig:${PKG_CONFIG_PATH:-}"

[[ -x "$FC" ]] || { echo "Error: MPI Fortran compiler not found: $FC" >&2; exit 1; }
[[ -x "$NETCDF_ROOT/bin/nc-config" ]] || { echo "Error: netCDF-C not found at $NETCDF_ROOT" >&2; exit 1; }

if [[ "$("$NETCDF_ROOT/bin/nc-config" --has-nc4)" != yes ]]; then
    echo "Error: netCDF-C at $NETCDF_ROOT does not have netCDF-4 support." >&2
    exit 1
fi

if [[ "$("$NETCDF_ROOT/bin/nc-config" --has-parallel4)" != yes ]]; then
    echo "Error: netCDF-C at $NETCDF_ROOT does not have parallel netCDF-4 support." >&2
    exit 1
fi

module list
"$CC" --showme:command 2>/dev/null
"$FC" --showme:command 2>/dev/null
"$NETCDF_ROOT/bin/nc-config" --version
echo "NETCDF_ROOT: $NETCDF_ROOT"
echo "CC: $CC"
echo "FC: $FC"
echo "CPPFLAGS: $CPPFLAGS"
echo "LDFLAGS: $LDFLAGS"
