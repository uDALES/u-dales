#!/usr/bin/env bash
set -euo pipefail

workdir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
export HDF5_ROOT="$workdir/hdf5-2.1.1/hdf5"
export NETCDF_ROOT="$workdir/netcdf-c-4.10.1/netcdfc"

# Modules and SDK paths for this machine (UDALES_DEP_SYSTEM: hx1 or common)
# come from tools/GPU_dependency/modules.sh, found by walking up from here.
d=$workdir; while [ ! -f "$d/modules.sh" ] && [ "$d" != / ]; do d=$(dirname "$d"); done
source "$d/modules.sh"

# OpenMPI 3.1.5 (NVHPC default, no UCX, no working InfiniBand path on HX1):
# export NVHPC_MPI="$DEP_MPI_DEFAULT/bin"
# OpenMPI 4.1.5 with UCX and CUDA support, also shipped with NVHPC 23.7:
export NVHPC_MPI="$DEP_MPI_OPENMPI4/bin"
export CC="$NVHPC_MPI/mpicc"
export CXX="$NVHPC_MPI/mpicxx"
export FC="$NVHPC_MPI/mpif90"
export F77="$FC"

export CFLAGS="-fPIC"
export CXXFLAGS="-fPIC"
export FFLAGS="-fPIC"
export FCFLAGS="-fPIC"

export CPPFLAGS="-I$NETCDF_ROOT/include -I$HDF5_ROOT/include -I$DEP_ZLIB_INCDIR"
export LDFLAGS="-L$NETCDF_ROOT/lib -L$HDF5_ROOT/lib -L$DEP_ZLIB_LIBDIR -Wl,-rpath,$NETCDF_ROOT/lib -Wl,-rpath,$HDF5_ROOT/lib -Wl,-rpath,$DEP_ZLIB_LIBDIR"
export LD_LIBRARY_PATH="$NETCDF_ROOT/lib:$HDF5_ROOT/lib:$DEP_ZLIB_LIBDIR:${LD_LIBRARY_PATH:-}"
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
