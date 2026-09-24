#!/usr/bin/env bash
set -euo pipefail

workdir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
export HDF5_ROOT="$workdir/hdf5-2.1.1/hdf5"

module purge
module load NVHPC/23.7-CUDA-12.2.0
module load zlib/1.2.13-GCCcore-12.3.0

export NVHPC_MPI="$EBROOTNVHPC/Linux_x86_64/23.7/comm_libs/mpi/bin"
export CC="$NVHPC_MPI/mpicc"
export CXX="$NVHPC_MPI/mpicxx"
export FC="$NVHPC_MPI/mpif90"

export CFLAGS="-fPIC"
export CXXFLAGS="-fPIC"
export FFLAGS="-fPIC"
export FCFLAGS="-fPIC"

export CPPFLAGS="-I$HDF5_ROOT/include -I$EBROOTZLIB/include"
export LDFLAGS="-L$HDF5_ROOT/lib -L$EBROOTZLIB/lib -Wl,-rpath,$HDF5_ROOT/lib -Wl,-rpath,$EBROOTZLIB/lib"
export LD_LIBRARY_PATH="$HDF5_ROOT/lib:$EBROOTZLIB/lib:${LD_LIBRARY_PATH:-}"
export PATH="$HDF5_ROOT/bin:$PATH"

[[ -x "$CC" ]] || { echo "Error: MPI C compiler not found: $CC" >&2; exit 1; }
[[ -f "$HDF5_ROOT/include/hdf5.h" ]] || { echo "Error: HDF5 not found at $HDF5_ROOT" >&2; exit 1; }

if ! "$HDF5_ROOT/bin/h5pcc" -showconfig | grep -qE 'Parallel HDF5:[[:space:]]+(ON|yes)'; then
    echo "Error: HDF5 at $HDF5_ROOT does not have parallel support." >&2
    exit 1
fi

if ! grep -qi 'I/O filters.*DEFLATE' "$HDF5_ROOT/lib/libhdf5.settings"; then
    echo "Error: HDF5 at $HDF5_ROOT does not have zlib/deflate support." >&2
    exit 1
fi

module list
"$CC" --showme:command
echo "HDF5_ROOT: $HDF5_ROOT"
echo "CC: $CC"
echo "CPPFLAGS: $CPPFLAGS"
echo "LDFLAGS: $LDFLAGS"
