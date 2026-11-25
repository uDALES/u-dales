#!/usr/bin/env bash
set -euo pipefail

workdir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
export HDF5_ROOT="$workdir/hdf5-2.1.1/hdf5"

# Modules and SDK paths for this machine (UDALES_DEP_SYSTEM: hx1 or common)
# come from tools/GPU_dependency/modules.sh, found by walking up from here.
d=$workdir; while [ ! -f "$d/modules.sh" ] && [ "$d" != / ]; do d=$(dirname "$d"); done
source "$d/modules.sh"

export NVHPC_MPI="$DEP_MPI_DEFAULT/bin"
export CC="$NVHPC_MPI/mpicc"
export CXX="$NVHPC_MPI/mpicxx"
export FC="$NVHPC_MPI/mpif90"

export CFLAGS="-fPIC"
export CXXFLAGS="-fPIC"
export FFLAGS="-fPIC"
export FCFLAGS="-fPIC"

export CPPFLAGS="-I$HDF5_ROOT/include -I$DEP_ZLIB_INCDIR"
export LDFLAGS="-L$HDF5_ROOT/lib -L$DEP_ZLIB_LIBDIR -Wl,-rpath,$HDF5_ROOT/lib -Wl,-rpath,$DEP_ZLIB_LIBDIR"
export LD_LIBRARY_PATH="$HDF5_ROOT/lib:$DEP_ZLIB_LIBDIR:${LD_LIBRARY_PATH:-}"
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
