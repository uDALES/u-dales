#!/usr/bin/env bash
set -euo pipefail

workdir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
cd "$workdir/hdf5-2.1.1"

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

export CFLAGS="-fPIC"
export CXXFLAGS="-fPIC"
export FFLAGS="-fPIC"
export FCFLAGS="-fPIC"

module list
echo "$DEP_ZLIB_ROOT"
ls "$DEP_ZLIB_INCDIR/zlib.h"

"$CC" --showme:command
"$CXX" --showme:command
"$FC" --showme:command

echo "CC: $CC"
echo "CXX: $CXX"
echo "FC: $FC"
echo "CFLAGS: $CFLAGS"
echo "CXXFLAGS: $CXXFLAGS"
echo "FFLAGS: $FFLAGS"
echo "FCFLAGS: $FCFLAGS"
