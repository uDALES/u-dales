#!/usr/bin/env bash
set -euo pipefail

workdir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
cd "$workdir/hdf5-2.1.1"

module purge
module load NVHPC/23.7-CUDA-12.2.0
module load CMake/3.26.3-GCCcore-12.3.0
module load zlib/1.2.13-GCCcore-12.3.0

export NVHPC_MPI="$EBROOTNVHPC/Linux_x86_64/23.7/comm_libs/mpi/bin"

export CC="$NVHPC_MPI/mpicc"
export CXX="$NVHPC_MPI/mpicxx"
export FC="$NVHPC_MPI/mpif90"

export CFLAGS="-fPIC"
export CXXFLAGS="-fPIC"
export FFLAGS="-fPIC"
export FCFLAGS="-fPIC"

module list
echo "$EBROOTZLIB"
ls "$EBROOTZLIB/include/zlib.h"

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
