#!/usr/bin/env bash
set -euo pipefail

workdir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
cd "$workdir/hdf5-2.1.1"

module purge
module load NVHPC/23.7-CUDA-12.2.0
module load CMake/3.26.3-GCCcore-12.3.0
module load zlib/1.2.13-GCCcore-12.3.0

# OpenMPI 3.1.5 (NVHPC default, no UCX, no working InfiniBand path on HX1):
# export NVHPC_MPI="$EBROOTNVHPC/Linux_x86_64/23.7/comm_libs/mpi/bin"
# OpenMPI 4.1.5 with UCX and CUDA support, also shipped with NVHPC 23.7:
# export NVHPC_MPI="$EBROOTNVHPC/Linux_x86_64/23.7/comm_libs/12.2/openmpi4/openmpi-4.1.5/bin"
# Stack 3 - OpenMPI 4.1.6 built for NVHPC 23.7 with IPv6, PBS (tm) launch and
# the site UCX; the NVHPC-bundled OpenMPIs cannot run across HX1 nodes.
export NVHPC_MPI="/gpfs/home/dmajumda/openmpi-4.1.6-NVHPC-23.7-CUDA-12.2.0/openmpi/bin"
# libmpi of this OpenMPI needs libpbs.so.0 (PBS tm launch); the RPATH covers it,
# this keeps configure-time test programs safe as well.
export LD_LIBRARY_PATH="/opt/pbs/lib:${LD_LIBRARY_PATH:-}"

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
