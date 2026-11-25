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
# export NVHPC_MPI="$DEP_MPI_OPENMPI4/bin"
# Stack 3 - OpenMPI 4.1.6 built for NVHPC 23.7 with IPv6, PBS (tm) launch and
# the site UCX; the NVHPC-bundled OpenMPIs cannot run across HX1 nodes.
# OMPI_PREFIX is set by install_dependencies.sh to the OpenMPI it built in
# build/openmpi/openmpi; without it the openmpi/openmpi next to this script
# (the same place when the scripts run inside build/) is used.
export NVHPC_MPI="${OMPI_PREFIX:-$workdir/openmpi/openmpi}/bin"
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
