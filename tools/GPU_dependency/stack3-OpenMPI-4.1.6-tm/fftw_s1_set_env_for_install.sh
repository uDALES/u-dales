#!/usr/bin/env bash
set -euo pipefail

workdir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)

module purge
module load NVHPC/23.7-CUDA-12.2.0

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
export NVHPC_OMP_LIB="$EBROOTNVHPC/Linux_x86_64/23.7/compilers/lib/libnvomp.so"

export CC="$NVHPC_MPI/mpicc"
export CXX="$NVHPC_MPI/mpicxx"
export FC="$NVHPC_MPI/mpif90"
export F77="$FC"
export MPICC="$CC"

export CFLAGS="-fPIC"
export CXXFLAGS="-fPIC"
export FFLAGS="-fPIC"
export FCFLAGS="-fPIC"

if [[ -n "${LD_PRELOAD:-}" ]]; then
    export LD_PRELOAD="$NVHPC_OMP_LIB:$LD_PRELOAD"
else
    export LD_PRELOAD="$NVHPC_OMP_LIB"
fi

[[ -x "$CC" ]] || { echo "Error: MPI C compiler not found: $CC" >&2; exit 1; }
[[ -x "$FC" ]] || { echo "Error: MPI Fortran compiler not found: $FC" >&2; exit 1; }
[[ -f "$NVHPC_OMP_LIB" ]] || { echo "Error: NVHPC OpenMP runtime not found: $NVHPC_OMP_LIB" >&2; exit 1; }

module list
"$CC" --showme:command 2>/dev/null
"$FC" --showme:command 2>/dev/null
echo "CC: $CC"
echo "FC: $FC"
echo "MPICC: $MPICC"
echo "NVHPC OpenMP runtime: $NVHPC_OMP_LIB"
