#!/usr/bin/env bash
set -euo pipefail

workdir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)

module purge
module load NVHPC/23.7-CUDA-12.2.0

export NVHPC_MPI="$EBROOTNVHPC/Linux_x86_64/23.7/comm_libs/mpi/bin"
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
