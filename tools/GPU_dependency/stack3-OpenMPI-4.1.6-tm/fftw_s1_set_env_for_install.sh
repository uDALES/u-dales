#!/usr/bin/env bash
set -euo pipefail

workdir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)

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
export NVHPC_OMP_LIB="$DEP_NVHPC_HOME/compilers/lib/libnvomp.so"

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
