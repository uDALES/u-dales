#!/usr/bin/env bash
# Environment for building OpenMPI 4.1.6 with the NVHPC compilers. On HX1 it
# mirrors the site's own OpenMPI/4.1.4-NVHPC-22.7-CUDA-11.7.0 recipe
# (--enable-ipv6, --with-tm, site UCX + UCX-CUDA, PMIx, hwloc, libevent) for
# NVHPC 23.7 / CUDA 12.2. Why not the OpenMPI bundled with NVHPC: it lacks IPv6
# (HX1's node networks are IPv6-only, so its launcher daemons cannot connect
# back to mpirun) and its UCX InfiniBand module was built against a different
# rdma-core than the system's and does not load.
# The module names and the UCX/PMIx/hwloc/libevent/PBS locations come from
# tools/GPU_dependency/modules.sh (UDALES_DEP_SYSTEM: hx1 or common). On a
# common machine whatever is not there is left out of the configure line and
# OpenMPI uses its internal copies.
set -euo pipefail
workdir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
d=$workdir; while [ ! -f "$d/modules.sh" ] && [ "$d" != / ]; do d=$(dirname "$d"); done
source "$d/modules.sh" openmpi
export NVHPC_BIN="$DEP_NVHPC_HOME/compilers/bin"
export CC="$NVHPC_BIN/nvc" CXX="$NVHPC_BIN/nvc++" FC="$NVHPC_BIN/nvfortran"
export OMPI_PREFIX="$workdir/openmpi"
export CUDA_ROOT_FOR_OMPI="$DEP_CUDA_HOME"
module list 2>&1 | grep -E '^ *[0-9]+\)' | tr -s ' ' | tr '\n' ' ' || true; echo
echo "CC=$CC FC=$FC CUDA=$CUDA_ROOT_FOR_OMPI"
echo "TM=${DEP_TM_ROOT:-none} UCX=${DEP_UCX_ROOT:-none} PMIX=${DEP_PMIX_ROOT:-internal} HWLOC=${DEP_HWLOC_ROOT:-internal} LIBEVENT=${DEP_LIBEVENT_ROOT:-internal}"
