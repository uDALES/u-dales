#!/usr/bin/env bash
# Environment for building OpenMPI 4.1.6 with the NVHPC 23.7 compilers for
# multi-node GPU runs on HX1. Mirrors the site's own
# OpenMPI/4.1.4-NVHPC-22.7-CUDA-11.7.0 recipe (--enable-ipv6, --with-tm,
# site UCX + UCX-CUDA, PMIx, hwloc, libevent) but for NVHPC 23.7 / CUDA 12.2.
# Why not the OpenMPI bundled with NVHPC: it lacks IPv6 (HX1's node networks
# are IPv6-only, so its launcher daemons cannot connect back to mpirun) and
# its UCX InfiniBand module was built against a different rdma-core than the
# system's and does not load.
set -euo pipefail
workdir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
module purge
module load NVHPC/23.7-CUDA-12.2.0
module load UCX-CUDA/1.14.1-GCCcore-12.3.0-CUDA-12.1.1   # pulls UCX/1.14.1, GDRCopy/2.3.1, CUDA/12.1.1
module load PMIx/4.2.4-GCCcore-12.3.0 hwloc/2.9.1-GCCcore-12.3.0 libevent/2.1.12-GCCcore-12.3.0
module load zlib/1.2.13-GCCcore-12.3.0
export NVHPC_BIN="$EBROOTNVHPC/Linux_x86_64/23.7/compilers/bin"
export CC="$NVHPC_BIN/nvc" CXX="$NVHPC_BIN/nvc++" FC="$NVHPC_BIN/nvfortran"
export OMPI_PREFIX="$workdir/openmpi"
export CUDA_ROOT_FOR_OMPI="$EBROOTNVHPC/Linux_x86_64/23.7/cuda/12.2"
module list 2>&1 | grep -E '^ *[0-9]+\)' | tr -s ' ' | tr '\n' ' '; echo
echo "CC=$CC FC=$FC UCX=$EBROOTUCX PMIX=$EBROOTPMIX HWLOC=$EBROOTHWLOC LIBEVENT=$EBROOTLIBEVENT CUDA=$CUDA_ROOT_FOR_OMPI"
