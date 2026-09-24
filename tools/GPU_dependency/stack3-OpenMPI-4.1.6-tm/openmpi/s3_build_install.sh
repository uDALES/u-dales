#!/usr/bin/env bash
set -euo pipefail
workdir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
source "$workdir/s1_set_env.sh"
cd "$workdir/openmpi-4.1.6"
make -j "${BUILD_JOBS:-8}" > "$workdir/logs/make.log" 2>&1
make install > "$workdir/logs/install.log" 2>&1
"$OMPI_PREFIX/bin/ompi_info" | grep -E '^ *Open MPI:|MPI extensions|Fort mpi_f08'
"$OMPI_PREFIX/bin/ompi_info" --parsable | grep -oE "enable-ipv6|with-tm[^ ']*|with-ucx[^ ']*|with-cuda[^ ']*" | sort -u | tr '\n' ' '; echo
"$OMPI_PREFIX/bin/ompi_info" --parsable | grep -E '^mca:(plm|ras):tm:|^mca:pml:ucx:|^mca:btl:smcuda:' | cut -d: -f1-3 | sort -u | tr '\n' ' '; echo
