#!/usr/bin/env bash
set -euo pipefail
workdir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
source "$workdir/s1_set_env.sh"
cd "$workdir/openmpi-4.1.6"
# -fPIC: nvfortran otherwise leaves one mpi_f08 object non-relocatable and the
# shared bindings library fails to link. LDFLAGS rpath: --with-tm makes libmpi
# and the launcher need libpbs.so.0, and /opt/pbs/lib is not on the system
# loader path, so record it in the RPATH rather than depend on LD_LIBRARY_PATH.
./configure --prefix="$OMPI_PREFIX" \
    CC="$CC" CXX="$CXX" FC="$FC" \
    CFLAGS="-fPIC" CXXFLAGS="-fPIC" FCFLAGS="-fPIC" FFLAGS="-fPIC" \
    LDFLAGS="-Wl,-rpath,/opt/pbs/lib" \
    --enable-ipv6 \
    --with-tm=/opt/pbs \
    --with-ucx="$EBROOTUCX" \
    --with-cuda="$CUDA_ROOT_FOR_OMPI" \
    --with-pmix="$EBROOTPMIX" \
    --with-hwloc="$EBROOTHWLOC" \
    --with-libevent="$EBROOTLIBEVENT" \
    --enable-mpirun-prefix-by-default \
    --enable-shared --enable-static \
    --without-verbs --disable-io-romio \
    > "$workdir/logs/configure.log" 2>&1
tail -3 "$workdir/logs/configure.log"
# head closing the pipe gives grep SIGPIPE, which pipefail would report as 141
grep -E 'IPv6|tm|ucx|CUDA support|PMIx|hwloc' "$workdir/logs/configure.log" | grep -iE 'yes|no|internal|external' | head -12 || true
