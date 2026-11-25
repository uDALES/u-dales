#!/usr/bin/env bash
set -euo pipefail
workdir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
source "$workdir/s1_set_env.sh"
cd "$workdir/openmpi-4.1.6"
# -fPIC: nvfortran otherwise leaves one mpi_f08 object non-relocatable and the
# shared bindings library fails to link. LDFLAGS rpath: --with-tm makes libmpi
# and the launcher need libpbs.so.0, and /opt/pbs/lib is not on the system
# loader path, so record it in the RPATH rather than depend on LD_LIBRARY_PATH.
# The --with-* below are added only where modules.sh found the package.
args=(--prefix="$OMPI_PREFIX"
      CC="$CC" CXX="$CXX" FC="$FC"
      CFLAGS="-fPIC" CXXFLAGS="-fPIC" FCFLAGS="-fPIC" FFLAGS="-fPIC")
[ -n "${DEP_TM_ROOT:-}" ]       && args+=(LDFLAGS="-Wl,-rpath,$DEP_TM_ROOT/lib" --with-tm="$DEP_TM_ROOT")
args+=(--enable-ipv6)
[ -n "${DEP_UCX_ROOT:-}" ]      && args+=(--with-ucx="$DEP_UCX_ROOT")
args+=(--with-cuda="$CUDA_ROOT_FOR_OMPI")
[ -n "${DEP_PMIX_ROOT:-}" ]     && args+=(--with-pmix="$DEP_PMIX_ROOT")
[ -n "${DEP_HWLOC_ROOT:-}" ]    && args+=(--with-hwloc="$DEP_HWLOC_ROOT")
[ -n "${DEP_LIBEVENT_ROOT:-}" ] && args+=(--with-libevent="$DEP_LIBEVENT_ROOT")
args+=(--enable-mpirun-prefix-by-default --enable-shared --enable-static
       --without-verbs --disable-io-romio)
echo "./configure ${args[*]}"
./configure "${args[@]}" > "$workdir/logs/configure.log" 2>&1
tail -3 "$workdir/logs/configure.log"
# head closing the pipe gives grep SIGPIPE, which pipefail would report as 141
grep -E 'IPv6|tm|ucx|CUDA support|PMIx|hwloc' "$workdir/logs/configure.log" | grep -iE 'yes|no|internal|external' | head -12 || true
