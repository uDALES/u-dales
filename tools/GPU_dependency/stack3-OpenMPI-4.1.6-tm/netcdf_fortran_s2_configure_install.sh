#!/usr/bin/env bash
set -euo pipefail

workdir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
source "$workdir/netcdf_fortran_s1_set_env_for_install.sh"

cd "$workdir/netcdf-fortran-4.6.4"

./configure \
    --prefix="$PWD/netcdff" \
    --enable-parallel-tests \
    --enable-shared \
    --enable-static \
    --disable-zstandard-plugin \
    CC="$CC" \
    FC="$FC" \
    F77="$F77" \
    CFLAGS="$CFLAGS" \
    FFLAGS="$FFLAGS" \
    FCFLAGS="$FCFLAGS" \
    CPPFLAGS="$CPPFLAGS" \
    LDFLAGS="$LDFLAGS"
